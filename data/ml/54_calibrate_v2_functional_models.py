"""Calibrate V2 functional baselines and quantify prediction uncertainty."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import pickle
from collections import Counter
from datetime import date
from pathlib import Path

from sklearn.calibration import CalibratedClassifierCV
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score

try:
    from sklearn.frozen import FrozenEstimator
except ImportError:
    FrozenEstimator = None


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
ML_RESULTS = ROOT / "data" / "ml_results"
ANALYSIS = ROOT / "data" / "analysis_results"
MODEL_DIR = ROOT / "data" / "model"


def load_baseline_helpers():
    path = Path(__file__).with_name("53_train_v2_functional_baselines.py")
    spec = importlib.util.spec_from_file_location("deepgene_v2_baselines", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load baseline feature helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def expected_calibration_error(y_true: list[int], probabilities: list[float], bins: int = 5) -> float:
    total = len(y_true)
    error = 0.0
    for index in range(bins):
        lower = index / bins
        upper = (index + 1) / bins
        selected = [i for i, probability in enumerate(probabilities) if lower <= probability < upper or (index == bins - 1 and probability == upper)]
        if not selected:
            continue
        confidence = sum(probabilities[i] for i in selected) / len(selected)
        observed = sum(y_true[i] for i in selected) / len(selected)
        error += len(selected) / total * abs(confidence - observed)
    return float(error)


def uncertainty(probability: float) -> tuple[str, float]:
    probability = min(max(probability, 1e-12), 1 - 1e-12)
    entropy = -(probability * math.log2(probability) + (1 - probability) * math.log2(1 - probability))
    confidence = max(probability, 1 - probability)
    if confidence >= 0.80:
        level = "high"
    elif confidence >= 0.65:
        level = "moderate"
    else:
        level = "uncertain"
    return level, entropy


def main() -> None:
    base = load_baseline_helpers()
    rows = read_csv(PROCESSED / "deepgene_v2_functional_evidence.csv")
    split_rows = read_csv(ML_RESULTS / "v2_grouped_split_ids.csv")
    split_by_id = {clean(row.get("row_id")): clean(row.get("split")) for row in split_rows if clean(row.get("product")) == "functional_evidence"}
    usable = []
    for row in rows:
        target = base.target(row)
        split = split_by_id.get(clean(row.get("source_row_id")), "")
        if target and split:
            item = dict(row)
            item["target"] = target
            item["split"] = split
            usable.append(item)
    x = [base.features(row) for row in usable]
    y = [row["target"] for row in usable]
    indexes = {name: [i for i, row in enumerate(usable) if row["split"] == name] for name in ("train", "validation", "test")}
    models = {
        "logistic_regression": LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42),
        "random_forest": RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=42, n_jobs=1),
        "gradient_boosting": GradientBoostingClassifier(random_state=42),
        "majority_baseline": DummyClassifier(strategy="most_frequent"),
    }
    report = {
        "created": date.today().isoformat(),
        "status": "RESEARCH_ONLY_CALIBRATION",
        "calibration_protocol": "Base model fit on train; sigmoid calibration fit on validation; final metrics reported on untouched test.",
        "rows": {name: len(values) for name, values in indexes.items()},
        "models": {},
        "limitations": [
            "The validation set is small and source-dependent.",
            "Calibration estimates are exploratory and require larger independent validation.",
            "Uncertainty describes model probability reliability, not disease risk or clinical certainty.",
        ],
    }
    calibrated_models = {}
    prediction_rows = []
    for name, model in models.items():
        train = indexes["train"]
        validation = indexes["validation"]
        test = indexes["test"]
        model.fit([x[i] for i in train], [y[i] for i in train])
        if FrozenEstimator is not None:
            calibrated = CalibratedClassifierCV(FrozenEstimator(model), method="sigmoid")
        else:
            try:
                calibrated = CalibratedClassifierCV(estimator=model, method="sigmoid", cv="prefit")
            except TypeError:
                calibrated = CalibratedClassifierCV(base_estimator=model, method="sigmoid", cv="prefit")
        calibrated.fit([x[i] for i in validation], [y[i] for i in validation])
        calibrated_models[name] = calibrated
        model_report = {}
        for split_name, split_indexes in (("validation", validation), ("test", test)):
            true_labels = [1 if y[i] == "GOF" else 0 for i in split_indexes]
            probabilities = [float(value[list(calibrated.classes_).index("GOF")]) for value in calibrated.predict_proba([x[i] for i in split_indexes])]
            predicted_labels = ["GOF" if probability >= 0.5 else "LOF" for probability in probabilities]
            level_counts = Counter()
            uncertainty_values = []
            for index, probability, predicted in zip(split_indexes, probabilities, predicted_labels):
                level, entropy = uncertainty(probability)
                level_counts[level] += 1
                uncertainty_values.append(entropy)
                prediction_rows.append({
                    "model": name,
                    "source_row_id": usable[index]["source_row_id"],
                    "variant_key": usable[index]["variant_key"],
                    "split": split_name,
                    "observed_target": y[index],
                    "predicted_target": predicted,
                    "gof_probability": f"{probability:.8f}",
                    "confidence": f"{max(probability, 1 - probability):.8f}",
                    "entropy_bits": f"{entropy:.8f}",
                    "uncertainty_level": level,
                })
            model_report[split_name] = {
                "rows": len(split_indexes),
                "class_counts": dict(Counter(y[i] for i in split_indexes)),
                "brier_score": float(brier_score_loss(true_labels, probabilities)),
                "log_loss": float(log_loss(true_labels, probabilities, labels=[0, 1])),
                "expected_calibration_error": expected_calibration_error(true_labels, probabilities),
                "roc_auc": float(roc_auc_score(true_labels, probabilities)) if len(set(true_labels)) == 2 else None,
                "mean_entropy_bits": sum(uncertainty_values) / len(uncertainty_values),
                "uncertainty_level_counts": dict(level_counts),
            }
        report["models"][name] = model_report

    fields = ["model", "source_row_id", "variant_key", "split", "observed_target", "predicted_target", "gof_probability", "confidence", "entropy_bits", "uncertainty_level"]
    write_csv(ML_RESULTS / "v2_uncertainty_predictions.csv", prediction_rows, fields)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    with (MODEL_DIR / "deepgene_v2_calibrated_functional_model.pkl").open("wb") as handle:
        pickle.dump({"models": calibrated_models, "feature_names": ["protein_position", "delta_hydrophobicity", "absolute_delta_hydrophobicity", "delta_charge", "absolute_delta_charge", "delta_volume", "absolute_delta_volume", "same_amino_acid"], "calibration": "sigmoid_on_validation_split"}, handle)
    (ML_RESULTS / "v2_calibration_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    (ANALYSIS / "v2_calibration_report.json").write_text(json.dumps({
        "status": report["status"],
        "models": report["models"],
        "interpretation": "Lower Brier score, log loss, and expected calibration error indicate better probability reliability within this exploratory test.",
    }, indent=2), encoding="utf-8")
    print(json.dumps({
        "rows": {name: len(values) for name, values in indexes.items()},
        "models": list(models),
        "outputs": [
            "data/ml_results/v2_calibration_report.json",
            "data/ml_results/v2_uncertainty_predictions.csv",
            "data/analysis_results/v2_calibration_report.json",
            "data/model/deepgene_v2_calibrated_functional_model.pkl",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
