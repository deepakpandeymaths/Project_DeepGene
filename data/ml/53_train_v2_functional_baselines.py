"""Train transparent V2 functional-effect baselines on frozen grouped splits."""

from __future__ import annotations

import csv
import json
import pickle
from collections import Counter
from datetime import date
from pathlib import Path

from sklearn.dummy import DummyClassifier
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
ML_RESULTS = ROOT / "data" / "ml_results"
ANALYSIS = ROOT / "data" / "analysis_results"
MODEL_DIR = ROOT / "data" / "model"

AA = {
    "A": (1.8, 0.0, 88.6), "R": (-4.5, 1.0, 173.4), "N": (-3.5, 0.0, 114.1),
    "D": (-3.5, -1.0, 111.1), "C": (2.5, 0.0, 108.5), "Q": (-3.5, 0.0, 143.8),
    "E": (-3.5, -1.0, 138.4), "G": (-0.4, 0.0, 60.1), "H": (-3.2, 0.1, 153.2),
    "I": (4.5, 0.0, 166.7), "L": (3.8, 0.0, 166.7), "K": (-3.9, 1.0, 168.6),
    "M": (1.9, 0.0, 162.9), "F": (2.8, 0.0, 189.9), "P": (-1.6, 0.0, 112.7),
    "S": (-0.8, 0.0, 89.0), "T": (-0.7, 0.0, 116.1), "W": (-0.9, 0.0, 227.8),
    "Y": (-1.3, 0.0, 193.6), "V": (4.2, 0.0, 140.0),
}
THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}


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


def aa_value(value: str) -> str:
    value = clean(value)
    return THREE_TO_ONE.get(value, value if len(value) == 1 else "")


def numeric(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def features(row: dict[str, str]) -> list[float]:
    ref = AA.get(aa_value(row.get("reference_aa", "")), (0.0, 0.0, 0.0))
    alt = AA.get(aa_value(row.get("alternate_aa", "")), (0.0, 0.0, 0.0))
    dh = alt[0] - ref[0]
    dc = alt[1] - ref[1]
    dv = alt[2] - ref[2]
    return [numeric(row.get("protein_position", "0")), dh, abs(dh), dc, abs(dc), dv, abs(dv), float(ref == alt)]


def target(row: dict[str, str]) -> str:
    label = clean(row.get("functional_label")).upper()
    if label in {"LOF", "LOSS"}:
        return "LOF"
    if label in {"GOF", "GAIN"}:
        return "GOF"
    return ""


def metrics(y_true: list[str], y_pred: list[str], probabilities: list[float]) -> dict[str, object]:
    labels = ["GOF", "LOF"]
    result: dict[str, object] = {
        "rows": len(y_true),
        "class_counts": dict(Counter(y_true)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, labels=labels, average="macro", zero_division=0)),
        "macro_precision": float(precision_score(y_true, y_pred, labels=labels, average="macro", zero_division=0)),
        "macro_recall": float(recall_score(y_true, y_pred, labels=labels, average="macro", zero_division=0)),
        "confusion_matrix_labels": labels,
        "confusion_matrix": confusion_matrix(y_true, y_pred, labels=labels).tolist(),
    }
    if len(set(y_true)) == 2:
        result["roc_auc"] = float(roc_auc_score([1 if value == "GOF" else 0 for value in y_true], probabilities))
    else:
        result["roc_auc"] = None
    return result


def main() -> None:
    rows = read_csv(PROCESSED / "deepgene_v2_functional_evidence.csv")
    split_rows = read_csv(ML_RESULTS / "v2_grouped_split_ids.csv")
    split_by_id = {
        clean(row.get("row_id")): clean(row.get("split"))
        for row in split_rows
        if clean(row.get("product")) == "functional_evidence"
    }
    usable = []
    for row in rows:
        row = dict(row)
        row["target"] = target(row)
        row["split"] = split_by_id.get(clean(row.get("source_row_id")), "")
        if row["target"] and row["split"]:
            usable.append(row)
    feature_names = ["protein_position", "delta_hydrophobicity", "absolute_delta_hydrophobicity", "delta_charge", "absolute_delta_charge", "delta_volume", "absolute_delta_volume", "same_amino_acid"]
    x = [features(row) for row in usable]
    y = [row["target"] for row in usable]
    train_idx = [i for i, row in enumerate(usable) if row["split"] == "train"]
    val_idx = [i for i, row in enumerate(usable) if row["split"] == "validation"]
    test_idx = [i for i, row in enumerate(usable) if row["split"] == "test"]

    models = {
        "majority_baseline": DummyClassifier(strategy="most_frequent"),
        "logistic_regression": LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42),
        "random_forest": RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=42, n_jobs=1),
        "gradient_boosting": GradientBoostingClassifier(random_state=42),
    }
    report: dict[str, object] = {
        "created": date.today().isoformat(),
        "status": "RESEARCH_ONLY_BASELINE",
        "target_definition": "LOF and Loss mapped to LOF; GOF and Gain mapped to GOF. Other labels are excluded.",
        "feature_names": feature_names,
        "source_rows": len(rows),
        "usable_rows": len(usable),
        "split_rows": {"train": len(train_idx), "validation": len(val_idx), "test": len(test_idx)},
        "models": {},
        "limitations": [
            "Functional assay labels are not clinical pathogenicity labels.",
            "The dataset is source-dependent and class-imbalanced.",
            "The grouped split is reproducible but is not an independent clinical benchmark.",
            "No final model is selected for deployment from this baseline comparison.",
        ],
    }
    prediction_rows: list[dict[str, object]] = []
    fitted_models = {}
    for name, model in models.items():
        model.fit([x[i] for i in train_idx], [y[i] for i in train_idx])
        fitted_models[name] = model
        model_report = {}
        for split_name, indexes in (("train", train_idx), ("validation", val_idx), ("test", test_idx)):
            true = [y[i] for i in indexes]
            predicted = list(model.predict([x[i] for i in indexes]))
            gof_index = list(model.classes_).index("GOF")
            probabilities = [float(value[gof_index]) for value in model.predict_proba([x[i] for i in indexes])]
            model_report[split_name] = metrics(true, predicted, probabilities)
            for index, prediction, probability in zip(indexes, predicted, probabilities):
                prediction_rows.append({
                    "model": name,
                    "source_row_id": usable[index]["source_row_id"],
                    "variant_key": usable[index]["variant_key"],
                    "split": usable[index]["split"],
                    "observed_target": y[index],
                    "predicted_target": prediction,
                    "predicted_gof_probability": f"{probability:.8f}",
                })
        report["models"][name] = model_report

    write_csv(ML_RESULTS / "v2_functional_predictions.csv", prediction_rows, ["model", "source_row_id", "variant_key", "split", "observed_target", "predicted_target", "predicted_gof_probability"])
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    with (MODEL_DIR / "deepgene_v2_functional_baseline.pkl").open("wb") as handle:
        pickle.dump({"models": fitted_models, "feature_names": feature_names, "target_mapping": {"LOF": ["LOF", "Loss"], "GOF": ["GOF", "Gain"]}}, handle)
    (ML_RESULTS / "v2_baseline_model_comparison.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    (ANALYSIS / "v2_baseline_model_report.json").write_text(json.dumps({
        "status": report["status"],
        "best_test_balanced_accuracy": max((details["test"]["balanced_accuracy"], name) for name, details in report["models"].items()),
        "models": report["models"],
        "interpretation": "Metrics are exploratory functional-effect results and must not be interpreted as clinical prediction or pathogenicity validation.",
    }, indent=2), encoding="utf-8")
    print(json.dumps({
        "usable_rows": len(usable),
        "split_rows": {"train": len(train_idx), "validation": len(val_idx), "test": len(test_idx)},
        "models": list(models),
        "outputs": [
            "data/ml_results/v2_baseline_model_comparison.json",
            "data/ml_results/v2_functional_predictions.csv",
            "data/model/deepgene_v2_functional_baseline.pkl",
            "data/analysis_results/v2_baseline_model_report.json",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
