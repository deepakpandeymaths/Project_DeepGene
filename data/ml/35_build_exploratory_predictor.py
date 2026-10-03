"""Build the auditable DeepGene Generation-1 exploratory predictor.

This script deliberately does not create a clinical model.  It trains the
three-class functional-effect candidates on the independent Brunklaus-derived
dataset, records repeated-CV performance, and saves a research-only artifact
for transparent inference.
"""

from pathlib import Path
import json
import pickle
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "processed" / "deepgene_ml_ready_v1.csv"
MODEL_DIR = ROOT / "data" / "model"
ARTIFACT = MODEL_DIR / "deepgene_functional_exploratory.pkl"
REPORT = MODEL_DIR / "deepgene_functional_exploratory_report.json"

TARGET = "functional_target"
CLASSES = ["Loss", "Mixed", "Gain"]
FEATURES = [
    "protein_position", "same_amino_acid", "reference_hydrophobicity",
    "alternate_hydrophobicity", "delta_hydrophobicity",
    "absolute_delta_hydrophobicity", "reference_charge", "alternate_charge",
    "delta_charge", "absolute_delta_charge", "reference_volume",
    "alternate_volume", "delta_volume", "absolute_delta_volume",
]
SEEDS = [42, 123, 2024, 7, 99]
N_SPLITS = 4


def candidates():
    return {
        "logistic_regression": Pipeline([
            ("scale", StandardScaler()),
            ("model", LogisticRegression(max_iter=5000, class_weight="balanced")),
        ]),
        "random_forest": RandomForestClassifier(
            n_estimators=200, class_weight="balanced", random_state=42,
            min_samples_leaf=2, n_jobs=1,
        ),
    }


def validate(df):
    missing = [c for c in [TARGET, *FEATURES] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if df[FEATURES].isna().any().any() or df[TARGET].isna().any():
        raise ValueError("Missing values found in predictors or target.")
    if sorted(df[TARGET].unique()) != sorted(CLASSES):
        raise ValueError(f"Expected target classes {CLASSES}, found {sorted(df[TARGET].unique())}")


def evaluate(df):
    X, y = df[FEATURES], df[TARGET]
    rows = []
    majority = y.value_counts().idxmax()
    for seed in SEEDS:
        cv = StratifiedKFold(N_SPLITS, shuffle=True, random_state=seed)
        for fold, (train_i, valid_i) in enumerate(cv.split(X, y), 1):
            baseline = np.repeat(majority, len(valid_i))
            rows.append({"model": "majority_baseline", "seed": seed, "fold": fold,
                         "balanced_accuracy": balanced_accuracy_score(y.iloc[valid_i], baseline),
                         "macro_f1": f1_score(y.iloc[valid_i], baseline, labels=CLASSES,
                                               average="macro", zero_division=0)})
            for name, model in candidates().items():
                model.fit(X.iloc[train_i], y.iloc[train_i])
                pred = model.predict(X.iloc[valid_i])
                rows.append({"model": name, "seed": seed, "fold": fold,
                             "balanced_accuracy": balanced_accuracy_score(y.iloc[valid_i], pred),
                             "macro_f1": f1_score(y.iloc[valid_i], pred, labels=CLASSES,
                                                   average="macro", zero_division=0)})
    scores = pd.DataFrame(rows)
    summary = scores.groupby("model")[["balanced_accuracy", "macro_f1"]].agg(["mean", "std"])
    return scores, summary


def main():
    df = pd.read_csv(INPUT)
    validate(df)
    scores, summary = evaluate(df)
    model_names = ["logistic_regression", "random_forest"]
    winner = max(model_names, key=lambda n: summary.loc[n, ("balanced_accuracy", "mean")])
    baseline_score = float(summary.loc["majority_baseline", ("balanced_accuracy", "mean")])
    winner_score = float(summary.loc[winner, ("balanced_accuracy", "mean")])
    status = "exploratory_only" if winner_score <= baseline_score else "exploratory_candidate"

    final_model = candidates()[winner]
    final_model.fit(df[FEATURES], df[TARGET])
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    bundle = {
        "model": final_model, "features": FEATURES, "classes": CLASSES,
        "selected_model": winner, "status": status,
        "training_rows": int(len(df)), "class_counts": df[TARGET].value_counts().to_dict(),
        "validation_design": {"folds": N_SPLITS, "repeats": len(SEEDS), "seeds": SEEDS},
        "summary": summary.to_dict(),
    }
    with ARTIFACT.open("wb") as handle:
        pickle.dump(bundle, handle)
    report = {
        "artifact": str(ARTIFACT), "status": status, "selected_model": winner,
        "classes": CLASSES, "features": FEATURES, "training_rows": len(df),
        "class_counts": df[TARGET].value_counts().to_dict(),
        "validation": scores.to_dict(orient="records"),
        "summary": json.loads(summary.to_json()),
        "scientific_boundary": "Functional-effect research output only; not pathogenicity, disease-risk, diagnosis, or treatment advice.",
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Saved: {ARTIFACT}")
    print(f"Selected candidate: {winner}")
    print(f"Status: {status}")
    print(summary.to_string())


if __name__ == "__main__":
    main()
