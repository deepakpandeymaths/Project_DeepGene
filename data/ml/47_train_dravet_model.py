"""Train the Dravet model only after reviewed clinical labels exist."""

from pathlib import Path
import pickle
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "processed" / "deepgene_dravet_ml_ready_v1.csv"
OUT = ROOT / "data" / "model" / "deepgene_dravet_research_candidate.pkl"

def main():
    if not INPUT.exists():
        print("STATUS: BLOCKED")
        print("Reason: no prepared clinically reviewed Dravet dataset exists.")
        print("Action: populate data/raw/dravet_clinical_labels.csv, then run 46_prepare_dravet_dataset.py.")
        return
    df = pd.read_csv(INPUT)
    metadata = {"variant_key", "dravet_label", "cohort_id", "label_source"}
    features = [c for c in df.columns if c not in metadata]
    if len(df) < 30 or df.dravet_label.nunique() != 2:
        raise ValueError("At least 30 labelled variants and both Dravet classes are required.")
    cv = StratifiedKFold(n_splits=min(5, int(df.dravet_label.value_counts().min())), shuffle=True, random_state=42)
    candidate = Pipeline([("scale", StandardScaler()), ("model", LogisticRegression(max_iter=5000, class_weight="balanced"))])
    scores = cross_validate(candidate, df[features], df.dravet_label, cv=cv, scoring=["balanced_accuracy", "roc_auc"])
    candidate.fit(df[features], df.dravet_label)
    bundle = {"model": candidate, "features": features, "training_rows": len(df), "cv_balanced_accuracy_mean": float(scores["test_balanced_accuracy"].mean()), "cv_roc_auc_mean": float(scores["test_roc_auc"].mean()), "status": "research_only_not_clinically_validated"}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("wb") as handle: pickle.dump(bundle, handle)
    print(f"Saved: {OUT}")
    print(f"CV balanced accuracy: {bundle['cv_balanced_accuracy_mean']:.4f}")
    print(f"CV ROC-AUC: {bundle['cv_roc_auc_mean']:.4f}")

if __name__ == "__main__":
    main()
