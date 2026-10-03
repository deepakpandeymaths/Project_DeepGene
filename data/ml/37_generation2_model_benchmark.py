"""Benchmark Generation-2 LOF/GOF models with gene-held-out validation."""

from pathlib import Path
import json
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "processed" / "deepgene_gen2_ml_ready_v1.csv"
REPORT = ROOT / "data" / "analysis_results" / "44_generation2_model_benchmark.json"


def main():
    df = pd.read_csv(INPUT)
    target = "functional_label"
    metadata = {"variant_key", "gene", "aa1", "aa2", target}
    features = [c for c in df.columns if c not in metadata]
    X, y, groups = df[features], df[target], df["gene"]
    models = {
        "logistic_regression": Pipeline([("scale", StandardScaler()), ("model", LogisticRegression(max_iter=5000, class_weight="balanced"))]),
        "svm_rbf": Pipeline([("scale", StandardScaler()), ("model", SVC(kernel="rbf", class_weight="balanced"))]),
        "gradient_boosting": GradientBoostingClassifier(random_state=42),
        "random_forest": RandomForestClassifier(n_estimators=300, class_weight="balanced", min_samples_leaf=2, random_state=42, n_jobs=1),
    }
    rows = []
    logo = LeaveOneGroupOut()
    for fold, (train_i, test_i) in enumerate(logo.split(X, y, groups), 1):
        for name, model in models.items():
            model.fit(X.iloc[train_i], y.iloc[train_i])
            pred = model.predict(X.iloc[test_i])
            rows.append({"fold": fold, "held_out_gene": groups.iloc[test_i].iloc[0], "model": name, "n_test": len(test_i), "balanced_accuracy": balanced_accuracy_score(y.iloc[test_i], pred), "macro_f1": f1_score(y.iloc[test_i], pred, labels=["LOF", "GOF"], average="macro", zero_division=0)})
    result = pd.DataFrame(rows)
    summary = result.groupby("model")[["balanced_accuracy", "macro_f1"]].agg(["mean", "std"]).round(4)
    payload = {"dataset": str(INPUT), "rows": len(df), "genes": sorted(groups.unique()), "features": features, "validation": "Leave-one-gene-out", "fold_results": rows, "summary": json.loads(summary.to_json()), "interpretation": "Generation-2 cross-channel benchmark; not a clinical or fully independent external validation.", "status": "PASS"}
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(summary.to_string())
    print(f"Saved: {REPORT}")


if __name__ == "__main__":
    main()
