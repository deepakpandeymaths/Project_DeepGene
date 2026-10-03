"""Fit the Generation-2 gradient-boosting research artifact."""

from pathlib import Path
import pickle
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "processed" / "deepgene_gen2_ml_ready_v1.csv"
OUT = ROOT / "data" / "model" / "deepgene_generation2_gbm.pkl"


def main():
    df = pd.read_csv(INPUT)
    metadata = {"variant_key", "gene", "aa1", "aa2", "functional_label"}
    features = [c for c in df.columns if c not in metadata]
    model = GradientBoostingClassifier(random_state=42)
    model.fit(df[features], df["functional_label"])
    bundle = {"model": model, "features": features, "classes": list(model.classes_), "training_rows": len(df), "genes": sorted(df.gene.unique()), "status": "generation2_research_candidate", "source": "SCION harmonized table with exact DeepGene V1 overlap excluded"}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("wb") as handle:
        pickle.dump(bundle, handle)
    print(f"Saved: {OUT}")
    print(f"Rows: {len(df)}; features: {len(features)}; classes: {list(model.classes_)}")


if __name__ == "__main__":
    main()
