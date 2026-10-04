"""Prepare and train a research-only SCN1A functional-effect candidate model.

Source: Brunklaus et al. (2020), UCL supplementary spreadsheet.
Only rows explicitly labelled Loss or Gain are used. Mixed, unclear,
insufficient-data, no-effect, and missing labels are excluded from this
binary LOF/GOF experiment rather than being forced into a class.
"""

from __future__ import annotations

import json
import pickle
import re
import sys
from datetime import date
from pathlib import Path

import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import balanced_accuracy_score, classification_report, roc_auc_score
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from app.server import feature_frame


RAW = ROOT / "data/raw/SCN1A_Brunklaus_2020_supplementary.xlsx"
NORMALIZED = ROOT / "data/processed/deepgene_scn1a_brunklaus_2020_functional.csv"
MODEL = ROOT / "data/model/deepgene_scn1a_supplement_candidate.pkl"
REPORT = ROOT / "data/analysis_results/49_scn1a_supplement_model.json"
SOURCE_URL = "https://discovery.ucl.ac.uk/id/eprint/10091164/10/Schorge_SCN1A%20variants%20from%20bench%20to%20bedside-improved%20clinical%20prediction%20from%20functional%20characterization_Spreadsheet.xlsx"

AA = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
VARIANT = re.compile(r"^(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2})$")


def prepare() -> pd.DataFrame:
    if NORMALIZED.exists():
        return pd.read_csv(NORMALIZED)
    source = pd.read_excel(RAW, sheet_name="Sheet1")
    keep = source[source["Function"].isin(["Loss", "Gain"])].copy()
    rows = []
    excluded = source.shape[0] - keep.shape[0]
    for _, row in keep.iterrows():
        match = VARIANT.match(str(row["Variant"]).strip())
        if not match or match.group("ref") not in AA or match.group("alt") not in AA:
            excluded += 1
            continue
        rows.append({
            "id": f"brunklaus2020_{len(rows) + 1:03d}",
            "gene": "SCN1A",
            "y": "LOF" if row["Function"] == "Loss" else "GOF",
            "pheno": str(row.get("Phenotype", "")),
            "aa1": AA[match.group("ref")],
            "aa2": AA[match.group("alt")],
            "pos": int(match.group("pos")),
            "variant_key": str(row["Variant"]).strip(),
            "functional_label": str(row["Function"]),
            "label_status": "published_experimental_function",
            "source_name": "Brunklaus et al. 2020 supplementary functional characterization",
            "source_url": SOURCE_URL,
        })
    result = pd.DataFrame(rows)
    result.to_csv(NORMALIZED, index=False)
    return result


def train(supplement: pd.DataFrame) -> dict:
    base = pd.read_csv(ROOT / "data/processed/deepgene_gen2_scion_harmonized_v1.csv")
    combined = pd.concat([base, supplement], ignore_index=True)
    combined = combined.drop_duplicates(["gene", "pos", "aa1", "aa2", "y"], keep="last")
    features_frame = feature_frame(combined)
    numeric = [c for c in features_frame.columns if c not in {"gene", "functional_label"}]
    genes = sorted(features_frame.gene.unique())
    for gene in genes:
        features_frame[f"gene_{gene}"] = (features_frame.gene == gene).astype(int)
    feature_names = numeric + [f"gene_{gene}" for gene in genes]
    X = features_frame[feature_names]
    y = features_frame["functional_label"]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42, stratify=y)
    evaluation_model = GradientBoostingClassifier(random_state=42).fit(X_train, y_train)
    probabilities = evaluation_model.predict_proba(X_test)[:, list(evaluation_model.classes_).index("GOF")]
    predictions = evaluation_model.predict(X_test)
    metrics = {
        "rows_in_source": int(len(supplement)),
        "source_label_counts": {str(k): int(v) for k, v in supplement["y"].value_counts().items()},
        "rows_in_combined_training_table": int(len(features_frame)),
        "combined_label_counts": {str(k): int(v) for k, v in y.value_counts().items()},
        "test_rows": int(len(y_test)),
        "balanced_accuracy": float(balanced_accuracy_score(y_test, predictions)),
        "roc_auc": float(roc_auc_score((y_test == "GOF").astype(int), probabilities)),
        "classification_report": classification_report(y_test, predictions, output_dict=True, zero_division=0),
    }
    final_model = GradientBoostingClassifier(random_state=42).fit(X, y)
    with MODEL.open("wb") as handle:
        pickle.dump({"model": final_model, "features": feature_names, "classes": list(final_model.classes_), "training_rows": len(features_frame), "source": "Brunklaus et al. 2020 supplementary + DeepGene Gen-2 harmonized data", "status": "research_candidate_not_clinical"}, handle)
    report = {"created": str(date.today()), "source_file": str(RAW), "source_url": SOURCE_URL, "excluded_source_rows": int(59 - len(supplement)), "model_file": str(MODEL), "metrics": metrics, "limitations": ["This is a binary LOF/GOF research experiment; mixed and uncertain effects were excluded.", "The random holdout is internal, not an independent external validation set.", "The model is not a pathogenicity classifier, diagnostic tool, or treatment recommendation system."]}
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


if __name__ == "__main__":
    prepared = prepare()
    result = train(prepared)
    print(json.dumps(result, indent=2))
