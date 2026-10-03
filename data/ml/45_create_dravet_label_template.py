"""Create the required clinically curated Dravet-label input template."""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "raw" / "dravet_clinical_labels.csv"
INSTRUCTIONS = ROOT / "data" / "raw" / "dravet_clinical_labels_README.txt"

def main():
    columns = ["variant_key", "dravet_label", "label_source", "source_accession", "cohort_id", "reviewer_status"]
    if not OUT.exists():
        pd.DataFrame(columns=columns).to_csv(OUT, index=False)
    INSTRUCTIONS.write_text("""DeepGene Dravet clinical-label schema

Required columns:
  variant_key       SCN1A:p.Arg1648His or SCN1A:p.R1648H
  dravet_label      1 = clinically confirmed Dravet syndrome; 0 = explicitly non-Dravet SCN1A-related phenotype
  label_source      publication, registry, or curated clinical cohort
  source_accession  PMID, database accession, or de-identified cohort identifier
  cohort_id         patient/cohort grouping used to prevent related samples crossing splits
  reviewer_status   independently reviewed / unresolved

Do not infer dravet_label from ClinVar pathogenicity, LOF/GOF, phenotype text alone,
or the prediction output of another model. Unknown and unresolved cases must be
excluded rather than assigned 0.
""", encoding="utf-8")
    print(f"Template: {OUT}")
    print("No clinical labels were invented. Add reviewed labels before training.")

if __name__ == "__main__":
    main()
