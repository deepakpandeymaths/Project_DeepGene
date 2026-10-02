import pandas as pd
from pathlib import Path


# ============================================================
# 15. FEATURE DISTRIBUTION ANALYSIS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_evidence_v1.csv"
)


print("=" * 60)
print("DEEPGENE FEATURE DISTRIBUTION ANALYSIS")
print("=" * 60)


# ------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------

df = pd.read_csv(INPUT_PATH)

print(f"\nVariants: {len(df)}")
print(f"Features: {len(df.columns)}")


# ------------------------------------------------------------
# 2. CLINICAL EVIDENCE DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("CLINICAL EVIDENCE DISTRIBUTION")
print("=" * 60)

clinical_columns = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
]

for column in clinical_columns:

    print(f"\n{column}")

    print(
        f"  Variants with value > 0: "
        f"{(df[column] > 0).sum()}"
    )

    print(
        f"  Maximum: "
        f"{df[column].max()}"
    )


# ------------------------------------------------------------
# 3. RELIABILITY DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("RELIABILITY DISTRIBUTION")
print("=" * 60)

reliability_columns = [
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
]

for column in reliability_columns:

    print(f"\n{column}")

    print(
        f"  Present: "
        f"{(df[column] == 1).sum()}"
    )

    print(
        f"  Absent: "
        f"{(df[column] == 0).sum()}"
    )


# ------------------------------------------------------------
# 4. PHENOTYPE / QC
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("PHENOTYPE / QC")
print("=" * 60)

qc_columns = [
    "phenotype_available",
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",
]

for column in qc_columns:

    print(f"\n{column}")

    print(
        f"  Present: "
        f"{(df[column] == 1).sum()}"
    )

    print(
        f"  Absent: "
        f"{(df[column] == 0).sum()}"
    )


# ------------------------------------------------------------
# 5. SUBMITTER DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("SUBMITTER DISTRIBUTION")
print("=" * 60)

print(
    df["submitter_count"].describe()
)


# ------------------------------------------------------------
# 6. CLINVAR RECORD DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("EVIDENCE DATASET SUMMARY")
print("=" * 60)

print(
    "Each row represents one unique SCN1A variant."
)

print(
    "ClinVar evidence was deduplicated using "
    "variant ID + clinical significance + review status "
    "+ submitter count + phenotype."
)

print(
    "\nThe previous total_clinvar_records feature was "
    "removed because it became constant (= 1) after "
    "evidence deduplication."
)


print("\nAnalysis complete.")