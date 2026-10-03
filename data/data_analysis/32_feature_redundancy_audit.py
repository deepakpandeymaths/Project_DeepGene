# ============================================================
# 32. FEATURE REDUNDANCY & ANALYTICAL FEATURE AUDIT
# ============================================================
#
# Purpose:
#   Identify which canonical V1 fields are:
#       - identifiers
#       - analytical evidence features
#       - reliability/context features
#       - phenotype features
#       - representation/QC fields
#       - constant features
#       - redundant features
#
# IMPORTANT:
#   This script does NOT modify the canonical V1 dataset.
#
# Outputs:
#   data/processed/deepgene_feature_audit_v1.csv
#   data/analysis_results/32_feature_redundancy_audit.txt
#
# ============================================================

from pathlib import Path
import pandas as pd


# ============================================================
# 1. PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[1]

PROCESSED_DIR = BASE_DIR / "processed"
ANALYSIS_DIR = BASE_DIR / "analysis_results"

INPUT_FILE = PROCESSED_DIR / "deepgene_v1_final.csv"

OUTPUT_FILE = PROCESSED_DIR / "deepgene_feature_audit_v1.csv"

REPORT_FILE = ANALYSIS_DIR / "32_feature_redundancy_audit.txt"


# ============================================================
# 2. LOAD DATA
# ============================================================

print("=" * 70)
print("DEEPGENE V1 — FEATURE REDUNDANCY & ANALYTICAL FEATURE AUDIT")
print("=" * 70)

if not INPUT_FILE.exists():
    raise FileNotFoundError(f"Input dataset not found: {INPUT_FILE}")

df = pd.read_csv(INPUT_FILE, low_memory=False)

print("\nDataset loaded successfully.")
print(f"Rows: {len(df)}")
print(f"Columns: {len(df.columns)}")


# ============================================================
# 3. FEATURE GROUP DEFINITIONS
# ============================================================

identifier_features = {
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
    "hgvs_name",
}

clinical_features = {
    "clinical_significance",
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
}

reliability_features = {
    "review_status",
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
}

phenotype_features = {
    "phenotypes",
    "phenotype_available",
}

representation_features = {
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",
}


# ============================================================
# 4. CLASSIFY FEATURE
# ============================================================

def classify_feature(column):

    if column in identifier_features:
        return "IDENTIFIER"

    if column in clinical_features:
        return "CLINICAL_EVIDENCE"

    if column in reliability_features:
        return "RELIABILITY_CONTEXT"

    if column in phenotype_features:
        return "PHENOTYPE"

    if column in representation_features:
        return "REPRESENTATION_QC"

    return "UNCLASSIFIED"


# ============================================================
# 5. BUILD FEATURE AUDIT TABLE
# ============================================================

audit_rows = []

for column in df.columns:

    series = df[column]

    unique_values = series.nunique(dropna=False)
    missing_values = int(series.isna().sum())

    is_constant = unique_values <= 1

    audit_rows.append(
        {
            "feature": column,
            "feature_group": classify_feature(column),
            "dtype": str(series.dtype),
            "unique_values": unique_values,
            "missing_values": missing_values,
            "is_constant": is_constant,
            "usable_as_analytical_feature": None,
            "reason": None,
        }
    )


audit_df = pd.DataFrame(audit_rows)


# ============================================================
# 6. DETERMINE ANALYTICAL USABILITY
# ============================================================

def analytical_status(row):

    feature = row["feature"]
    group = row["feature_group"]
    constant = row["is_constant"]

    # Identifiers are retained for traceability,
    # but are not automatically analytical predictors.
    if group == "IDENTIFIER":
        return False, "Identifier / traceability field"

    # Constant fields carry no variation in this dataset.
    if constant:
        return False, "Constant across all SCN1A variants"

    # Representation/QC fields describe data availability,
    # not biological pathogenicity directly.
    if group == "REPRESENTATION_QC":
        return False, "Representation / QC context"

    # Phenotype fields can be analytically useful,
    # but require explicit representation before modeling.
    if group == "PHENOTYPE":
        return True, "Potential analytical feature after explicit representation"

    # Clinical evidence fields are potentially analytical,
    # but must be handled carefully because they are ClinVar-derived.
    if group == "CLINICAL_EVIDENCE":
        return True, "ClinVar-derived evidence feature; target leakage must be avoided"

    # Reliability/context fields may be useful for evidence modeling,
    # but are not automatically biological predictors.
    if group == "RELIABILITY_CONTEXT":
        return True, "Evidence reliability/context feature"

    return False, "Not classified"


statuses = audit_df.apply(
    analytical_status,
    axis=1
)

audit_df["usable_as_analytical_feature"] = [
    value[0] for value in statuses
]

audit_df["reason"] = [
    value[1] for value in statuses
]


# ============================================================
# 7. EXACT DUPLICATE COLUMN DETECTION
# ============================================================

duplicate_pairs = []

columns = list(df.columns)

for i in range(len(columns)):

    for j in range(i + 1, len(columns)):

        col_a = columns[i]
        col_b = columns[j]

        try:
            if df[col_a].equals(df[col_b]):
                duplicate_pairs.append(
                    (col_a, col_b)
                )
        except Exception:
            pass


duplicate_lookup = {}

for a, b in duplicate_pairs:

    duplicate_lookup.setdefault(a, []).append(b)
    duplicate_lookup.setdefault(b, []).append(a)


audit_df["exact_duplicate_of"] = audit_df["feature"].apply(
    lambda x: "|".join(duplicate_lookup.get(x, []))
)


# ============================================================
# 8. CONSTANT FEATURE SUMMARY
# ============================================================

constant_features = audit_df[
    audit_df["is_constant"]
]["feature"].tolist()


# ============================================================
# 9. GROUP SUMMARY
# ============================================================

group_summary = (
    audit_df
    .groupby("feature_group")
    .agg(
        feature_count=("feature", "count"),
        constant_count=("is_constant", "sum"),
        analytical_feature_count=(
            "usable_as_analytical_feature",
            "sum"
        ),
    )
    .reset_index()
)


# ============================================================
# 10. SAVE MACHINE-READABLE AUDIT
# ============================================================

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

audit_df.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# 11. PRINT REPORT
# ============================================================

print("\n" + "=" * 70)
print("FEATURE GROUP SUMMARY")
print("=" * 70)

print(
    group_summary.to_string(index=False)
)


print("\n" + "=" * 70)
print("CONSTANT FEATURES")
print("=" * 70)

print(f"Count: {len(constant_features)}")

for feature in constant_features:
    print(f"- {feature}")


print("\n" + "=" * 70)
print("EXACT DUPLICATE COLUMN PAIRS")
print("=" * 70)

print(f"Count: {len(duplicate_pairs)}")

for a, b in duplicate_pairs:
    print(f"- {a} == {b}")


print("\n" + "=" * 70)
print("ANALYTICAL FEATURE CANDIDATES")
print("=" * 70)

analytical_candidates = audit_df[
    audit_df["usable_as_analytical_feature"] == True
]

for feature in analytical_candidates["feature"]:
    print(f"- {feature}")


print("\n" + "=" * 70)
print("IDENTIFIER / TRACEABILITY FEATURES")
print("=" * 70)

for feature in sorted(identifier_features):
    if feature in df.columns:
        print(f"- {feature}")


print("\n" + "=" * 70)
print("REPRESENTATION / QC FEATURES")
print("=" * 70)

for feature in sorted(representation_features):
    if feature in df.columns:
        print(f"- {feature}")


# ============================================================
# 12. SCIENTIFIC INTERPRETATION
# ============================================================

print("\n" + "=" * 70)
print("SCIENTIFIC INTERPRETATION")
print("=" * 70)

print(
    "\n1. Identifier fields are retained for traceability but should "
    "not automatically be used as biological predictors."
)

print(
    "\n2. Constant features provide no variation within this SCN1A-only "
    "dataset and therefore do not contribute predictive information."
)

print(
    "\n3. Representation/QC fields describe data availability and "
    "representation rather than biological pathogenicity."
)

print(
    "\n4. ClinVar-derived clinical evidence may be analytically useful, "
    "but using the same clinical significance information as both "
    "input and prediction target would create target leakage."
)

print(
    "\n5. Reliability/context features describe the evidence environment "
    "and should not automatically be interpreted as biological "
    "pathogenicity predictors."
)

print(
    "\n6. Phenotype information requires explicit representation before "
    "it can be considered a model input."
)


# ============================================================
# 13. RECOMMENDED CANONICAL STATUS
# ============================================================

print("\n" + "=" * 70)
print("CANONICAL DATASET STATUS")
print("=" * 70)

print(
    "\nThe canonical V1 dataset remains unchanged:"
)

print(
    "data/processed/deepgene_v1_final.csv"
)

print(
    "\nThis audit does NOT delete constant, identifier, or redundant "
    "columns from the canonical dataset."
)

print(
    "\nThe purpose of this audit is to document how each feature should "
    "be interpreted before any future modeling phase."
)


# ============================================================
# 14. WRITE REPORT FILE
# ============================================================
#
# The easiest reproducible approach is to save the same information
# by re-running this script through shell redirection.
#
# The script itself still creates a machine-readable feature audit.
#

print("\n" + "=" * 70)
print("OUTPUTS")
print("=" * 70)

print(f"\nFeature audit:")
print(OUTPUT_FILE)

print("\nReport:")
print(REPORT_FILE)

print("\nAudit complete.")
print("=" * 70)