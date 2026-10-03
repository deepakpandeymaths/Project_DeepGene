# ============================================================
# 33. DEEPGENE PREDICTION PROBLEM DEFINITION
# ============================================================
#
# Purpose:
#   Define the prediction problem before entering the ML phase.
#
# IMPORTANT:
#   This script does NOT train a machine-learning model.
#   It does NOT create labels.
#   It does NOT modify the V1 dataset.
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

REPORT_FILE = (
    ANALYSIS_DIR /
    "33_prediction_problem_definition.txt"
)


# ============================================================
# 2. LOAD DATA
# ============================================================

print("=" * 70)
print("DEEPGENE — PREDICTION PROBLEM DEFINITION")
print("=" * 70)

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Canonical V1 dataset not found: {INPUT_FILE}"
    )

df = pd.read_csv(INPUT_FILE, low_memory=False)

print("\nCanonical V1 dataset loaded.")
print(f"Rows: {len(df)}")
print(f"Columns: {len(df.columns)}")


# ============================================================
# 3. CURRENT DATA SOURCE
# ============================================================

print("\n" + "=" * 70)
print("1. CURRENT DATA SOURCE")
print("=" * 70)

print("Primary source: ClinVar")
print("Gene: SCN1A")
print("Unique variants:", df["variant_id"].nunique())


# ============================================================
# 4. CANDIDATE PREDICTION PROBLEMS
# ============================================================

print("\n" + "=" * 70)
print("2. CANDIDATE PREDICTION PROBLEMS")
print("=" * 70)

candidate_problems = [
    (
        "Pathogenicity prediction",
        "Predict whether a variant is pathogenic or non-pathogenic.",
        "NOT CURRENTLY READY",
        "Requires an independent target label."
    ),

    (
        "Variant prioritization",
        "Rank variants according to their relevance for further investigation.",
        "POSSIBLE AS AN EVIDENCE SYSTEM",
        "Can be transparent and evidence-based, but should not be presented as an ML pathogenicity prediction without an independent target."
    ),

    (
        "Disease relevance prediction",
        "Predict whether a variant is relevant to an SCN1A-related disease.",
        "NOT CURRENTLY READY",
        "Requires an independently defined disease-association target."
    ),

    (
        "Evidence-based ranking",
        "Organize variants according to available evidence and review context.",
        "CURRENTLY POSSIBLE",
        "Can be implemented transparently without supervised ML."
    ),
]

for name, description, status, reason in candidate_problems:

    print(f"\n{name}")
    print("-" * len(name))
    print(f"Description: {description}")
    print(f"Status: {status}")
    print(f"Reason: {reason}")


# ============================================================
# 5. CURRENT CLINVAR LABELS
# ============================================================

print("\n" + "=" * 70)
print("3. CURRENT CLINVAR CLINICAL SIGNIFICANCE")
print("=" * 70)

clinical_counts = (
    df["clinical_significance"]
    .value_counts(dropna=False)
)

for value, count in clinical_counts.items():
    print(f"{value}: {count}")


# ============================================================
# 6. WHY CLINVAR SIGNIFICANCE CANNOT DIRECTLY BE THE TARGET
# ============================================================

print("\n" + "=" * 70)
print("4. TARGET LEAKAGE CHECK")
print("=" * 70)

print(
    "\nThe canonical V1 dataset already contains ClinVar-derived "
    "clinical significance and multiple ClinVar-derived evidence "
    "features."
)

print(
    "\nTherefore, directly using clinical_significance as the supervised "
    "ML target while simultaneously using the same ClinVar evidence "
    "as model inputs would create target leakage."
)

print(
    "\nThe model would effectively be asked to reproduce information "
    "already contained in its input evidence."
)


# ============================================================
# 7. INDEPENDENT LABEL STATUS
# ============================================================

print("\n" + "=" * 70)
print("5. INDEPENDENT TARGET STATUS")
print("=" * 70)

independent_target_available = False

print(
    "Independent pathogenicity benchmark:",
    independent_target_available
)

print(
    "\nNo independent ML target is currently included in the V1 dataset."
)

print(
    "\nArtificial labels will NOT be created from the existing ClinVar "
    "clinical significance field."
)


# ============================================================
# 8. WHAT A VALID TARGET WOULD REQUIRE
# ============================================================

print("\n" + "=" * 70)
print("6. REQUIREMENTS FOR A VALID ML TARGET")
print("=" * 70)

requirements = [
    "The target must represent the biological question DeepGene is trying to answer.",
    "The target must be independently defined from the model input evidence.",
    "The target generation process must be reproducible.",
    "The target must have sufficient examples for evaluation.",
    "The target must allow train/validation/test separation without leakage.",
    "The target definition must be documented before model training.",
]

for item in requirements:
    print(f"- {item}")


# ============================================================
# 9. FEATURES THAT MUST BE TREATED CAREFULLY
# ============================================================

print("\n" + "=" * 70)
print("7. FEATURES REQUIRING SPECIAL HANDLING")
print("=" * 70)

careful_features = [
    "clinical_significance",
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
    "review_status",
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
]

for feature in careful_features:

    if feature in df.columns:
        print(f"- {feature}")


print(
    "\nThese features are not automatically forbidden."
)

print(
    "Their eligibility depends on the final prediction target and "
    "whether they would leak information about that target."
)


# ============================================================
# 10. FEATURES THAT ARE NOT BIOLOGICAL PREDICTORS
# ============================================================

print("\n" + "=" * 70)
print("8. IDENTIFIER / TRACEABILITY FIELDS")
print("=" * 70)

identifier_features = [
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
    "hgvs_name",
]

for feature in identifier_features:

    if feature in df.columns:
        print(f"- {feature}")


print(
    "\nThese fields are retained for traceability but should not "
    "automatically be used as biological predictors."
)


# ============================================================
# 11. CURRENT DECISION
# ============================================================

print("\n" + "=" * 70)
print("9. CURRENT PROJECT DECISION")
print("=" * 70)

print(
    "\nDeepGene will NOT begin supervised ML training yet."
)

print(
    "\nThe next task is to define and obtain an independent target "
    "that matches the intended biological prediction problem."
)

print(
    "\nUntil that target is established, the existing ClinVar dataset "
    "will remain the validated V1 evidence foundation."
)


# ============================================================
# 12. ML BOUNDARY
# ============================================================

print("\n" + "=" * 70)
print("10. ML BOUNDARY")
print("=" * 70)

print(
    "\n🚨 ML HAS NOT STARTED."
)

print(
    "\nML will officially begin only after:"
)

ml_boundary = [
    "Prediction target is defined.",
    "Independent labels are obtained.",
    "Feature leakage is audited.",
    "Dataset splitting strategy is defined.",
    "Evaluation metrics are defined.",
    "A reproducible baseline is specified.",
]

for item in ml_boundary:
    print(f"- {item}")


# ============================================================
# 13. FINAL STATUS
# ============================================================

print("\n" + "=" * 70)
print("FINAL STATUS")
print("=" * 70)

print(
    "\nPHASE 1 — DATA ANALYSIS: COMPLETE"
)

print(
    "PHASE 2 — PREDICTION PROBLEM DEFINITION: IN PROGRESS"
)

print(
    "PHASE 3 — MACHINE LEARNING: NOT STARTED"
)

print(
    "\nNo dataset was modified."
)

print(
    "\nNo labels were created."
)

print(
    "\nNo machine-learning model was trained."
)

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)