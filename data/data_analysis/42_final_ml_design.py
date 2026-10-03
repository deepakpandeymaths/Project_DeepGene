# ============================================================
# 42. FINAL ML TARGET + DESIGN DEFINITION
# ============================================================
#
# Purpose:
#   Finalize the prediction target, feature policy,
#   leakage policy, split strategy and evaluation metrics
#   BEFORE any machine-learning model is trained.
#
# IMPORTANT:
#   This script does NOT train a model.
#
# ============================================================

from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATASET_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_independent_ml_dataset_v1.csv"
)

REPORT_FILE = (
    PROJECT_ROOT
    / "data"
    / "analysis_results"
    / "42_final_ml_design.txt"
)


# ============================================================
# LOAD DATASET
# ============================================================

print("=" * 70)
print("42. FINAL ML TARGET + DESIGN DEFINITION")
print("=" * 70)

if not DATASET_FILE.exists():
    raise FileNotFoundError(
        f"Dataset not found:\n{DATASET_FILE}"
    )

df = pd.read_csv(DATASET_FILE)

print()
print("Dataset:")
print("Rows:", len(df))
print("Columns:", len(df.columns))


# ============================================================
# BASIC INTEGRITY
# ============================================================

assert len(df) > 0

assert df["variant_id"].nunique() == len(df)

print()
print("Unique variants:", df["variant_id"].nunique())

print()
print("Target distribution:")

target_distribution = (
    df["functional_target"]
    .value_counts()
)

print(target_distribution)


# ============================================================
# TARGET DEFINITION
# ============================================================

print()
print("=" * 70)
print("TARGET DEFINITION")
print("=" * 70)

raw_target = df["functional_target"].astype(str).str.strip()

print()
print("Raw experimental functional categories:")

for category, count in raw_target.value_counts().items():

    print(
        f" - {category}: {count}"
    )


# ============================================================
# TARGET ELIGIBILITY
# ============================================================
#
# We do NOT silently merge biological categories.
#
# Categories with extremely small sample counts cannot support
# reliable standalone supervised classes in a 56-variant dataset.
#
# Therefore:
#
#   Loss 35
#   Mixed 12
#   Gain 4
#
# are the only categories with enough observations to be
# considered for a multiclass exploratory model.
#
# "Insufficient data", "No effect", "Unclear effect" and
# "Loss/Gain" are retained in the source dataset but excluded
# from the primary supervised target because their sample sizes
# are too small and/or their meaning is not equivalent to the
# three main functional-effect classes.
#
# ============================================================

PRIMARY_CLASSES = [
    "Loss",
    "Mixed",
    "Gain",
]

EXCLUDED_CLASSES = [
    "Insufficient data",
    "No effect",
    "Unclear effect",
    "Loss/Gain",
]


df["ml_target"] = raw_target.where(
    raw_target.isin(PRIMARY_CLASSES),
    np.nan
)


print()
print("Primary ML classes:")

for category in PRIMARY_CLASSES:

    count = (
        df["ml_target"]
        .eq(category)
        .sum()
    )

    print(
        f" - {category}: {count}"
    )


print()
print("Excluded from primary supervised target:")

for category in EXCLUDED_CLASSES:

    count = (
        raw_target
        .eq(category)
        .sum()
    )

    print(
        f" - {category}: {count}"
    )


# ============================================================
# PRIMARY ML DATASET
# ============================================================

ml_df = df[
    df["ml_target"].notna()
].copy()

ml_df = ml_df.reset_index(drop=True)

print()
print(
    "Primary supervised dataset:",
    len(ml_df),
    "variants"
)

print()
print("Primary target distribution:")

print(
    ml_df["ml_target"]
    .value_counts()
)


# ============================================================
# CLASS BALANCE CHECK
# ============================================================

class_counts = (
    ml_df["ml_target"]
    .value_counts()
)

minimum_class_count = class_counts.min()

print()
print(
    "Smallest primary class:",
    minimum_class_count
)

if minimum_class_count < 3:

    raise ValueError(
        "A primary class has fewer than 3 variants."
    )


# ============================================================
# FEATURE POLICY
# ============================================================
#
# Only variant-derived biological features are allowed.
#
# Metadata:
#   variant_id
#   functional_protein
#   protein_variant
#   reference_aa
#   alternate_aa
#   amino_acid_change_type
#   functional_effect
#   functional_target
#   phenotype
#   inheritance
#
# are NOT model predictors.
#
# ============================================================

FEATURE_COLUMNS = [

    "protein_position",

    "same_amino_acid",

    "reference_hydrophobicity",
    "alternate_hydrophobicity",
    "delta_hydrophobicity",
    "absolute_delta_hydrophobicity",

    "reference_charge",
    "alternate_charge",
    "delta_charge",
    "absolute_delta_charge",

    "reference_volume",
    "alternate_volume",
    "delta_volume",
    "absolute_delta_volume",

]


print()
print("=" * 70)
print("FEATURE POLICY")
print("=" * 70)

print()
print("Allowed biological features:")

for feature in FEATURE_COLUMNS:

    print(
        f" - {feature}"
    )


# ============================================================
# FEATURE VALIDATION
# ============================================================

missing_features = [
    feature
    for feature in FEATURE_COLUMNS
    if feature not in df.columns
]

if missing_features:

    raise ValueError(
        "Missing required features:\n"
        + "\n".join(missing_features)
    )


feature_missingness = (
    ml_df[FEATURE_COLUMNS]
    .isna()
    .sum()
)

print()
print("Feature missingness:")

print(feature_missingness)


if feature_missingness.sum() != 0:

    raise ValueError(
        "Missing feature values detected."
    )


# ============================================================
# EXCLUDED INFORMATION
# ============================================================

EXCLUDED_FROM_MODEL = [

    "variant_id",
    "functional_protein",
    "protein_variant",
    "reference_aa",
    "alternate_aa",
    "amino_acid_change_type",

    "functional_effect",
    "functional_target",

    "phenotype",
    "inheritance",

    # ClinVar-derived information
    "clinical_significance",
    "review_status",
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",

]


print()
print("Excluded from model:")

for feature in EXCLUDED_FROM_MODEL:

    print(
        f" - {feature}"
    )


# ============================================================
# LEAKAGE CHECK
# ============================================================

print()
print("=" * 70)
print("LEAKAGE CHECK")
print("=" * 70)

leakage_columns = [

    "functional_effect",
    "functional_target",

]

leakage_present = [
    column
    for column in leakage_columns
    if column in FEATURE_COLUMNS
]

if leakage_present:

    raise ValueError(
        "Target leakage detected:\n"
        + "\n".join(leakage_present)
    )

print()
print("Functional target used as predictor: NO")

print(
    "ClinVar pathogenicity used as predictor: NO"
)

print(
    "Phenotype used as predictor: NO"
)

print(
    "Identifier used as predictor: NO"
)

print()
print("LEAKAGE CHECK: PASS")


# ============================================================
# SPLIT STRATEGY
# ============================================================
#
# Only 51 usable variants exist after removing sparse classes.
#
# Because the dataset is very small, a single arbitrary
# train/test split would produce unstable estimates.
#
# Therefore the primary evaluation strategy will be:
#
#   Stratified repeated cross-validation
#
# with an additional held-out test set only if the class counts
# permit it without producing unusably small classes.
#
# The exact random seed and folds will be locked in the
# dataset-preparation stage.
#
# ============================================================

N_VARIANTS = len(ml_df)
N_CLASSES = len(PRIMARY_CLASSES)

print()
print("=" * 70)
print("SPLIT / VALIDATION STRATEGY")
print("=" * 70)

print()
print("Usable variants:", N_VARIANTS)

print("Classes:", N_CLASSES)

print()
print(
    "Primary validation strategy:"
)

print(
    "Stratified cross-validation"
)

print(
    "Reason: very small dataset"
)

print(
    "Class stratification: REQUIRED"
)

print(
    "Random seed: WILL BE FIXED BEFORE TRAINING"
)

print(
    "Test-set strategy: determined during dataset preparation"
)


# ============================================================
# METRICS
# ============================================================

METRICS = [

    "balanced_accuracy",
    "macro_precision",
    "macro_recall",
    "macro_f1",
    "confusion_matrix",

]


print()
print("=" * 70)
print("EVALUATION METRICS")
print("=" * 70)

for metric in METRICS:

    print(
        f" - {metric}"
    )


# ============================================================
# BASELINE REQUIREMENT
# ============================================================

print()
print("=" * 70)
print("BASELINE")
print("=" * 70)

print()
print(
    "A simple baseline classifier is REQUIRED before"
)

print(
    "more complex models are evaluated."
)

print()
print(
    "Primary baseline:"
)

print(
    "Majority-class classifier"
)


# ============================================================
# MODELING SCOPE
# ============================================================

print()
print("=" * 70)
print("MODELING SCOPE")
print("=" * 70)

print()
print(
    "Prediction task:"
)

print(
    "Experimental SCN1A functional-effect classification"
)

print()
print(
    "Primary classes:"
)

for category in PRIMARY_CLASSES:

    print(
        f" - {category}"
    )

print()
print(
    "This is NOT:"
)

print(
    " - a clinical diagnostic model"
)

print(
    " - a pathogenic/benign ClinVar classifier"
)

print(
    " - a disease-risk predictor"
)

print(
    " - a clinical decision system"
)


# ============================================================
# FINAL DESIGN STATUS
# ============================================================

print()
print("=" * 70)
print("FINAL DESIGN STATUS")
print("=" * 70)

print()
print(
    "Target defined: YES"
)

print(
    "Independent experimental target: YES"
)

print(
    "Biological features defined: YES"
)

print(
    "ClinVar leakage excluded: YES"
)

print(
    "Identifier leakage excluded: YES"
)

print(
    "Phenotype leakage excluded: YES"
)

print(
    "Evaluation metrics defined: YES"
)

print(
    "Baseline defined: YES"
)

print()
print(
    "ML training started: NO"
)


# ============================================================
# WRITE REPORT
# ============================================================

report = []

report.append(
    "DEEPGENE STEP 42 - FINAL ML DESIGN"
)

report.append("=" * 70)

report.append("")

report.append(
    "DATASET"
)

report.append(
    f"Total Step 41 variants: {len(df)}"
)

report.append(
    f"Primary ML variants: {len(ml_df)}"
)

report.append(
    f"Primary classes: {len(PRIMARY_CLASSES)}"
)

report.append("")

report.append(
    "TARGET"
)

report.append(
    "Independent experimental functional-effect classification"
)

report.append(
    "NOT pathogenic/benign"
)

report.append("")

report.append(
    "PRIMARY CLASSES"
)

for category in PRIMARY_CLASSES:

    report.append(
        f" - {category}: "
        f"{class_counts[category]}"
    )

report.append("")

report.append(
    "EXCLUDED CLASSES"
)

for category in EXCLUDED_CLASSES:

    report.append(
        f" - {category}: "
        f"{raw_target.eq(category).sum()}"
    )

report.append("")

report.append(
    "FEATURES"
)

for feature in FEATURE_COLUMNS:

    report.append(
        f" - {feature}"
    )

report.append("")

report.append(
    "LEAKAGE POLICY"
)

report.append(
    "ClinVar clinical significance: EXCLUDED"
)

report.append(
    "ClinVar pathogenicity counts: EXCLUDED"
)

report.append(
    "ClinVar review/reliability fields: EXCLUDED"
)

report.append(
    "Phenotype: EXCLUDED"
)

report.append(
    "Inheritance: EXCLUDED"
)

report.append(
    "Functional target: EXCLUDED from predictors"
)

report.append(
    "Identifiers: EXCLUDED from predictors"
)

report.append("")

report.append(
    "VALIDATION"
)

report.append(
    "Stratified cross-validation"
)

report.append(
    "Fixed random seed before training"
)

report.append("")

report.append(
    "METRICS"
)

for metric in METRICS:

    report.append(
        f" - {metric}"
    )

report.append("")

report.append(
    "BASELINE"
)

report.append(
    "Majority-class classifier"
)

report.append("")

report.append(
    "ML TRAINING"
)

report.append(
    "NOT STARTED"
)

report.append("")

report.append(
    "STATUS"
)

report.append(
    "FINAL DESIGN DEFINED"
)


REPORT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

REPORT_FILE.write_text(
    "\n".join(report),
    encoding="utf-8"
)


# ============================================================
# SAVE PRIMARY ML DATASET
# ============================================================

PRIMARY_DATASET_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_ready_v1.csv"
)

ml_df.to_csv(
    PRIMARY_DATASET_FILE,
    index=False
)


print()
print("=" * 70)
print("STEP 42 COMPLETE")
print("=" * 70)

print()
print(
    "Primary ML dataset:"
)

print(
    PRIMARY_DATASET_FILE
)

print()
print(
    "Design report:"
)

print(
    REPORT_FILE
)

print()
print(
    "STATUS: PASS"
)

print()
print(
    "ML TRAINING: NOT STARTED"
)

print()
print(
    "NEXT: FINAL VALIDATION OF DESIGN"
)