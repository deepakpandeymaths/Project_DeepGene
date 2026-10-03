from pathlib import Path
import pandas as pd


# ============================================================
# Project DeepGene
# Step 02: Feature Selection
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_ready_v1.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
)

RESULTS_DIR = (
    PROJECT_ROOT
    / "data"
    / "ml_results"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

FEATURE_DATASET = (
    OUTPUT_DIR
    / "deepgene_ml_features_v1.csv"
)

FEATURE_MANIFEST = (
    RESULTS_DIR
    / "02_feature_manifest_v1.csv"
)

REPORT_FILE = (
    RESULTS_DIR
    / "02_feature_selection.txt"
)


# ============================================================
# Load dataset
# ============================================================

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Input dataset not found:\n{INPUT_FILE}"
    )

df = pd.read_csv(INPUT_FILE)


# ============================================================
# Expected feature groups
# ============================================================

# ------------------------------------------------------------
# Evidence-count features
# ------------------------------------------------------------

evidence_features = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
]


# ------------------------------------------------------------
# Reliability features
# ------------------------------------------------------------

reliability_features = [
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
]


# ------------------------------------------------------------
# Conflict feature
# ------------------------------------------------------------

conflict_features = [
    "conflict_flag",
]


# ------------------------------------------------------------
# Phenotype feature
# ------------------------------------------------------------

phenotype_features = [
    "phenotype_available",
]


# ------------------------------------------------------------
# Representation / QC features
# ------------------------------------------------------------

representation_features = [
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",
]


# ============================================================
# Combine candidate features
# ============================================================

candidate_features = (
    evidence_features
    + reliability_features
    + conflict_features
    + phenotype_features
    + representation_features
)


# ============================================================
# Validate feature existence
# ============================================================

missing_features = [
    feature
    for feature in candidate_features
    if feature not in df.columns
]

if missing_features:
    raise ValueError(
        "The following expected features are missing:\n"
        + "\n".join(missing_features)
    )


# ============================================================
# Features requiring leakage assessment
# ============================================================

leakage_risk_features = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
]


# ============================================================
# Lower-risk initial features
# ============================================================

lower_risk_features = [
    "submitter_count",
    "phenotype_available",
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",
]


# ============================================================
# Build feature manifest
# ============================================================

manifest_rows = []

for feature in candidate_features:

    if feature in leakage_risk_features:

        category = "LEAKAGE_RISK"
        initial_use = "AUDIT_BEFORE_MODELING"

        reason = (
            "Feature is derived from ClinVar evidence or "
            "review information and may encode information "
            "closely related to the eventual target."
        )

    elif feature in lower_risk_features:

        category = "CANDIDATE_FEATURE"
        initial_use = "ELIGIBLE_AFTER_AUDIT"

        reason = (
            "Feature represents phenotype availability, "
            "variant representation, or basic evidence "
            "metadata rather than a direct clinical label."
        )

    else:

        category = "REVIEW"
        initial_use = "REQUIRES_REVIEW"

        reason = (
            "Feature requires additional assessment before "
            "being used in supervised learning."
        )

    manifest_rows.append(
        {
            "feature": feature,
            "category": category,
            "initial_use": initial_use,
            "reason": reason,
        }
    )


manifest = pd.DataFrame(manifest_rows)


# ============================================================
# Verify no duplicate features
# ============================================================

if manifest["feature"].duplicated().any():
    raise ValueError(
        "Duplicate feature names detected in feature manifest."
    )


# ============================================================
# Create ML feature dataset
# ============================================================

feature_df = df[candidate_features].copy()


# ============================================================
# Validate missing values
# ============================================================

missing_total = int(
    feature_df.isnull().sum().sum()
)

if missing_total != 0:
    raise ValueError(
        f"Missing values detected in selected features: "
        f"{missing_total}"
    )


# ============================================================
# Validate binary features
# ============================================================

binary_features = [
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
    "phenotype_available",
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",
]

for feature in binary_features:

    if feature not in feature_df.columns:
        continue

    values = set(
        feature_df[feature]
        .dropna()
        .unique()
    )

    if not values.issubset({0, 1}):
        raise ValueError(
            f"Binary feature '{feature}' contains "
            f"unexpected values: {values}"
        )


# ============================================================
# Validate numeric features
# ============================================================

numeric_features = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
    "submitter_count",
]

for feature in numeric_features:

    if not pd.api.types.is_numeric_dtype(
        feature_df[feature]
    ):
        raise TypeError(
            f"Feature '{feature}' is not numeric."
        )


# ============================================================
# Save outputs
# ============================================================

feature_df.to_csv(
    FEATURE_DATASET,
    index=False
)

manifest.to_csv(
    FEATURE_MANIFEST,
    index=False
)


# ============================================================
# Create report
# ============================================================

report = []

report.append("=" * 70)
report.append("PROJECT DEEPGENE")
report.append("STEP 02 — FEATURE SELECTION")
report.append("=" * 70)

report.append("")
report.append("INPUT DATASET")
report.append("-" * 70)
report.append(f"Input: {INPUT_FILE}")
report.append(f"Variants: {len(df)}")
report.append(f"Original ML columns: {len(df.columns)}")

report.append("")
report.append("SELECTED FEATURE GROUPS")
report.append("-" * 70)

report.append(
    f"Evidence features: {len(evidence_features)}"
)

for feature in evidence_features:
    report.append(f"  - {feature}")

report.append(
    f"\nReliability features: {len(reliability_features)}"
)

for feature in reliability_features:
    report.append(f"  - {feature}")

report.append(
    f"\nConflict features: {len(conflict_features)}"
)

for feature in conflict_features:
    report.append(f"  - {feature}")

report.append(
    f"\nPhenotype features: {len(phenotype_features)}"
)

for feature in phenotype_features:
    report.append(f"  - {feature}")

report.append(
    f"\nRepresentation/QC features: "
    f"{len(representation_features)}"
)

for feature in representation_features:
    report.append(f"  - {feature}")

report.append("")
report.append("FEATURE SUMMARY")
report.append("-" * 70)

report.append(
    f"Total candidate features: "
    f"{len(candidate_features)}"
)

report.append(
    f"Leakage-risk features: "
    f"{len(leakage_risk_features)}"
)

report.append(
    f"Lower-risk candidate features: "
    f"{len(lower_risk_features)}"
)

report.append(
    f"Missing values: {missing_total}"
)

report.append("")
report.append("IMPORTANT ML CONSTRAINT")
report.append("-" * 70)

report.append(
    "No independent target is currently available."
)

report.append(
    "Clinical significance must not be used as the "
    "independent target."
)

report.append(
    "Leakage-risk features must be assessed before "
    "supervised model training."
)

report.append("")
report.append("OUTPUTS")
report.append("-" * 70)

report.append(
    f"Feature dataset: {FEATURE_DATASET}"
)

report.append(
    f"Feature manifest: {FEATURE_MANIFEST}"
)

report.append("")
report.append("STATUS")
report.append("-" * 70)

report.append(
    "PASS — candidate feature dataset created."
)

report.append(
    "PASS — feature types validated."
)

report.append(
    "PASS — binary features validated."
)

report.append(
    "PASS — no missing values detected."
)

report.append(
    "NEXT — perform train/test strategy and "
    "leakage assessment before modeling."
)

report.append("")
report.append("=" * 70)


REPORT_FILE.write_text(
    "\n".join(report),
    encoding="utf-8"
)


# ============================================================
# Console output
# ============================================================

print("=" * 70)
print("PROJECT DEEPGENE")
print("STEP 02 — FEATURE SELECTION")
print("=" * 70)

print(f"Input variants: {len(df)}")
print(f"Candidate features: {len(candidate_features)}")
print(
    f"Leakage-risk features: "
    f"{len(leakage_risk_features)}"
)
print(
    f"Lower-risk candidate features: "
    f"{len(lower_risk_features)}"
)
print(f"Missing values: {missing_total}")

print("")
print("Feature dataset:")
print(FEATURE_DATASET)

print("")
print("Feature manifest:")
print(FEATURE_MANIFEST)

print("")
print("Report:")
print(REPORT_FILE)

print("")
print("STATUS: PASS")
print("=" * 70)