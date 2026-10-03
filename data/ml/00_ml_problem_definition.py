from pathlib import Path
import pandas as pd


# ============================================================
# Project DeepGene
# Step 00: ML Problem Definition
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "data" / "processed" / "deepgene_v1_final.csv"
OUTPUT_DIR = PROJECT_ROOT / "data" / "ml_results"
REPORT_FILE = OUTPUT_DIR / "00_ml_problem_definition.txt"


# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# ML problem definition
# ------------------------------------------------------------

ML_PROBLEM = {
    "project": "Project_DeepGene",
    "gene": "SCN1A",

    "primary_ml_goal": (
        "Develop a reproducible machine-learning framework for "
        "SCN1A variant prediction using appropriately defined "
        "independent labels and non-leaking features."
    ),

    "current_dataset": "deepgene_v1_final.csv",

    "current_data_source": "ClinVar-derived variant evidence",

    "current_target_status": (
        "NOT AVAILABLE — an independent ML benchmark label has "
        "not yet been established."
    ),

    "clinical_significance_status": (
        "Clinical significance is ClinVar-derived evidence and "
        "must not be used as an independent target for the current "
        "ML benchmark."
    ),

    "training_status": (
        "MODEL TRAINING MUST NOT START until an independent target "
        "definition is available."
    ),

    "feature_policy": (
        "Candidate features may be prepared and audited, but "
        "features derived directly from the target source must "
        "be assessed for leakage before supervised training."
    ),

    "evaluation_policy": (
        "Final model evaluation must use data that are independent "
        "of model fitting and target construction."
    ),
}


# ------------------------------------------------------------
# Load and validate current dataset
# ------------------------------------------------------------

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Input dataset not found:\n{INPUT_FILE}"
    )

df = pd.read_csv(INPUT_FILE)


# ------------------------------------------------------------
# Basic dataset validation
# ------------------------------------------------------------

required_columns = [
    "variant_id",
    "clinical_significance",
    "review_status",
]

missing_columns = [
    column for column in required_columns
    if column not in df.columns
]

if missing_columns:
    raise ValueError(
        f"Required columns are missing: {missing_columns}"
    )


if df["variant_id"].duplicated().any():
    raise ValueError(
        "variant_id is not unique. "
        "ML dataset preparation requires unique variants."
    )


# ------------------------------------------------------------
# Identify current target-related information
# ------------------------------------------------------------

clinical_significance_values = (
    df["clinical_significance"]
    .fillna("MISSING")
    .value_counts()
    .to_dict()
)

review_status_values = (
    df["review_status"]
    .fillna("MISSING")
    .value_counts()
    .to_dict()
)


# ------------------------------------------------------------
# Feature policy
# ------------------------------------------------------------

identifier_columns = [
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
]

raw_text_columns = [
    "hgvs_name",
    "phenotypes",
    "review_status",
]

target_candidate_columns = [
    "clinical_significance",
]

candidate_feature_columns = [
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
    "phenotype_available",
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",
]


# ------------------------------------------------------------
# Leakage-risk features
# ------------------------------------------------------------

leakage_risk_columns = [
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
    "review_status",
    "clinical_significance",
]


# ------------------------------------------------------------
# Build feature policy table
# ------------------------------------------------------------

feature_policy_rows = []

for column in df.columns:

    if column in candidate_feature_columns:
        role = "CANDIDATE_FEATURE"

        if column in leakage_risk_columns:
            leakage_status = (
                "REQUIRES_LEAKAGE_ASSESSMENT"
            )
        else:
            leakage_status = "LOWER_INITIAL_RISK"

    elif column in identifier_columns:
        role = "IDENTIFIER"
        leakage_status = "EXCLUDE_FROM_MODEL"

    elif column in raw_text_columns:
        role = "RAW_TEXT_OR_REVIEW_FIELD"
        leakage_status = "EXCLUDE_FROM_INITIAL_MODEL"

    elif column in target_candidate_columns:
        role = "TARGET_CANDIDATE"
        leakage_status = (
            "DO_NOT_USE_AS_INDEPENDENT_TARGET"
        )

    else:
        role = "UNCLASSIFIED"
        leakage_status = "REQUIRES_REVIEW"

    feature_policy_rows.append(
        {
            "column": column,
            "role": role,
            "leakage_status": leakage_status,
        }
    )

feature_policy = pd.DataFrame(feature_policy_rows)


# ------------------------------------------------------------
# Write ML problem-definition report
# ------------------------------------------------------------

report_lines = []

report_lines.append("=" * 70)
report_lines.append("PROJECT DEEPGENE")
report_lines.append("STEP 00 — ML PROBLEM DEFINITION")
report_lines.append("=" * 70)

report_lines.append("")
report_lines.append("PROJECT")
report_lines.append("-" * 70)
report_lines.append(f"Project: {ML_PROBLEM['project']}")
report_lines.append(f"Gene: {ML_PROBLEM['gene']}")

report_lines.append("")
report_lines.append("CURRENT DATASET")
report_lines.append("-" * 70)
report_lines.append(f"Input file: {INPUT_FILE}")
report_lines.append(f"Rows: {len(df)}")
report_lines.append(f"Columns: {len(df.columns)}")
report_lines.append(
    f"Unique variants: {df['variant_id'].nunique()}"
)
report_lines.append(
    f"Data source: {ML_PROBLEM['current_data_source']}"
)

report_lines.append("")
report_lines.append("ML OBJECTIVE")
report_lines.append("-" * 70)
report_lines.append(ML_PROBLEM["primary_ml_goal"])

report_lines.append("")
report_lines.append("TARGET STATUS")
report_lines.append("-" * 70)
report_lines.append(ML_PROBLEM["current_target_status"])

report_lines.append("")
report_lines.append("IMPORTANT TARGET-LEAKAGE RULE")
report_lines.append("-" * 70)
report_lines.append(
    ML_PROBLEM["clinical_significance_status"]
)

report_lines.append("")
report_lines.append("MODEL TRAINING STATUS")
report_lines.append("-" * 70)
report_lines.append(
    ML_PROBLEM["training_status"]
)

report_lines.append("")
report_lines.append("FEATURE POLICY")
report_lines.append("-" * 70)
report_lines.append(ML_PROBLEM["feature_policy"])

report_lines.append("")
report_lines.append("EVALUATION POLICY")
report_lines.append("-" * 70)
report_lines.append(ML_PROBLEM["evaluation_policy"])

report_lines.append("")
report_lines.append("IDENTIFIER COLUMNS")
report_lines.append("-" * 70)

for column in identifier_columns:
    report_lines.append(f"- {column}")

report_lines.append("")
report_lines.append("RAW TEXT / REVIEW COLUMNS")
report_lines.append("-" * 70)

for column in raw_text_columns:
    report_lines.append(f"- {column}")

report_lines.append("")
report_lines.append("TARGET CANDIDATE COLUMNS")
report_lines.append("-" * 70)

for column in target_candidate_columns:
    report_lines.append(
        f"- {column} → NOT an independent ML target"
    )

report_lines.append("")
report_lines.append("CANDIDATE ML FEATURES")
report_lines.append("-" * 70)

for column in candidate_feature_columns:
    report_lines.append(f"- {column}")

report_lines.append("")
report_lines.append("FEATURES REQUIRING LEAKAGE ASSESSMENT")
report_lines.append("-" * 70)

for column in leakage_risk_columns:
    report_lines.append(f"- {column}")

report_lines.append("")
report_lines.append("CLINICAL SIGNIFICANCE DISTRIBUTION")
report_lines.append("-" * 70)

for category, count in clinical_significance_values.items():
    report_lines.append(
        f"- {category}: {count}"
    )

report_lines.append("")
report_lines.append("REVIEW STATUS DISTRIBUTION")
report_lines.append("-" * 70)

for category, count in review_status_values.items():
    report_lines.append(
        f"- {category}: {count}"
    )

report_lines.append("")
report_lines.append("FINAL STATUS")
report_lines.append("-" * 70)
report_lines.append(
    "PASS — ML problem definition documented."
)
report_lines.append(
    "PASS — No independent target was invented."
)
report_lines.append(
    "PASS — ClinVar clinical significance is not treated "
    "as an independent benchmark target."
)
report_lines.append(
    "PASS — Candidate features have been identified."
)
report_lines.append(
    "NEXT — Proceed to ML dataset preparation and feature audit."
)

report_lines.append("")
report_lines.append("=" * 70)

REPORT_FILE.write_text(
    "\n".join(report_lines),
    encoding="utf-8"
)


# ------------------------------------------------------------
# Save feature policy
# ------------------------------------------------------------

FEATURE_POLICY_FILE = (
    OUTPUT_DIR / "00_feature_policy_v1.csv"
)

feature_policy.to_csv(
    FEATURE_POLICY_FILE,
    index=False
)


# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

print("=" * 70)
print("PROJECT DEEPGENE")
print("STEP 00 — ML PROBLEM DEFINITION")
print("=" * 70)

print(f"Input dataset: {INPUT_FILE.name}")
print(f"Variants: {len(df)}")
print(f"Columns: {len(df.columns)}")

print("")
print("ML TARGET STATUS:")
print("NOT AVAILABLE — independent target not established.")

print("")
print("ClinVar clinical_significance:")
print("NOT used as an independent ML target.")

print("")
print(
    f"Candidate ML features: "
    f"{len(candidate_feature_columns)}"
)

print(
    f"Leakage-risk fields requiring assessment: "
    f"{len(leakage_risk_columns)}"
)

print("")
print(f"Report written to:")
print(REPORT_FILE)

print("")
print(f"Feature policy written to:")
print(FEATURE_POLICY_FILE)

print("")
print("STATUS: PASS")
print("=" * 70)