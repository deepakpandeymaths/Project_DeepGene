from pathlib import Path
import pandas as pd


# ============================================================
# ML STEP 12 — STEP-42 HANDOFF CHECK
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "data" / "processed" / "deepgene_ml_ready_v1.csv"
REPORT_FILE = PROJECT_ROOT / "data" / "analysis_results" / "12_ml_handoff_check.txt"


# ------------------------------------------------------------
# Expected target
# ------------------------------------------------------------

TARGET = "functional_target"


# ------------------------------------------------------------
# Expected biological features from Step 42
# ------------------------------------------------------------

EXPECTED_FEATURES = [
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


# ------------------------------------------------------------
# Expected primary target classes
# ------------------------------------------------------------

EXPECTED_CLASSES = {
    "Loss",
    "Mixed",
    "Gain",
}


# ------------------------------------------------------------
# Load dataset
# ------------------------------------------------------------

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Step-42 ML dataset not found:\n{INPUT_FILE}"
    )

df = pd.read_csv(INPUT_FILE)


# ------------------------------------------------------------
# Basic dataset information
# ------------------------------------------------------------

rows, columns = df.shape


# ------------------------------------------------------------
# Validate target
# ------------------------------------------------------------

if TARGET not in df.columns:
    raise ValueError(
        f"Required target column '{TARGET}' not found."
    )

target_values = set(df[TARGET].dropna().unique())

unexpected_classes = target_values - EXPECTED_CLASSES

if unexpected_classes:
    raise ValueError(
        f"Unexpected target classes found: {sorted(unexpected_classes)}"
    )


# ------------------------------------------------------------
# Validate features
# ------------------------------------------------------------

missing_features = [
    feature
    for feature in EXPECTED_FEATURES
    if feature not in df.columns
]

if missing_features:
    raise ValueError(
        "Missing required Step-42 features:\n"
        + "\n".join(missing_features)
    )


# ------------------------------------------------------------
# Check missing values
# ------------------------------------------------------------

feature_missing = df[EXPECTED_FEATURES].isna().sum()

target_missing = df[TARGET].isna().sum()


# ------------------------------------------------------------
# Check finite numeric values
# ------------------------------------------------------------

numeric_features = df[EXPECTED_FEATURES].apply(
    pd.to_numeric,
    errors="coerce"
)

non_numeric_values = numeric_features.isna().sum().sum()


# ------------------------------------------------------------
# Target distribution
# ------------------------------------------------------------

target_distribution = df[TARGET].value_counts()


# ------------------------------------------------------------
# Write report
# ------------------------------------------------------------

report_lines = []

report_lines.append("PROJECT DEEPGENE — ML STEP 12")
report_lines.append("STEP-42 ML HANDOFF CHECK")
report_lines.append("=" * 60)
report_lines.append("")

report_lines.append(f"Input file: {INPUT_FILE}")
report_lines.append(f"Rows: {rows}")
report_lines.append(f"Columns: {columns}")
report_lines.append("")

report_lines.append("TARGET")
report_lines.append("-" * 60)
report_lines.append(f"Target column: {TARGET}")
report_lines.append(
    f"Target missing values: {target_missing}"
)
report_lines.append(
    f"Target classes: {sorted(target_values)}"
)
report_lines.append("")

report_lines.append("TARGET DISTRIBUTION")
report_lines.append("-" * 60)

for class_name, count in target_distribution.items():
    report_lines.append(f"{class_name}: {count}")

report_lines.append("")

report_lines.append("EXPECTED BIOLOGICAL FEATURES")
report_lines.append("-" * 60)

for feature in EXPECTED_FEATURES:
    report_lines.append(feature)

report_lines.append("")

report_lines.append("FEATURE VALIDATION")
report_lines.append("-" * 60)
report_lines.append(
    f"Missing required features: {len(missing_features)}"
)
report_lines.append(
    f"Missing feature values: {int(feature_missing.sum())}"
)
report_lines.append(
    f"Non-numeric/invalid feature values: {int(non_numeric_values)}"
)
report_lines.append("")

report_lines.append("ML HANDOFF DECISION")
report_lines.append("-" * 60)

if (
    not missing_features
    and target_missing == 0
    and feature_missing.sum() == 0
    and non_numeric_values == 0
    and target_values == EXPECTED_CLASSES
):
    report_lines.append("STATUS: PASS")
    report_lines.append(
        "Step-42 dataset is ready for baseline and stratified CV."
    )
else:
    report_lines.append("STATUS: FAIL")
    report_lines.append(
        "Step-42 dataset requires correction before ML training."
    )

REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)
REPORT_FILE.write_text(
    "\n".join(report_lines),
    encoding="utf-8"
)

print("\n".join(report_lines))