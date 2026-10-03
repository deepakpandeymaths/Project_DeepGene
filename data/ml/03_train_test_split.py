from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split


# ============================================================
# Project DeepGene
# Step 03: Train/Test Split
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_features_v1.csv"
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


# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------

TRAIN_FILE = (
    OUTPUT_DIR
    / "deepgene_ml_train_v1.csv"
)

TEST_FILE = (
    OUTPUT_DIR
    / "deepgene_ml_test_v1.csv"
)

SPLIT_IDS_FILE = (
    RESULTS_DIR
    / "03_train_test_split_ids_v1.csv"
)

REPORT_FILE = (
    RESULTS_DIR
    / "03_train_test_split.txt"
)


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

RANDOM_STATE = 42
TEST_SIZE = 0.20


# ------------------------------------------------------------
# Load dataset
# ------------------------------------------------------------

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Input dataset not found:\n{INPUT_FILE}"
    )

df = pd.read_csv(INPUT_FILE)


# ------------------------------------------------------------
# Basic validation
# ------------------------------------------------------------

if len(df) == 0:
    raise ValueError("Input dataset is empty.")

if df.isnull().sum().sum() != 0:
    raise ValueError(
        "Missing values detected in the ML feature dataset."
    )


# ------------------------------------------------------------
# IMPORTANT:
# The feature-selection dataset does not contain variant_id.
#
# Therefore, create a stable row-level identifier ONLY for
# partition tracking.
#
# This identifier is NOT a biological feature and will not be
# used by models.
# ------------------------------------------------------------

df = df.copy()

df.insert(
    0,
    "ml_row_id",
    range(1, len(df) + 1)
)


# ------------------------------------------------------------
# Create reproducible train/test split
# ------------------------------------------------------------

train_df, test_df = train_test_split(
    df,
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
    shuffle=True
)


# ------------------------------------------------------------
# Sort by row ID for reproducibility/readability
# ------------------------------------------------------------

train_df = train_df.sort_values(
    "ml_row_id"
).reset_index(drop=True)

test_df = test_df.sort_values(
    "ml_row_id"
).reset_index(drop=True)


# ------------------------------------------------------------
# Validate partition
# ------------------------------------------------------------

train_ids = set(train_df["ml_row_id"])
test_ids = set(test_df["ml_row_id"])

overlap = train_ids.intersection(test_ids)

if overlap:
    raise ValueError(
        "Train/test overlap detected."
    )


if len(train_df) + len(test_df) != len(df):
    raise ValueError(
        "Train + test rows do not equal original dataset."
    )


# ------------------------------------------------------------
# Create split manifest
# ------------------------------------------------------------

split_ids = pd.concat(
    [
        train_df[["ml_row_id"]].assign(
            split="train"
        ),
        test_df[["ml_row_id"]].assign(
            split="test"
        ),
    ],
    ignore_index=True
)

split_ids = split_ids.sort_values(
    "ml_row_id"
).reset_index(drop=True)


# ------------------------------------------------------------
# Save datasets
# ------------------------------------------------------------

train_df.to_csv(
    TRAIN_FILE,
    index=False
)

test_df.to_csv(
    TEST_FILE,
    index=False
)

split_ids.to_csv(
    SPLIT_IDS_FILE,
    index=False
)


# ------------------------------------------------------------
# Create report
# ------------------------------------------------------------

report = []

report.append("=" * 70)
report.append("PROJECT DEEPGENE")
report.append("STEP 03 — TRAIN/TEST SPLIT")
report.append("=" * 70)

report.append("")
report.append("INPUT")
report.append("-" * 70)
report.append(f"Input file: {INPUT_FILE}")
report.append(f"Total variants/rows: {len(df)}")
report.append(f"Total features: {len(df.columns) - 1}")

report.append("")
report.append("SPLIT CONFIGURATION")
report.append("-" * 70)
report.append(f"Random state: {RANDOM_STATE}")
report.append(f"Test size: {TEST_SIZE}")
report.append("Training fraction: 0.80")
report.append("Testing fraction: 0.20")

report.append("")
report.append("SPLIT RESULTS")
report.append("-" * 70)
report.append(f"Training rows: {len(train_df)}")
report.append(f"Testing rows: {len(test_df)}")

report.append("")
report.append("OVERLAP CHECK")
report.append("-" * 70)
report.append(f"Train/test overlapping rows: {len(overlap)}")

report.append("")
report.append("TARGET STATUS")
report.append("-" * 70)
report.append(
    "No independent ML target is currently available."
)
report.append(
    "This split is therefore a reproducible partition "
    "for future supervised learning and dataset auditing."
)

report.append("")
report.append("IMPORTANT")
report.append("-" * 70)
report.append(
    "The ml_row_id column is a partition-tracking identifier "
    "and must not be used as a model feature."
)
report.append(
    "Clinical significance is not included in this dataset "
    "and is not used as a target."
)

report.append("")
report.append("OUTPUTS")
report.append("-" * 70)
report.append(f"Training dataset: {TRAIN_FILE}")
report.append(f"Testing dataset: {TEST_FILE}")
report.append(f"Split manifest: {SPLIT_IDS_FILE}")

report.append("")
report.append("STATUS")
report.append("-" * 70)
report.append(
    "PASS — reproducible train/test partition created."
)
report.append(
    "PASS — no train/test overlap."
)
report.append(
    "PASS — all original rows accounted for."
)
report.append(
    "PASS — no target leakage introduced by this step."
)

report.append("")
report.append("=" * 70)


REPORT_FILE.write_text(
    "\n".join(report),
    encoding="utf-8"
)


# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

print("=" * 70)
print("PROJECT DEEPGENE")
print("STEP 03 — TRAIN/TEST SPLIT")
print("=" * 70)

print(f"Total rows: {len(df)}")
print(f"Training rows: {len(train_df)}")
print(f"Testing rows: {len(test_df)}")

print(
    f"Training percentage: "
    f"{len(train_df) / len(df) * 100:.2f}%"
)

print(
    f"Testing percentage: "
    f"{len(test_df) / len(df) * 100:.2f}%"
)

print(f"Train/test overlap: {len(overlap)}")

print("")
print("Training dataset:")
print(TRAIN_FILE)

print("")
print("Testing dataset:")
print(TEST_FILE)

print("")
print("Split manifest:")
print(SPLIT_IDS_FILE)

print("")
print("Report:")
print(REPORT_FILE)

print("")
print("STATUS: PASS")
print("=" * 70)