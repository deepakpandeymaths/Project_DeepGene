from pathlib import Path

import pandas as pd
from sklearn.model_selection import StratifiedKFold


# ============================================================
# ML STEP 14 — STRATIFIED CROSS-VALIDATION SETUP
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_ready_v1.csv"
)

REPORT_FILE = (
    PROJECT_ROOT
    / "data"
    / "analysis_results"
    / "14_stratified_cv_setup.txt"
)


TARGET = "functional_target"

N_SPLITS = 4
RANDOM_STATE = 42

FEATURES = [
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
# Load dataset
# ------------------------------------------------------------

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Dataset not found:\n{INPUT_FILE}"
    )

df = pd.read_csv(INPUT_FILE)


# ------------------------------------------------------------
# Validate dataset
# ------------------------------------------------------------

if TARGET not in df.columns:
    raise ValueError(
        f"Target column '{TARGET}' not found."
    )

missing_features = [
    feature for feature in FEATURES
    if feature not in df.columns
]

if missing_features:
    raise ValueError(
        f"Missing required features: {missing_features}"
    )

if df[TARGET].isna().any():
    raise ValueError(
        "Target contains missing values."
    )


# ------------------------------------------------------------
# Target distribution
# ------------------------------------------------------------

class_counts = df[TARGET].value_counts()

smallest_class_count = class_counts.min()

if N_SPLITS > smallest_class_count:
    raise ValueError(
        f"N_SPLITS={N_SPLITS} exceeds the smallest "
        f"class count ({smallest_class_count})."
    )


# ------------------------------------------------------------
# Create stratified folds
# ------------------------------------------------------------

X = df[FEATURES]
y = df[TARGET]

skf = StratifiedKFold(
    n_splits=N_SPLITS,
    shuffle=True,
    random_state=RANDOM_STATE,
)


# ------------------------------------------------------------
# Validate fold distributions
# ------------------------------------------------------------

fold_records = []

for fold_number, (train_index, validation_index) in enumerate(
    skf.split(X, y),
    start=1,
):

    train_targets = y.iloc[train_index]
    validation_targets = y.iloc[validation_index]

    fold_records.append(
        {
            "fold": fold_number,
            "train_size": len(train_index),
            "validation_size": len(validation_index),
            "train_loss": int(
                (train_targets == "Loss").sum()
            ),
            "train_mixed": int(
                (train_targets == "Mixed").sum()
            ),
            "train_gain": int(
                (train_targets == "Gain").sum()
            ),
            "validation_loss": int(
                (validation_targets == "Loss").sum()
            ),
            "validation_mixed": int(
                (validation_targets == "Mixed").sum()
            ),
            "validation_gain": int(
                (validation_targets == "Gain").sum()
            ),
        }
    )


# ------------------------------------------------------------
# Write report
# ------------------------------------------------------------

lines = []

lines.append("PROJECT DEEPGENE — ML STEP 14")
lines.append("STRATIFIED CROSS-VALIDATION SETUP")
lines.append("=" * 60)
lines.append("")

lines.append(f"Input dataset: {INPUT_FILE}")
lines.append(f"Number of variants: {len(df)}")
lines.append(f"Target: {TARGET}")
lines.append("")

lines.append("CROSS-VALIDATION DESIGN")
lines.append("-" * 60)
lines.append(
    f"Number of folds: {N_SPLITS}"
)
lines.append(
    f"Shuffle: True"
)
lines.append(
    f"Random state: {RANDOM_STATE}"
)
lines.append(
    "Validation strategy: StratifiedKFold"
)
lines.append("")

lines.append("TARGET DISTRIBUTION")
lines.append("-" * 60)

for class_name in ["Loss", "Mixed", "Gain"]:
    count = class_counts.get(class_name, 0)

    lines.append(
        f"{class_name}: {count}"
    )

lines.append(
    f"Smallest class count: {smallest_class_count}"
)
lines.append("")

lines.append("FOLD DISTRIBUTIONS")
lines.append("-" * 60)

for record in fold_records:

    lines.append(
        f"Fold {record['fold']}: "
        f"train={record['train_size']}, "
        f"validation={record['validation_size']}"
    )

    lines.append(
        "  Train: "
        f"Loss={record['train_loss']}, "
        f"Mixed={record['train_mixed']}, "
        f"Gain={record['train_gain']}"
    )

    lines.append(
        "  Validation: "
        f"Loss={record['validation_loss']}, "
        f"Mixed={record['validation_mixed']}, "
        f"Gain={record['validation_gain']}"
    )

    lines.append("")

lines.append("FEATURE SET")
lines.append("-" * 60)

for feature in FEATURES:
    lines.append(feature)

lines.append("")

lines.append("INTERPRETATION")
lines.append("-" * 60)
lines.append(
    "Four-fold stratified cross-validation is used because "
    "the smallest target class contains four variants."
)
lines.append(
    "Stratification preserves class representation across "
    "training and validation folds as far as the sample size permits."
)
lines.append(
    "The fixed random state makes the fold assignment reproducible."
)
lines.append(
    "This step defines the validation framework; no supervised "
    "model is trained in this step."
)

lines.append("")
lines.append("STATUS: PASS")


REPORT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)

REPORT_FILE.write_text(
    "\n".join(lines),
    encoding="utf-8",
)


print("\n".join(lines))