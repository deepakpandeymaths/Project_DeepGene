from pathlib import Path

import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    balanced_accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix,
)
from sklearn.model_selection import StratifiedKFold


# ============================================================
# ML STEP 16 — RANDOM FOREST WITH STRATIFIED CV
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
    / "16_random_forest_cv.txt"
)


TARGET = "functional_target"

CLASSES = [
    "Loss",
    "Mixed",
    "Gain",
]

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

N_SPLITS = 4
RANDOM_STATE = 42


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

if df[FEATURES].isna().any().any():
    raise ValueError(
        "Feature matrix contains missing values."
    )


# ------------------------------------------------------------
# Prepare X and y
# ------------------------------------------------------------

X = df[FEATURES].copy()
y = df[TARGET].copy()


# ------------------------------------------------------------
# Stratified cross-validation
# ------------------------------------------------------------

cv = StratifiedKFold(
    n_splits=N_SPLITS,
    shuffle=True,
    random_state=RANDOM_STATE,
)


# ------------------------------------------------------------
# Random Forest
# ------------------------------------------------------------
# Fixed parameters are used here.
# Hyperparameter tuning is intentionally NOT performed
# in this step.

model = RandomForestClassifier(
    n_estimators=300,
    random_state=RANDOM_STATE,
    class_weight="balanced",
    n_jobs=-1,
)


# ------------------------------------------------------------
# Storage
# ------------------------------------------------------------

fold_results = []

all_true = []
all_pred = []


# ------------------------------------------------------------
# Cross-validation
# ------------------------------------------------------------

for fold_number, (train_index, validation_index) in enumerate(
    cv.split(X, y),
    start=1,
):

    X_train = X.iloc[train_index]
    X_validation = X.iloc[validation_index]

    y_train = y.iloc[train_index]
    y_validation = y.iloc[validation_index]

    # Train only on the training fold
    model.fit(
        X_train,
        y_train,
    )

    # Predict validation fold
    y_prediction = model.predict(
        X_validation
    )

    # Fold metrics
    balanced_accuracy = balanced_accuracy_score(
        y_validation,
        y_prediction,
    )

    macro_precision = precision_score(
        y_validation,
        y_prediction,
        labels=CLASSES,
        average="macro",
        zero_division=0,
    )

    macro_recall = recall_score(
        y_validation,
        y_prediction,
        labels=CLASSES,
        average="macro",
        zero_division=0,
    )

    macro_f1 = f1_score(
        y_validation,
        y_prediction,
        labels=CLASSES,
        average="macro",
        zero_division=0,
    )

    fold_results.append(
        {
            "fold": fold_number,
            "balanced_accuracy": balanced_accuracy,
            "macro_precision": macro_precision,
            "macro_recall": macro_recall,
            "macro_f1": macro_f1,
        }
    )

    all_true.extend(
        y_validation.tolist()
    )

    all_pred.extend(
        y_prediction.tolist()
    )


# ------------------------------------------------------------
# Aggregate metrics
# ------------------------------------------------------------

results_df = pd.DataFrame(
    fold_results
)

mean_balanced_accuracy = (
    results_df["balanced_accuracy"].mean()
)

std_balanced_accuracy = (
    results_df["balanced_accuracy"].std(
        ddof=1
    )
)

mean_macro_precision = (
    results_df["macro_precision"].mean()
)

std_macro_precision = (
    results_df["macro_precision"].std(
        ddof=1
    )
)

mean_macro_recall = (
    results_df["macro_recall"].mean()
)

std_macro_recall = (
    results_df["macro_recall"].std(
        ddof=1
    )
)

mean_macro_f1 = (
    results_df["macro_f1"].mean()
)

std_macro_f1 = (
    results_df["macro_f1"].std(
        ddof=1
    )
)


# ------------------------------------------------------------
# Aggregate confusion matrix
# ------------------------------------------------------------

aggregate_confusion = confusion_matrix(
    all_true,
    all_pred,
    labels=CLASSES,
)


# ------------------------------------------------------------
# Write report
# ------------------------------------------------------------

lines = []

lines.append("PROJECT DEEPGENE — ML STEP 16")
lines.append("RANDOM FOREST WITH STRATIFIED CV")
lines.append("=" * 60)
lines.append("")

lines.append(
    f"Input dataset: {INPUT_FILE}"
)
lines.append(
    f"Number of variants: {len(df)}"
)
lines.append(
    f"Target: {TARGET}"
)
lines.append("")

lines.append("MODEL")
lines.append("-" * 60)
lines.append(
    "Algorithm: Random Forest Classifier"
)
lines.append(
    "Number of trees: 300"
)
lines.append(
    "Class weighting: balanced"
)
lines.append(
    "Hyperparameter tuning: No"
)
lines.append(
    f"Random state: {RANDOM_STATE}"
)
lines.append("")

lines.append("CROSS-VALIDATION")
lines.append("-" * 60)
lines.append(
    f"Folds: {N_SPLITS}"
)
lines.append(
    "Method: StratifiedKFold"
)
lines.append(
    "Shuffle: True"
)
lines.append(
    f"Random state: {RANDOM_STATE}"
)
lines.append("")

lines.append("FEATURES")
lines.append("-" * 60)

for feature in FEATURES:
    lines.append(feature)

lines.append("")

lines.append("FOLD RESULTS")
lines.append("-" * 60)

for result in fold_results:

    lines.append(
        f"Fold {result['fold']}: "
        f"Balanced accuracy="
        f"{result['balanced_accuracy']:.4f}, "
        f"Macro precision="
        f"{result['macro_precision']:.4f}, "
        f"Macro recall="
        f"{result['macro_recall']:.4f}, "
        f"Macro F1="
        f"{result['macro_f1']:.4f}"
    )

lines.append("")

lines.append("CROSS-VALIDATED SUMMARY")
lines.append("-" * 60)

lines.append(
    f"Balanced accuracy: "
    f"{mean_balanced_accuracy:.4f} "
    f"+/- {std_balanced_accuracy:.4f}"
)

lines.append(
    f"Macro precision:   "
    f"{mean_macro_precision:.4f} "
    f"+/- {std_macro_precision:.4f}"
)

lines.append(
    f"Macro recall:      "
    f"{mean_macro_recall:.4f} "
    f"+/- {std_macro_recall:.4f}"
)

lines.append(
    f"Macro F1:          "
    f"{mean_macro_f1:.4f} "
    f"+/- {std_macro_f1:.4f}"
)

lines.append("")

lines.append("AGGREGATE CONFUSION MATRIX")
lines.append("-" * 60)
lines.append(
    "Rows = actual classes"
)
lines.append(
    "Columns = predicted classes"
)
lines.append("")

lines.append(
    "             " + "  ".join(
        f"{c:>8}" for c in CLASSES
    )
)

for i, class_name in enumerate(CLASSES):

    row = aggregate_confusion[i]

    lines.append(
        f"{class_name:>10} "
        + "  ".join(
            f"{value:8d}"
            for value in row
        )
    )

lines.append("")

lines.append("INTERPRETATION")
lines.append("-" * 60)
lines.append(
    "This is the first tree-based supervised model."
)
lines.append(
    "The model uses only the 14 Step-42 biological features."
)
lines.append(
    "Class weighting is used because the target classes are imbalanced."
)
lines.append(
    "No hyperparameter tuning is performed in this step."
)
lines.append(
    "Performance should be compared with both the Step-13 "
    "majority baseline and Step-15 Logistic Regression."
)
lines.append(
    "Because only four Gain variants exist, Gain-specific "
    "performance has high uncertainty."
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