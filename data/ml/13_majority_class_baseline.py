from pathlib import Path

import pandas as pd
from sklearn.metrics import (
    balanced_accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix,
)


# ============================================================
# ML STEP 13 — MAJORITY-CLASS BASELINE
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
    / "13_majority_class_baseline.txt"
)


TARGET = "functional_target"

EXPECTED_CLASSES = [
    "Loss",
    "Mixed",
    "Gain",
]


# ------------------------------------------------------------
# Load dataset
# ------------------------------------------------------------

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"Dataset not found:\n{INPUT_FILE}"
    )

df = pd.read_csv(INPUT_FILE)


if TARGET not in df.columns:
    raise ValueError(
        f"Target column '{TARGET}' not found."
    )


# ------------------------------------------------------------
# Target distribution
# ------------------------------------------------------------

class_counts = df[TARGET].value_counts()

majority_class = class_counts.idxmax()
majority_count = class_counts.max()

n_samples = len(df)


# ------------------------------------------------------------
# Create majority-class predictions
# ------------------------------------------------------------

y_true = df[TARGET]

y_pred = [majority_class] * n_samples


# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------

balanced_accuracy = balanced_accuracy_score(
    y_true,
    y_pred,
)

macro_precision = precision_score(
    y_true,
    y_pred,
    labels=EXPECTED_CLASSES,
    average="macro",
    zero_division=0,
)

macro_recall = recall_score(
    y_true,
    y_pred,
    labels=EXPECTED_CLASSES,
    average="macro",
    zero_division=0,
)

macro_f1 = f1_score(
    y_true,
    y_pred,
    labels=EXPECTED_CLASSES,
    average="macro",
    zero_division=0,
)

confusion = confusion_matrix(
    y_true,
    y_pred,
    labels=EXPECTED_CLASSES,
)


# ------------------------------------------------------------
# Write report
# ------------------------------------------------------------

lines = []

lines.append("PROJECT DEEPGENE — ML STEP 13")
lines.append("MAJORITY-CLASS BASELINE")
lines.append("=" * 60)
lines.append("")

lines.append(f"Input dataset: {INPUT_FILE}")
lines.append(f"Number of variants: {n_samples}")
lines.append(f"Target: {TARGET}")
lines.append("")

lines.append("TARGET DISTRIBUTION")
lines.append("-" * 60)

for class_name in EXPECTED_CLASSES:
    count = class_counts.get(class_name, 0)
    percentage = (count / n_samples) * 100

    lines.append(
        f"{class_name}: {count} "
        f"({percentage:.2f}%)"
    )

lines.append("")

lines.append("BASELINE")
lines.append("-" * 60)
lines.append(
    f"Majority class: {majority_class}"
)
lines.append(
    f"Majority count: {majority_count}"
)
lines.append(
    f"Majority-class accuracy: "
    f"{majority_count / n_samples:.4f}"
)
lines.append("")

lines.append("BASELINE METRICS")
lines.append("-" * 60)
lines.append(
    f"Balanced accuracy: {balanced_accuracy:.4f}"
)
lines.append(
    f"Macro precision:   {macro_precision:.4f}"
)
lines.append(
    f"Macro recall:      {macro_recall:.4f}"
)
lines.append(
    f"Macro F1:          {macro_f1:.4f}"
)
lines.append("")

lines.append("CONFUSION MATRIX")
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
        f"{c:>8}" for c in EXPECTED_CLASSES
    )
)

for i, class_name in enumerate(EXPECTED_CLASSES):
    row = confusion[i]

    lines.append(
        f"{class_name:>10} "
        + "  ".join(
            f"{value:8d}" for value in row
        )
    )

lines.append("")

lines.append("INTERPRETATION")
lines.append("-" * 60)
lines.append(
    "This baseline predicts the majority class "
    "for every variant."
)
lines.append(
    "It provides the reference point for subsequent "
    "supervised models."
)
lines.append(
    "No biological features are used by this baseline."
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