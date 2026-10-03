from pathlib import Path

import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    precision_score,
    recall_score,
    f1_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# ============================================================
# ML STEP 18 — REPEATED STRATIFIED CROSS-VALIDATION
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
    / "18_repeated_stratified_cv.txt"
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

REPEAT_SEEDS = [
    42,
    123,
    2024,
    7,
    99,
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

if df[FEATURES].isna().any().any():
    raise ValueError(
        "Feature matrix contains missing values."
    )


# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------

X = df[FEATURES].copy()
y = df[TARGET].copy()


# ------------------------------------------------------------
# Model definitions
# ------------------------------------------------------------

models = {
    "Logistic Regression": Pipeline(
        steps=[
            (
                "scaler",
                StandardScaler(),
            ),
            (
                "classifier",
                LogisticRegression(
                    max_iter=5000,
                    random_state=42,
                ),
            ),
        ]
    ),

    "Random Forest": RandomForestClassifier(
        n_estimators=300,
        class_weight="balanced",
        random_state=42,
        n_jobs=-1,
    ),
}


# ------------------------------------------------------------
# Repeated stratified CV
# ------------------------------------------------------------

results = []

for repeat_number, seed in enumerate(
    REPEAT_SEEDS,
    start=1,
):

    cv = StratifiedKFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=seed,
    )

    for fold_number, (
        train_index,
        validation_index,
    ) in enumerate(
        cv.split(X, y),
        start=1,
    ):

        X_train = X.iloc[train_index]
        X_validation = X.iloc[validation_index]

        y_train = y.iloc[train_index]
        y_validation = y.iloc[validation_index]

        for model_name, model in models.items():

            # Train only on current training fold
            model.fit(
                X_train,
                y_train,
            )

            # Predict validation fold
            y_prediction = model.predict(
                X_validation
            )

            balanced_accuracy = (
                balanced_accuracy_score(
                    y_validation,
                    y_prediction,
                )
            )

            macro_precision = (
                precision_score(
                    y_validation,
                    y_prediction,
                    labels=CLASSES,
                    average="macro",
                    zero_division=0,
                )
            )

            macro_recall = (
                recall_score(
                    y_validation,
                    y_prediction,
                    labels=CLASSES,
                    average="macro",
                    zero_division=0,
                )
            )

            macro_f1 = (
                f1_score(
                    y_validation,
                    y_prediction,
                    labels=CLASSES,
                    average="macro",
                    zero_division=0,
                )
            )

            results.append(
                {
                    "repeat": repeat_number,
                    "seed": seed,
                    "fold": fold_number,
                    "model": model_name,
                    "balanced_accuracy": balanced_accuracy,
                    "macro_precision": macro_precision,
                    "macro_recall": macro_recall,
                    "macro_f1": macro_f1,
                }
            )


# ------------------------------------------------------------
# Results dataframe
# ------------------------------------------------------------

results_df = pd.DataFrame(results)


# ------------------------------------------------------------
# Overall summaries
# ------------------------------------------------------------

summary = (
    results_df
    .groupby("model")
    [
        [
            "balanced_accuracy",
            "macro_precision",
            "macro_recall",
            "macro_f1",
        ]
    ]
    .agg(["mean", "std"])
)


# ------------------------------------------------------------
# Per-repeat summaries
# ------------------------------------------------------------

repeat_summary = (
    results_df
    .groupby(
        [
            "model",
            "repeat",
            "seed",
        ]
    )
    [
        [
            "balanced_accuracy",
            "macro_precision",
            "macro_recall",
            "macro_f1",
        ]
    ]
    .mean()
)


# ------------------------------------------------------------
# Write report
# ------------------------------------------------------------

lines = []

lines.append("PROJECT DEEPGENE — ML STEP 18")
lines.append("REPEATED STRATIFIED CROSS-VALIDATION")
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

lines.append("VALIDATION DESIGN")
lines.append("-" * 60)
lines.append(
    f"Stratified folds per repeat: {N_SPLITS}"
)
lines.append(
    f"Number of repeats: {len(REPEAT_SEEDS)}"
)
lines.append(
    f"Total validation evaluations per model: "
    f"{N_SPLITS * len(REPEAT_SEEDS)}"
)
lines.append(
    f"Seeds: {REPEAT_SEEDS}"
)
lines.append("")

lines.append("MODELS")
lines.append("-" * 60)
lines.append(
    "Logistic Regression"
)
lines.append(
    "Random Forest"
)
lines.append("")

lines.append("PER-REPEAT RESULTS")
lines.append("-" * 60)

for (
    model_name,
    repeat_number,
    seed,
), row in repeat_summary.iterrows():

    lines.append(
        f"{model_name} | "
        f"Repeat {repeat_number} | "
        f"Seed {seed}"
    )

    lines.append(
        f"  Balanced accuracy: "
        f"{row['balanced_accuracy']:.4f}"
    )

    lines.append(
        f"  Macro precision:   "
        f"{row['macro_precision']:.4f}"
    )

    lines.append(
        f"  Macro recall:      "
        f"{row['macro_recall']:.4f}"
    )

    lines.append(
        f"  Macro F1:          "
        f"{row['macro_f1']:.4f}"
    )

    lines.append("")


lines.append("OVERALL REPEATED-CV SUMMARY")
lines.append("-" * 60)

for model_name in models.keys():

    model_summary = summary.loc[
        model_name
    ]

    lines.append(
        f"{model_name}:"
    )

    lines.append(
        f"  Balanced accuracy: "
        f"{model_summary[('balanced_accuracy', 'mean')]:.4f} "
        f"+/- "
        f"{model_summary[('balanced_accuracy', 'std')]:.4f}"
    )

    lines.append(
        f"  Macro precision:   "
        f"{model_summary[('macro_precision', 'mean')]:.4f} "
        f"+/- "
        f"{model_summary[('macro_precision', 'std')]:.4f}"
    )

    lines.append(
        f"  Macro recall:      "
        f"{model_summary[('macro_recall', 'mean')]:.4f} "
        f"+/- "
        f"{model_summary[('macro_recall', 'std')]:.4f}"
    )

    lines.append(
        f"  Macro F1:          "
        f"{model_summary[('macro_f1', 'mean')]:.4f} "
        f"+/- "
        f"{model_summary[('macro_f1', 'std')]:.4f}"
    )

    lines.append("")


lines.append("INTERPRETATION")
lines.append("-" * 60)
lines.append(
    "Repeated stratified cross-validation evaluates model "
    "stability across multiple fold assignments."
)
lines.append(
    "The same 14 Step-42 biological features are used "
    "throughout."
)
lines.append(
    "Logistic Regression scaling is performed inside each "
    "training fold."
)
lines.append(
    "No hyperparameter tuning is performed."
)
lines.append(
    "No model is selected as a final model in this step."
)
lines.append(
    "The four Gain variants remain a major limitation on "
    "performance-estimate stability."
)
lines.append(
    "Results are exploratory and do not establish clinical "
    "validity or clinical predictive performance."
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