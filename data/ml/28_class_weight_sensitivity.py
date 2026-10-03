from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, precision_score, recall_score, f1_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# ============================================================
# ML STEP 28
# Class-Weight Sensitivity Audit
# ============================================================

INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
REPORT_FILE = Path("data/analysis_results/28_class_weight_sensitivity.txt")

RANDOM_SEEDS = [42, 123, 2024, 7, 99]
N_SPLITS = 4

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

TARGET = "functional_target"


def evaluate_model(model, X, y):
    results = []

    for seed in RANDOM_SEEDS:
        cv = StratifiedKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=seed,
        )

        for fold_number, (train_idx, val_idx) in enumerate(
            cv.split(X, y), start=1
        ):
            X_train = X.iloc[train_idx]
            X_val = X.iloc[val_idx]
            y_train = y.iloc[train_idx]
            y_val = y.iloc[val_idx]

            model.fit(X_train, y_train)

            predictions = model.predict(X_val)

            results.append(
                {
                    "seed": seed,
                    "fold": fold_number,
                    "balanced_accuracy": balanced_accuracy_score(
                        y_val, predictions
                    ),
                    "macro_precision": precision_score(
                        y_val,
                        predictions,
                        average="macro",
                        zero_division=0,
                    ),
                    "macro_recall": recall_score(
                        y_val,
                        predictions,
                        average="macro",
                        zero_division=0,
                    ),
                    "macro_f1": f1_score(
                        y_val,
                        predictions,
                        average="macro",
                        zero_division=0,
                    ),
                }
            )

    return pd.DataFrame(results)


def summarize(results):
    metrics = [
        "balanced_accuracy",
        "macro_precision",
        "macro_recall",
        "macro_f1",
    ]

    summary = {}

    for metric in metrics:
        values = results[metric].to_numpy()

        summary[metric] = {
            "mean": float(np.mean(values)),
            "sd": float(np.std(values, ddof=1)),
            "median": float(np.median(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
        }

    return summary


def main():

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT_FILE)

    X = df[FEATURES].copy()
    y = df[TARGET].copy()

    class_counts = y.value_counts().sort_index()

    models = {
        "Logistic Regression - unweighted": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=5000,
                        class_weight=None,
                        random_state=42,
                    ),
                ),
            ]
        ),
        "Logistic Regression - balanced": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=42,
                    ),
                ),
            ]
        ),
        "Random Forest - unweighted": RandomForestClassifier(
            n_estimators=300,
            class_weight=None,
            random_state=42,
            n_jobs=-1,
        ),
        "Random Forest - balanced": RandomForestClassifier(
            n_estimators=300,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
    }

    all_results = {}
    all_summaries = {}

    for model_name, model in models.items():
        results = evaluate_model(model, X, y)

        all_results[model_name] = results
        all_summaries[model_name] = summarize(results)

    lines = []

    lines.append("=" * 70)
    lines.append("DEEPGENE ML STEP 28")
    lines.append("CLASS-WEIGHT SENSITIVITY AUDIT")
    lines.append("=" * 70)
    lines.append("")

    lines.append("PURPOSE")
    lines.append("-" * 70)
    lines.append(
        "Assess whether model performance is sensitive to class weighting."
    )
    lines.append(
        "This is a sensitivity analysis, not model selection."
    )
    lines.append(
        "No hyperparameter tuning was performed."
    )
    lines.append("")

    lines.append("INPUT")
    lines.append("-" * 70)
    lines.append(f"Dataset: {INPUT_FILE}")
    lines.append(f"Rows: {len(df)}")
    lines.append(f"Features: {len(FEATURES)}")
    lines.append(f"Target: {TARGET}")
    lines.append("")

    lines.append("TARGET DISTRIBUTION")
    lines.append("-" * 70)

    for cls, count in class_counts.items():
        lines.append(
            f"{cls}: {count} "
            f"({count / len(y):.4f})"
        )

    lines.append("")

    lines.append("CROSS-VALIDATION DESIGN")
    lines.append("-" * 70)
    lines.append(
        f"StratifiedKFold: {N_SPLITS} folds"
    )
    lines.append(
        f"Seeds: {RANDOM_SEEDS}"
    )
    lines.append(
        f"Total validation evaluations per model: "
        f"{len(RANDOM_SEEDS) * N_SPLITS}"
    )
    lines.append(
        "Logistic Regression: StandardScaler fitted inside each training fold."
    )
    lines.append(
        "Random Forest: no feature standardization."
    )
    lines.append("")

    lines.append("MODEL CONFIGURATIONS")
    lines.append("-" * 70)
    lines.append(
        "Logistic Regression: max_iter=5000, random_state=42"
    )
    lines.append(
        "Random Forest: n_estimators=300, random_state=42, n_jobs=-1"
    )
    lines.append(
        "Compared class_weight values: None vs balanced"
    )
    lines.append("")

    for model_name in models:
        lines.append("=" * 70)
        lines.append(model_name.upper())
        lines.append("=" * 70)

        summary = all_summaries[model_name]

        for metric, values in summary.items():
            lines.append(
                f"{metric}: "
                f"mean={values['mean']:.6f}, "
                f"SD={values['sd']:.6f}, "
                f"median={values['median']:.6f}, "
                f"min={values['min']:.6f}, "
                f"max={values['max']:.6f}"
            )

        lines.append("")

    # --------------------------------------------------------
    # Direct sensitivity comparison
    # --------------------------------------------------------

    lines.append("=" * 70)
    lines.append("CLASS-WEIGHT SENSITIVITY COMPARISON")
    lines.append("=" * 70)
    lines.append("")

    pairs = [
        (
            "Logistic Regression",
            "Logistic Regression - unweighted",
            "Logistic Regression - balanced",
        ),
        (
            "Random Forest",
            "Random Forest - unweighted",
            "Random Forest - balanced",
        ),
    ]

    for family, unweighted_name, balanced_name in pairs:

        lines.append(f"{family}")
        lines.append("-" * 70)

        unweighted = all_summaries[unweighted_name]
        balanced = all_summaries[balanced_name]

        for metric in [
            "balanced_accuracy",
            "macro_precision",
            "macro_recall",
            "macro_f1",
        ]:
            difference = (
                balanced[metric]["mean"]
                - unweighted[metric]["mean"]
            )

            lines.append(
                f"{metric}: "
                f"unweighted={unweighted[metric]['mean']:.6f}, "
                f"balanced={balanced[metric]['mean']:.6f}, "
                f"balanced-minus-unweighted={difference:.6f}"
            )

        lines.append("")

    # --------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------

    lines.append("=" * 70)
    lines.append("INTERPRETATION")
    lines.append("=" * 70)
    lines.append(
        "Class weighting changes the training objective by giving greater "
        "weight to minority classes."
    )
    lines.append(
        "This audit determines whether conclusions are sensitive to that "
        "choice."
    )
    lines.append(
        "Because the Gain class contains only 4 variants, class-weight "
        "sensitivity is particularly important."
    )
    lines.append(
        "Differences between weighted and unweighted models are descriptive "
        "and should not be interpreted as evidence that one configuration "
        "is the final model."
    )
    lines.append(
        "Repeated cross-validation variability remains important because "
        "the dataset contains only 51 variants."
    )
    lines.append("")

    lines.append("DECISION")
    lines.append("-" * 70)
    lines.append(
        "No model selected in Step 28."
    )
    lines.append(
        "No feature selection performed."
    )
    lines.append(
        "No hyperparameter tuning performed."
    )
    lines.append(
        "No clinical interpretation made."
    )
    lines.append("")

    lines.append("=" * 70)
    lines.append("STATUS: PASS")
    lines.append("=" * 70)

    REPORT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print("\n".join(lines))


if __name__ == "__main__":
    main()