from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# ============================================================
# DEEPGENE ML STEP 24
# Error-pattern analysis by functional-effect class
# ============================================================

INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
OUTPUT_FILE = Path("data/analysis_results/24_error_pattern_analysis.txt")

TARGET = "functional_target"

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

SEEDS = [42, 123, 2024, 7, 99]
N_SPLITS = 4

CLASSES = ["Loss", "Mixed", "Gain"]


def make_models():
    """Create the two controlled models used in Steps 18–23."""
    return {
        "Logistic Regression": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
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


def safe_mean(series):
    """Return a numeric mean or NaN for empty input."""
    if len(series) == 0:
        return np.nan
    return float(series.mean())


def format_number(value):
    """Format numeric values consistently."""
    if pd.isna(value):
        return "NA"
    return f"{float(value):.4f}"


def main():
    # --------------------------------------------------------
    # Load and validate dataset
    # --------------------------------------------------------
    df = pd.read_csv(INPUT_FILE)

    required_columns = FEATURES + [TARGET]
    missing = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    if df[required_columns].isnull().any().any():
        raise ValueError(
            "Missing values detected in target/features."
        )

    X = df[FEATURES].copy()
    y = df[TARGET].copy()

    observed_classes = sorted(y.unique())

    if observed_classes != sorted(CLASSES):
        raise ValueError(
            f"Unexpected target classes: {observed_classes}"
        )

    # --------------------------------------------------------
    # Collect out-of-fold predictions
    # --------------------------------------------------------
    prediction_records = []

    for seed in SEEDS:
        cv = StratifiedKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=seed,
        )

        for fold_number, (train_idx, valid_idx) in enumerate(
            cv.split(X, y),
            start=1,
        ):
            X_train = X.iloc[train_idx]
            X_valid = X.iloc[valid_idx]

            y_train = y.iloc[train_idx]

            models = make_models()

            for model_name, model in models.items():
                model.fit(X_train, y_train)

                predictions = model.predict(X_valid)

                for row_index, prediction in zip(
                    valid_idx,
                    predictions,
                ):
                    prediction_records.append(
                        {
                            "row_index": int(row_index),
                            "seed": seed,
                            "fold": fold_number,
                            "model": model_name,
                            "actual": y.iloc[row_index],
                            "predicted": prediction,
                        }
                    )

    pred_df = pd.DataFrame(prediction_records)

    # --------------------------------------------------------
    # Attach biological features to prediction records
    #
    # Each variant appears 20 times per model because of the
    # 5 seeds × 4 folds design.
    # --------------------------------------------------------
    feature_df = df[FEATURES].copy()
    feature_df["row_index"] = feature_df.index

    pred_df = pred_df.merge(
        feature_df,
        on="row_index",
        how="left",
        validate="many_to_one",
    )

    # --------------------------------------------------------
    # Define error categories
    # --------------------------------------------------------
    pred_df["correct"] = (
        pred_df["actual"] == pred_df["predicted"]
    )

    pred_df["error_type"] = np.where(
        pred_df["correct"],
        "Correct",
        pred_df["actual"] + " -> " + pred_df["predicted"],
    )

    # --------------------------------------------------------
    # Output report
    # --------------------------------------------------------
    lines = []

    lines.append(
        "DEEPGENE ML STEP 24 — ERROR-PATTERN ANALYSIS"
    )
    lines.append("=" * 72)
    lines.append("")
    lines.append("Purpose:")
    lines.append(
        "Characterize repeated out-of-fold classification errors "
        "by functional-effect class and examine whether error "
        "patterns correspond to different biological feature "
        "distributions."
    )
    lines.append(
        "This is an exploratory error analysis and does not "
        "establish clinical validity."
    )
    lines.append("")

    # --------------------------------------------------------
    # Dataset information
    # --------------------------------------------------------
    lines.append("INPUT")
    lines.append("-" * 72)
    lines.append(f"File: {INPUT_FILE}")
    lines.append(f"Rows: {len(df)}")
    lines.append(f"Features: {len(FEATURES)}")
    lines.append(f"Target: {TARGET}")
    lines.append(f"Seeds: {SEEDS}")
    lines.append(f"Folds per seed: {N_SPLITS}")
    lines.append(
        f"Validation evaluations per model: "
        f"{len(SEEDS) * N_SPLITS}"
    )
    lines.append("")

    lines.append("TARGET DISTRIBUTION")
    lines.append("-" * 72)

    target_counts = y.value_counts()

    for class_name in CLASSES:
        count = int(target_counts.get(class_name, 0))
        lines.append(
            f"{class_name}: {count} "
            f"({count / len(y):.4f})"
        )

    lines.append("")

    # ========================================================
    # MODEL-BY-MODEL ERROR ANALYSIS
    # ========================================================
    for model_name in [
        "Logistic Regression",
        "Random Forest",
    ]:
        model_df = pred_df[
            pred_df["model"] == model_name
        ].copy()

        lines.append("")
        lines.append(model_name.upper())
        lines.append("=" * 72)

        # ----------------------------------------------------
        # Overall counts
        # ----------------------------------------------------
        lines.append("")
        lines.append("OVERALL PREDICTION COUNTS")
        lines.append("-" * 72)

        actual_counts = model_df["actual"].value_counts()
        predicted_counts = model_df["predicted"].value_counts()

        for class_name in CLASSES:
            lines.append(
                f"{class_name}: "
                f"actual={int(actual_counts.get(class_name, 0))}, "
                f"predicted={int(predicted_counts.get(class_name, 0))}"
            )

        lines.append("")

        # ----------------------------------------------------
        # Confusion matrix
        # ----------------------------------------------------
        cm = confusion_matrix(
            model_df["actual"],
            model_df["predicted"],
            labels=CLASSES,
        )

        lines.append("AGGREGATE CONFUSION MATRIX")
        lines.append("-" * 72)
        lines.append(
            "Rows = actual; columns = predicted"
        )
        lines.append(
            "             "
            + "  ".join(
                f"{class_name:>7}"
                for class_name in CLASSES
            )
        )

        for class_name, row in zip(CLASSES, cm):
            lines.append(
                f"{class_name:>12} "
                + "  ".join(
                    f"{int(value):7d}"
                    for value in row
                )
            )

        lines.append("")

        # ----------------------------------------------------
        # Error counts by direction
        # ----------------------------------------------------
        lines.append("ERROR DIRECTIONS")
        lines.append("-" * 72)

        errors = model_df[
            model_df["actual"] != model_df["predicted"]
        ]

        error_counts = (
            errors.groupby(
                ["actual", "predicted"]
            )
            .size()
            .sort_values(ascending=False)
        )

        if len(error_counts) == 0:
            lines.append("No errors detected.")
        else:
            for (actual, predicted), count in error_counts.items():
                lines.append(
                    f"{actual} -> {predicted}: {int(count)}"
                )

        lines.append("")

        # ----------------------------------------------------
        # Per-class recall from repeated OOF predictions
        # ----------------------------------------------------
        lines.append(
            "CLASS-WISE REPEATED OUT-OF-FOLD RECALL"
        )
        lines.append("-" * 72)

        for class_name in CLASSES:
            actual_class = model_df[
                model_df["actual"] == class_name
            ]

            correct_class = actual_class[
                actual_class["predicted"] == class_name
            ]

            recall = (
                len(correct_class) / len(actual_class)
                if len(actual_class) > 0
                else np.nan
            )

            lines.append(
                f"{class_name}: "
                f"{len(correct_class)}/{len(actual_class)} "
                f"= {format_number(recall)}"
            )

        lines.append("")

        # ----------------------------------------------------
        # Biological feature distributions by actual class
        # ----------------------------------------------------
        lines.append(
            "BIOLOGICAL FEATURE MEANS BY ACTUAL CLASS"
        )
        lines.append("-" * 72)

        class_feature_means = (
            df.groupby(TARGET)[FEATURES]
            .mean()
            .reindex(CLASSES)
        )

        for feature in FEATURES:
            values = []

            for class_name in CLASSES:
                value = class_feature_means.loc[
                    class_name,
                    feature,
                ]

                values.append(
                    f"{class_name}={format_number(value)}"
                )

            lines.append(
                f"{feature}: " + ", ".join(values)
            )

        lines.append("")

        # ----------------------------------------------------
        # Feature means for correct vs incorrect predictions
        # ----------------------------------------------------
        lines.append(
            "FEATURE MEANS: CORRECT VS INCORRECT PREDICTIONS"
        )
        lines.append("-" * 72)

        correct_df = model_df[
            model_df["correct"]
        ]

        incorrect_df = model_df[
            ~model_df["correct"]
        ]

        for feature in FEATURES:
            correct_mean = safe_mean(
                correct_df[feature]
            )

            incorrect_mean = safe_mean(
                incorrect_df[feature]
            )

            lines.append(
                f"{feature}: "
                f"correct={format_number(correct_mean)}, "
                f"incorrect={format_number(incorrect_mean)}"
            )

        lines.append("")

        # ----------------------------------------------------
        # Feature means by error direction
        # ----------------------------------------------------
        lines.append(
            "FEATURE MEANS BY ERROR DIRECTION"
        )
        lines.append("-" * 72)

        if len(errors) == 0:
            lines.append(
                "No error directions available."
            )
        else:
            unique_error_types = (
                errors["error_type"]
                .value_counts()
                .index
                .tolist()
            )

            for error_type in unique_error_types:
                error_subset = errors[
                    errors["error_type"] == error_type
                ]

                lines.append("")
                lines.append(
                    f"{error_type} "
                    f"(n={len(error_subset)})"
                )

                for feature in FEATURES:
                    mean_value = safe_mean(
                        error_subset[feature]
                    )

                    lines.append(
                        f"  {feature}: "
                        f"{format_number(mean_value)}"
                    )

        lines.append("")

        # ----------------------------------------------------
        # Variant-level error persistence
        # ----------------------------------------------------
        lines.append(
            "VARIANT-LEVEL ERROR PERSISTENCE"
        )
        lines.append("-" * 72)

        variant_summary = (
            model_df
            .groupby(["row_index", "actual"])
            .agg(
                total_predictions=("correct", "size"),
                correct_predictions=("correct", "sum"),
            )
            .reset_index()
        )

        variant_summary["error_predictions"] = (
            variant_summary["total_predictions"]
            - variant_summary["correct_predictions"]
        )

        variant_summary["error_fraction"] = (
            variant_summary["error_predictions"]
            / variant_summary["total_predictions"]
        )

        for class_name in CLASSES:
            subset = variant_summary[
                variant_summary["actual"] == class_name
            ]

            if len(subset) == 0:
                continue

            lines.append(
                f"{class_name}: "
                f"variants={len(subset)}, "
                f"mean_error_fraction="
                f"{subset['error_fraction'].mean():.4f}, "
                f"median_error_fraction="
                f"{subset['error_fraction'].median():.4f}, "
                f"always_wrong="
                f"{int((subset['error_fraction'] == 1.0).sum())}"
            )

        lines.append("")

    # ========================================================
    # Cross-model comparison
    # ========================================================
    lines.append("")
    lines.append("CROSS-MODEL ERROR COMPARISON")
    lines.append("=" * 72)
    lines.append("")

    model_summaries = []

    for model_name in [
        "Logistic Regression",
        "Random Forest",
    ]:
        model_df = pred_df[
            pred_df["model"] == model_name
        ]

        errors = model_df[
            ~model_df["correct"]
        ]

        summary = {
            "model": model_name,
            "total_predictions": len(model_df),
            "correct_predictions": int(
                model_df["correct"].sum()
            ),
            "incorrect_predictions": len(errors),
        }

        model_summaries.append(summary)

    for summary in model_summaries:
        accuracy = (
            summary["correct_predictions"]
            / summary["total_predictions"]
        )

        lines.append(
            f"{summary['model']}: "
            f"correct={summary['correct_predictions']}, "
            f"incorrect={summary['incorrect_predictions']}, "
            f"overall_repeated_OOF_accuracy="
            f"{accuracy:.4f}"
        )

    lines.append("")

    # --------------------------------------------------------
    # Shared error directions
    # --------------------------------------------------------
    lines.append("SHARED ERROR DIRECTIONS")
    lines.append("-" * 72)

    error_direction_tables = {}

    for model_name in [
        "Logistic Regression",
        "Random Forest",
    ]:
        model_df = pred_df[
            pred_df["model"] == model_name
        ]

        errors = model_df[
            model_df["actual"] != model_df["predicted"]
        ]

        error_direction_tables[model_name] = set(
            errors["error_type"].unique()
        )

    shared_errors = (
        error_direction_tables["Logistic Regression"]
        & error_direction_tables["Random Forest"]
    )

    if shared_errors:
        for error_type in sorted(shared_errors):
            lines.append(error_type)
    else:
        lines.append(
            "No shared error directions."
        )

    lines.append("")

    # ========================================================
    # Interpretation
    # ========================================================
    lines.append("INTERPRETATION")
    lines.append("=" * 72)
    lines.append("")
    lines.append(
        "This analysis identifies the directions and persistence "
        "of classification errors across repeated out-of-fold "
        "validation."
    )
    lines.append(
        "Feature means are descriptive only. Differences between "
        "correct and incorrect groups are not treated as evidence "
        "of predictive significance."
    )
    lines.append(
        "Repeated error patterns involving Gain or Mixed must be "
        "interpreted cautiously because Gain contains only 4 "
        "variants and Mixed contains only 12."
    )
    lines.append(
        "The analysis does not perform feature selection, "
        "hyperparameter tuning, or model selection."
    )
    lines.append(
        "No causal biological interpretation is assigned to "
        "observed feature differences."
    )
    lines.append(
        "No clinical interpretation or pathogenicity/benignity "
        "conclusion is made."
    )
    lines.append(
        "No final model is selected."
    )
    lines.append("")

    # --------------------------------------------------------
    # Status
    # --------------------------------------------------------
    lines.append("STATUS")
    lines.append("-" * 72)
    lines.append("PASS")
    lines.append(
        "Error-pattern analysis completed successfully."
    )

    # --------------------------------------------------------
    # Save report
    # --------------------------------------------------------
    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print("\n".join(lines))


if __name__ == "__main__":
    main()