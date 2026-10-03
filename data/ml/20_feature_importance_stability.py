from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold


# ============================================================
# ML STEP 20
# FEATURE-IMPORTANCE STABILITY AUDIT
# ============================================================

INPUT_FILE = Path(
    "data/processed/deepgene_ml_ready_v1.csv"
)

REPORT_FILE = Path(
    "data/analysis_results/20_feature_importance_stability.txt"
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

SEEDS = [
    42,
    123,
    2024,
    7,
    99,
]

N_SPLITS = 4

N_PERMUTATIONS = 20


def main():

    print("=" * 70)
    print("ML STEP 20 - FEATURE-IMPORTANCE STABILITY AUDIT")
    print("=" * 70)

    # --------------------------------------------------------
    # 1. Load data
    # --------------------------------------------------------

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE}"
        )

    df = pd.read_csv(INPUT_FILE)

    required_columns = [TARGET] + FEATURES

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required columns:\n"
            + "\n".join(missing_columns)
        )

    # --------------------------------------------------------
    # 2. Restrict to primary target classes
    # --------------------------------------------------------

    df = df[
        df[TARGET].isin(CLASSES)
    ].copy()

    if len(df) != 51:
        raise ValueError(
            f"Expected 51 primary ML variants, "
            f"found {len(df)}."
        )

    # --------------------------------------------------------
    # 3. Prepare X and y
    # --------------------------------------------------------

    X = df[FEATURES].apply(
        pd.to_numeric,
        errors="coerce"
    )

    if X.isna().any().any():
        raise ValueError(
            "Missing or non-numeric feature values detected."
        )

    y = df[TARGET].astype(str)

    # --------------------------------------------------------
    # 4. Storage for importance results
    # --------------------------------------------------------

    logistic_results = []
    random_forest_results = []

    # --------------------------------------------------------
    # 5. Repeated stratified CV
    # --------------------------------------------------------

    for seed in SEEDS:

        print()
        print(f"Processing seed: {seed}")

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
            y_valid = y.iloc[valid_idx]

            # ------------------------------------------------
            # Logistic Regression
            # ------------------------------------------------

            logistic_model = Pipeline(
                [
                    (
                        "scaler",
                        StandardScaler(),
                    ),
                    (
                        "model",
                        LogisticRegression(
                            max_iter=5000,
                            random_state=seed,
                        ),
                    ),
                ]
            )

            logistic_model.fit(
                X_train,
                y_train,
            )

            logistic_permutation = permutation_importance(
                logistic_model,
                X_valid,
                y_valid,
                scoring="balanced_accuracy",
                n_repeats=N_PERMUTATIONS,
                random_state=seed,
                n_jobs=-1,
            )

            for feature, mean_importance, std_importance in zip(
                FEATURES,
                logistic_permutation.importances_mean,
                logistic_permutation.importances_std,
            ):

                logistic_results.append(
                    {
                        "seed": seed,
                        "fold": fold_number,
                        "feature": feature,
                        "importance_mean": mean_importance,
                        "importance_std": std_importance,
                    }
                )

            # ------------------------------------------------
            # Random Forest
            # ------------------------------------------------

            random_forest_model = RandomForestClassifier(
                n_estimators=300,
                class_weight="balanced",
                random_state=seed,
                n_jobs=-1,
            )

            random_forest_model.fit(
                X_train,
                y_train,
            )

            rf_permutation = permutation_importance(
                random_forest_model,
                X_valid,
                y_valid,
                scoring="balanced_accuracy",
                n_repeats=N_PERMUTATIONS,
                random_state=seed,
                n_jobs=-1,
            )

            for feature, mean_importance, std_importance in zip(
                FEATURES,
                rf_permutation.importances_mean,
                rf_permutation.importances_std,
            ):

                random_forest_results.append(
                    {
                        "seed": seed,
                        "fold": fold_number,
                        "feature": feature,
                        "importance_mean": mean_importance,
                        "importance_std": std_importance,
                    }
                )

            print(
                f"  Fold {fold_number}: completed"
            )

    # --------------------------------------------------------
    # 6. Convert to dataframes
    # --------------------------------------------------------

    logistic_df = pd.DataFrame(
        logistic_results
    )

    random_forest_df = pd.DataFrame(
        random_forest_results
    )

    expected_rows = (
        len(SEEDS)
        * N_SPLITS
        * len(FEATURES)
    )

    if len(logistic_df) != expected_rows:
        raise RuntimeError(
            "Unexpected Logistic Regression "
            "importance result count."
        )

    if len(random_forest_df) != expected_rows:
        raise RuntimeError(
            "Unexpected Random Forest "
            "importance result count."
        )

    # --------------------------------------------------------
    # 7. Aggregate Logistic Regression importance
    # --------------------------------------------------------

    logistic_summary = (
        logistic_df
        .groupby("feature")
        .agg(
            mean_importance=(
                "importance_mean",
                "mean",
            ),
            std_across_evaluations=(
                "importance_mean",
                "std",
            ),
            median_importance=(
                "importance_mean",
                "median",
            ),
            positive_fraction=(
                "importance_mean",
                lambda values: np.mean(
                    values > 0
                ),
            ),
        )
        .reset_index()
    )

    # --------------------------------------------------------
    # 8. Aggregate Random Forest importance
    # --------------------------------------------------------

    rf_summary = (
        random_forest_df
        .groupby("feature")
        .agg(
            mean_importance=(
                "importance_mean",
                "mean",
            ),
            std_across_evaluations=(
                "importance_mean",
                "std",
            ),
            median_importance=(
                "importance_mean",
                "median",
            ),
            positive_fraction=(
                "importance_mean",
                lambda values: np.mean(
                    values > 0
                ),
            ),
        )
        .reset_index()
    )

    # --------------------------------------------------------
    # 9. Create report
    # --------------------------------------------------------

    report_lines = []

    def add(line=""):
        report_lines.append(line)

    add("=" * 70)
    add("ML STEP 20 - FEATURE-IMPORTANCE STABILITY AUDIT")
    add("=" * 70)
    add("")

    add("DATASET")
    add("-" * 70)
    add(f"Input file: {INPUT_FILE}")
    add(f"Variants: {len(df)}")
    add(f"Features: {len(FEATURES)}")
    add(f"Classes: {', '.join(CLASSES)}")
    add("")

    add("VALIDATION DESIGN")
    add("-" * 70)
    add(
        f"Repeated stratified CV: "
        f"{len(SEEDS)} repeats × {N_SPLITS} folds"
    )
    add(
        f"Total validation folds: "
        f"{len(SEEDS) * N_SPLITS}"
    )
    add(
        f"Permutation repeats per fold: "
        f"{N_PERMUTATIONS}"
    )
    add(
        "Permutation importance metric: "
        "balanced accuracy"
    )
    add("")

    # --------------------------------------------------------
    # 10. Logistic Regression report
    # --------------------------------------------------------

    add("=" * 70)
    add("LOGISTIC REGRESSION")
    add("=" * 70)
    add("")

    add(
        "Mean permutation importance across all "
        "20 validation folds:"
    )
    add("")

    logistic_summary = logistic_summary.sort_values(
        "mean_importance",
        ascending=False,
    )

    for _, row in logistic_summary.iterrows():

        add(
            f"{row['feature']}: "
            f"mean={row['mean_importance']:.6f}, "
            f"std={row['std_across_evaluations']:.6f}, "
            f"median={row['median_importance']:.6f}, "
            f"positive_fraction="
            f"{row['positive_fraction']:.3f}"
        )

    # --------------------------------------------------------
    # 11. Random Forest report
    # --------------------------------------------------------

    add("")
    add("=" * 70)
    add("RANDOM FOREST")
    add("=" * 70)
    add("")

    add(
        "Mean permutation importance across all "
        "20 validation folds:"
    )
    add("")

    rf_summary = rf_summary.sort_values(
        "mean_importance",
        ascending=False,
    )

    for _, row in rf_summary.iterrows():

        add(
            f"{row['feature']}: "
            f"mean={row['mean_importance']:.6f}, "
            f"std={row['std_across_evaluations']:.6f}, "
            f"median={row['median_importance']:.6f}, "
            f"positive_fraction="
            f"{row['positive_fraction']:.3f}"
        )

    # --------------------------------------------------------
    # 12. Stability interpretation
    # --------------------------------------------------------

    add("")
    add("=" * 70)
    add("STABILITY INTERPRETATION")
    add("=" * 70)
    add("")

    add(
        "Positive permutation importance means that "
        "permuting the feature tended to reduce "
        "balanced accuracy on validation data."
    )

    add(
        "Importance values near zero indicate little "
        "measurable contribution under this validation setup."
    )

    add(
        "Negative permutation importance can occur when "
        "permuting a feature happens to improve validation "
        "performance or due to sampling variability."
    )

    add(
        "Because the dataset contains only 51 variants and "
        "the Gain class contains 4 variants, importance "
        "estimates may be unstable."
    )

    add(
        "Feature importance from this audit is descriptive "
        "and is not sufficient by itself to justify feature "
        "removal or a final model."
    )

    add(
        "No feature was removed from the ML feature set "
        "in Step 20."
    )

    # --------------------------------------------------------
    # 13. Save report
    # --------------------------------------------------------

    add("")
    add("=" * 70)
    add("STATUS: PASS")
    add("=" * 70)

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(
        REPORT_FILE,
        "w",
        encoding="utf-8",
    ) as report_file:

        report_file.write(
            "\n".join(report_lines)
        )

    # --------------------------------------------------------
    # 14. Terminal summary
    # --------------------------------------------------------

    print()
    print(
        "Logistic Regression permutation importance:"
    )

    print(
        logistic_summary[
            [
                "feature",
                "mean_importance",
                "std_across_evaluations",
                "positive_fraction",
            ]
        ].to_string(index=False)
    )

    print()
    print(
        "Random Forest permutation importance:"
    )

    print(
        rf_summary[
            [
                "feature",
                "mean_importance",
                "std_across_evaluations",
                "positive_fraction",
            ]
        ].to_string(index=False)
    )

    print()
    print(
        f"Report saved to: {REPORT_FILE}"
    )

    print("STATUS: PASS")


if __name__ == "__main__":
    main()