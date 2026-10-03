from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score


# ============================================================
# ML STEP 21
# LEARNING-CURVE / SAMPLE-SIZE STABILITY ANALYSIS
# ============================================================

INPUT_FILE = Path(
    "data/processed/deepgene_ml_ready_v1.csv"
)

REPORT_FILE = Path(
    "data/analysis_results/21_learning_curve_analysis.txt"
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

# Fractions of the available training data to use.
TRAIN_FRACTIONS = [
    0.50,
    0.75,
    1.00,
]


def make_logistic_model(seed):
    return Pipeline(
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


def make_random_forest_model(seed):
    return RandomForestClassifier(
        n_estimators=300,
        class_weight="balanced",
        random_state=seed,
        n_jobs=-1,
    )


def main():

    print("=" * 70)
    print("ML STEP 21 - LEARNING-CURVE / SAMPLE-SIZE STABILITY ANALYSIS")
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
    # 2. Keep primary classes
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
        errors="coerce",
    )

    if X.isna().any().any():
        raise ValueError(
            "Missing or non-numeric feature values detected."
        )

    y = df[TARGET].astype(str)

    # --------------------------------------------------------
    # 4. Class counts
    # --------------------------------------------------------

    class_counts = (
        y.value_counts()
        .reindex(CLASSES)
        .fillna(0)
        .astype(int)
    )

    # --------------------------------------------------------
    # 5. Results storage
    # --------------------------------------------------------

    results = []

    # --------------------------------------------------------
    # 6. Repeated stratified CV
    # --------------------------------------------------------

    for seed in SEEDS:

        print()
        print(f"Processing seed: {seed}")

        cv = StratifiedKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=seed,
        )

        for fold_number, (
            train_idx,
            valid_idx,
        ) in enumerate(
            cv.split(X, y),
            start=1,
        ):

            X_train_full = X.iloc[train_idx].copy()
            y_train_full = y.iloc[train_idx].copy()

            X_valid = X.iloc[valid_idx].copy()
            y_valid = y.iloc[valid_idx].copy()

            # ------------------------------------------------
            # For each requested training fraction
            # ------------------------------------------------

            for fraction in TRAIN_FRACTIONS:

                if fraction == 1.00:

                    subset_indices = np.arange(
                        len(X_train_full)
                    )

                else:

                    # Stratified subset:
                    # preserve class proportions as much
                    # as possible within each training fold.
                    subset_indices = []

                    rng = np.random.RandomState(
                        seed + fold_number
                    )

                    for class_name in CLASSES:

                        class_positions = np.where(
                            y_train_full.to_numpy()
                            == class_name
                        )[0]

                        desired_count = int(
                            np.floor(
                                len(class_positions)
                                * fraction
                            )
                        )

                        # Always keep at least one sample
                        # when the class exists.
                        desired_count = max(
                            1,
                            desired_count,
                        )

                        desired_count = min(
                            desired_count,
                            len(class_positions),
                        )

                        selected = rng.choice(
                            class_positions,
                            size=desired_count,
                            replace=False,
                        )

                        subset_indices.extend(
                            selected.tolist()
                        )

                    subset_indices = np.array(
                        subset_indices,
                        dtype=int,
                    )

                X_train = X_train_full.iloc[
                    subset_indices
                ]

                y_train = y_train_full.iloc[
                    subset_indices
                ]

                # --------------------------------------------
                # Logistic Regression
                # --------------------------------------------

                logistic_model = make_logistic_model(
                    seed
                )

                logistic_model.fit(
                    X_train,
                    y_train,
                )

                logistic_predictions = (
                    logistic_model.predict(
                        X_valid
                    )
                )

                logistic_balanced_accuracy = (
                    balanced_accuracy_score(
                        y_valid,
                        logistic_predictions,
                    )
                )

                results.append(
                    {
                        "model": "Logistic Regression",
                        "seed": seed,
                        "fold": fold_number,
                        "fraction": fraction,
                        "training_size": len(X_train),
                        "balanced_accuracy":
                            logistic_balanced_accuracy,
                    }
                )

                # --------------------------------------------
                # Random Forest
                # --------------------------------------------

                random_forest_model = (
                    make_random_forest_model(seed)
                )

                random_forest_model.fit(
                    X_train,
                    y_train,
                )

                rf_predictions = (
                    random_forest_model.predict(
                        X_valid
                    )
                )

                rf_balanced_accuracy = (
                    balanced_accuracy_score(
                        y_valid,
                        rf_predictions,
                    )
                )

                results.append(
                    {
                        "model": "Random Forest",
                        "seed": seed,
                        "fold": fold_number,
                        "fraction": fraction,
                        "training_size": len(X_train),
                        "balanced_accuracy":
                            rf_balanced_accuracy,
                    }
                )

            print(
                f"  Fold {fold_number}: completed"
            )

    # --------------------------------------------------------
    # 7. Results dataframe
    # --------------------------------------------------------

    results_df = pd.DataFrame(results)

    expected_rows = (
        len(SEEDS)
        * N_SPLITS
        * len(TRAIN_FRACTIONS)
        * 2
    )

    if len(results_df) != expected_rows:
        raise RuntimeError(
            "Unexpected number of learning-curve results."
        )

    # --------------------------------------------------------
    # 8. Aggregate results
    # --------------------------------------------------------

    summary = (
        results_df
        .groupby(
            [
                "model",
                "fraction",
                "training_size",
            ]
        )
        .agg(
            mean_balanced_accuracy=(
                "balanced_accuracy",
                "mean",
            ),
            std_balanced_accuracy=(
                "balanced_accuracy",
                "std",
            ),
            median_balanced_accuracy=(
                "balanced_accuracy",
                "median",
            ),
            min_balanced_accuracy=(
                "balanced_accuracy",
                "min",
            ),
            max_balanced_accuracy=(
                "balanced_accuracy",
                "max",
            ),
        )
        .reset_index()
    )

    # --------------------------------------------------------
    # 9. Report
    # --------------------------------------------------------

    report_lines = []

    def add(line=""):
        report_lines.append(line)

    add("=" * 70)
    add(
        "ML STEP 21 - LEARNING-CURVE / "
        "SAMPLE-SIZE STABILITY ANALYSIS"
    )
    add("=" * 70)
    add("")

    add("DATASET")
    add("-" * 70)
    add(f"Input file: {INPUT_FILE}")
    add(f"Variants: {len(df)}")
    add(f"Features: {len(FEATURES)}")
    add("")

    add("CLASS COUNTS")
    add("-" * 70)

    for class_name in CLASSES:
        add(
            f"{class_name}: "
            f"{class_counts[class_name]}"
        )

    add("")

    add("VALIDATION DESIGN")
    add("-" * 70)
    add(
        f"Repeated stratified CV: "
        f"{len(SEEDS)} repeats × {N_SPLITS} folds"
    )
    add(
        f"Training fractions: "
        f"{TRAIN_FRACTIONS}"
    )
    add(
        "Metric: balanced accuracy"
    )
    add(
        "Validation sets remain untouched while "
        "training subsets are reduced."
    )
    add("")

    # --------------------------------------------------------
    # 10. Detailed summary
    # --------------------------------------------------------

    for model_name in [
        "Logistic Regression",
        "Random Forest",
    ]:

        add("=" * 70)
        add(model_name)
        add("=" * 70)
        add("")

        model_summary = summary[
            summary["model"] == model_name
        ]

        for _, row in model_summary.iterrows():

            add(
                f"Training fraction: "
                f"{row['fraction']:.2f}"
            )

            add(
                f"Training size: "
                f"{int(row['training_size'])}"
            )

            add(
                f"Mean balanced accuracy: "
                f"{row['mean_balanced_accuracy']:.4f}"
            )

            add(
                f"SD balanced accuracy: "
                f"{row['std_balanced_accuracy']:.4f}"
            )

            add(
                f"Median balanced accuracy: "
                f"{row['median_balanced_accuracy']:.4f}"
            )

            add(
                f"Minimum balanced accuracy: "
                f"{row['min_balanced_accuracy']:.4f}"
            )

            add(
                f"Maximum balanced accuracy: "
                f"{row['max_balanced_accuracy']:.4f}"
            )

            add("")

    # --------------------------------------------------------
    # 11. Interpretation
    # --------------------------------------------------------

    add("=" * 70)
    add("INTERPRETATION")
    add("=" * 70)
    add("")

    add(
        "This analysis examines whether model performance "
        "changes systematically with the amount of training "
        "data available."
    )

    add(
        "It is not a substitute for independent validation."
    )

    add(
        "Because the complete dataset contains only 51 "
        "variants and the Gain class contains 4 variants, "
        "all learning-curve estimates are subject to substantial "
        "sampling variability."
    )

    add(
        "Training subsets were constructed within each "
        "training fold, so validation observations remain "
        "outside the training subset."
    )

    add(
        "No feature selection or hyperparameter tuning "
        "was performed in Step 21."
    )

    # --------------------------------------------------------
    # 12. Final status
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
    # 13. Terminal output
    # --------------------------------------------------------

    print()
    print("Learning-curve summary:")
    print()

    print(
        summary[
            [
                "model",
                "fraction",
                "training_size",
                "mean_balanced_accuracy",
                "std_balanced_accuracy",
                "median_balanced_accuracy",
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