from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal


# ============================================================
# ML STEP 19
# FEATURE DISTRIBUTION / CLASS-SEPARABILITY AUDIT
# ============================================================

INPUT_FILE = Path(
    "data/processed/deepgene_ml_ready_v1.csv"
)

REPORT_FILE = Path(
    "data/analysis_results/19_feature_class_separability.txt"
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


def main():

    print("=" * 70)
    print("ML STEP 19 - FEATURE CLASS-SEPARABILITY AUDIT")
    print("=" * 70)

    # --------------------------------------------------------
    # 1. Check input file
    # --------------------------------------------------------

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE}"
        )

    # --------------------------------------------------------
    # 2. Load dataset
    # --------------------------------------------------------

    df = pd.read_csv(INPUT_FILE)

    print(f"Input file: {INPUT_FILE}")
    print(f"Original rows: {len(df)}")
    print()

    # --------------------------------------------------------
    # 3. Check required columns
    # --------------------------------------------------------

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
    # 4. Keep only primary ML classes
    # --------------------------------------------------------

    df = df[
        df[TARGET].isin(CLASSES)
    ].copy()

    expected_rows = 51

    if len(df) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} primary ML variants, "
            f"found {len(df)}."
        )

    # --------------------------------------------------------
    # 5. Convert features to numeric
    # --------------------------------------------------------

    X = df[FEATURES].apply(
        pd.to_numeric,
        errors="coerce"
    )

    if X.isna().any().any():

        bad_columns = X.columns[
            X.isna().any()
        ].tolist()

        raise ValueError(
            "Missing or non-numeric values detected in:\n"
            + "\n".join(bad_columns)
        )

    # --------------------------------------------------------
    # 6. Class counts
    # --------------------------------------------------------

    class_counts = (
        df[TARGET]
        .value_counts()
        .reindex(CLASSES)
        .fillna(0)
        .astype(int)
    )

    # --------------------------------------------------------
    # 7. Prepare report
    # --------------------------------------------------------

    report_lines = []

    def add(line=""):
        report_lines.append(line)

    add("=" * 70)
    add("ML STEP 19 - FEATURE CLASS-SEPARABILITY AUDIT")
    add("=" * 70)
    add("")

    add("INPUT")
    add("-" * 70)
    add(f"Input file: {INPUT_FILE}")
    add(f"Rows analyzed: {len(df)}")
    add(f"Features analyzed: {len(FEATURES)}")
    add("")

    add("PURPOSE")
    add("-" * 70)
    add(
        "Audit the distributions of the 14 biological features "
        "across the Loss, Mixed, and Gain functional-effect classes."
    )
    add(
        "This is an exploratory descriptive analysis."
    )
    add(
        "No features are removed or selected based on this audit."
    )
    add(
        "No model tuning is performed in this step."
    )
    add("")

    # --------------------------------------------------------
    # 8. Class counts
    # --------------------------------------------------------

    add("CLASS COUNTS")
    add("-" * 70)

    for class_name in CLASSES:
        add(
            f"{class_name}: "
            f"{class_counts[class_name]}"
        )

    add("")

    # --------------------------------------------------------
    # 9. Feature analysis
    # --------------------------------------------------------

    results = []

    for feature in FEATURES:

        add("=" * 70)
        add(f"FEATURE: {feature}")
        add("=" * 70)

        values_by_class = []

        class_statistics = {}

        # ----------------------------------------------------
        # Calculate descriptive statistics
        # ----------------------------------------------------

        for class_name in CLASSES:

            values = X.loc[
                df[TARGET] == class_name,
                feature
            ].to_numpy(dtype=float)

            values_by_class.append(values)

            q1 = np.percentile(values, 25)
            median = np.percentile(values, 50)
            q3 = np.percentile(values, 75)

            statistics = {
                "n": len(values),
                "mean": np.mean(values),
                "std": np.std(values, ddof=1),
                "median": median,
                "q1": q1,
                "q3": q3,
                "iqr": q3 - q1,
                "min": np.min(values),
                "max": np.max(values),
            }

            class_statistics[class_name] = statistics

            add(
                f"{class_name}: "
                f"n={statistics['n']}, "
                f"mean={statistics['mean']:.4f}, "
                f"std={statistics['std']:.4f}, "
                f"median={statistics['median']:.4f}, "
                f"IQR={statistics['iqr']:.4f}, "
                f"min={statistics['min']:.4f}, "
                f"max={statistics['max']:.4f}"
            )

        # ----------------------------------------------------
        # Kruskal-Wallis test
        # ----------------------------------------------------
        #
        # Rank-based exploratory comparison across the
        # three functional classes.
        #
        # Constant features are handled explicitly.
        # ----------------------------------------------------

        all_values = np.concatenate(values_by_class)

        if np.all(
            all_values == all_values[0]
        ):
            h_statistic = 0.0
            p_value = 1.0

        else:
            try:

                h_statistic, p_value = kruskal(
                    *values_by_class
                )

            except (ValueError, RuntimeWarning):

                h_statistic = np.nan
                p_value = np.nan

        # ----------------------------------------------------
        # Exploratory epsilon-squared effect size
        # ----------------------------------------------------

        number_of_classes = len(CLASSES)
        number_of_samples = len(df)

        if (
            np.isfinite(h_statistic)
            and number_of_samples > number_of_classes
        ):

            epsilon_squared = (
                h_statistic
                - number_of_classes
                + 1
            ) / (
                number_of_samples
                - number_of_classes
            )

            epsilon_squared = max(
                0.0,
                epsilon_squared
            )

        else:

            epsilon_squared = np.nan

        # ----------------------------------------------------
        # Store results
        # ----------------------------------------------------

        results.append(
            {
                "feature": feature,
                "kruskal_h": h_statistic,
                "p_value": p_value,
                "epsilon_squared": epsilon_squared,
            }
        )

        # ----------------------------------------------------
        # Add statistical results to report
        # ----------------------------------------------------

        if np.isfinite(h_statistic):

            add(
                f"Kruskal-Wallis H: "
                f"{h_statistic:.4f}"
            )

        else:

            add(
                "Kruskal-Wallis H: NA"
            )

        if np.isfinite(p_value):

            add(
                f"Kruskal-Wallis p-value: "
                f"{p_value:.6f}"
            )

        else:

            add(
                "Kruskal-Wallis p-value: NA"
            )

        if np.isfinite(epsilon_squared):

            add(
                f"Epsilon-squared: "
                f"{epsilon_squared:.4f}"
            )

        else:

            add(
                "Epsilon-squared: NA"
            )

        add("")

    # --------------------------------------------------------
    # 10. Results dataframe
    # --------------------------------------------------------

    results_df = pd.DataFrame(results)

    if len(results_df) != len(FEATURES):

        raise RuntimeError(
            "Number of feature results does not match "
            "number of features."
        )

    # --------------------------------------------------------
    # 11. Identify undefined statistical tests
    # --------------------------------------------------------

    invalid_tests = results_df[
        ~np.isfinite(
            results_df["p_value"].to_numpy(
                dtype=float
            )
        )
    ]

    # --------------------------------------------------------
    # 12. Summary
    # --------------------------------------------------------

    add("=" * 70)
    add("SUMMARY")
    add("=" * 70)
    add("")

    add(
        "Kruskal-Wallis statistics are exploratory "
        "rank-based comparisons across the three "
        "functional-effect classes."
    )

    add(
        "The Gain class contains only 4 variants, "
        "so statistical estimates involving Gain "
        "are highly uncertain."
    )

    add(
        "Several features are mathematically related, "
        "including delta and absolute-delta features. "
        "This is consistent with the earlier redundancy audit."
    )

    add(
        "No feature was removed or selected based on "
        "this audit."
    )

    add(
        "No model was trained or tuned in Step 19."
    )

    if len(invalid_tests) > 0:

        add("")
        add("UNDEFINED STATISTICAL TESTS")
        add("-" * 70)

        for feature in invalid_tests["feature"]:

            add(
                f"{feature}: "
                "Kruskal-Wallis statistic undefined"
            )

        add(
            "These features were retained for descriptive analysis."
        )

    # --------------------------------------------------------
    # 13. Compact results table
    # --------------------------------------------------------

    add("")
    add("COMPACT STATISTICAL TABLE")
    add("-" * 70)

    for _, row in results_df.iterrows():

        h_text = (
            f"{row['kruskal_h']:.4f}"
            if np.isfinite(row["kruskal_h"])
            else "NA"
        )

        p_text = (
            f"{row['p_value']:.6f}"
            if np.isfinite(row["p_value"])
            else "NA"
        )

        effect_text = (
            f"{row['epsilon_squared']:.4f}"
            if np.isfinite(
                row["epsilon_squared"]
            )
            else "NA"
        )

        add(
            f"{row['feature']}: "
            f"H={h_text}, "
            f"p={p_text}, "
            f"epsilon²={effect_text}"
        )

    # --------------------------------------------------------
    # 14. Final status
    # --------------------------------------------------------

    add("")
    add("=" * 70)
    add("STATUS: PASS")
    add("=" * 70)

    # --------------------------------------------------------
    # 15. Save report
    # --------------------------------------------------------

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(
        REPORT_FILE,
        "w",
        encoding="utf-8"
    ) as report_file:

        report_file.write(
            "\n".join(report_lines)
        )

    # --------------------------------------------------------
    # 16. Terminal output
    # --------------------------------------------------------

    print()
    print("Rows analyzed:", len(df))
    print(
        "Features analyzed:",
        len(FEATURES)
    )

    print()
    print("Class counts:")

    for class_name in CLASSES:

        print(
            f"  {class_name}: "
            f"{class_counts[class_name]}"
        )

    print()
    print(
        "Top-level Kruskal-Wallis results:"
    )

    print(
        results_df[
            [
                "feature",
                "kruskal_h",
                "p_value",
                "epsilon_squared",
            ]
        ].to_string(index=False)
    )

    if len(invalid_tests) > 0:

        print()
        print(
            "Undefined statistical tests:"
        )

        for feature in invalid_tests[
            "feature"
        ]:

            print(
                f"  {feature}"
            )

    print()
    print(
        f"Report saved to: {REPORT_FILE}"
    )

    print("STATUS: PASS")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    main()