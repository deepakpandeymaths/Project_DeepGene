from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# ML STEP 26 — FEATURE VARIANCE AUDIT
# ============================================================

INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
OUTPUT_FILE = Path("data/analysis_results/26_ml_feature_variance_audit.txt")

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

# Thresholds used for descriptive screening only.
# These do NOT automatically remove features.
ZERO_VARIANCE_THRESHOLD = 0.0
NEAR_ZERO_VARIANCE_THRESHOLD = 1e-6


def main():

    print("=" * 70)
    print("ML STEP 26 — FEATURE VARIANCE AUDIT")
    print("=" * 70)

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE}"
        )

    df = pd.read_csv(INPUT_FILE)

    missing_features = [
        feature for feature in FEATURES
        if feature not in df.columns
    ]

    if missing_features:
        raise ValueError(
            f"Missing expected features: {missing_features}"
        )

    X = df[FEATURES].copy()

    # --------------------------------------------------------
    # Basic validation
    # --------------------------------------------------------

    if X.isna().any().any():
        missing_counts = X.isna().sum()
        raise ValueError(
            f"Missing feature values detected:\n{missing_counts}"
        )

    non_numeric = [
        feature
        for feature in FEATURES
        if not pd.api.types.is_numeric_dtype(X[feature])
    ]

    if non_numeric:
        raise ValueError(
            f"Non-numeric features detected: {non_numeric}"
        )

    # --------------------------------------------------------
    # Variance statistics
    # --------------------------------------------------------

    records = []

    for feature in FEATURES:

        values = X[feature]

        variance = values.var(ddof=1)
        standard_deviation = values.std(ddof=1)
        minimum = values.min()
        maximum = values.max()
        mean = values.mean()
        unique_count = values.nunique(dropna=False)

        value_range = maximum - minimum

        if variance == 0:
            variance_status = "ZERO_VARIANCE"
        elif variance <= NEAR_ZERO_VARIANCE_THRESHOLD:
            variance_status = "NEAR_ZERO_VARIANCE"
        else:
            variance_status = "VARIABLE"

        records.append(
            {
                "feature": feature,
                "unique_values": unique_count,
                "mean": mean,
                "variance": variance,
                "standard_deviation": standard_deviation,
                "minimum": minimum,
                "maximum": maximum,
                "range": value_range,
                "variance_status": variance_status,
            }
        )

    variance_df = pd.DataFrame(records)

    # --------------------------------------------------------
    # Identify zero / near-zero variance features
    # --------------------------------------------------------

    zero_variance = variance_df[
        variance_df["variance"] <= ZERO_VARIANCE_THRESHOLD
    ].copy()

    near_zero_variance = variance_df[
        (variance_df["variance"] > ZERO_VARIANCE_THRESHOLD)
        & (
            variance_df["variance"]
            <= NEAR_ZERO_VARIANCE_THRESHOLD
        )
    ].copy()

    variable_features = variance_df[
        variance_df["variance"] > NEAR_ZERO_VARIANCE_THRESHOLD
    ].copy()

    # --------------------------------------------------------
    # Sort by variance
    # --------------------------------------------------------

    variance_sorted = variance_df.sort_values(
        "variance",
        ascending=True
    )

    # --------------------------------------------------------
    # Build report
    # --------------------------------------------------------

    lines = []

    lines.append("ML STEP 26 — FEATURE VARIANCE AUDIT")
    lines.append("=" * 70)
    lines.append("")

    lines.append(
        "Purpose: audit feature variance and identify constant or "
        "near-zero-variance biological ML features."
    )

    lines.append(
        "This is a feature-structure audit only. It does not establish "
        "predictive importance, biological significance, or clinical validity."
    )

    lines.append("")

    lines.append("INPUT")
    lines.append("-" * 70)
    lines.append(f"File: {INPUT_FILE}")
    lines.append(f"Rows: {len(df)}")
    lines.append(f"Features audited: {len(FEATURES)}")
    lines.append("")

    lines.append("VARIANCE THRESHOLDS")
    lines.append("-" * 70)
    lines.append(
        f"Zero variance threshold: {ZERO_VARIANCE_THRESHOLD}"
    )
    lines.append(
        f"Near-zero variance threshold: "
        f"{NEAR_ZERO_VARIANCE_THRESHOLD}"
    )
    lines.append("")

    lines.append("FEATURE VARIANCE SUMMARY")
    lines.append("-" * 70)

    for _, row in variance_sorted.iterrows():

        lines.append(
            f"{row['feature']}: "
            f"unique={int(row['unique_values'])}, "
            f"mean={row['mean']:.6f}, "
            f"variance={row['variance']:.12f}, "
            f"SD={row['standard_deviation']:.6f}, "
            f"min={row['minimum']:.6f}, "
            f"max={row['maximum']:.6f}, "
            f"range={row['range']:.6f}, "
            f"status={row['variance_status']}"
        )

    lines.append("")

    lines.append("ZERO-VARIANCE FEATURES")
    lines.append("-" * 70)

    if len(zero_variance) == 0:
        lines.append("None")
    else:
        for _, row in zero_variance.iterrows():
            lines.append(
                f"- {row['feature']} "
                f"(unique values={int(row['unique_values'])})"
            )

    lines.append("")

    lines.append("NEAR-ZERO-VARIANCE FEATURES")
    lines.append("-" * 70)

    if len(near_zero_variance) == 0:
        lines.append("None")
    else:
        for _, row in near_zero_variance.iterrows():
            lines.append(
                f"- {row['feature']} "
                f"(variance={row['variance']:.12f})"
            )

    lines.append("")

    lines.append("VARIABLE FEATURES")
    lines.append("-" * 70)

    if len(variable_features) == 0:
        lines.append("None")
    else:
        for _, row in variable_features.iterrows():
            lines.append(
                f"- {row['feature']} "
                f"(variance={row['variance']:.12f})"
            )

    lines.append("")

    lines.append("INTERPRETATION")
    lines.append("-" * 70)

    if len(zero_variance) > 0:
        lines.append(
            f"{len(zero_variance)} feature(s) have zero variance and "
            "therefore contain no within-dataset variation."
        )
    else:
        lines.append(
            "No zero-variance features were detected."
        )

    if len(near_zero_variance) > 0:
        lines.append(
            f"{len(near_zero_variance)} feature(s) fall within the "
            "near-zero-variance screening threshold."
        )
    else:
        lines.append(
            "No additional near-zero-variance features were detected."
        )

    lines.append(
        "A constant feature is not informative for distinguishing "
        "observations within this dataset."
    )

    lines.append(
        "Variance magnitude alone does not determine predictive value; "
        "features on different biological scales can legitimately have "
        "different variances."
    )

    lines.append(
        "No feature was automatically removed in this step."
    )

    lines.append(
        "No model was trained, tuned, or selected in this step."
    )

    lines.append(
        "No clinical interpretation is made."
    )

    lines.append("")

    lines.append("STATUS: PASS")
    lines.append("=" * 70)

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    OUTPUT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8"
    )

    print("\n".join(lines))

    print("")
    print(f"Report saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()