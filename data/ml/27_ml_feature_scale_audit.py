from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# ML STEP 27 — FEATURE SCALE AND STANDARDIZATION AUDIT
# ============================================================

INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
OUTPUT_FILE = Path("data/analysis_results/27_ml_feature_scale_audit.txt")

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
    print("ML STEP 27 — FEATURE SCALE AND STANDARDIZATION AUDIT")
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
    # Validation
    # --------------------------------------------------------

    if X.isna().any().any():
        raise ValueError(
            "Missing values detected in ML features."
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
    # Raw feature scale statistics
    # --------------------------------------------------------

    records = []

    for feature in FEATURES:

        values = X[feature]

        mean = values.mean()
        std = values.std(ddof=1)
        minimum = values.min()
        maximum = values.max()
        feature_range = maximum - minimum

        records.append(
            {
                "feature": feature,
                "mean": mean,
                "std": std,
                "minimum": minimum,
                "maximum": maximum,
                "range": feature_range,
            }
        )

    scale_df = pd.DataFrame(records)

    # --------------------------------------------------------
    # Standardized statistics for variable features
    # --------------------------------------------------------

    standardized_records = []

    for feature in FEATURES:

        values = X[feature]
        std = values.std(ddof=1)

        if std == 0:
            standardized_records.append(
                {
                    "feature": feature,
                    "standardized_mean": np.nan,
                    "standardized_std": np.nan,
                    "status": "CONSTANT_FEATURE"
                }
            )
            continue

        z = (values - values.mean()) / std

        standardized_records.append(
            {
                "feature": feature,
                "standardized_mean": z.mean(),
                "standardized_std": z.std(ddof=1),
                "status": "VARIABLE_FEATURE"
            }
        )

    standardized_df = pd.DataFrame(
        standardized_records
    )

    # --------------------------------------------------------
    # Scale ratio
    # --------------------------------------------------------

    variable_scale_df = scale_df[
        scale_df["std"] > 0
    ].copy()

    if len(variable_scale_df) > 0:
        largest_std = variable_scale_df["std"].max()
        smallest_std = variable_scale_df["std"].min()
        scale_ratio = largest_std / smallest_std
    else:
        scale_ratio = np.nan

    # --------------------------------------------------------
    # Identify broad scale differences
    # --------------------------------------------------------

    high_scale_features = scale_df[
        scale_df["std"] > 10
    ]

    low_scale_features = scale_df[
        (scale_df["std"] > 0)
        & (scale_df["std"] < 1)
    ]

    constant_features = scale_df[
        scale_df["std"] == 0
    ]

    # --------------------------------------------------------
    # Build report
    # --------------------------------------------------------

    lines = []

    lines.append(
        "ML STEP 27 — FEATURE SCALE AND STANDARDIZATION AUDIT"
    )
    lines.append("=" * 70)
    lines.append("")

    lines.append(
        "Purpose: audit the numerical scale of the 14 biological ML "
        "features and document whether standardization is relevant "
        "for scale-sensitive models."
    )

    lines.append(
        "This is a preprocessing/design audit only. It does not "
        "establish predictive importance, biological significance, "
        "or clinical validity."
    )

    lines.append("")

    lines.append("INPUT")
    lines.append("-" * 70)
    lines.append(f"File: {INPUT_FILE}")
    lines.append(f"Rows: {len(df)}")
    lines.append(f"Features audited: {len(FEATURES)}")
    lines.append("")

    lines.append("RAW FEATURE SCALE")
    lines.append("-" * 70)

    for _, row in scale_df.iterrows():

        lines.append(
            f"{row['feature']}: "
            f"mean={row['mean']:.6f}, "
            f"SD={row['std']:.6f}, "
            f"min={row['minimum']:.6f}, "
            f"max={row['maximum']:.6f}, "
            f"range={row['range']:.6f}"
        )

    lines.append("")

    lines.append("STANDARDIZATION CHECK")
    lines.append("-" * 70)

    for _, row in standardized_df.iterrows():

        if row["status"] == "CONSTANT_FEATURE":
            lines.append(
                f"{row['feature']}: CONSTANT_FEATURE "
                "(standardized statistics not defined)"
            )
        else:
            lines.append(
                f"{row['feature']}: "
                f"standardized_mean={row['standardized_mean']:.12f}, "
                f"standardized_SD={row['standardized_std']:.12f}"
            )

    lines.append("")

    lines.append("FEATURE SCALE SPREAD")
    lines.append("-" * 70)

    if np.isnan(scale_ratio):
        lines.append(
            "Scale ratio could not be calculated because no "
            "variable features were present."
        )
    else:
        lines.append(
            f"Largest variable-feature SD: "
            f"{largest_std:.6f}"
        )
        lines.append(
            f"Smallest variable-feature SD: "
            f"{smallest_std:.6f}"
        )
        lines.append(
            f"Largest/smallest SD ratio: "
            f"{scale_ratio:.6f}"
        )

    lines.append("")

    lines.append("CONSTANT FEATURES")
    lines.append("-" * 70)

    if len(constant_features) == 0:
        lines.append("None")
    else:
        for _, row in constant_features.iterrows():
            lines.append(
                f"- {row['feature']}"
            )

    lines.append("")

    lines.append("FEATURES WITH SD > 10")
    lines.append("-" * 70)

    if len(high_scale_features) == 0:
        lines.append("None")
    else:
        for _, row in high_scale_features.iterrows():
            lines.append(
                f"- {row['feature']}: SD={row['std']:.6f}"
            )

    lines.append("")

    lines.append("VARIABLE FEATURES WITH SD < 1")
    lines.append("-" * 70)

    if len(low_scale_features) == 0:
        lines.append("None")
    else:
        for _, row in low_scale_features.iterrows():
            lines.append(
                f"- {row['feature']}: SD={row['std']:.6f}"
            )

    lines.append("")

    lines.append("INTERPRETATION")
    lines.append("-" * 70)

    lines.append(
        "The biological features are measured on different numerical "
        "scales, including amino-acid physicochemical values and "
        "protein-position/volume variables."
    )

    lines.append(
        "Scale differences are relevant for scale-sensitive models "
        "such as logistic regression."
    )

    lines.append(
        "The existing logistic-regression pipeline already applies "
        "StandardScaler inside each cross-validation training fold, "
        "which prevents validation-data leakage during scaling."
    )

    lines.append(
        "Random Forest does not require feature standardization for "
        "tree split construction."
    )

    if len(constant_features) > 0:
        lines.append(
            "The constant feature(s) contain no within-dataset "
            "variation and do not contribute useful scale information."
        )

    lines.append(
        "No feature was removed or transformed in this step."
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