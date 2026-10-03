from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# ML STEP 25 — FEATURE CORRELATION AUDIT
# ============================================================

INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
OUTPUT_FILE = Path("data/analysis_results/25_ml_feature_correlation_audit.txt")

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

CORRELATION_THRESHOLD = 0.80


def main():
    print("=" * 70)
    print("ML STEP 25 — FEATURE CORRELATION AUDIT")
    print("=" * 70)

    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE)

    missing_features = [f for f in FEATURES if f not in df.columns]

    if missing_features:
        raise ValueError(
            f"Missing expected ML features: {missing_features}"
        )

    X = df[FEATURES].copy()

    # --------------------------------------------------------
    # Basic validation
    # --------------------------------------------------------
    if X.isna().any().any():
        raise ValueError("Missing values detected in ML features.")

    non_numeric = [
        col for col in FEATURES
        if not pd.api.types.is_numeric_dtype(X[col])
    ]

    if non_numeric:
        raise ValueError(
            f"Non-numeric ML features detected: {non_numeric}"
        )

    # --------------------------------------------------------
    # Pearson correlation
    # --------------------------------------------------------
    corr = X.corr(method="pearson")

    # Collect unique feature pairs
    pairs = []

    for i in range(len(FEATURES)):
        for j in range(i + 1, len(FEATURES)):
            f1 = FEATURES[i]
            f2 = FEATURES[j]

            value = corr.loc[f1, f2]

            pairs.append(
                {
                    "feature_1": f1,
                    "feature_2": f2,
                    "correlation": value,
                    "absolute_correlation": abs(value),
                }
            )

    pairs_df = pd.DataFrame(pairs)

    pairs_df = pairs_df.sort_values(
        "absolute_correlation",
        ascending=False
    ).reset_index(drop=True)

    high_corr = pairs_df[
        pairs_df["absolute_correlation"] >= CORRELATION_THRESHOLD
    ].copy()

    # --------------------------------------------------------
    # Constant features
    # --------------------------------------------------------
    constant_features = [
        col for col in FEATURES
        if X[col].nunique(dropna=False) <= 1
    ]

    # --------------------------------------------------------
    # Output report
    # --------------------------------------------------------
    lines = []

    lines.append("ML STEP 25 — FEATURE CORRELATION AUDIT")
    lines.append("=" * 70)
    lines.append("")
    lines.append(
        "Purpose: audit redundancy and pairwise correlation among the "
        "14 biological ML features."
    )
    lines.append(
        "This is a feature-structure audit only. It does not establish "
        "predictive importance, causality, or biological significance."
    )
    lines.append("")

    lines.append("INPUT")
    lines.append("-" * 70)
    lines.append(f"File: {INPUT_FILE}")
    lines.append(f"Rows: {len(df)}")
    lines.append(f"Features audited: {len(FEATURES)}")
    lines.append("")

    lines.append("FEATURES")
    lines.append("-" * 70)

    for feature in FEATURES:
        lines.append(f"- {feature}")

    lines.append("")

    lines.append("CONSTANT FEATURES")
    lines.append("-" * 70)

    if constant_features:
        for feature in constant_features:
            lines.append(f"- {feature}")
    else:
        lines.append("None")

    lines.append("")

    lines.append("PAIRWISE PEARSON CORRELATIONS")
    lines.append("-" * 70)

    for _, row in pairs_df.iterrows():
        lines.append(
            f"{row['feature_1']} <-> {row['feature_2']}: "
            f"{row['correlation']:.6f}"
        )

    lines.append("")

    lines.append(
        f"HIGH-CORRELATION PAIRS "
        f"(|r| >= {CORRELATION_THRESHOLD:.2f})"
    )
    lines.append("-" * 70)

    if len(high_corr) == 0:
        lines.append("None")
    else:
        for _, row in high_corr.iterrows():
            lines.append(
                f"{row['feature_1']} <-> {row['feature_2']}: "
                f"r={row['correlation']:.6f}"
            )

    lines.append("")

    # --------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------
    lines.append("INTERPRETATION")
    lines.append("-" * 70)

    if constant_features:
        lines.append(
            "At least one feature is constant and therefore contains "
            "no within-dataset variation."
        )
    else:
        lines.append(
            "No constant ML features were detected."
        )

    if len(high_corr) == 0:
        lines.append(
            "No feature pair reached the predefined high-correlation "
            "threshold of |r| >= 0.80."
        )
        lines.append(
            "The 14-feature representation therefore does not show "
            "strong pairwise linear redundancy under this threshold."
        )
    else:
        lines.append(
            f"{len(high_corr)} feature pair(s) reached "
            f"|r| >= {CORRELATION_THRESHOLD:.2f}."
        )
        lines.append(
            "These pairs represent potentially redundant feature "
            "representations and should be considered when interpreting "
            "model behavior."
        )

    lines.append("")
    lines.append(
        "Important: correlation does not imply causation and does not "
        "demonstrate predictive usefulness."
    )
    lines.append(
        "No feature was removed or selected in this step."
    )
    lines.append(
        "No model was tuned or selected in this step."
    )
    lines.append(
        "No clinical interpretation is made."
    )

    lines.append("")
    lines.append("STATUS: PASS")
    lines.append("=" * 70)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8"
    )

    print("\n".join(lines))

    print("")
    print(f"Report saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()