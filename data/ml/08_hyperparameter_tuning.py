from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# PROJECT DEEPGENE
# STEP 08 — FEATURE REDUNDANCY AUDIT
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

TRAIN_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_train_v1.csv"
)

RESULTS_DIR = PROJECT_ROOT / "data" / "ml_results"

CORRELATION_FILE = (
    RESULTS_DIR / "08_feature_correlations_v1.csv"
)

REPORT_FILE = (
    RESULTS_DIR / "08_hyperparameter_tuning.txt"
)


FEATURES = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
    "phenotype_available",
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",
]


def main():

    print("=" * 70)
    print("PROJECT DEEPGENE")
    print("STEP 08 — FEATURE REDUNDANCY AUDIT")
    print("=" * 70)

    if not TRAIN_FILE.exists():
        raise FileNotFoundError(
            f"Training dataset not found:\n{TRAIN_FILE}"
        )

    train = pd.read_csv(TRAIN_FILE)

    missing = [
        feature
        for feature in FEATURES
        if feature not in train.columns
    ]

    if missing:
        raise ValueError(
            "Missing expected features:\n"
            + "\n".join(f"  - {x}" for x in missing)
        )

    X = train[FEATURES].copy()

    if X.isna().sum().sum() != 0:
        raise ValueError(
            "Missing values detected in feature matrix."
        )

    # --------------------------------------------------------
    # Correlation matrix
    # --------------------------------------------------------
    correlation = X.corr(method="pearson")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    correlation.to_csv(CORRELATION_FILE)

    # --------------------------------------------------------
    # Extract unique feature pairs
    # --------------------------------------------------------
    pairs = []

    for i in range(len(FEATURES)):
        for j in range(i + 1, len(FEATURES)):

            feature_a = FEATURES[i]
            feature_b = FEATURES[j]

            value = correlation.loc[
                feature_a,
                feature_b
            ]

            if pd.notna(value):
                pairs.append(
                    {
                        "feature_a": feature_a,
                        "feature_b": feature_b,
                        "pearson_correlation": float(value),
                        "absolute_correlation": abs(float(value)),
                    }
                )

    pair_df = pd.DataFrame(pairs)

    pair_df = pair_df.sort_values(
        "absolute_correlation",
        ascending=False
    ).reset_index(drop=True)

    # --------------------------------------------------------
    # Highly correlated pairs
    # --------------------------------------------------------
    high_correlation = pair_df[
        pair_df["absolute_correlation"] >= 0.90
    ].copy()

    # --------------------------------------------------------
    # Constant features
    # --------------------------------------------------------
    constant_features = [
        feature
        for feature in FEATURES
        if X[feature].nunique(dropna=False) <= 1
    ]

    # --------------------------------------------------------
    # Feature summary
    # --------------------------------------------------------
    feature_summary = []

    for feature in FEATURES:

        feature_summary.append(
            {
                "feature": feature,
                "unique_values": int(X[feature].nunique()),
                "mean": float(X[feature].mean()),
                "std": float(X[feature].std()),
                "minimum": float(X[feature].min()),
                "maximum": float(X[feature].max()),
            }
        )

    summary_df = pd.DataFrame(feature_summary)

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------
    report = [
        "PROJECT DEEPGENE",
        "STEP 08 — FEATURE REDUNDANCY AUDIT",
        "",
        f"Training rows: {len(X)}",
        f"Feature count: {len(FEATURES)}",
        "",
        "PURPOSE:",
        "Assess numerical redundancy among candidate ML features",
        "before an independently labeled target is introduced.",
        "",
        f"Feature pairs evaluated: {len(pair_df)}",
        f"Pairs with |Pearson r| >= 0.90: {len(high_correlation)}",
        f"Constant features: {len(constant_features)}",
        "",
        "HIGH-CORRELATION PAIRS:",
    ]

    if len(high_correlation) > 0:

        for _, row in high_correlation.iterrows():

            report.append(
                f"  - {row['feature_a']} <-> "
                f"{row['feature_b']}: "
                f"r={row['pearson_correlation']:.4f}"
            )

    else:

        report.append(
            "  None at threshold |r| >= 0.90"
        )

    report.extend([
        "",
        "CONSTANT FEATURES:",
    ])

    if constant_features:

        report.extend(
            f"  - {feature}"
            for feature in constant_features
        )

    else:

        report.append("  None")

    report.extend([
        "",
        "TARGET STATUS:",
        "  Independent target: NOT AVAILABLE",
        "  Supervised model tuning: NOT PERFORMED",
        "",
        "INTERPRETATION:",
        "  Correlation analysis is an unsupervised feature audit.",
        "  It does not establish predictive performance.",
        "",
        "STATUS: PASS",
    ])

    REPORT_FILE.write_text(
        "\n".join(report),
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # Terminal output
    # --------------------------------------------------------
    print(f"Training rows: {len(X)}")
    print(f"Feature count: {len(FEATURES)}")
    print(f"Feature pairs evaluated: {len(pair_df)}")
    print(
        f"Pairs with |Pearson r| >= 0.90: "
        f"{len(high_correlation)}"
    )
    print(
        f"Constant features: "
        f"{len(constant_features)}"
    )
    print()

    if len(high_correlation) > 0:
        print("Highly correlated feature pairs:")

        for _, row in high_correlation.iterrows():
            print(
                f"  {row['feature_a']} <-> "
                f"{row['feature_b']}: "
                f"r={row['pearson_correlation']:.4f}"
            )

        print()

    print("Independent target: NOT AVAILABLE")
    print("Supervised model tuning: NOT PERFORMED")
    print()
    print("Correlation matrix:")
    print(CORRELATION_FILE)
    print()
    print("Report:")
    print(REPORT_FILE)
    print()
    print("STATUS: PASS")
    print("=" * 70)


if __name__ == "__main__":
    main()