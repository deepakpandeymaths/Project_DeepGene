from pathlib import Path
import pandas as pd


# ============================================================
# PROJECT DEEPGENE
# STEP 10 — FINAL FEATURE-SET AUDIT
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

TRAIN_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_train_v1.csv"
)

TEST_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_test_v1.csv"
)

RESULTS_DIR = PROJECT_ROOT / "data" / "ml_results"

FINAL_FEATURE_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_features_v2.csv"
)

MANIFEST_FILE = (
    RESULTS_DIR
    / "10_final_feature_manifest_v2.csv"
)

REPORT_FILE = (
    RESULTS_DIR
    / "10_error_analysis.txt"
)


# Current 18-feature candidate set.
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


# Constant features identified in Step 09.
CONSTANT_FEATURES = [
    "phenotype_available",
    "genomic_coordinates_present",
    "hgvs_present",
]


# Keep conflicting_count as the primary quantitative feature.
# conflict_flag is a deterministic duplicate because r = 1.0.
REDUNDANT_FEATURE = "conflict_flag"


def main():

    print("=" * 70)
    print("PROJECT DEEPGENE")
    print("STEP 10 — FINAL FEATURE-SET AUDIT")
    print("=" * 70)

    if not TRAIN_FILE.exists():
        raise FileNotFoundError(
            f"Training dataset not found:\n{TRAIN_FILE}"
        )

    if not TEST_FILE.exists():
        raise FileNotFoundError(
            f"Testing dataset not found:\n{TEST_FILE}"
        )

    train = pd.read_csv(TRAIN_FILE)
    test = pd.read_csv(TEST_FILE)

    # --------------------------------------------------------
    # Validate current feature set
    # --------------------------------------------------------
    missing_train = [
        feature
        for feature in FEATURES
        if feature not in train.columns
    ]

    missing_test = [
        feature
        for feature in FEATURES
        if feature not in test.columns
    ]

    if missing_train or missing_test:
        raise ValueError(
            f"Missing training features: {missing_train}\n"
            f"Missing testing features: {missing_test}"
        )

    # --------------------------------------------------------
    # Construct final candidate feature set
    # --------------------------------------------------------
    final_features = [
        feature
        for feature in FEATURES
        if feature not in CONSTANT_FEATURES
        and feature != REDUNDANT_FEATURE
    ]

    # --------------------------------------------------------
    # Verify no duplicates
    # --------------------------------------------------------
    if len(final_features) != len(set(final_features)):
        raise ValueError(
            "Duplicate feature names detected."
        )

    # --------------------------------------------------------
    # Create final feature dataset from full source
    # --------------------------------------------------------
    source = pd.read_csv(
        PROJECT_ROOT
        / "data"
        / "processed"
        / "deepgene_v1_final.csv"
    )

    final_dataset = source[
        [
            "variant_id",
            *final_features
        ]
    ].copy()

    # --------------------------------------------------------
    # Validate final matrix
    # --------------------------------------------------------
    missing_values = int(
        final_dataset[final_features]
        .isna()
        .sum()
        .sum()
    )

    if missing_values != 0:
        raise ValueError(
            "Missing values detected in final feature dataset."
        )

    # --------------------------------------------------------
    # Manifest
    # --------------------------------------------------------
    manifest_rows = []

    for feature in FEATURES:

        if feature in CONSTANT_FEATURES:

            decision = "REMOVED_CONSTANT"

            reason = (
                "Constant in both training and testing datasets."
            )

        elif feature == REDUNDANT_FEATURE:

            decision = "REMOVED_REDUNDANT"

            reason = (
                "Perfectly redundant with conflicting_count "
                "(Pearson r = 1.0000)."
            )

        else:

            decision = "RETAIN"

            reason = (
                "Retained as a non-constant candidate feature."
            )

        manifest_rows.append({
            "feature": feature,
            "decision": decision,
            "reason": reason,
        })

    manifest = pd.DataFrame(manifest_rows)

    # --------------------------------------------------------
    # Save outputs
    # --------------------------------------------------------
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    final_dataset.to_csv(
        FINAL_FEATURE_FILE,
        index=False
    )

    manifest.to_csv(
        MANIFEST_FILE,
        index=False
    )

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------
    report = [
        "PROJECT DEEPGENE",
        "STEP 10 — FINAL FEATURE-SET AUDIT",
        "",
        f"Original candidate features: {len(FEATURES)}",
        f"Constant features removed: {len(CONSTANT_FEATURES)}",
        "Perfectly redundant feature removed: 1",
        f"Final candidate features: {len(final_features)}",
        "",
        "CONSTANT FEATURES REMOVED:",
    ]

    report.extend(
        f"  - {feature}"
        for feature in CONSTANT_FEATURES
    )

    report.extend([
        "",
        "REDUNDANT FEATURE REMOVED:",
        f"  - {REDUNDANT_FEATURE}",
        "    Reason: perfectly correlated with conflicting_count.",
        "",
        "RETAINED FEATURES:",
    ])

    report.extend(
        f"  - {feature}"
        for feature in final_features
    )

    report.extend([
        "",
        f"Final dataset rows: {len(final_dataset)}",
        f"Final dataset columns: {len(final_dataset.columns)}",
        f"Final feature missing values: {missing_values}",
        "",
        "TARGET STATUS:",
        "  Independent target: NOT AVAILABLE",
        "  Error analysis: NOT PERFORMED",
        "",
        "IMPORTANT:",
        "  Feature removal is based only on feature structure.",
        "  No predictive performance claim is made.",
        "  Supervised modeling remains blocked until an independent",
        "  target is established.",
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
    print(
        f"Original candidate features: "
        f"{len(FEATURES)}"
    )

    print(
        f"Constant features removed: "
        f"{len(CONSTANT_FEATURES)}"
    )

    print("Redundant feature removed: 1")
    print(
        f"Final candidate features: "
        f"{len(final_features)}"
    )

    print()
    print("Removed constant features:")

    for feature in CONSTANT_FEATURES:
        print(f"  - {feature}")

    print()
    print(
        "Removed redundant feature: "
        f"{REDUNDANT_FEATURE}"
    )

    print()
    print(
        "Final feature dataset:"
    )
    print(FINAL_FEATURE_FILE)

    print()
    print("Feature manifest:")
    print(MANIFEST_FILE)

    print()
    print("Independent target: NOT AVAILABLE")
    print("Supervised modeling: BLOCKED")
    print()
    print("Report:")
    print(REPORT_FILE)
    print()
    print("STATUS: PASS")
    print("=" * 70)


if __name__ == "__main__":
    main()