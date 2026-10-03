from pathlib import Path
import pandas as pd


# ============================================================
# PROJECT DEEPGENE
# STEP 04 — BASELINE MODEL READINESS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

TRAIN_FILE = PROJECT_ROOT / "data" / "processed" / "deepgene_ml_train_v1.csv"
TEST_FILE = PROJECT_ROOT / "data" / "processed" / "deepgene_ml_test_v1.csv"

RESULTS_DIR = PROJECT_ROOT / "data" / "ml_results"
REPORT_FILE = RESULTS_DIR / "04_baseline_models.txt"


def main():

    print("=" * 70)
    print("PROJECT DEEPGENE")
    print("STEP 04 — BASELINE MODEL READINESS")
    print("=" * 70)

    # --------------------------------------------------------
    # Check input files
    # --------------------------------------------------------
    if not TRAIN_FILE.exists():
        raise FileNotFoundError(f"Training file not found:\n{TRAIN_FILE}")

    if not TEST_FILE.exists():
        raise FileNotFoundError(f"Testing file not found:\n{TEST_FILE}")

    train = pd.read_csv(TRAIN_FILE)
    test = pd.read_csv(TEST_FILE)

    print(f"Training rows: {len(train)}")
    print(f"Testing rows: {len(test)}")
    print(f"Training columns: {len(train.columns)}")
    print(f"Testing columns: {len(test.columns)}")

    # --------------------------------------------------------
    # Expected feature structure
    # --------------------------------------------------------
    expected_features = [
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

    missing_train_features = [
        col for col in expected_features
        if col not in train.columns
    ]

    missing_test_features = [
        col for col in expected_features
        if col not in test.columns
    ]

    if missing_train_features:
        raise ValueError(
            "Missing expected training features:\n"
            + "\n".join(f"  - {x}" for x in missing_train_features)
        )

    if missing_test_features:
        raise ValueError(
            "Missing expected testing features:\n"
            + "\n".join(f"  - {x}" for x in missing_test_features)
        )

    # --------------------------------------------------------
    # Check that no clinical significance target is present
    # --------------------------------------------------------
    forbidden_target_fields = [
        "clinical_significance",
        "review_status",
    ]

    leakage_columns = [
        col for col in forbidden_target_fields
        if col in train.columns or col in test.columns
    ]

    # --------------------------------------------------------
    # Check missing values
    # --------------------------------------------------------
    train_missing = int(train[expected_features].isna().sum().sum())
    test_missing = int(test[expected_features].isna().sum().sum())

    # --------------------------------------------------------
    # Determine whether a valid independent target exists
    # --------------------------------------------------------
    independent_target_available = False

    # The current dataset does not contain an independent target.
    # We explicitly do not derive one from ClinVar clinical significance.
    target_status = "NOT AVAILABLE"

    # --------------------------------------------------------
    # Build report
    # --------------------------------------------------------
    report_lines = [
        "PROJECT DEEPGENE",
        "STEP 04 — BASELINE MODEL READINESS",
        "",
        f"Training rows: {len(train)}",
        f"Testing rows: {len(test)}",
        f"Feature count: {len(expected_features)}",
        "",
        "EXPECTED FEATURES:",
    ]

    report_lines.extend(
        f"  - {feature}" for feature in expected_features
    )

    report_lines.extend([
        "",
        f"Training feature missing values: {train_missing}",
        f"Testing feature missing values: {test_missing}",
        "",
        "INDEPENDENT TARGET:",
        f"  Status: {target_status}",
        "",
        "CLINVAR CLINICAL SIGNIFICANCE:",
        "  NOT used as an independent ML target.",
        "  Using it as the prediction target would create target leakage",
        "  because the current feature set contains ClinVar-derived evidence.",
        "",
        "LEAKAGE-CHECK FIELDS PRESENT:",
    ])

    if leakage_columns:
        report_lines.extend(
            f"  - {col}" for col in leakage_columns
        )
    else:
        report_lines.append("  None")

    report_lines.extend([
        "",
        "BASELINE MODEL STATUS:",
        "  SUPERVISED MODEL TRAINING BLOCKED",
        "",
        "Reason:",
        "  No independent ground-truth target is currently available.",
        "",
        "VALID NEXT INPUT:",
        "  An independently established target label for each variant,",
        "  with documented provenance and labeling rules.",
        "",
        "STATUS: PASS",
    ])

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    REPORT_FILE.write_text(
        "\n".join(report_lines),
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # Terminal output
    # --------------------------------------------------------
    print()
    print("Independent ML target: NOT AVAILABLE")
    print("ClinVar clinical_significance: NOT used as target")
    print(f"Training feature missing values: {train_missing}")
    print(f"Testing feature missing values: {test_missing}")
    print()
    print("Supervised baseline training: BLOCKED")
    print("Reason: independent ground-truth target is not available.")
    print()
    print("Report:")
    print(REPORT_FILE)
    print()
    print("STATUS: PASS")
    print("=" * 70)


if __name__ == "__main__":
    main()