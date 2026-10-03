from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# PROJECT DEEPGENE
# STEP 05 — MODEL INPUT MATRIX AUDIT
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

TRAIN_FILE = PROJECT_ROOT / "data" / "processed" / "deepgene_ml_train_v1.csv"
TEST_FILE = PROJECT_ROOT / "data" / "processed" / "deepgene_ml_test_v1.csv"

RESULTS_DIR = PROJECT_ROOT / "data" / "ml_results"
REPORT_FILE = RESULTS_DIR / "05_model_training.txt"


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
    print("STEP 05 — MODEL INPUT MATRIX AUDIT")
    print("=" * 70)

    if not TRAIN_FILE.exists():
        raise FileNotFoundError(f"Training dataset not found:\n{TRAIN_FILE}")

    if not TEST_FILE.exists():
        raise FileNotFoundError(f"Testing dataset not found:\n{TEST_FILE}")

    train = pd.read_csv(TRAIN_FILE)
    test = pd.read_csv(TEST_FILE)

    # --------------------------------------------------------
    # Validate feature columns
    # --------------------------------------------------------
    for name, df in [("training", train), ("testing", test)]:

        missing_features = [
            feature for feature in FEATURES
            if feature not in df.columns
        ]

        if missing_features:
            raise ValueError(
                f"Missing features in {name} dataset:\n"
                + "\n".join(f"  - {x}" for x in missing_features)
            )

    X_train = train[FEATURES].copy()
    X_test = test[FEATURES].copy()

    # --------------------------------------------------------
    # Check missing values
    # --------------------------------------------------------
    train_missing = int(X_train.isna().sum().sum())
    test_missing = int(X_test.isna().sum().sum())

    if train_missing != 0 or test_missing != 0:
        raise ValueError("Missing values detected in model features.")

    # --------------------------------------------------------
    # Check numeric compatibility
    # --------------------------------------------------------
    non_numeric_train = [
        col for col in FEATURES
        if not pd.api.types.is_numeric_dtype(X_train[col])
    ]

    non_numeric_test = [
        col for col in FEATURES
        if not pd.api.types.is_numeric_dtype(X_test[col])
    ]

    if non_numeric_train or non_numeric_test:
        raise ValueError(
            "Non-numeric feature detected.\n"
            f"Training: {non_numeric_train}\n"
            f"Testing: {non_numeric_test}"
        )

    # --------------------------------------------------------
    # Check finite values
    # --------------------------------------------------------
    train_finite = np.isfinite(X_train.to_numpy()).all()
    test_finite = np.isfinite(X_test.to_numpy()).all()

    if not train_finite or not test_finite:
        raise ValueError("Infinite or invalid numerical values detected.")

    # --------------------------------------------------------
    # Check feature ranges
    # --------------------------------------------------------
    binary_features = [
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

    invalid_binary = {}

    for feature in binary_features:
        values = sorted(X_train[feature].dropna().unique().tolist())

        if not set(values).issubset({0, 1}):
            invalid_binary[feature] = values

    if invalid_binary:
        raise ValueError(
            f"Invalid binary feature values: {invalid_binary}"
        )

    # --------------------------------------------------------
    # Basic matrix statistics
    # --------------------------------------------------------
    train_shape = X_train.shape
    test_shape = X_test.shape

    train_min = float(X_train.min().min())
    train_max = float(X_train.max().max())

    # --------------------------------------------------------
    # Create report
    # --------------------------------------------------------
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    report = [
        "PROJECT DEEPGENE",
        "STEP 05 — MODEL INPUT MATRIX AUDIT",
        "",
        f"Training matrix shape: {train_shape}",
        f"Testing matrix shape: {test_shape}",
        f"Feature count: {len(FEATURES)}",
        "",
        f"Training missing values: {train_missing}",
        f"Testing missing values: {test_missing}",
        "",
        f"Training matrix finite: {train_finite}",
        f"Testing matrix finite: {test_finite}",
        "",
        f"Training minimum value: {train_min}",
        f"Training maximum value: {train_max}",
        "",
        "FEATURES:",
    ]

    report.extend(f"  - {feature}" for feature in FEATURES)

    report.extend([
        "",
        "TARGET STATUS:",
        "  Independent target: NOT AVAILABLE",
        "  Supervised model training: NOT PERFORMED",
        "",
        "INTERPRETATION:",
        "  The feature matrix is numerically ready for supervised ML.",
        "  Model fitting remains blocked until an independent target",
        "  is supplied and validated.",
        "",
        "STATUS: PASS",
    ])

    REPORT_FILE.write_text(
        "\n".join(report),
        encoding="utf-8"
    )

    print(f"Training matrix: {train_shape}")
    print(f"Testing matrix: {test_shape}")
    print(f"Feature count: {len(FEATURES)}")
    print(f"Training missing values: {train_missing}")
    print(f"Testing missing values: {test_missing}")
    print(f"Training matrix finite: {train_finite}")
    print(f"Testing matrix finite: {test_finite}")
    print()
    print("Independent target: NOT AVAILABLE")
    print("Supervised model fitting: NOT PERFORMED")
    print()
    print("Report:")
    print(REPORT_FILE)
    print()
    print("STATUS: PASS")
    print("=" * 70)


if __name__ == "__main__":
    main()