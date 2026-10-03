from pathlib import Path
import pandas as pd


# ============================================================
# PROJECT DEEPGENE
# STEP 09 — FEATURE STABILITY AND DISTRIBUTION AUDIT
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

SUMMARY_FILE = (
    RESULTS_DIR / "09_feature_distribution_v1.csv"
)

REPORT_FILE = (
    RESULTS_DIR / "09_model_interpretability.txt"
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
    print("STEP 09 — FEATURE STABILITY AND DISTRIBUTION AUDIT")
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

    X_train = train[FEATURES].copy()
    X_test = test[FEATURES].copy()

    # --------------------------------------------------------
    # Build feature summary
    # --------------------------------------------------------
    rows = []

    for feature in FEATURES:

        train_values = X_train[feature]
        test_values = X_test[feature]

        train_unique = int(train_values.nunique())
        test_unique = int(test_values.nunique())

        train_mean = float(train_values.mean())
        test_mean = float(test_values.mean())

        train_std = float(train_values.std())
        test_std = float(test_values.std())

        train_min = float(train_values.min())
        test_min = float(test_values.min())

        train_max = float(train_values.max())
        test_max = float(test_values.max())

        rows.append({
            "feature": feature,
            "train_unique_values": train_unique,
            "test_unique_values": test_unique,
            "train_mean": train_mean,
            "test_mean": test_mean,
            "train_std": train_std,
            "test_std": test_std,
            "train_min": train_min,
            "test_min": test_min,
            "train_max": train_max,
            "test_max": test_max,
            "train_constant": train_unique <= 1,
            "test_constant": test_unique <= 1,
        })

    summary = pd.DataFrame(rows)

    # --------------------------------------------------------
    # Identify constant features
    # --------------------------------------------------------
    train_constant = summary[
        summary["train_constant"]
    ]["feature"].tolist()

    test_constant = summary[
        summary["test_constant"]
    ]["feature"].tolist()

    # --------------------------------------------------------
    # Identify train/test distribution shifts
    # --------------------------------------------------------
    distribution_notes = []

    for _, row in summary.iterrows():

        feature = row["feature"]

        train_mean = row["train_mean"]
        test_mean = row["test_mean"]

        if train_mean == 0 and test_mean == 0:
            relative_difference = 0.0

        elif train_mean != 0:
            relative_difference = abs(
                test_mean - train_mean
            ) / abs(train_mean)

        else:
            relative_difference = float("inf")

        distribution_notes.append({
            "feature": feature,
            "absolute_mean_difference": abs(
                test_mean - train_mean
            ),
            "relative_mean_difference": relative_difference,
        })

    distribution_df = pd.DataFrame(distribution_notes)

    summary = summary.merge(
        distribution_df,
        on="feature",
        how="left"
    )

    # --------------------------------------------------------
    # Save summary
    # --------------------------------------------------------
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    summary.to_csv(
        SUMMARY_FILE,
        index=False
    )

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------
    report = [
        "PROJECT DEEPGENE",
        "STEP 09 — FEATURE STABILITY AND DISTRIBUTION AUDIT",
        "",
        f"Training rows: {len(train)}",
        f"Testing rows: {len(test)}",
        f"Feature count: {len(FEATURES)}",
        "",
        "TRAINING-SET CONSTANT FEATURES:",
    ]

    if train_constant:
        report.extend(
            f"  - {feature}"
            for feature in train_constant
        )
    else:
        report.append("  None")

    report.extend([
        "",
        "TEST-SET CONSTANT FEATURES:",
    ])

    if test_constant:
        report.extend(
            f"  - {feature}"
            for feature in test_constant
        )
    else:
        report.append("  None")

    report.extend([
        "",
        "FEATURE DISTRIBUTION SUMMARY:",
    ])

    for _, row in summary.iterrows():

        report.append(
            f"  {row['feature']}: "
            f"train_unique={int(row['train_unique_values'])}, "
            f"test_unique={int(row['test_unique_values'])}, "
            f"train_mean={row['train_mean']:.6f}, "
            f"test_mean={row['test_mean']:.6f}"
        )

    report.extend([
        "",
        "TARGET STATUS:",
        "  Independent target: NOT AVAILABLE",
        "  Model interpretability analysis: NOT PERFORMED",
        "",
        "NOTE:",
        "  This step evaluates feature stability only.",
        "  It does not estimate model performance or feature importance.",
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
    print(f"Training rows: {len(train)}")
    print(f"Testing rows: {len(test)}")
    print(f"Feature count: {len(FEATURES)}")
    print()

    print("Training constant features:")

    if train_constant:
        for feature in train_constant:
            print(f"  - {feature}")
    else:
        print("  None")

    print()
    print("Testing constant features:")

    if test_constant:
        for feature in test_constant:
            print(f"  - {feature}")
    else:
        print("  None")

    print()
    print("Independent target: NOT AVAILABLE")
    print("Model interpretability: NOT PERFORMED")
    print()
    print("Distribution summary:")
    print(SUMMARY_FILE)
    print()
    print("Report:")
    print(REPORT_FILE)
    print()
    print("STATUS: PASS")
    print("=" * 70)


if __name__ == "__main__":
    main()