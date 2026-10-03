from pathlib import Path
import pandas as pd


# ============================================================
# PROJECT DEEPGENE
# STEP 06 — INDEPENDENT TARGET AUDIT
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "data" / "ml_results"

REPORT_FILE = RESULTS_DIR / "06_model_evaluation.txt"


def main():

    print("=" * 70)
    print("PROJECT DEEPGENE")
    print("STEP 06 — INDEPENDENT TARGET AUDIT")
    print("=" * 70)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(PROCESSED_DIR.glob("*.csv"))

    evidence_dataset = (
        PROCESSED_DIR / "deepgene_v1_final.csv"
    )

    feature_dataset = (
        PROCESSED_DIR / "deepgene_ml_features_v1.csv"
    )

    train_dataset = (
        PROCESSED_DIR / "deepgene_ml_train_v1.csv"
    )

    test_dataset = (
        PROCESSED_DIR / "deepgene_ml_test_v1.csv"
    )

    known_ml_files = {
        evidence_dataset.name,
        feature_dataset.name,
        train_dataset.name,
        test_dataset.name,
    }

    candidate_external_files = [
        path for path in csv_files
        if path.name not in known_ml_files
    ]

    print(f"Processed CSV files found: {len(csv_files)}")
    print()

    print("Candidate external target files:")

    if candidate_external_files:
        for path in candidate_external_files:
            print(f"  - {path.name}")
    else:
        print("  None")

    # --------------------------------------------------------
    # Inspect candidate files for possible target columns
    # --------------------------------------------------------
    target_keywords = [
        "target",
        "label",
        "class",
        "outcome",
        "pathogenicity",
        "pathogenic",
        "benign",
        "ground_truth",
        "groundtruth",
    ]

    candidates_with_target = []

    for path in candidate_external_files:

        try:
            df = pd.read_csv(path, nrows=5)
        except Exception:
            continue

        matching_columns = []

        for column in df.columns:
            column_lower = str(column).lower()

            if any(
                keyword in column_lower
                for keyword in target_keywords
            ):
                matching_columns.append(column)

        if matching_columns:
            candidates_with_target.append(
                (path.name, matching_columns)
            )

    # --------------------------------------------------------
    # Important: ClinVar evidence is not an independent target
    # --------------------------------------------------------
    target_available = len(candidates_with_target) > 0

    report_lines = [
        "PROJECT DEEPGENE",
        "STEP 06 — INDEPENDENT TARGET AUDIT",
        "",
        f"Processed CSV files found: {len(csv_files)}",
        "",
        "Known ML/evidence datasets excluded from target search:",
    ]

    report_lines.extend(
        f"  - {name}" for name in sorted(known_ml_files)
    )

    report_lines.extend([
        "",
        "Candidate external files:",
    ])

    if candidate_external_files:
        report_lines.extend(
            f"  - {path.name}"
            for path in candidate_external_files
        )
    else:
        report_lines.append("  None")

    report_lines.extend([
        "",
        "Potential target-bearing files:",
    ])

    if candidates_with_target:
        for filename, columns in candidates_with_target:
            report_lines.append(
                f"  - {filename}: {', '.join(map(str, columns))}"
            )
    else:
        report_lines.append("  None identified")

    report_lines.extend([
        "",
        "TARGET DECISION:",
    ])

    if target_available:
        report_lines.extend([
            "  Potential target source detected.",
            "  It requires manual provenance and leakage validation",
            "  before it can be used as an ML target.",
            "",
            "  Supervised modeling remains BLOCKED until validated.",
        ])
    else:
        report_lines.extend([
            "  No independent target identified.",
            "",
            "  Supervised modeling remains BLOCKED.",
            "",
            "  ClinVar clinical_significance must not be converted",
            "  into an ML target because the current feature set",
            "  contains ClinVar-derived evidence.",
        ])

    report_lines.extend([
        "",
        "STATUS: PASS",
    ])

    REPORT_FILE.write_text(
        "\n".join(report_lines),
        encoding="utf-8"
    )

    print()

    if candidate_external_files:
        print(
            f"Candidate external CSV files: "
            f"{len(candidate_external_files)}"
        )
    else:
        print("Candidate external target files: 0")

    if candidates_with_target:
        print(
            f"Potential target-bearing files: "
            f"{len(candidates_with_target)}"
        )

        for filename, columns in candidates_with_target:
            print(
                f"  {filename}: "
                f"{', '.join(map(str, columns))}"
            )
    else:
        print("Potential target-bearing files: 0")

    print()
    print("Independent target: NOT YET VALIDATED")
    print("Supervised modeling: BLOCKED")
    print()
    print("Report:")
    print(REPORT_FILE)
    print()
    print("STATUS: PASS")
    print("=" * 70)


if __name__ == "__main__":
    main()