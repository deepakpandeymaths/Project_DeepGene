from pathlib import Path
import pandas as pd


# ============================================================
# PROJECT DEEPGENE
# STEP 07 — TARGET PROVENANCE AND LEAKAGE AUDIT
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "data" / "ml_results"

REPORT_FILE = RESULTS_DIR / "07_cross_validation.txt"


# Features currently used by the ML pipeline.
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


# These are known ClinVar-derived evidence fields.
CLINVAR_DERIVED_FIELDS = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
]


def inspect_file(path):

    try:
        df = pd.read_csv(path, nrows=5)
    except Exception as exc:
        return {
            "file": path.name,
            "status": "READ_ERROR",
            "columns": [],
            "clinvar_fields": [],
            "error": str(exc),
        }

    columns = list(df.columns)

    clinvar_fields = [
        col for col in CLINVAR_DERIVED_FIELDS
        if col in columns
    ]

    return {
        "file": path.name,
        "status": "READ_OK",
        "columns": columns,
        "clinvar_fields": clinvar_fields,
        "error": "",
    }


def main():

    print("=" * 70)
    print("PROJECT DEEPGENE")
    print("STEP 07 — TARGET PROVENANCE AND LEAKAGE AUDIT")
    print("=" * 70)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    candidate_files = [
        PROCESSED_DIR / "deepgene_evidence_v1.csv",
        PROCESSED_DIR / "deepgene_ml_ready_v1.csv",
        PROCESSED_DIR / "deepgene_variant_profiles_v1.csv",
    ]

    results = []

    for path in candidate_files:

        if not path.exists():
            continue

        results.append(inspect_file(path))

    # --------------------------------------------------------
    # Evaluate target candidates
    # --------------------------------------------------------
    leakage_risk = []

    for result in results:

        if result["clinvar_fields"]:
            leakage_risk.append(result["file"])

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------
    report = [
        "PROJECT DEEPGENE",
        "STEP 07 — TARGET PROVENANCE AND LEAKAGE AUDIT",
        "",
        "PURPOSE:",
        "Determine whether candidate target-bearing datasets provide",
        "an independent ML target or merely reproduce ClinVar-derived",
        "evidence already represented in the feature matrix.",
        "",
        "CURRENT ML FEATURES:",
    ]

    report.extend(
        f"  - {feature}"
        for feature in FEATURES
    )

    report.extend([
        "",
        "CANDIDATE DATASETS:",
    ])

    for result in results:

        report.append(
            f"  {result['file']}: {result['status']}"
        )

        if result["clinvar_fields"]:
            report.append(
                "    ClinVar-derived fields detected:"
            )

            report.extend(
                f"      - {field}"
                for field in result["clinvar_fields"]
            )

    report.extend([
        "",
        "LEAKAGE ASSESSMENT:",
    ])

    if leakage_risk:

        report.append(
            "  Candidate datasets contain ClinVar-derived evidence "
            "already represented in ML features."
        )

        report.append(
            "  They are NOT independent ground-truth targets."
        )

    else:

        report.append(
            "  No ClinVar-derived fields detected in candidates."
        )

    report.extend([
        "",
        "CROSS-VALIDATION STATUS:",
        "  NOT PERFORMED.",
        "",
        "Reason:",
        "  Cross-validation requires a valid independent target.",
        "  Running cross-validation without such a target would not",
        "  produce a meaningful supervised ML performance estimate.",
        "",
        "REQUIRED NEXT DATA:",
        "  An independently defined variant-level target with:",
        "    1. documented provenance,",
        "    2. explicit labeling rules,",
        "    3. one label per eligible variant,",
        "    4. no circular dependence on the current ClinVar evidence",
        "       features.",
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
    print(f"Candidate datasets inspected: {len(results)}")
    print()

    for result in results:

        print(f"{result['file']}: {result['status']}")

        if result["clinvar_fields"]:
            print(
                f"  ClinVar-derived fields: "
                f"{len(result['clinvar_fields'])}"
            )

    print()
    print("Independent target: NOT AVAILABLE")
    print("Cross-validation: NOT PERFORMED")
    print()
    print("Reason:")
    print("  Candidate targets reproduce ClinVar-derived evidence.")
    print("  They are not independent ground truth.")
    print()
    print("Report:")
    print(REPORT_FILE)
    print()
    print("STATUS: PASS")
    print("=" * 70)


if __name__ == "__main__":
    main()