# ============================================================
# 19. VALIDATE DEEPGENE EVIDENCE DATASET
# ============================================================

import os
import sys
import pandas as pd


# ============================================================
# PATHS
# ============================================================

BASE_DIR = os.path.dirname(
    os.path.dirname(
        os.path.abspath(__file__)
    )
)

PROCESSED_DIR = os.path.join(BASE_DIR, "processed")
RESULTS_DIR = os.path.join(BASE_DIR, "analysis_results")

INPUT_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_evidence_v1.csv"
)

PRIORITIZED_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_prioritized_v1.csv"
)

OUTPUT_FILE = os.path.join(
    RESULTS_DIR,
    "19_validate_evidence_dataset.txt"
)


# ============================================================
# REQUIRED COLUMNS
# ============================================================

REQUIRED_COLUMNS = [
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
    "hgvs_name",
    "clinical_significance",
    "review_status",
    "phenotypes",
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
    "hgvs_present"
]


# ============================================================
# HELPER
# ============================================================

results = []


def check(name, passed, detail):
    status = "PASS" if passed else "FAIL"

    results.append({
        "check": name,
        "status": status,
        "detail": detail
    })

    symbol = "✓" if passed else "✗"

    print(f"{symbol} {name}")
    print(f"  {detail}")


# ============================================================
# HEADER
# ============================================================

print("=" * 70)
print("DEEPGENE EVIDENCE DATASET VALIDATION")
print("=" * 70)

print()
print(f"Input: {INPUT_FILE}")


# ============================================================
# LOAD DATA
# ============================================================

if not os.path.exists(INPUT_FILE):
    print()
    print("ERROR: Input file not found.")
    print(INPUT_FILE)
    sys.exit(1)

df = pd.read_csv(INPUT_FILE)

print()
print(f"Rows loaded: {len(df)}")
print(f"Columns loaded: {len(df.columns)}")


# ============================================================
# 1. ROW COUNT
# ============================================================

print()
print("=" * 70)
print("1. DATASET SIZE")
print("=" * 70)

check(
    "Expected variant count",
    len(df) == 5381,
    f"Found {len(df)} rows; expected 5381."
)


# ============================================================
# 2. REQUIRED COLUMNS
# ============================================================

print()
print("=" * 70)
print("2. REQUIRED COLUMNS")
print("=" * 70)

missing_columns = [
    col for col in REQUIRED_COLUMNS
    if col not in df.columns
]

check(
    "Required columns present",
    len(missing_columns) == 0,
    (
        "All required columns are present."
        if not missing_columns
        else f"Missing columns: {missing_columns}"
    )
)


# ============================================================
# 3. UNIQUE VARIANT IDS
# ============================================================

print()
print("=" * 70)
print("3. VARIANT UNIQUENESS")
print("=" * 70)

duplicate_variant_ids = df[
    df["variant_id"].duplicated(keep=False)
]["variant_id"].unique()

check(
    "One row per variant",
    len(duplicate_variant_ids) == 0,
    (
        "No duplicate variant_id values found."
        if len(duplicate_variant_ids) == 0
        else (
            f"Found {len(duplicate_variant_ids)} duplicated "
            "variant_id values."
        )
    )
)

check(
    "Unique variant count",
    df["variant_id"].nunique() == len(df),
    (
        f"{df['variant_id'].nunique()} unique variant IDs "
        f"across {len(df)} rows."
    )
)


# ============================================================
# 4. ESSENTIAL IDENTIFIERS
# ============================================================

print()
print("=" * 70)
print("4. ESSENTIAL IDENTIFIERS")
print("=" * 70)

for column in [
    "variant_id",
    "hgvs_name",
    "gene_id"
]:

    missing = df[column].isna().sum()

    check(
        f"No missing {column}",
        missing == 0,
        f"{missing} missing values."
    )


# ============================================================
# 5. CLINICAL SIGNIFICANCE
# ============================================================

print()
print("=" * 70)
print("5. CLINICAL SIGNIFICANCE")
print("=" * 70)

clinical_counts = df["clinical_significance"].value_counts(
    dropna=False
)

print()
print(clinical_counts.to_string())

allowed_clinical_significance = {
    "Uncertain significance",
    "Pathogenic",
    "Likely benign",
    "Likely pathogenic",
    "Conflicting classifications of pathogenicity",
    "Pathogenic/Likely pathogenic",
    "Benign",
    "Benign/Likely benign",
    "not provided",
    "-",
    "drug response"
}

unexpected_clinical = set(
    clinical_counts.index
) - allowed_clinical_significance

check(
    "Clinical significance categories",
    len(unexpected_clinical) == 0,
    (
        "All categories are expected."
        if not unexpected_clinical
        else f"Unexpected categories: {unexpected_clinical}"
    )
)


# ============================================================
# 6. REVIEW STATUS
# ============================================================

print()
print("=" * 70)
print("6. REVIEW STATUS")
print("=" * 70)

review_counts = df["review_status"].value_counts(
    dropna=False
)

print()
print(review_counts.to_string())

allowed_review_status = {
    "criteria provided, single submitter",
    "criteria provided, multiple submitters, no conflicts",
    "criteria provided, conflicting classifications",
    "no assertion criteria provided",
    "no classification provided",
    "reviewed by expert panel",
    "-"
}

unexpected_review = set(
    review_counts.index
) - allowed_review_status

check(
    "Review status categories",
    len(unexpected_review) == 0,
    (
        "All categories are expected."
        if not unexpected_review
        else f"Unexpected categories: {unexpected_review}"
    )
)


# ============================================================
# 7. SUBMITTER COUNT
# ============================================================

print()
print("=" * 70)
print("7. SUBMITTER COUNTS")
print("=" * 70)

missing_submitters = df["submitter_count"].isna().sum()
negative_submitters = (
    df["submitter_count"].dropna() < 0
).sum()
zero_submitters = (
    df["submitter_count"].dropna() == 0
).sum()

check(
    "Submitter count has no missing values",
    missing_submitters == 0,
    f"{missing_submitters} missing values."
)

check(
    "Submitter count is non-negative",
    negative_submitters == 0,
    f"{negative_submitters} negative values."
)

check(
    "Submitter count is positive",
    zero_submitters == 0,
    f"{zero_submitters} variants have zero submitters."
)

print()
print("Submitter count summary:")
print(df["submitter_count"].describe().to_string())


# ============================================================
# 8. MULTIPLE SUBMITTER FLAG
# ============================================================

print()
print("=" * 70)
print("8. MULTIPLE-SUBMITTER FLAG")
print("=" * 70)

expected_multiple = (
    df["submitter_count"] > 1
).astype(int)

multiple_consistent = (
    df["multiple_submitters"] == expected_multiple
).all()

inconsistent_multiple = (
    df["multiple_submitters"] != expected_multiple
).sum()

check(
    "multiple_submitters flag is consistent",
    multiple_consistent,
    (
        "Flag matches submitter_count for all variants."
        if multiple_consistent
        else (
            f"{inconsistent_multiple} variants have "
            "inconsistent multiple_submitters flags."
        )
    )
)


# ============================================================
# 9. EXPERT PANEL FLAG
# ============================================================

print()
print("=" * 70)
print("9. EXPERT PANEL FLAG")
print("=" * 70)

expected_expert = (
    df["review_status"]
    == "reviewed by expert panel"
).astype(int)

expert_consistent = (
    df["expert_panel_review"] == expected_expert
).all()

inconsistent_expert = (
    df["expert_panel_review"] != expected_expert
).sum()

check(
    "expert_panel_review flag is consistent",
    expert_consistent,
    (
        "Flag matches review status for all variants."
        if expert_consistent
        else (
            f"{inconsistent_expert} variants have "
            "inconsistent expert-panel flags."
        )
    )
)


# ============================================================
# 10. CONFLICT FLAG
# ============================================================

print()
print("=" * 70)
print("10. CONFLICT FLAG")
print("=" * 70)

expected_conflict = (
    df["clinical_significance"]
    == "Conflicting classifications of pathogenicity"
).astype(int)

conflict_consistent = (
    df["conflict_flag"] == expected_conflict
).all()

inconsistent_conflict = (
    df["conflict_flag"] != expected_conflict
).sum()

check(
    "conflict_flag is consistent",
    conflict_consistent,
    (
        "Flag matches clinical significance for all variants."
        if conflict_consistent
        else (
            f"{inconsistent_conflict} variants have "
            "inconsistent conflict flags."
        )
    )
)


# ============================================================
# 11. EXPERT PANEL / SUBMITTER RELATIONSHIP
# ============================================================

print()
print("=" * 70)
print("11. EXPERT PANEL CONSISTENCY")
print("=" * 70)

expert_variants = df[
    df["expert_panel_review"] == 1
]

check(
    "Expert-panel variants have expert review status",
    (
        expert_variants["review_status"]
        == "reviewed by expert panel"
    ).all(),
    f"{len(expert_variants)} expert-panel variants checked."
)


# ============================================================
# 12. PHENOTYPE FLAG
# ============================================================

print()
print("=" * 70)
print("12. PHENOTYPE AVAILABILITY")
print("=" * 70)

phenotype_values = set(
    df["phenotype_available"].dropna().unique()
)

check(
    "Phenotype flag uses binary values",
    phenotype_values.issubset({0, 1}),
    f"Observed values: {sorted(phenotype_values)}"
)

check(
    "Phenotype information available",
    (df["phenotype_available"] == 1).sum() == len(df),
    (
        f"{(df['phenotype_available'] == 1).sum()} / "
        f"{len(df)} variants have phenotype information."
    )
)


# ============================================================
# 13. GENOMIC REPRESENTATION FLAGS
# ============================================================

print()
print("=" * 70)
print("13. GENOMIC REPRESENTATION")
print("=" * 70)

for column in [
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present"
]:

    values = set(
        df[column].dropna().unique()
    )

    check(
        f"{column} uses binary values",
        values.issubset({0, 1}),
        f"Observed values: {sorted(values)}"
    )


# ============================================================
# 14. REQUIRED GENOMIC DATA
# ============================================================

check(
    "All variants have genomic coordinates",
    (df["genomic_coordinates_present"] == 1).all(),
    (
        f"{(df['genomic_coordinates_present'] == 1).sum()} / "
        f"{len(df)} variants have genomic coordinates."
    )
)

check(
    "All variants have HGVS representation",
    (df["hgvs_present"] == 1).all(),
    (
        f"{(df['hgvs_present'] == 1).sum()} / "
        f"{len(df)} variants have HGVS."
    )
)


# ============================================================
# 15. EVIDENCE COUNT CONSISTENCY
# ============================================================

print()
print("=" * 70)
print("14. EVIDENCE COUNTS")
print("=" * 70)

evidence_columns = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count"
]

for column in evidence_columns:

    missing = df[column].isna().sum()
    negative = (df[column].dropna() < 0).sum()

    check(
        f"{column} is valid",
        missing == 0 and negative == 0,
        f"{missing} missing, {negative} negative values."
    )


# ============================================================
# 16. CONFLICTING COUNT CONSISTENCY
# ============================================================

conflict_count_consistent = (
    (
        df["conflicting_count"] > 0
    ).astype(int)
    == df["conflict_flag"]
).all()

check(
    "conflicting_count agrees with conflict_flag",
    conflict_count_consistent,
    (
        "All variants are internally consistent."
        if conflict_count_consistent
        else "Some variants have inconsistent conflict evidence."
    )
)


# ============================================================
# 17. PRIORITIZED DATASET
# ============================================================

print()
print("=" * 70)
print("15. PRIORITIZED DATASET")
print("=" * 70)

if os.path.exists(PRIORITIZED_FILE):

    prioritized = pd.read_csv(PRIORITIZED_FILE)

    check(
        "Prioritized dataset exists",
        True,
        f"{len(prioritized)} rows loaded."
    )

    check(
        "Prioritized dataset preserves variants",
        (
            len(prioritized) == len(df)
            and prioritized["variant_id"].nunique()
            == df["variant_id"].nunique()
        ),
        (
            f"Evidence dataset: {len(df)} rows / "
            f"{df['variant_id'].nunique()} variants; "
            f"prioritized dataset: {len(prioritized)} rows / "
            f"{prioritized['variant_id'].nunique()} variants."
        )
    )

    priority_values = set(
        prioritized["review_priority"]
        .dropna()
        .unique()
    )

    expected_priorities = {
        "HIGH REVIEW PRIORITY",
        "MEDIUM REVIEW PRIORITY",
        "LOW REVIEW PRIORITY",
        "REVIEW CLASSIFICATION"
    }

    check(
        "Prioritization categories are valid",
        priority_values.issubset(expected_priorities),
        f"Observed categories: {sorted(priority_values)}"
    )

else:

    check(
        "Prioritized dataset exists",
        False,
        "deepgene_prioritized_v1.csv was not found."
    )


# ============================================================
# 18. MISSING VALUES
# ============================================================

print()
print("=" * 70)
print("16. MISSING VALUES")
print("=" * 70)

missing_summary = df.isna().sum()

missing_columns = missing_summary[
    missing_summary > 0
]

if len(missing_columns) == 0:

    check(
        "No missing values",
        True,
        "No missing values detected."
    )

else:

    print()
    print(missing_columns.to_string())

    # Missing values are reported rather than automatically
    # treated as failure because some biological fields may
    # legitimately be unavailable.

    check(
        "Missing values documented",
        True,
        (
            f"{len(missing_columns)} columns contain missing "
            "values; see summary above."
        )
    )


# ============================================================
# 19. FINAL DATASET SUMMARY
# ============================================================

print()
print("=" * 70)
print("FINAL DATASET SUMMARY")
print("=" * 70)

print()
print(f"Rows:                    {len(df)}")
print(f"Unique variants:         {df['variant_id'].nunique()}")
print(f"Features:                {len(df.columns)}")
print(
    f"Conflicting variants:    "
    f"{df['conflict_flag'].sum()}"
)
print(
    f"Expert-panel variants:   "
    f"{df['expert_panel_review'].sum()}"
)
print(
    f"Multiple-submitter:      "
    f"{df['multiple_submitters'].sum()}"
)
print(
    f"GRCh37 present:          "
    f"{df['grch37_present'].sum()}"
)
print(
    f"GRCh38 present:          "
    f"{df['grch38_present'].sum()}"
)
print(
    f"dbSNP present:           "
    f"{df['dbsnp_present'].sum()}"
)


# ============================================================
# 20. OVERALL RESULT
# ============================================================

failed_checks = [
    r for r in results
    if r["status"] == "FAIL"
]

print()
print("=" * 70)
print("VALIDATION RESULT")
print("=" * 70)

if len(failed_checks) == 0:

    overall_status = "PASS"

    print()
    print("STATUS: PASS")
    print()
    print(
        "The DeepGene V1 evidence dataset passed "
        "all validation checks."
    )

else:

    overall_status = "FAIL"

    print()
    print("STATUS: FAIL")
    print()
    print(
        f"{len(failed_checks)} validation check(s) failed."
    )

    print()
    print("FAILED CHECKS:")

    for item in failed_checks:
        print(
            f"- {item['check']}: "
            f"{item['detail']}"
        )


# ============================================================
# SAVE REPORT
# ============================================================

report_lines = []

report_lines.append(
    "DEEPGENE EVIDENCE DATASET VALIDATION"
)

report_lines.append("=" * 70)

report_lines.append(
    f"Rows: {len(df)}"
)

report_lines.append(
    f"Unique variants: {df['variant_id'].nunique()}"
)

report_lines.append(
    f"Features: {len(df.columns)}"
)

report_lines.append(
    f"Overall status: {overall_status}"
)

report_lines.append("")

for item in results:

    report_lines.append(
        f"[{item['status']}] "
        f"{item['check']}: "
        f"{item['detail']}"
    )

report_lines.append("")
report_lines.append("=" * 70)
report_lines.append("END OF VALIDATION REPORT")


with open(
    OUTPUT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "\n".join(report_lines)
    )


# ============================================================
# END
# ============================================================

print()
print("=" * 70)
print("OUTPUT")
print("=" * 70)

print()
print("Validation report:")
print(OUTPUT_FILE)

print()
print("=" * 70)
print("DATASET VALIDATION COMPLETE")
print("=" * 70)