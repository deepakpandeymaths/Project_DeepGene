# ============================================================
# 30. VALIDATE PHENOTYPE NORMALIZATION
# ============================================================
#
# Purpose:
# Validate the normalized phenotype dataset created in Step 29.
#
# We compare:
#
#   Step 26 phenotype profiles
#       VS
#   Step 29 normalized phenotype profiles
#
# Validation checks:
#
#   1. Row count
#   2. Variant ID uniqueness
#   3. Variant ID preservation
#   4. Original phenotype preservation
#   5. Normalized phenotype completeness
#   6. Semicolon removal
#   7. Duplicate removal
#   8. Empty-term check
#   9. Normalized term-count consistency
#  10. Non-phenotype column preservation
#  11. Expected normalization changes
#
# IMPORTANT:
# This script DOES NOT modify any dataset.
#
# ============================================================

import os
import pandas as pd


# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

BASE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

PROCESSED_DIR = os.path.join(
    BASE_DIR,
    "processed"
)

RESULTS_DIR = os.path.join(
    BASE_DIR,
    "analysis_results"
)

os.makedirs(RESULTS_DIR, exist_ok=True)

ORIGINAL_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_phenotype_profiles_v1.csv"
)

NORMALIZED_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_phenotype_normalized_v1.csv"
)

REPORT_FILE = os.path.join(
    RESULTS_DIR,
    "30_validate_phenotype_normalization.txt"
)


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 30")
print("VALIDATE PHENOTYPE NORMALIZATION")
print("=" * 70)


# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

print("\nLoading original phenotype dataset...")

original = pd.read_csv(
    ORIGINAL_FILE
)

print(
    f"Original rows: {len(original)}"
)

print(
    f"Original columns: {len(original.columns)}"
)


print("\nLoading normalized phenotype dataset...")

normalized = pd.read_csv(
    NORMALIZED_FILE
)

print(
    f"Normalized rows: {len(normalized)}"
)

print(
    f"Normalized columns: {len(normalized.columns)}"
)


# ------------------------------------------------------------
# VALIDATION STATUS
# ------------------------------------------------------------

checks = []


def record_check(name, passed, detail):

    checks.append(
        {
            "name": name,
            "passed": passed,
            "detail": detail
        }
    )

    status = "PASS" if passed else "FAIL"

    print(
        f"[{status}] {name}: {detail}"
    )


# ------------------------------------------------------------
# 1. ROW COUNT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. ROW COUNT")
print("=" * 70)

row_match = (
    len(original)
    ==
    len(normalized)
)

record_check(
    "Row count",
    row_match,
    f"{len(original)} -> {len(normalized)}"
)


# ------------------------------------------------------------
# 2. VARIANT ID UNIQUENESS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("2. VARIANT ID UNIQUENESS")
print("=" * 70)

original_duplicate_ids = int(
    original["variant_id"]
    .duplicated()
    .sum()
)

normalized_duplicate_ids = int(
    normalized["variant_id"]
    .duplicated()
    .sum()
)

record_check(
    "Original variant uniqueness",
    original_duplicate_ids == 0,
    f"duplicates = {original_duplicate_ids}"
)

record_check(
    "Normalized variant uniqueness",
    normalized_duplicate_ids == 0,
    f"duplicates = {normalized_duplicate_ids}"
)


# ------------------------------------------------------------
# 3. VARIANT ID PRESERVATION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. VARIANT ID PRESERVATION")
print("=" * 70)

original_ids = set(
    original["variant_id"]
)

normalized_ids = set(
    normalized["variant_id"]
)

missing_ids = (
    original_ids
    -
    normalized_ids
)

extra_ids = (
    normalized_ids
    -
    original_ids
)

record_check(
    "No original variants missing",
    len(missing_ids) == 0,
    f"missing = {len(missing_ids)}"
)

record_check(
    "No unexpected variants",
    len(extra_ids) == 0,
    f"extra = {len(extra_ids)}"
)


# ------------------------------------------------------------
# 4. ORIGINAL PHENOTYPE PRESERVATION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. ORIGINAL PHENOTYPE PRESERVATION")
print("=" * 70)

original_compare = (
    original[
        [
            "variant_id",
            "meaningful_phenotype_terms"
        ]
    ]
    .copy()
)

normalized_compare = (
    normalized[
        [
            "variant_id",
            "meaningful_phenotype_terms"
        ]
    ]
    .copy()
)

original_compare = (
    original_compare
    .sort_values("variant_id")
    .reset_index(drop=True)
)

normalized_compare = (
    normalized_compare
    .sort_values("variant_id")
    .reset_index(drop=True)
)

original_phenotype_equal = (
    original_compare[
        "meaningful_phenotype_terms"
    ]
    .fillna("")
    .astype(str)
    ==
    normalized_compare[
        "meaningful_phenotype_terms"
    ]
    .fillna("")
    .astype(str)
)

original_phenotype_mismatches = int(
    (~original_phenotype_equal).sum()
)

record_check(
    "Original phenotype field preserved",
    original_phenotype_mismatches == 0,
    f"mismatches = {original_phenotype_mismatches}"
)


# ------------------------------------------------------------
# 5. NORMALIZED COLUMN EXISTS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. NORMALIZED FIELD")
print("=" * 70)

required_columns = [
    "normalized_phenotype_terms",
    "original_phenotype_term_count",
    "normalized_phenotype_term_count",
    "phenotype_normalization_changed"
]

missing_required_columns = [
    column
    for column in required_columns
    if column not in normalized.columns
]

record_check(
    "Required normalized columns",
    len(missing_required_columns) == 0,
    f"missing = {missing_required_columns}"
)


# ------------------------------------------------------------
# 6. SEMICOLON CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("6. SEMICOLON REMOVAL")
print("=" * 70)

remaining_semicolons = int(
    normalized[
        "normalized_phenotype_terms"
    ]
    .fillna("")
    .astype(str)
    .str.contains(
        ";",
        regex=False
    )
    .sum()
)

record_check(
    "No semicolons remain",
    remaining_semicolons == 0,
    f"remaining = {remaining_semicolons}"
)


# ------------------------------------------------------------
# 7. EMPTY TERM CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("7. EMPTY TERM CHECK")
print("=" * 70)


def has_empty_term(value):

    if pd.isna(value):
        return False

    value = str(value).strip()

    if not value:
        return False

    terms = value.split("|")

    return any(
        not term.strip()
        for term in terms
    )


empty_term_count = int(
    normalized[
        "normalized_phenotype_terms"
    ]
    .apply(has_empty_term)
    .sum()
)

record_check(
    "No empty normalized terms",
    empty_term_count == 0,
    f"affected variants = {empty_term_count}"
)


# ------------------------------------------------------------
# 8. DUPLICATE TERM CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("8. DUPLICATE TERM CHECK")
print("=" * 70)


def has_duplicate_terms(value):

    if pd.isna(value):
        return False

    value = str(value).strip()

    if not value:
        return False

    terms = [
        term.strip().lower()
        for term in value.split("|")
        if term.strip()
    ]

    return len(terms) != len(
        set(terms)
    )


duplicate_term_count = int(
    normalized[
        "normalized_phenotype_terms"
    ]
    .apply(has_duplicate_terms)
    .sum()
)

record_check(
    "No duplicate normalized terms",
    duplicate_term_count == 0,
    f"affected variants = {duplicate_term_count}"
)


# ------------------------------------------------------------
# 9. NORMALIZED TERM COUNT CONSISTENCY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("9. NORMALIZED TERM COUNT CONSISTENCY")
print("=" * 70)


def count_terms(value):

    if pd.isna(value):
        return 0

    value = str(value).strip()

    if not value:
        return 0

    return len(
        [
            term
            for term in value.split("|")
            if term.strip()
        ]
    )


calculated_counts = (
    normalized[
        "normalized_phenotype_terms"
    ]
    .apply(count_terms)
)

stored_counts = (
    normalized[
        "normalized_phenotype_term_count"
    ]
    .fillna(0)
    .astype(int)
)

count_mismatches = int(
    (
        calculated_counts
        !=
        stored_counts
    ).sum()
)

record_check(
    "Normalized term counts are consistent",
    count_mismatches == 0,
    f"mismatches = {count_mismatches}"
)


# ------------------------------------------------------------
# 10. NORMALIZATION FLAG CONSISTENCY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("10. NORMALIZATION FLAG CONSISTENCY")
print("=" * 70)

calculated_change_flag = (
    normalized[
        "normalized_phenotype_terms"
    ]
    .fillna("")
    .astype(str)
    !=
    normalized[
        "meaningful_phenotype_terms"
    ]
    .fillna("")
    .astype(str)
)

stored_change_flag = (
    normalized[
        "phenotype_normalization_changed"
    ]
    .fillna(False)
    .astype(bool)
)

flag_mismatches = int(
    (
        calculated_change_flag
        !=
        stored_change_flag
    ).sum()
)

record_check(
    "Normalization flags are consistent",
    flag_mismatches == 0,
    f"mismatches = {flag_mismatches}"
)


# ------------------------------------------------------------
# 11. EXPECTED NUMBER OF CHANGED VARIANTS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("11. EXPECTED NORMALIZATION CHANGES")
print("=" * 70)

semicolon_variants_original = int(
    original[
        "meaningful_phenotype_terms"
    ]
    .fillna("")
    .astype(str)
    .str.contains(
        ";",
        regex=False
    )
    .sum()
)

changed_variants = int(
    normalized[
        "phenotype_normalization_changed"
    ]
    .sum()
)

print(
    f"Original semicolon-containing variants: "
    f"{semicolon_variants_original}"
)

print(
    f"Normalized changed variants: "
    f"{changed_variants}"
)

record_check(
    "Changed variants match semicolon-containing variants",
    changed_variants
    ==
    semicolon_variants_original,
    f"{changed_variants} vs {semicolon_variants_original}"
)


# ------------------------------------------------------------
# 12. NON-PHENOTYPE COLUMN PRESERVATION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("12. NON-PHENOTYPE COLUMN PRESERVATION")
print("=" * 70)

phenotype_specific_columns = {
    "normalized_phenotype_terms",
    "original_phenotype_term_count",
    "normalized_phenotype_term_count",
    "phenotype_normalization_changed"
}

original_columns = set(
    original.columns
)

normalized_columns = set(
    normalized.columns
)

expected_original_columns = (
    normalized_columns
    -
    phenotype_specific_columns
)

missing_original_columns = (
    expected_original_columns
    -
    original_columns
)

record_check(
    "Original columns preserved",
    len(missing_original_columns) == 0,
    f"missing = {len(missing_original_columns)}"
)


# ------------------------------------------------------------
# 13. FINAL VALIDATION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("13. FINAL VALIDATION")
print("=" * 70)

failed_checks = [
    check
    for check in checks
    if not check["passed"]
]

passed_checks = [
    check
    for check in checks
    if check["passed"]
]

print(
    f"Checks passed: {len(passed_checks)}"
)

print(
    f"Checks failed: {len(failed_checks)}"
)


if len(failed_checks) == 0:

    overall_status = "PASS"

    print(
        "\nOVERALL STATUS: PASS"
    )

    print(
        "\nThe normalized phenotype dataset passed "
        "all validation checks."
    )

else:

    overall_status = "FAIL"

    print(
        "\nOVERALL STATUS: FAIL"
    )

    print(
        "\nFailed checks:"
    )

    for check in failed_checks:

        print(
            f"- {check['name']}: "
            f"{check['detail']}"
        )


# ------------------------------------------------------------
# 14. SAVE REPORT
# ------------------------------------------------------------

with open(
    REPORT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "DEEPGENE STEP 30\n"
        "VALIDATE PHENOTYPE NORMALIZATION\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Original rows: {len(original)}\n"
    )

    f.write(
        f"Normalized rows: {len(normalized)}\n"
    )

    f.write(
        f"Original semicolon variants: "
        f"{semicolon_variants_original}\n"
    )

    f.write(
        f"Normalized changed variants: "
        f"{changed_variants}\n"
    )

    f.write(
        f"Remaining semicolons: "
        f"{remaining_semicolons}\n"
    )

    f.write(
        f"Duplicate normalized variants: "
        f"{duplicate_term_count}\n"
    )

    f.write(
        f"Empty normalized variants: "
        f"{empty_term_count}\n"
    )

    f.write(
        f"Term count mismatches: "
        f"{count_mismatches}\n"
    )

    f.write(
        f"Normalization flag mismatches: "
        f"{flag_mismatches}\n"
    )

    f.write("\n")

    f.write(
        f"Checks passed: {len(passed_checks)}\n"
    )

    f.write(
        f"Checks failed: {len(failed_checks)}\n"
    )

    f.write(
        f"\nOVERALL STATUS: {overall_status}\n"
    )


# ------------------------------------------------------------
# FINISH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("REPORT SAVED")
print("=" * 70)

print(REPORT_FILE)

print("\nStep 30 complete.")