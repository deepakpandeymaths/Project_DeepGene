# ============================================================
# 29. PHENOTYPE NORMALIZATION
# ============================================================
#
# Purpose:
# Create a normalized phenotype representation from the
# phenotype profiles generated in Step 26.
#
# Normalization rules:
#
#   1. Preserve the original phenotype fields.
#   2. Treat "|" and ";" as phenotype separators.
#   3. Remove leading/trailing whitespace.
#   4. Remove empty terms.
#   5. Remove exact duplicate phenotype terms within a variant
#      using case-insensitive comparison.
#   6. Preserve the spelling of the first occurrence.
#
# IMPORTANT:
# This creates a NEW derived dataset.
# It does NOT overwrite Step 26 data.
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

INPUT_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_phenotype_profiles_v1.csv"
)

OUTPUT_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_phenotype_normalized_v1.csv"
)

REPORT_FILE = os.path.join(
    RESULTS_DIR,
    "29_phenotype_normalization.txt"
)


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 29")
print("PHENOTYPE NORMALIZATION")
print("=" * 70)


# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

print("\nLoading phenotype profile dataset...")

df = pd.read_csv(INPUT_FILE)

print(f"Rows loaded: {len(df)}")
print(f"Columns loaded: {len(df.columns)}")


# ------------------------------------------------------------
# BASIC TEXT HANDLING
# ------------------------------------------------------------

df["meaningful_phenotype_terms"] = (
    df["meaningful_phenotype_terms"]
    .fillna("")
    .astype(str)
    .str.strip()
)


# ------------------------------------------------------------
# NORMALIZATION FUNCTION
# ------------------------------------------------------------

def normalize_phenotype_terms(value):

    if pd.isna(value):
        return []

    value = str(value).strip()

    if not value:
        return []

    # Treat both separators as delimiters.
    value = value.replace(";", "|")

    raw_terms = value.split("|")

    normalized_terms = []
    seen = set()

    for term in raw_terms:

        term = term.strip()

        if not term:
            continue

        key = term.lower()

        if key in seen:
            continue

        seen.add(key)

        normalized_terms.append(term)

    return normalized_terms


# ------------------------------------------------------------
# CREATE NORMALIZED TERMS
# ------------------------------------------------------------

print("\nNormalizing phenotype terms...")

df["normalized_phenotype_terms_list"] = (
    df["meaningful_phenotype_terms"]
    .apply(normalize_phenotype_terms)
)


# ------------------------------------------------------------
# CREATE STRING REPRESENTATION
# ------------------------------------------------------------

df["normalized_phenotype_terms"] = (
    df["normalized_phenotype_terms_list"]
    .apply(
        lambda terms: "|".join(terms)
    )
)


# ------------------------------------------------------------
# REMOVE TEMPORARY LIST COLUMN
# ------------------------------------------------------------

df.drop(
    columns=["normalized_phenotype_terms_list"],
    inplace=True
)


# ------------------------------------------------------------
# TERM COUNT FUNCTIONS
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# ORIGINAL TERM COUNT
# ------------------------------------------------------------

df["original_phenotype_term_count"] = (
    df["meaningful_phenotype_terms"]
    .apply(count_terms)
)


# ------------------------------------------------------------
# NORMALIZED TERM COUNT
# ------------------------------------------------------------

df["normalized_phenotype_term_count"] = (
    df["normalized_phenotype_terms"]
    .apply(count_terms)
)


# ------------------------------------------------------------
# NORMALIZATION CHANGE FLAG
# ------------------------------------------------------------

df["phenotype_normalization_changed"] = (
    df["normalized_phenotype_terms"]
    !=
    df["meaningful_phenotype_terms"]
)


# ------------------------------------------------------------
# 1. BASIC SUMMARY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. NORMALIZATION SUMMARY")
print("=" * 70)

changed_count = int(
    df["phenotype_normalization_changed"].sum()
)

unchanged_count = (
    len(df) - changed_count
)

print(
    f"Variants processed: "
    f"{len(df)}"
)

print(
    f"Variants changed by normalization: "
    f"{changed_count}"
)

print(
    f"Variants unchanged: "
    f"{unchanged_count}"
)


# ------------------------------------------------------------
# 2. TERM COUNTS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("2. TERM COUNT COMPARISON")
print("=" * 70)

original_total = int(
    df["original_phenotype_term_count"].sum()
)

normalized_total = int(
    df["normalized_phenotype_term_count"].sum()
)

removed_total = (
    original_total - normalized_total
)

print(
    f"Original term occurrences: "
    f"{original_total}"
)

print(
    f"Normalized term occurrences: "
    f"{normalized_total}"
)

print(
    f"Terms removed during normalization: "
    f"{removed_total}"
)


# ------------------------------------------------------------
# 3. NORMALIZATION CHANGE DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. NORMALIZATION CHANGE DISTRIBUTION")
print("=" * 70)

changed_df = df[
    df["phenotype_normalization_changed"]
].copy()

print(
    f"Changed variants: "
    f"{len(changed_df)}"
)

if len(changed_df) > 0:

    changed_df["terms_removed"] = (
        changed_df["original_phenotype_term_count"]
        -
        changed_df["normalized_phenotype_term_count"]
    )

    print(
        "\nTerms removed per changed variant:"
    )

    print(
        changed_df["terms_removed"]
        .value_counts()
        .sort_index()
        .to_string()
    )


# ------------------------------------------------------------
# 4. EXAMPLES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. NORMALIZATION EXAMPLES")
print("=" * 70)

if len(changed_df) > 0:

    for _, row in changed_df.head(20).iterrows():

        print(
            "\nVariant:",
            row["variant_id"]
        )

        print("\nOriginal:")

        print(
            row["meaningful_phenotype_terms"]
        )

        print("\nNormalized:")

        print(
            row["normalized_phenotype_terms"]
        )

        print(
            "\nOriginal term count:",
            row["original_phenotype_term_count"]
        )

        print(
            "Normalized term count:",
            row["normalized_phenotype_term_count"]
        )


# ------------------------------------------------------------
# 5. DUPLICATE CHECK AFTER NORMALIZATION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. DUPLICATE CHECK AFTER NORMALIZATION")
print("=" * 70)


def has_duplicate_normalized_terms(value):

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

    return len(terms) != len(set(terms))


duplicate_mask = (
    df["normalized_phenotype_terms"]
    .apply(
        has_duplicate_normalized_terms
    )
)

duplicate_count = int(
    duplicate_mask.sum()
)

print(
    f"Variants with duplicate normalized terms: "
    f"{duplicate_count}"
)


# ------------------------------------------------------------
# 6. SEMICOLON CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("6. SEMICOLON CHECK")
print("=" * 70)

remaining_semicolons = int(
    df["normalized_phenotype_terms"]
    .str.contains(
        ";",
        regex=False,
        na=False
    )
    .sum()
)

print(
    f"Normalized entries still containing semicolons: "
    f"{remaining_semicolons}"
)


# ------------------------------------------------------------
# 7. EMPTY TERM CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("7. EMPTY TERM CHECK")
print("=" * 70)


def has_empty_normalized_terms(value):

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


empty_normalized = int(
    df["normalized_phenotype_terms"]
    .apply(
        has_empty_normalized_terms
    )
    .sum()
)

print(
    f"Variants containing empty normalized terms: "
    f"{empty_normalized}"
)


# ------------------------------------------------------------
# 8. UNIQUE TERM COUNT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("8. UNIQUE NORMALIZED PHENOTYPE TERMS")
print("=" * 70)

all_normalized_terms = []

for value in df[
    "normalized_phenotype_terms"
]:

    if pd.isna(value):
        continue

    value = str(value).strip()

    if not value:
        continue

    all_normalized_terms.extend(
        [
            term.strip()
            for term in value.split("|")
            if term.strip()
        ]
    )


normalized_term_series = pd.Series(
    all_normalized_terms,
    dtype="object"
)

unique_normalized_terms = (
    normalized_term_series
    .str.lower()
    .nunique()
)

print(
    f"Total normalized term occurrences: "
    f"{len(normalized_term_series)}"
)

print(
    f"Unique normalized phenotype terms: "
    f"{unique_normalized_terms}"
)

print(
    "\nTop 30 normalized phenotype terms:"
)

print(
    normalized_term_series
    .value_counts()
    .head(30)
    .to_string()
)


# ------------------------------------------------------------
# 9. DATASET INTEGRITY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("9. DATASET INTEGRITY")
print("=" * 70)

print(
    f"Rows: {len(df)}"
)

print(
    f"Unique variant IDs: "
    f"{df['variant_id'].nunique()}"
)

print(
    f"Duplicate variant IDs: "
    f"{df['variant_id'].duplicated().sum()}"
)

print(
    f"Columns: {len(df.columns)}"
)


# ------------------------------------------------------------
# 10. SAVE NORMALIZED DATASET
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("10. SAVING NORMALIZED DATASET")
print("=" * 70)

df.to_csv(
    OUTPUT_FILE,
    index=False
)

print(
    "Normalized dataset saved:"
)

print(OUTPUT_FILE)


# ------------------------------------------------------------
# 11. SAVE REPORT
# ------------------------------------------------------------

with open(
    REPORT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "DEEPGENE STEP 29\n"
        "PHENOTYPE NORMALIZATION\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Variants processed: {len(df)}\n"
    )

    f.write(
        f"Variants changed: {changed_count}\n"
    )

    f.write(
        f"Variants unchanged: {unchanged_count}\n"
    )

    f.write(
        f"Original term occurrences: "
        f"{original_total}\n"
    )

    f.write(
        f"Normalized term occurrences: "
        f"{normalized_total}\n"
    )

    f.write(
        f"Terms removed: "
        f"{removed_total}\n"
    )

    f.write(
        f"Duplicate normalized variants: "
        f"{duplicate_count}\n"
    )

    f.write(
        f"Remaining semicolon entries: "
        f"{remaining_semicolons}\n"
    )

    f.write(
        f"Empty normalized entries: "
        f"{empty_normalized}\n"
    )

    f.write(
        f"Unique normalized phenotype terms: "
        f"{unique_normalized_terms}\n"
    )

    f.write("\n")

    f.write(
        "IMPORTANT:\n"
    )

    f.write(
        "The original phenotype representation remains preserved.\n"
    )

    f.write(
        "The normalized representation is a separate derived field.\n"
    )


# ------------------------------------------------------------
# FINISH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("STEP 29 COMPLETE")
print("=" * 70)

print(
    "\nNormalized dataset:"
)

print(OUTPUT_FILE)

print(
    "\nReport:"
)

print(REPORT_FILE)