# ============================================================
# 26. PHENOTYPE REPRESENTATION
# ============================================================
#
# Purpose:
# Convert the raw phenotype string into a structured
# phenotype representation while preserving the original
# phenotype field.
#
# This step does NOT assign pathogenicity.
# It does NOT create an ML target.
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
    "deepgene_v1_final.csv"
)

OUTPUT_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_phenotype_profiles_v1.csv"
)

REPORT_FILE = os.path.join(
    RESULTS_DIR,
    "26_phenotype_representation.txt"
)


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 26")
print("PHENOTYPE REPRESENTATION")
print("=" * 70)


# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

print("\nLoading final V1 dataset...")

df = pd.read_csv(INPUT_FILE)

print(f"Rows loaded: {len(df)}")
print(f"Columns loaded: {len(df.columns)}")


# ------------------------------------------------------------
# NORMALIZE ORIGINAL PHENOTYPE FIELD
# ------------------------------------------------------------

df["phenotypes_clean"] = (
    df["phenotypes"]
    .fillna("")
    .astype(str)
    .str.strip()
)


# ------------------------------------------------------------
# NON-INFORMATIVE TERMS
# ------------------------------------------------------------

NON_INFORMATIVE_TERMS = {
    "",
    "not provided",
    "not specified",
    "unknown",
    "unspecified",
    "none",
    "n/a",
    "na",
    "-"
}


# ------------------------------------------------------------
# SPLIT PHENOTYPE TERMS
# ------------------------------------------------------------

def split_terms(value):

    if not value:
        return []

    terms = value.split("|")

    cleaned = []

    for term in terms:

        term = term.strip()

        if term:
            cleaned.append(term)

    return cleaned


df["phenotype_terms"] = (
    df["phenotypes_clean"]
    .apply(split_terms)
)


# ------------------------------------------------------------
# MEANINGFUL TERMS
# ------------------------------------------------------------

def get_meaningful_terms(terms):

    meaningful = []

    for term in terms:

        if term.lower() not in NON_INFORMATIVE_TERMS:

            meaningful.append(term)

    return meaningful


df["meaningful_phenotype_terms_list"] = (
    df["phenotype_terms"]
    .apply(get_meaningful_terms)
)


# ------------------------------------------------------------
# NON-INFORMATIVE TERMS
# ------------------------------------------------------------

def get_non_informative_terms(terms):

    non_informative = []

    for term in terms:

        if term.lower() in NON_INFORMATIVE_TERMS:

            non_informative.append(term)

    return non_informative


df["non_informative_phenotype_terms_list"] = (
    df["phenotype_terms"]
    .apply(get_non_informative_terms)
)


# ------------------------------------------------------------
# COUNTS
# ------------------------------------------------------------

df["phenotype_term_count"] = (
    df["phenotype_terms"]
    .apply(len)
)

df["meaningful_phenotype_term_count"] = (
    df["meaningful_phenotype_terms_list"]
    .apply(len)
)

df["non_informative_term_count"] = (
    df["non_informative_phenotype_terms_list"]
    .apply(len)
)


# ------------------------------------------------------------
# STRUCTURED TEXT REPRESENTATION
# ------------------------------------------------------------

def join_terms(terms):

    return "|".join(terms)


df["meaningful_phenotype_terms"] = (
    df["meaningful_phenotype_terms_list"]
    .apply(join_terms)
)

df["non_informative_phenotype_terms"] = (
    df["non_informative_phenotype_terms_list"]
    .apply(join_terms)
)


# ------------------------------------------------------------
# INFORMATION CATEGORY
# ------------------------------------------------------------

def classify_information(row):

    meaningful = row[
        "meaningful_phenotype_term_count"
    ]

    total = row[
        "phenotype_term_count"
    ]

    if total == 0:

        return "NO PHENOTYPE INFORMATION"

    if meaningful == 0:

        return "NON-INFORMATIVE PHENOTYPE ONLY"

    if meaningful == 1:

        return "ONE MEANINGFUL PHENOTYPE TERM"

    return "MULTIPLE MEANINGFUL PHENOTYPE TERMS"


df["phenotype_information_category"] = (
    df.apply(
        classify_information,
        axis=1
    )
)


# ------------------------------------------------------------
# INFORMATION FLAG
# ------------------------------------------------------------

df["meaningful_phenotype_available"] = (
    df["meaningful_phenotype_term_count"] > 0
).astype(int)


# ------------------------------------------------------------
# CREATE OUTPUT DATASET
# ------------------------------------------------------------
#
# Keep the original fields.
# Add structured phenotype representation.
#
# ------------------------------------------------------------

output_columns = list(
    pd.read_csv(INPUT_FILE, nrows=1).columns
)

new_columns = [
    "phenotype_term_count",
    "meaningful_phenotype_term_count",
    "non_informative_term_count",
    "meaningful_phenotype_available",
    "phenotype_information_category",
    "meaningful_phenotype_terms",
    "non_informative_phenotype_terms"
]

output_columns.extend(new_columns)

output_df = df[output_columns].copy()


# ------------------------------------------------------------
# VALIDATION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. REPRESENTATION VALIDATION")
print("=" * 70)

print(
    f"Output rows: "
    f"{len(output_df)}"
)

print(
    f"Output columns: "
    f"{len(output_df.columns)}"
)

print(
    f"Unique variant IDs: "
    f"{output_df['variant_id'].nunique()}"
)

print(
    f"Duplicate variant IDs: "
    f"{output_df['variant_id'].duplicated().sum()}"
)


# ------------------------------------------------------------
# CATEGORY DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("2. PHENOTYPE INFORMATION CATEGORIES")
print("=" * 70)

category_counts = (
    output_df[
        "phenotype_information_category"
    ]
    .value_counts()
)

print(
    category_counts.to_string()
)


# ------------------------------------------------------------
# MEANINGFUL PHENOTYPE COVERAGE
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. MEANINGFUL PHENOTYPE COVERAGE")
print("=" * 70)

meaningful_count = (
    output_df[
        "meaningful_phenotype_available"
    ] == 1
).sum()

no_meaningful_count = (
    output_df[
        "meaningful_phenotype_available"
    ] == 0
).sum()

total = len(output_df)

print(
    f"Meaningful phenotype available: "
    f"{meaningful_count}"
)

print(
    f"No meaningful phenotype: "
    f"{no_meaningful_count}"
)

print(
    f"Meaningful phenotype coverage: "
    f"{meaningful_count / total * 100:.2f}%"
)


# ------------------------------------------------------------
# TERM COUNT DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. MEANINGFUL PHENOTYPE TERM DISTRIBUTION")
print("=" * 70)

print(
    output_df[
        "meaningful_phenotype_term_count"
    ]
    .value_counts()
    .sort_index()
    .to_string()
)


# ------------------------------------------------------------
# EXAMPLES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. EXAMPLE PHENOTYPE PROFILES")
print("=" * 70)

print(
    output_df[
        [
            "variant_id",
            "phenotypes",
            "meaningful_phenotype_terms",
            "non_informative_phenotype_terms",
            "meaningful_phenotype_term_count",
            "phenotype_information_category"
        ]
    ]
    .head(15)
    .to_string(index=False)
)


# ------------------------------------------------------------
# SAVE DATASET
# ------------------------------------------------------------

output_df.to_csv(
    OUTPUT_FILE,
    index=False
)


# ------------------------------------------------------------
# SAVE REPORT
# ------------------------------------------------------------

with open(
    REPORT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "DEEPGENE STEP 26\n"
        "PHENOTYPE REPRESENTATION\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Input rows: {len(df)}\n"
    )

    f.write(
        f"Output rows: {len(output_df)}\n"
    )

    f.write(
        f"Output columns: {len(output_df.columns)}\n"
    )

    f.write(
        f"Unique variants: "
        f"{output_df['variant_id'].nunique()}\n"
    )

    f.write(
        f"Duplicate variant IDs: "
        f"{output_df['variant_id'].duplicated().sum()}\n\n"
    )

    f.write(
        "PHENOTYPE INFORMATION CATEGORIES\n"
    )

    f.write("-" * 70 + "\n")

    for category, count in category_counts.items():

        f.write(
            f"{category}: {count}\n"
        )

    f.write("\n")

    f.write(
        f"Meaningful phenotype coverage: "
        f"{meaningful_count / total * 100:.2f}%\n"
    )

    f.write(
        f"No meaningful phenotype: "
        f"{no_meaningful_count / total * 100:.2f}%\n"
    )


# ------------------------------------------------------------
# FINISH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("OUTPUT SAVED")
print("=" * 70)

print(
    f"Dataset: {OUTPUT_FILE}"
)

print(
    f"Report: {REPORT_FILE}"
)

print("\nStep 26 complete.")