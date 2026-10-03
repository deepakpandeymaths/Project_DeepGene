# ============================================================
# 25. PHENOTYPE INFORMATION QUALITY
# ============================================================
#
# Purpose:
# Characterize how informative the phenotype field is for
# each SCN1A variant.
#
# This step does NOT modify the V1 dataset.
#
# It measures:
#   - number of phenotype terms
#   - missing-style terms
#   - meaningful phenotype terms
#   - phenotype information categories
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
    RESULTS_DIR,
    "25_phenotype_information_quality.txt"
)


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 25")
print("PHENOTYPE INFORMATION QUALITY")
print("=" * 70)


# ------------------------------------------------------------
# LOAD DATASET
# ------------------------------------------------------------

print("\nLoading final V1 dataset...")

df = pd.read_csv(INPUT_FILE)

print(f"Rows loaded: {len(df)}")
print(f"Columns loaded: {len(df.columns)}")


# ------------------------------------------------------------
# NORMALIZE PHENOTYPE FIELD
# ------------------------------------------------------------

df["phenotypes_clean"] = (
    df["phenotypes"]
    .fillna("")
    .astype(str)
    .str.strip()
)


# ------------------------------------------------------------
# MISSING / NON-INFORMATIVE TERMS
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

    if value == "":
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
# COUNT TOTAL TERMS
# ------------------------------------------------------------

df["phenotype_term_count"] = (
    df["phenotype_terms"]
    .apply(len)
)


# ------------------------------------------------------------
# COUNT NON-INFORMATIVE TERMS
# ------------------------------------------------------------

def count_non_informative(terms):

    count = 0

    for term in terms:

        if term.lower() in NON_INFORMATIVE_TERMS:
            count += 1

    return count


df["non_informative_term_count"] = (
    df["phenotype_terms"]
    .apply(count_non_informative)
)


# ------------------------------------------------------------
# COUNT MEANINGFUL TERMS
# ------------------------------------------------------------

df["meaningful_phenotype_term_count"] = (
    df["phenotype_term_count"]
    - df["non_informative_term_count"]
)


# ------------------------------------------------------------
# INFORMATION CATEGORY
# ------------------------------------------------------------

def classify_information(row):

    total = row["phenotype_term_count"]

    meaningful = (
        row["meaningful_phenotype_term_count"]
    )

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
# 1. BASIC PHENOTYPE INFORMATION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. PHENOTYPE INFORMATION SUMMARY")
print("=" * 70)

print(
    f"Variants: {len(df)}"
)

print(
    f"Variants with phenotype field: "
    f"{(df['phenotype_term_count'] > 0).sum()}"
)

print(
    f"Variants with at least one meaningful "
    f"phenotype term: "
    f"{(df['meaningful_phenotype_term_count'] > 0).sum()}"
)

print(
    f"Variants with zero meaningful phenotype terms: "
    f"{(df['meaningful_phenotype_term_count'] == 0).sum()}"
)


# ------------------------------------------------------------
# 2. TOTAL PHENOTYPE TERM DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("2. TOTAL PHENOTYPE TERM COUNT")
print("=" * 70)

print(
    df["phenotype_term_count"]
    .value_counts()
    .sort_index()
    .to_string()
)


# ------------------------------------------------------------
# 3. MEANINGFUL PHENOTYPE TERM DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. MEANINGFUL PHENOTYPE TERM COUNT")
print("=" * 70)

print(
    df["meaningful_phenotype_term_count"]
    .value_counts()
    .sort_index()
    .to_string()
)


# ------------------------------------------------------------
# 4. NON-INFORMATIVE TERM DISTRIBUTION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. NON-INFORMATIVE PHENOTYPE TERMS")
print("=" * 70)

print(
    df["non_informative_term_count"]
    .value_counts()
    .sort_index()
    .to_string()
)


# ------------------------------------------------------------
# 5. INFORMATION CATEGORY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. PHENOTYPE INFORMATION CATEGORY")
print("=" * 70)

category_counts = (
    df["phenotype_information_category"]
    .value_counts()
)

print(
    category_counts.to_string()
)


# ------------------------------------------------------------
# 6. PHENOTYPE COVERAGE PERCENTAGES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("6. PHENOTYPE INFORMATION COVERAGE")
print("=" * 70)

total = len(df)

meaningful_count = (
    df["meaningful_phenotype_term_count"] > 0
).sum()

no_meaningful_count = (
    df["meaningful_phenotype_term_count"] == 0
).sum()

print(
    f"At least one meaningful phenotype term: "
    f"{meaningful_count} "
    f"({meaningful_count / total * 100:.2f}%)"
)

print(
    f"No meaningful phenotype term: "
    f"{no_meaningful_count} "
    f"({no_meaningful_count / total * 100:.2f}%)"
)


# ------------------------------------------------------------
# 7. COMMON MEANINGFUL PHENOTYPE TERMS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("7. MOST COMMON MEANINGFUL PHENOTYPE TERMS")
print("=" * 70)

meaningful_terms = []

for terms in df["phenotype_terms"]:

    for term in terms:

        if term.lower() not in NON_INFORMATIVE_TERMS:

            meaningful_terms.append(term)


meaningful_term_series = pd.Series(
    meaningful_terms
)

print(
    meaningful_term_series
    .value_counts()
    .head(30)
    .to_string()
)


# ------------------------------------------------------------
# 8. VARIANTS WITH ONLY NON-INFORMATIVE PHENOTYPES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("8. NON-INFORMATIVE-ONLY VARIANTS")
print("=" * 70)

non_informative_only = df[
    df["meaningful_phenotype_term_count"] == 0
]

print(
    f"Variants with no meaningful phenotype term: "
    f"{len(non_informative_only)}"
)

if len(non_informative_only) > 0:

    print("\nExamples:")

    print(
        non_informative_only[
            [
                "variant_id",
                "phenotypes"
            ]
        ]
        .head(20)
        .to_string(index=False)
    )


# ------------------------------------------------------------
# 9. EXAMPLES WITH MEANINGFUL PHENOTYPES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("9. EXAMPLES WITH MEANINGFUL PHENOTYPES")
print("=" * 70)

meaningful_examples = df[
    df["meaningful_phenotype_term_count"] > 0
]

print(
    meaningful_examples[
        [
            "variant_id",
            "phenotypes",
            "meaningful_phenotype_term_count"
        ]
    ]
    .head(20)
    .to_string(index=False)
)


# ------------------------------------------------------------
# 10. QUALITY CONCLUSION
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("10. QUALITY CONCLUSION")
print("=" * 70)

if meaningful_count == total:

    print(
        "All variants contain at least one meaningful "
        "phenotype term."
    )

elif meaningful_count > 0:

    print(
        "Phenotype information is partially informative."
    )

    print(
        "Some variants contain only missing-style or "
        "non-informative phenotype terms."
    )

else:

    print(
        "No variants contain meaningful phenotype terms."
    )


# ------------------------------------------------------------
# SAVE REPORT
# ------------------------------------------------------------

with open(
    OUTPUT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "DEEPGENE STEP 25\n"
        "PHENOTYPE INFORMATION QUALITY\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Total variants: {total}\n"
    )

    f.write(
        f"Variants with at least one meaningful "
        f"phenotype term: {meaningful_count}\n"
    )

    f.write(
        f"Variants with no meaningful phenotype "
        f"term: {no_meaningful_count}\n"
    )

    f.write(
        f"Unique phenotype strings: "
        f"{df['phenotypes_clean'].nunique()}\n"
    )

    f.write("\nPHENOTYPE INFORMATION CATEGORIES\n")
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
        f"No meaningful phenotype coverage: "
        f"{no_meaningful_count / total * 100:.2f}%\n"
    )


# ------------------------------------------------------------
# FINISH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("REPORT SAVED")
print("=" * 70)

print(OUTPUT_FILE)

print("\nStep 25 complete.")