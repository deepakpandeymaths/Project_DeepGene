# ============================================================
# 28. PHENOTYPE PARSING / REPRESENTATION ANALYSIS
# ============================================================
#
# Purpose:
# Analyze phenotype strings containing semicolons and determine
# whether they represent multiple phenotype terms bundled inside
# a single field.
#
# IMPORTANT:
# This step DOES NOT modify the phenotype dataset.
#
# We are only measuring:
#   - semicolon-containing entries
#   - number of sub-terms created by semicolon splitting
#   - overlap with existing "|" separated terms
#   - duplicate terms after considering both separators
#   - common compound structures
#   - examples requiring careful review
#
# We will use this information before deciding whether a future
# normalization step is scientifically justified.
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
    RESULTS_DIR,
    "28_phenotype_parsing_analysis.txt"
)


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 28")
print("PHENOTYPE PARSING / REPRESENTATION ANALYSIS")
print("=" * 70)


# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

print("\nLoading phenotype profile dataset...")

df = pd.read_csv(INPUT_FILE)

print(f"Rows loaded: {len(df)}")
print(f"Columns loaded: {len(df.columns)}")


# ------------------------------------------------------------
# BASIC CLEANING FOR ANALYSIS ONLY
# ------------------------------------------------------------

df["meaningful_phenotype_terms"] = (
    df["meaningful_phenotype_terms"]
    .fillna("")
    .astype(str)
    .str.strip()
)


# ------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------

def split_existing_terms(value):
    """
    Split the current representation using '|'.
    This represents the structure created in Step 26.
    """

    if not value:
        return []

    return [
        term.strip()
        for term in value.split("|")
        if term.strip()
    ]


def split_semicolon_terms(value):
    """
    Split terms using ';'.

    This is ONLY for analysis.
    The original dataset is not changed.
    """

    if not value:
        return []

    return [
        term.strip()
        for term in value.split(";")
        if term.strip()
    ]


def split_both_separators(value):
    """
    Treat both '|' and ';' as separators.

    This lets us estimate how many individual phenotype
    terms are hidden inside compound entries.
    """

    if not value:
        return []

    value = value.replace(";", "|")

    return [
        term.strip()
        for term in value.split("|")
        if term.strip()
    ]


# ------------------------------------------------------------
# 1. IDENTIFY SEMICOLON-CONTAINING VARIANTS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. SEMICOLON-CONTAINING VARIANTS")
print("=" * 70)

semicolon_mask = (
    df["meaningful_phenotype_terms"]
    .str.contains(";", regex=False)
)

semicolon_df = df[
    semicolon_mask
].copy()

print(
    f"Variants containing semicolon: "
    f"{len(semicolon_df)}"
)

print(
    f"Percentage of all variants: "
    f"{len(semicolon_df) / len(df) * 100:.2f}%"
)


# ------------------------------------------------------------
# 2. COUNT "|" TERMS
# ------------------------------------------------------------

df["pipe_term_count"] = (
    df["meaningful_phenotype_terms"]
    .apply(
        lambda x: len(
            split_existing_terms(x)
        )
    )
)


# ------------------------------------------------------------
# 3. COUNT SEMICOLON-SPLIT TERMS
# ------------------------------------------------------------

df["semicolon_term_count"] = (
    df["meaningful_phenotype_terms"]
    .apply(
        lambda x: len(
            split_semicolon_terms(x)
        )
    )
)


# ------------------------------------------------------------
# 4. COUNT TERMS USING BOTH SEPARATORS
# ------------------------------------------------------------

df["combined_term_count"] = (
    df["meaningful_phenotype_terms"]
    .apply(
        lambda x: len(
            split_both_separators(x)
        )
    )
)


print("\n" + "=" * 70)
print("2. TERM COUNTS AFTER INTERPRETING SEPARATORS")
print("=" * 70)

print(
    "\nCurrent '|' representation:"
)

print(
    df["pipe_term_count"]
    .value_counts()
    .sort_index()
    .to_string()
)

print(
    "\nCombined '|' + ';' representation:"
)

print(
    df["combined_term_count"]
    .value_counts()
    .sort_index()
    .to_string()
)


# ------------------------------------------------------------
# 5. HOW MANY TERMS ARE HIDDEN BY SEMICOLONS?
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. HIDDEN TERMS INSIDE COMPOUND ENTRIES")
print("=" * 70)

df["additional_terms_from_semicolon"] = (
    df["combined_term_count"]
    - df["pipe_term_count"]
)

affected_variants = (
    df["additional_terms_from_semicolon"] > 0
).sum()

total_additional_terms = (
    df["additional_terms_from_semicolon"]
    .sum()
)

print(
    f"Variants gaining additional terms after "
    f"semicolon splitting: "
    f"{affected_variants}"
)

print(
    f"Total additional term occurrences exposed: "
    f"{total_additional_terms}"
)


# ------------------------------------------------------------
# 6. DISTRIBUTION OF ADDITIONAL TERMS
# ------------------------------------------------------------

print("\nAdditional terms revealed per variant:")

print(
    df.loc[
        df["additional_terms_from_semicolon"] > 0,
        "additional_terms_from_semicolon"
    ]
    .value_counts()
    .sort_index()
    .to_string()
)


# ------------------------------------------------------------
# 7. CHECK FOR DUPLICATES AFTER SPLITTING BOTH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. DUPLICATES AFTER BOTH-SEPARATOR PARSING")
print("=" * 70)


def duplicate_after_combined_split(value):

    terms = split_both_separators(value)

    normalized = [
        term.lower()
        for term in terms
    ]

    return len(normalized) != len(
        set(normalized)
    )


combined_duplicate_mask = (
    df["meaningful_phenotype_terms"]
    .apply(
        duplicate_after_combined_split
    )
)

combined_duplicate_variants = df[
    combined_duplicate_mask
]

print(
    f"Variants containing duplicate phenotype "
    f"terms after combined parsing: "
    f"{len(combined_duplicate_variants)}"
)


if len(combined_duplicate_variants) > 0:

    print("\nExamples:")

    print(
        combined_duplicate_variants[
            [
                "variant_id",
                "meaningful_phenotype_terms"
            ]
        ]
        .head(20)
        .to_string(index=False)
    )


# ------------------------------------------------------------
# 8. UNIQUE TERMS BEFORE AND AFTER PARSING
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. UNIQUE TERM COMPARISON")
print("=" * 70)

pipe_terms = []

combined_terms = []

for value in df[
    "meaningful_phenotype_terms"
]:

    pipe_terms.extend(
        split_existing_terms(value)
    )

    combined_terms.extend(
        split_both_separators(value)
    )


pipe_series = pd.Series(
    pipe_terms,
    dtype="object"
)

combined_series = pd.Series(
    combined_terms,
    dtype="object"
)

print(
    f"Unique terms using '|': "
    f"{pipe_series.str.lower().nunique()}"
)

print(
    f"Unique terms using '|' + ';': "
    f"{combined_series.str.lower().nunique()}"
)

print(
    f"Total occurrences using '|': "
    f"{len(pipe_series)}"
)

print(
    f"Total occurrences using '|' + ';': "
    f"{len(combined_series)}"
)


# ------------------------------------------------------------
# 9. MOST COMMON TERMS AFTER COMBINED PARSING
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("6. TOP TERMS AFTER COMBINED PARSING")
print("=" * 70)

print(
    combined_series
    .value_counts()
    .head(30)
    .to_string()
)


# ------------------------------------------------------------
# 10. COMPOUND ENTRY PATTERNS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("7. COMMON COMPOUND ENTRIES")
print("=" * 70)


compound_entries = (
    semicolon_df[
        "meaningful_phenotype_terms"
    ]
    .value_counts()
)

print(
    f"Unique semicolon-containing entries: "
    f"{len(compound_entries)}"
)

print("\nTop 20 compound entries:")

print(
    compound_entries
    .head(20)
    .to_string()
)


# ------------------------------------------------------------
# 11. EXAMPLE PARSING
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("8. EXAMPLE PARSING")
print("=" * 70)

examples = semicolon_df.head(10)

for _, row in examples.iterrows():

    original = row[
        "meaningful_phenotype_terms"
    ]

    parsed = split_both_separators(
        original
    )

    print("\nVariant:", row["variant_id"])

    print("Original:")
    print(original)

    print("\nParsed terms:")

    for number, term in enumerate(
        parsed,
        start=1
    ):

        print(
            f"  {number}. {term}"
        )


# ------------------------------------------------------------
# 12. CHECK WHETHER SEMICOLONS CREATE
#     TRUE INDIVIDUAL PHENOTYPE TERMS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("9. STRUCTURAL INTERPRETATION")
print("=" * 70)

print(
    """
The analysis treats both '|' and ';' as possible
phenotype separators ONLY for auditing.

This does not establish that semicolon splitting is
scientifically correct.

The important question is whether the text appearing
between semicolons behaves like individual phenotype
terms rather than a single biological phenotype label.
"""
)

semicolon_subterms = []

for value in semicolon_df[
    "meaningful_phenotype_terms"
]:

    parts = [
        term.strip()
        for term in value.split(";")
        if term.strip()
    ]

    semicolon_subterms.extend(
        parts
    )

semicolon_subterm_series = pd.Series(
    semicolon_subterms,
    dtype="object"
)

print(
    f"Total semicolon-separated components: "
    f"{len(semicolon_subterm_series)}"
)

print(
    f"Unique semicolon-separated components: "
    f"{semicolon_subterm_series.str.lower().nunique()}"
)

print("\nMost common semicolon-separated components:")

print(
    semicolon_subterm_series
    .value_counts()
    .head(30)
    .to_string()
)


# ------------------------------------------------------------
# 13. FINAL ASSESSMENT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("10. FINAL ASSESSMENT")
print("=" * 70)

print(
    f"Total variants: {len(df)}"
)

print(
    f"Variants containing semicolons: "
    f"{len(semicolon_df)}"
)

print(
    f"Variants affected by semicolon parsing: "
    f"{affected_variants}"
)

print(
    f"Additional term occurrences exposed: "
    f"{total_additional_terms}"
)

print(
    f"Duplicates after combined parsing: "
    f"{len(combined_duplicate_variants)}"
)

if len(combined_duplicate_variants) == 0:

    print(
        "\nNo duplicate phenotype terms were introduced "
        "by treating ';' as an additional separator."
    )

else:

    print(
        "\nSome duplicate phenotype terms appear after "
        "combined parsing and require review."
    )


print(
    """
\nIMPORTANT:
No phenotype values have been modified by Step 28.

The purpose of this step is to determine whether a future
normalization step can safely represent compound phenotype
entries without losing information.
"""
)


# ------------------------------------------------------------
# 14. SAVE REPORT
# ------------------------------------------------------------

with open(
    OUTPUT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "DEEPGENE STEP 28\n"
        "PHENOTYPE PARSING / REPRESENTATION ANALYSIS\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Total variants: {len(df)}\n"
    )

    f.write(
        f"Variants containing semicolons: "
        f"{len(semicolon_df)}\n"
    )

    f.write(
        f"Percentage containing semicolons: "
        f"{len(semicolon_df) / len(df) * 100:.2f}%\n"
    )

    f.write(
        f"Variants affected by semicolon parsing: "
        f"{affected_variants}\n"
    )

    f.write(
        f"Additional term occurrences exposed: "
        f"{total_additional_terms}\n"
    )

    f.write(
        f"Unique terms using '|': "
        f"{pipe_series.str.lower().nunique()}\n"
    )

    f.write(
        f"Unique terms using '|' + ';': "
        f"{combined_series.str.lower().nunique()}\n"
    )

    f.write(
        f"Duplicates after combined parsing: "
        f"{len(combined_duplicate_variants)}\n"
    )

    f.write(
        f"Unique semicolon-separated components: "
        f"{semicolon_subterm_series.str.lower().nunique()}\n"
    )

    f.write("\n")

    f.write(
        "IMPORTANT: No phenotype values were modified.\n"
    )

    f.write(
        "STATUS: AUDIT COMPLETE\n"
    )


# ------------------------------------------------------------
# FINISH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("REPORT SAVED")
print("=" * 70)

print(OUTPUT_FILE)

print("\nStep 28 complete.")