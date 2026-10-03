# ============================================================
# 27. PHENOTYPE TERM NORMALIZATION AUDIT
# ============================================================
#
# Purpose:
# Audit the structured phenotype terms created in Step 26.
#
# We specifically check:
#   - terms containing semicolons
#   - duplicate phenotype terms within a variant
#   - repeated terms across variants
#   - whitespace inconsistencies
#   - possible compound phenotype entries
#
# IMPORTANT:
# This step does NOT modify phenotype data.
# It only identifies normalization issues.
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
    "27_phenotype_term_normalization_audit.txt"
)


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 27")
print("PHENOTYPE TERM NORMALIZATION AUDIT")
print("=" * 70)


# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

print("\nLoading phenotype profile dataset...")

df = pd.read_csv(INPUT_FILE)

print(f"Rows loaded: {len(df)}")
print(f"Columns loaded: {len(df.columns)}")


# ------------------------------------------------------------
# NORMALIZE BASIC TEXT
# ------------------------------------------------------------

df["meaningful_phenotype_terms"] = (
    df["meaningful_phenotype_terms"]
    .fillna("")
    .astype(str)
    .str.strip()
)


# ------------------------------------------------------------
# 1. SEMICOLON-CONTAINING TERMS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. TERMS CONTAINING SEMICOLONS")
print("=" * 70)

semicolon_mask = (
    df["meaningful_phenotype_terms"]
    .str.contains(";", regex=False)
)

semicolon_variants = df[
    semicolon_mask
]

print(
    f"Variants containing semicolon characters: "
    f"{len(semicolon_variants)}"
)

if len(semicolon_variants) > 0:

    print("\nExamples:")

    print(
        semicolon_variants[
            [
                "variant_id",
                "meaningful_phenotype_terms"
            ]
        ]
        .head(20)
        .to_string(index=False)
    )


# ------------------------------------------------------------
# 2. DUPLICATE TERMS WITHIN EACH VARIANT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("2. DUPLICATE TERMS WITHIN VARIANT")
print("=" * 70)


def has_duplicate_terms(value):

    if not value:
        return False

    terms = [
        term.strip()
        for term in value.split("|")
        if term.strip()
    ]

    normalized = [
        term.lower()
        for term in terms
    ]

    return len(normalized) != len(set(normalized))


duplicate_term_mask = (
    df["meaningful_phenotype_terms"]
    .apply(has_duplicate_terms)
)

duplicate_term_variants = df[
    duplicate_term_mask
]

print(
    f"Variants containing repeated phenotype terms: "
    f"{len(duplicate_term_variants)}"
)

if len(duplicate_term_variants) > 0:

    print("\nExamples:")

    print(
        duplicate_term_variants[
            [
                "variant_id",
                "meaningful_phenotype_terms"
            ]
        ]
        .head(20)
        .to_string(index=False)
    )


# ------------------------------------------------------------
# 3. WHITESPACE INCONSISTENCY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. WHITESPACE CHECK")
print("=" * 70)


def has_whitespace_issue(value):

    if not value:
        return False

    terms = value.split("|")

    for term in terms:

        if term != term.strip():
            return True

    return False


whitespace_mask = (
    df["meaningful_phenotype_terms"]
    .apply(has_whitespace_issue)
)

whitespace_variants = df[
    whitespace_mask
]

print(
    f"Variants with leading/trailing whitespace "
    f"inside phenotype terms: "
    f"{len(whitespace_variants)}"
)


# ------------------------------------------------------------
# 4. EMPTY TERMS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. EMPTY TERM CHECK")
print("=" * 70)


def count_empty_terms(value):

    if not value:
        return 0

    parts = value.split("|")

    return sum(
        1
        for part in parts
        if not part.strip()
    )


df["empty_term_count"] = (
    df["meaningful_phenotype_terms"]
    .apply(count_empty_terms)
)

empty_term_variants = (
    df["empty_term_count"] > 0
).sum()

print(
    f"Variants containing empty phenotype terms: "
    f"{empty_term_variants}"
)


# ------------------------------------------------------------
# 5. TERM FREQUENCY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. PHENOTYPE TERM FREQUENCY")
print("=" * 70)

all_terms = []

for value in df[
    "meaningful_phenotype_terms"
]:

    if not value:
        continue

    terms = [
        term.strip()
        for term in value.split("|")
        if term.strip()
    ]

    all_terms.extend(terms)


term_series = pd.Series(
    all_terms,
    dtype="object"
)

print(
    f"Total phenotype term occurrences: "
    f"{len(term_series)}"
)

print(
    f"Unique phenotype terms: "
    f"{term_series.str.lower().nunique()}"
)

print("\nTop 30 phenotype terms:")

print(
    term_series
    .value_counts()
    .head(30)
    .to_string()
)


# ------------------------------------------------------------
# 6. CASE-INSENSITIVE DUPLICATE TERM CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("6. CASE-INSENSITIVE TERM DUPLICATION")
print("=" * 70)

normalized_term_series = (
    term_series
    .str.strip()
    .str.lower()
)

case_insensitive_duplicates = (
    normalized_term_series
    .value_counts()
)

duplicate_normalized_terms = (
    case_insensitive_duplicates[
        case_insensitive_duplicates > 1
    ]
)

print(
    f"Unique normalized phenotype terms: "
    f"{len(case_insensitive_duplicates)}"
)

print(
    f"Normalized terms appearing more than once: "
    f"{len(duplicate_normalized_terms)}"
)


# ------------------------------------------------------------
# 7. COMPOUND ENTRY CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("7. POSSIBLE COMPOUND PHENOTYPE ENTRIES")
print("=" * 70)

compound_patterns = [
    ";",
    ","
]

compound_candidates = set()

for value in term_series:

    for pattern in compound_patterns:

        if pattern in value:

            compound_candidates.add(value)

            break


print(
    f"Unique phenotype terms containing comma or "
    f"semicolon: "
    f"{len(compound_candidates)}"
)

print("\nExamples:")

for index, term in enumerate(
    list(compound_candidates)[:20],
    start=1
):

    print(
        f"{index}. {term}"
    )


# ------------------------------------------------------------
# 8. NORMALIZATION ASSESSMENT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("8. NORMALIZATION ASSESSMENT")
print("=" * 70)

issues = 0

if len(semicolon_variants) > 0:
    issues += 1

if len(duplicate_term_variants) > 0:
    issues += 1

if len(whitespace_variants) > 0:
    issues += 1

if empty_term_variants > 0:
    issues += 1


if issues == 0:

    print(
        "No structural phenotype normalization issues "
        "were detected."
    )

else:

    print(
        f"{issues} type(s) of phenotype representation "
        f"issue(s) were detected."
    )

    print(
        "\nThese issues should be reviewed before applying "
        "aggressive phenotype normalization."
    )


# ------------------------------------------------------------
# 9. SAVE REPORT
# ------------------------------------------------------------

with open(
    OUTPUT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "DEEPGENE STEP 27\n"
        "PHENOTYPE TERM NORMALIZATION AUDIT\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Rows checked: {len(df)}\n"
    )

    f.write(
        f"Variants containing semicolons: "
        f"{len(semicolon_variants)}\n"
    )

    f.write(
        f"Variants with duplicate terms: "
        f"{len(duplicate_term_variants)}\n"
    )

    f.write(
        f"Variants with whitespace issues: "
        f"{len(whitespace_variants)}\n"
    )

    f.write(
        f"Variants with empty terms: "
        f"{empty_term_variants}\n"
    )

    f.write(
        f"Total phenotype term occurrences: "
        f"{len(term_series)}\n"
    )

    f.write(
        f"Unique phenotype terms: "
        f"{term_series.str.lower().nunique()}\n"
    )

    f.write(
        f"Compound phenotype candidates: "
        f"{len(compound_candidates)}\n"
    )

    f.write("\n")

    if issues == 0:

        f.write(
            "NORMALIZATION STATUS: "
            "NO STRUCTURAL ISSUES DETECTED\n"
        )

    else:

        f.write(
            "NORMALIZATION STATUS: "
            "REVIEW REQUIRED\n"
        )


# ------------------------------------------------------------
# FINISH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("REPORT SAVED")
print("=" * 70)

print(OUTPUT_FILE)

print("\nStep 27 complete.")