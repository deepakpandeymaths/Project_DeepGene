# ============================================================
# 37. FUNCTIONAL VARIANT MATCHING AUDIT
# ============================================================
#
# Purpose:
#   Determine whether the 58 functional variants in the
#   Brunklaus 2020 dataset can be matched to DeepGene V1.
#
# Matching levels:
#   1. Exact text match
#   2. Normalized cDNA match
#   3. Normalized protein-change match
#   4. Protein-position + amino-acid change inspection
#
# IMPORTANT:
#   This step does NOT create ML labels.
#   This step does NOT train a model.
#
# ============================================================

import os
import re
import pandas as pd


# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)

FUNCTIONAL_FILE = os.path.join(
    PROJECT_ROOT,
    "data",
    "raw",
    "SCN1A_functional_Brunklaus_2020.xlsx"
)

V1_FILE = os.path.join(
    PROJECT_ROOT,
    "data",
    "processed",
    "deepgene_v1_final.csv"
)

RESULTS_DIR = os.path.join(
    PROJECT_ROOT,
    "data",
    "analysis_results"
)

OUTPUT_FILE = os.path.join(
    RESULTS_DIR,
    "37_functional_variant_matching.txt"
)


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

def section(title):
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def clean(value):
    if pd.isna(value):
        return ""

    return str(value).strip()


def normalize_text(value):
    """
    Conservative text normalization.

    Does NOT attempt to infer biological equivalence.
    """

    text = clean(value)

    if not text:
        return ""

    text = text.upper()
    text = re.sub(r"\s+", "", text)

    return text


def normalize_coding(value):
    """
    Normalize common cDNA HGVS strings.

    Examples:
        c.434T>C
        C.434T>C
        c. 434 T > C

    become:

        C.434T>C
    """

    text = normalize_text(value)

    if not text:
        return ""

    text = text.replace(">", ">")
    text = text.replace(" ", "")

    return text


def normalize_protein(value):
    """
    Conservative protein normalization.

    Examples:
        Met145Thr
        M145T

    are normalized only for formatting.

    We do NOT assume that all three-letter and
    one-letter amino-acid representations are
    equivalent unless explicitly converted.
    """

    text = normalize_text(value)

    if not text:
        return ""

    text = text.replace("P.", "")

    return text


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

os.makedirs(RESULTS_DIR, exist_ok=True)

section("DEEPGENE — SCN1A FUNCTIONAL VARIANT MATCHING AUDIT")


# ------------------------------------------------------------
# LOAD FUNCTIONAL DATA
# ------------------------------------------------------------

section("1. LOAD FUNCTIONAL DATASET")

if not os.path.exists(FUNCTIONAL_FILE):
    print("ERROR: Functional dataset not found.")
    print(FUNCTIONAL_FILE)
    raise SystemExit(1)

functional = pd.read_excel(
    FUNCTIONAL_FILE,
    sheet_name="Sheet1"
)

functional = functional.dropna(
    how="all"
).reset_index(drop=True)

print(f"Functional rows: {len(functional)}")


# ------------------------------------------------------------
# LOAD DEEPGENE V1
# ------------------------------------------------------------

section("2. LOAD DEEPGENE V1")

if not os.path.exists(V1_FILE):
    print("ERROR: DeepGene V1 not found.")
    print(V1_FILE)
    raise SystemExit(1)

v1 = pd.read_csv(V1_FILE)

print(f"DeepGene V1 rows: {len(v1)}")
print(
    f"DeepGene V1 unique variants: "
    f"{v1['variant_id'].nunique()}"
)


# ------------------------------------------------------------
# FUNCTIONAL VARIANT FIELDS
# ------------------------------------------------------------

section("3. FUNCTIONAL VARIANT FIELDS")

functional_coding = []

if "cDNA" in functional.columns:

    functional_coding = (
        functional["cDNA"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    functional_coding = functional_coding[
        functional_coding != ""
    ]

print(
    f"Functional cDNA descriptions: "
    f"{len(functional_coding)}"
)

print(
    f"Unique functional cDNA descriptions: "
    f"{functional_coding.nunique()}"
)


functional_protein = []

if "Variant" in functional.columns:

    functional_protein = (
        functional["Variant"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    functional_protein = functional_protein[
        functional_protein != ""
    ]

print(
    f"Functional protein descriptions: "
    f"{len(functional_protein)}"
)

print(
    f"Unique functional protein descriptions: "
    f"{functional_protein.nunique()}"
)


# ------------------------------------------------------------
# DEEPGENE HGVS FIELDS
# ------------------------------------------------------------

section("4. DEEPGENE HGVS FIELDS")

v1_hgvs = (
    v1["hgvs_name"]
    .dropna()
    .astype(str)
    .str.strip()
)

v1_hgvs = v1_hgvs[
    v1_hgvs != ""
]

print(
    f"DeepGene HGVS descriptions: "
    f"{len(v1_hgvs)}"
)

print(
    f"Unique DeepGene HGVS descriptions: "
    f"{v1_hgvs.nunique()}"
)


# ------------------------------------------------------------
# EXACT MATCH
# ------------------------------------------------------------

section("5. EXACT TEXT MATCH")

functional_coding_exact = set(functional_coding)
v1_hgvs_exact = set(v1_hgvs)

exact_matches = (
    functional_coding_exact
    & v1_hgvs_exact
)

print(
    f"Exact cDNA/HGVS matches: "
    f"{len(exact_matches)}"
)

if exact_matches:

    print()
    print("Examples:")

    for value in sorted(exact_matches)[:20]:
        print(f"- {value}")

else:

    print("No exact matches.")


# ------------------------------------------------------------
# NORMALIZED CODING MATCH
# ------------------------------------------------------------

section("6. NORMALIZED CODING MATCH")

functional_coding_normalized = {
    normalize_coding(value): value
    for value in functional_coding
    if normalize_coding(value)
}

v1_hgvs_normalized = {
    normalize_coding(value): value
    for value in v1_hgvs
    if normalize_coding(value)
}

normalized_coding_matches = (
    set(functional_coding_normalized.keys())
    &
    set(v1_hgvs_normalized.keys())
)

print(
    f"Normalized coding matches: "
    f"{len(normalized_coding_matches)}"
)

if normalized_coding_matches:

    print()
    print("Examples:")

    for normalized in sorted(
        normalized_coding_matches
    )[:20]:

        print(
            f"- Functional: "
            f"{functional_coding_normalized[normalized]}"
        )

        print(
            f"  DeepGene: "
            f"{v1_hgvs_normalized[normalized]}"
        )

else:

    print("No normalized coding matches.")


# ------------------------------------------------------------
# PROTEIN REPRESENTATION ANALYSIS
# ------------------------------------------------------------

section("7. PROTEIN REPRESENTATION ANALYSIS")

print()
print(
    "Functional dataset uses protein-level descriptions "
    "such as Met145Thr."
)

print()
print(
    "DeepGene V1 primarily stores HGVS descriptions."
)

print()
print(
    "Therefore protein-level matching cannot be considered "
    "a validated biological match yet."
)


# ------------------------------------------------------------
# EXTRACT SIMPLE PROTEIN POSITION
# ------------------------------------------------------------

def extract_protein_position(value):

    text = clean(value)

    if not text:
        return None

    match = re.search(
        r"(?:[A-Za-z]{1,3})?(\d+)",
        text
    )

    if match:
        return int(match.group(1))

    return None


functional_positions = []

if "Variant" in functional.columns:

    for value in functional["Variant"]:

        position = extract_protein_position(value)

        if position is not None:
            functional_positions.append(
                position
            )

print(
    f"Functional variants with detectable "
    f"protein position: {len(functional_positions)}"
)

print(
    f"Unique protein positions: "
    f"{len(set(functional_positions))}"
)


# ------------------------------------------------------------
# BUILD MATCHING TABLE
# ------------------------------------------------------------

section("8. FUNCTIONAL VARIANT MATCH TABLE")

records = []

for _, row in functional.iterrows():

    coding = clean(
        row.get("cDNA", "")
    )

    protein = clean(
        row.get("Variant", "")
    )

    function = clean(
        row.get("Function", "")
    )

    position = extract_protein_position(
        protein
    )

    records.append(
        {
            "functional_cDNA": coding,
            "functional_protein": protein,
            "function": function,
            "protein_position": position,
            "normalized_cDNA":
                normalize_coding(coding),
            "exact_cDNA_match":
                coding in v1_hgvs_exact,
            "normalized_cDNA_match":
                normalize_coding(coding)
                in v1_hgvs_normalized
        }
    )

match_df = pd.DataFrame(records)

print(
    f"Matching table rows: {len(match_df)}"
)

print()
print(
    "Exact cDNA matches: "
    f"{match_df['exact_cDNA_match'].sum()}"
)

print(
    "Normalized cDNA matches: "
    f"{match_df['normalized_cDNA_match'].sum()}"
)


# ------------------------------------------------------------
# FUNCTIONAL CATEGORY BY MATCH STATUS
# ------------------------------------------------------------

section("9. FUNCTION CATEGORY BY MATCH STATUS")

if "Function" in functional.columns:

    category_table = (
        match_df
        .groupby(
            [
                "function",
                "normalized_cDNA_match"
            ]
        )
        .size()
        .reset_index(
            name="count"
        )
    )

    for _, row in category_table.iterrows():

        status = (
            "MATCH"
            if row["normalized_cDNA_match"]
            else "NO MATCH"
        )

        print(
            f"- {row['function']} | "
            f"{status}: {row['count']}"
        )


# ------------------------------------------------------------
# IMPORTANT BIOLOGICAL CHECK
# ------------------------------------------------------------

section("10. BIOLOGICAL MATCHING STATUS")

print()
print(
    "A protein-level description such as Met145Thr "
    "does not by itself prove identity with a DeepGene "
    "HGVS variant."
)

print()
print(
    "Transcript, reference sequence, amino-acid "
    "representation and genomic representation must "
    "be reconciled before constructing an ML dataset."
)

print()
print(
    "Therefore protein-position similarity is treated "
    "as an inspection signal only."
)


# ------------------------------------------------------------
# TARGET STATUS
# ------------------------------------------------------------

section("11. TARGET STATUS")

print()
print(
    "Independent functional measurements: YES"
)

print(
    "Validated biological overlap: NOT YET"
)

print(
    "Independent ML target: NOT YET DEFINED"
)

print(
    "Pathogenic/Benign labels created: NO"
)

print(
    "ML training started: NO"
)


# ------------------------------------------------------------
# SAVE MATCH TABLE
# ------------------------------------------------------------

MATCH_FILE = os.path.join(
    PROJECT_ROOT,
    "data",
    "processed",
    "deepgene_functional_variant_matching_v1.csv"
)

match_df.to_csv(
    MATCH_FILE,
    index=False
)

print()
print(
    "Matching table saved:"
)

print(MATCH_FILE)


# ------------------------------------------------------------
# SAVE REPORT
# ------------------------------------------------------------

section("12. SAVE REPORT")

with open(
    OUTPUT_FILE,
    "w",
    encoding="utf-8"
) as report:

    report.write(
        "DEEPGENE — SCN1A FUNCTIONAL VARIANT MATCHING AUDIT\n"
    )

    report.write("=" * 70 + "\n\n")

    report.write(
        f"Functional rows: {len(functional)}\n"
    )

    report.write(
        f"DeepGene V1 rows: {len(v1)}\n"
    )

    report.write(
        f"Exact coding matches: {len(exact_matches)}\n"
    )

    report.write(
        f"Normalized coding matches: "
        f"{len(normalized_coding_matches)}\n"
    )

    report.write("\n")
    report.write(
        "TARGET STATUS\n"
    )

    report.write(
        "Validated biological overlap: NOT YET\n"
    )

    report.write(
        "Independent ML target: NOT YET DEFINED\n"
    )

    report.write(
        "Pathogenic/Benign labels: NOT CREATED\n"
    )

    report.write(
        "ML training: NOT STARTED\n"
    )


# ------------------------------------------------------------
# FINAL STATUS
# ------------------------------------------------------------

section("13. STEP 37 STATUS")

print()
print(
    "Functional variants inspected: "
    f"{len(functional)}"
)

print(
    "Exact cDNA matches: "
    f"{len(exact_matches)}"
)

print(
    "Normalized cDNA matches: "
    f"{len(normalized_coding_matches)}"
)

print(
    "Protein-level inspection: YES"
)

print(
    "Validated biological overlap: NOT YET"
)

print(
    "Independent ML target: NOT YET"

)

print()
print("🚨 ML HAS NOT STARTED.")

print()
print("STEP 37 COMPLETE")