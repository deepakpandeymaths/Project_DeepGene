# ============================================================
# 38–39. FUNCTIONAL OVERLAP + TARGET ASSESSMENT
# ============================================================
#
# DeepGene — SCN1A
#
# PURPOSE
# -------
# Combine:
#
#   Step 38 -> deeper functional variant matching
#   Step 39 -> independent target/benchmark assessment
#
# This script investigates whether the Brunklaus 2020
# functional dataset can provide a defensible independent
# target for DeepGene.
#
# IMPORTANT:
#   - No pathogenic/benign labels are created.
#   - No ClinVar labels are copied into the target.
#   - No ML model is trained.
#   - DeepGene V1 is not modified.
#
# ============================================================

import os
import re
import pandas as pd


# ============================================================
# 1. PATHS
# ============================================================

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

PROCESSED_DIR = os.path.join(
    PROJECT_ROOT,
    "data",
    "processed"
)

REPORT_FILE = os.path.join(
    RESULTS_DIR,
    "38_39_functional_overlap_and_target_assessment.txt"
)

MATCH_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_functional_overlap_v1.csv"
)


# ============================================================
# 2. HELPERS
# ============================================================

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

    text = clean(value)

    if not text:
        return ""

    text = text.upper()

    # Remove whitespace
    text = re.sub(r"\s+", "", text)

    return text


# ------------------------------------------------------------
# Amino-acid conversion
# ------------------------------------------------------------

AA_THREE_TO_ONE = {

    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V"
}


def protein_to_one_letter(value):

    text = normalize_text(value)

    if not text:
        return ""

    # Remove HGVS protein prefix
    text = text.replace("P.", "")

    # Three-letter amino-acid notation
    match = re.fullmatch(
        r"([A-Z]{3})(\d+)([A-Z]{3})",
        text
    )

    if match:

        ref = match.group(1)
        pos = match.group(2)
        alt = match.group(3)

        if ref in AA_THREE_TO_ONE and alt in AA_THREE_TO_ONE:

            return (
                AA_THREE_TO_ONE[ref]
                + pos
                + AA_THREE_TO_ONE[alt]
            )

    # Already one-letter notation
    match = re.fullmatch(
        r"([A-Z])(\d+)([A-Z])",
        text
    )

    if match:

        return text

    return ""


def extract_protein_position(value):

    text = clean(value)

    if not text:
        return None

    match = re.search(
        r"[A-Za-z]{1,3}(\d+)",
        text
    )

    if match:
        return int(match.group(1))

    return None


def extract_protein_change(value):

    text = normalize_text(value)

    if not text:
        return None

    text = text.replace("P.", "")

    # Three-letter form
    match = re.fullmatch(
        r"([A-Z]{3})(\d+)([A-Z]{3})",
        text
    )

    if match:

        ref = match.group(1)
        pos = match.group(2)
        alt = match.group(3)

        if (
            ref in AA_THREE_TO_ONE
            and alt in AA_THREE_TO_ONE
        ):

            return (
                AA_THREE_TO_ONE[ref]
                + pos
                + AA_THREE_TO_ONE[alt]
            )

    # One-letter form
    match = re.fullmatch(
        r"([A-Z])(\d+)([A-Z])",
        text
    )

    if match:

        return text

    return None


# ============================================================
# 3. START
# ============================================================

os.makedirs(
    RESULTS_DIR,
    exist_ok=True
)

os.makedirs(
    PROCESSED_DIR,
    exist_ok=True
)

section(
    "DEEPGENE — STEPS 38–39 FUNCTIONAL OVERLAP "
    "AND TARGET ASSESSMENT"
)


# ============================================================
# 4. LOAD FUNCTIONAL DATA
# ============================================================

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
).reset_index(
    drop=True
)

print(
    f"Functional dataset rows: "
    f"{len(functional)}"
)

print(
    f"Functional dataset columns: "
    f"{len(functional.columns)}"
)


# ============================================================
# 5. LOAD DEEPGENE V1
# ============================================================

section("2. LOAD DEEPGENE V1")

if not os.path.exists(V1_FILE):

    print("ERROR: DeepGene V1 not found.")

    print(V1_FILE)

    raise SystemExit(1)


v1 = pd.read_csv(
    V1_FILE
)

print(
    f"DeepGene V1 rows: "
    f"{len(v1)}"
)

print(
    f"DeepGene V1 unique variants: "
    f"{v1['variant_id'].nunique()}"
)


# ============================================================
# 6. FUNCTIONAL VARIANTS
# ============================================================

section("3. FUNCTIONAL VARIANT INVENTORY")

functional = functional[
    functional["Variant"].notna()
].copy()

functional = functional[
    functional["Variant"].astype(str).str.strip() != ""
].copy()

print(
    f"Functional variants with protein description: "
    f"{len(functional)}"
)

print(
    f"Unique functional protein descriptions: "
    f"{functional['Variant'].nunique()}"
)


# ============================================================
# 7. FUNCTIONAL PROTEIN NORMALIZATION
# ============================================================

section("4. NORMALIZE FUNCTIONAL PROTEIN VARIANTS")

functional["protein_original"] = (
    functional["Variant"]
    .astype(str)
    .str.strip()
)

functional["protein_normalized"] = (
    functional["Variant"]
    .apply(protein_to_one_letter)
)

functional["protein_position"] = (
    functional["Variant"]
    .apply(extract_protein_position)
)

print(
    "Protein variants converted to one-letter notation: "
    f"{(functional['protein_normalized'] != '').sum()}"
)

print(
    "Protein positions detected: "
    f"{functional['protein_position'].notna().sum()}"
)

print()
print("Examples:")

shown = 0

for _, row in functional.iterrows():

    if row["protein_normalized"]:

        print(
            f"- {row['protein_original']} "
            f"→ {row['protein_normalized']}"
        )

        shown += 1

        if shown >= 10:
            break


# ============================================================
# 8. EXTRACT PROTEIN REPRESENTATIONS FROM DEEPGENE HGVS
# ============================================================

section(
    "5. EXTRACT PROTEIN REPRESENTATIONS FROM DEEPGENE"
)

v1["hgvs_clean"] = (
    v1["hgvs_name"]
    .fillna("")
    .astype(str)
    .str.strip()
)

# Search for p.X123Y-like representations
# and common three-letter representations.

def extract_hgvs_protein(value):

    text = clean(value)

    if not text:
        return ""

    # Search for p.A123B
    match = re.search(
        r"p\.([A-Za-z])(\d+)([A-Za-z])",
        text,
        flags=re.IGNORECASE
    )

    if match:

        return (
            match.group(1).upper()
            + match.group(2)
            + match.group(3).upper()
        )

    # Search for p.Met145Thr
    match = re.search(
        r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})",
        text,
        flags=re.IGNORECASE
    )

    if match:

        ref = match.group(1).upper()
        pos = match.group(2)
        alt = match.group(3).upper()

        if (
            ref in AA_THREE_TO_ONE
            and alt in AA_THREE_TO_ONE
        ):

            return (
                AA_THREE_TO_ONE[ref]
                + pos
                + AA_THREE_TO_ONE[alt]
            )

    return ""


v1["protein_normalized"] = (
    v1["hgvs_clean"]
    .apply(extract_hgvs_protein)
)

protein_available = (
    v1["protein_normalized"] != ""
)

print(
    "DeepGene variants with detectable protein HGVS: "
    f"{protein_available.sum()}"
)


# ============================================================
# 9. PROTEIN-LEVEL MATCHING
# ============================================================

section("6. PROTEIN-LEVEL MATCHING")

v1_protein_set = set(
    v1.loc[
        protein_available,
        "protein_normalized"
    ]
)

functional_protein_set = set(
    functional.loc[
        functional["protein_normalized"] != "",
        "protein_normalized"
    ]
)

protein_matches = (
    functional_protein_set
    &
    v1_protein_set
)

print(
    f"Unique functional protein variants: "
    f"{len(functional_protein_set)}"
)

print(
    f"Unique DeepGene protein variants: "
    f"{len(v1_protein_set)}"
)

print(
    f"Protein-level exact normalized matches: "
    f"{len(protein_matches)}"
)


# ============================================================
# 10. MATCH TABLE
# ============================================================

section("7. BUILD FUNCTIONAL OVERLAP TABLE")

v1_lookup = {}

for _, row in v1.iterrows():

    protein = row["protein_normalized"]

    if not protein:
        continue

    if protein not in v1_lookup:

        v1_lookup[protein] = []

    v1_lookup[protein].append(
        row["variant_id"]
    )


records = []

for _, row in functional.iterrows():

    protein = row["protein_normalized"]

    matched_ids = v1_lookup.get(
        protein,
        []
    )

    records.append(
        {
            "functional_cDNA":
                clean(row.get("cDNA", "")),

            "functional_protein":
                clean(row.get("Variant", "")),

            "normalized_protein":
                protein,

            "protein_position":
                row["protein_position"],

            "functional_class":
                clean(row.get("Function", "")),

            "phenotype":
                clean(row.get("Phenotype", "")),

            "inheritance":
                clean(row.get("Inheritance", "")),

            "deepgene_match":
                len(matched_ids) > 0,

            "deepgene_match_count":
                len(matched_ids),

            "deepgene_variant_ids":
                "|".join(map(str, matched_ids))
        }
    )


overlap = pd.DataFrame(
    records
)

print(
    f"Functional variants in overlap table: "
    f"{len(overlap)}"
)

print(
    "Functional variants with protein-level "
    "DeepGene matches: "
    f"{overlap['deepgene_match'].sum()}"
)


# ============================================================
# 11. MATCHED FUNCTIONAL CLASSES
# ============================================================

section(
    "8. FUNCTIONAL CLASSES AMONG MATCHED VARIANTS"
)

matched = overlap[
    overlap["deepgene_match"]
].copy()

if len(matched) == 0:

    print(
        "No validated protein-level matches found."
    )

else:

    print(
        f"Matched functional variants: "
        f"{len(matched)}"
    )

    print()

    print(
        matched[
            "functional_class"
        ].value_counts()
    )


# ============================================================
# 12. FUNCTIONAL DATASET CLASS DISTRIBUTION
# ============================================================

section(
    "9. COMPLETE FUNCTION CLASS DISTRIBUTION"
)

function_distribution = (
    overlap[
        "functional_class"
    ]
    .value_counts()
)

for category, count in function_distribution.items():

    print(
        f"- {category}: {count}"
    )


# ============================================================
# 13. TARGET CANDIDATE ASSESSMENT
# ============================================================

section(
    "10. INDEPENDENT TARGET CANDIDATE ASSESSMENT"
)

print()
print(
    "The functional dataset contains experimentally "
    "derived information."
)

print()
print(
    "The qualitative Function field contains categories "
    "such as Loss, Mixed, Gain and No effect."
)

print()
print(
    "However, these categories must be treated as "
    "experimental functional classifications, not "
    "automatically as pathogenic/benign labels."
)

print()
print(
    "Continuous measurements such as current density, "
    "activation shift and inactivation shift are also "
    "available for subsets of variants."
)

print()
print(
    "The target must therefore be defined from the "
    "experimental measurement itself."
)


# ============================================================
# 14. TARGET LEAKAGE CHECK
# ============================================================

section("11. TARGET LEAKAGE CHECK")

print()
print(
    "The following fields from the functional workbook "
    "must NOT automatically become model inputs if they "
    "directly encode the target:"
)

print(
    "- Function"
)

print(
    "- ACMG score with functional data"
)

print(
    "- Any post-hoc interpretation derived from "
    "the experimental result"
)

print()
print(
    "The target should be kept separate from the "
    "predictor feature set."
)


# ============================================================
# 15. CLINVAR DEPENDENCY CHECK
# ============================================================

section("12. CLINVAR DEPENDENCY CHECK")

print()
print(
    "DeepGene V1 contains ClinVar-derived clinical "
    "significance and evidence features."
)

print()
print(
    "These must be audited before ML because a model "
    "trained to predict functional measurements should "
    "not simply reproduce ClinVar classification."
)

print()
print(
    "Therefore the functional measurement remains "
    "conceptually independent, but the feature set "
    "must still be filtered for leakage."
)


# ============================================================
# 16. USABLE SAMPLE SIZE
# ============================================================

section("13. USABLE SAMPLE SIZE")

print(
    f"Total functional variants: "
    f"{len(overlap)}"
)

print(
    f"Protein-level DeepGene matches: "
    f"{len(matched)}"
)

if len(overlap) > 0:

    coverage = (
        100
        * len(matched)
        / len(overlap)
    )

else:

    coverage = 0

print(
    f"Functional-to-DeepGene overlap coverage: "
    f"{coverage:.2f}%"
)


# ============================================================
# 17. TARGET READINESS
# ============================================================

section("14. TARGET READINESS DECISION")

if len(matched) == 0:

    decision = (
        "NOT READY — no validated biological "
        "overlap established"
    )

elif len(matched) < 10:

    decision = (
        "NOT READY — overlap too small for "
        "a defensible benchmark without further "
        "independent data"
    )

else:

    decision = (
        "PROMISING — overlap exists; further "
        "measurement harmonization and leakage "
        "audit required before ML"
    )

print()
print(
    f"Decision: {decision}"
)


# ============================================================
# 18. FINAL SCIENTIFIC STATUS
# ============================================================

section(
    "15. FINAL SCIENTIFIC STATUS"
)

print()
print(
    "Functional dataset available: YES"
)

print(
    "Experimental measurements available: YES"
)

print(
    "Protein-level normalization attempted: YES"
)

print(
    "Biological overlap investigated: YES"
)

print(
    f"Validated protein-level overlap: "
    f"{len(matched)}"
)

print(
    "Independent pathogenic/benign labels: NO"
)

print(
    "Final ML target: NOT YET ACCEPTED"
)

print(
    "ML training: NO"
)

print(
    "DeepGene V1 modified: NO"
)


# ============================================================
# 19. SAVE OVERLAP DATASET
# ============================================================

section(
    "16. SAVE OVERLAP DATASET"
)

overlap.to_csv(
    MATCH_FILE,
    index=False
)

print(
    "Saved:"
)

print(
    MATCH_FILE
)


# ============================================================
# 20. SAVE REPORT
# ============================================================

section(
    "17. SAVE REPORT"
)

with open(
    REPORT_FILE,
    "w",
    encoding="utf-8"
) as report:

    report.write(
        "DEEPGENE — STEPS 38–39 FUNCTIONAL "
        "OVERLAP AND TARGET ASSESSMENT\n"
    )

    report.write("=" * 70 + "\n\n")

    report.write(
        f"Functional variants: {len(overlap)}\n"
    )

    report.write(
        f"Protein-level DeepGene matches: "
        f"{len(matched)}\n"
    )

    report.write(
        f"Overlap coverage: {coverage:.2f}%\n"
    )

    report.write(
        f"Decision: {decision}\n"
    )

    report.write("\n")
    report.write(
        "TARGET STATUS\n"
    )

    report.write(
        "Independent pathogenic/benign labels: NO\n"
    )

    report.write(
        "Final ML target accepted: NO\n"
    )

    report.write(
        "ML training started: NO\n"
    )

    report.write(
        "DeepGene V1 modified: NO\n"
    )


# ============================================================
# 21. FINAL STATUS
# ============================================================

section(
    "18. STEPS 38–39 STATUS"
)

print()
print(
    "Functional dataset inspected: YES"
)

print(
    "Protein-level matching performed: YES"
)

print(
    f"Validated protein-level overlap: "
    f"{len(matched)}"
)

print(
    f"Overlap coverage: "
    f"{coverage:.2f}%"
)

print(
    f"Target assessment: "
    f"{decision}"
)

print()
print(
    "🚨 ML HAS NOT STARTED."
)

print()
print(
    "STEPS 38–39 COMPLETE"
)