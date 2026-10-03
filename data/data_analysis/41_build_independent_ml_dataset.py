# ============================================================
# 41. BUILD INDEPENDENT FUNCTIONAL ML DATASET
# ============================================================
#
# Purpose:
# Build a clean dataset using:
#
#   INPUT FEATURES:
#       Variant-derived biological properties
#
#   TARGET:
#       Independent experimental functional classification
#
# IMPORTANT:
#   ClinVar pathogenicity/evidence is NOT used as a feature.
#   This is NOT a pathogenic-vs-benign model.
#   ML training has NOT started.
#
# ============================================================

from pathlib import Path
import pandas as pd
import numpy as np
import re


# ============================================================
# PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

OVERLAP_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_functional_overlap_v1.csv"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_independent_ml_dataset_v1.csv"
)

REPORT_FILE = (
    PROJECT_ROOT
    / "data"
    / "analysis_results"
    / "41_independent_ml_dataset.txt"
)


# ============================================================
# AMINO ACID PROPERTIES
# ============================================================

AA_PROPERTIES = {

    "A": {"hydrophobicity": 1.8,  "charge": 0,  "volume": 88.6},
    "R": {"hydrophobicity": -4.5, "charge": 1,  "volume": 173.4},
    "N": {"hydrophobicity": -3.5, "charge": 0,  "volume": 114.1},
    "D": {"hydrophobicity": -3.5, "charge": -1, "volume": 111.1},
    "C": {"hydrophobicity": 2.5,  "charge": 0,  "volume": 108.5},
    "Q": {"hydrophobicity": -3.5, "charge": 0,  "volume": 143.8},
    "E": {"hydrophobicity": -3.5, "charge": -1, "volume": 138.4},
    "G": {"hydrophobicity": -0.4, "charge": 0,  "volume": 60.1},
    "H": {"hydrophobicity": -3.2, "charge": 1,  "volume": 153.2},
    "I": {"hydrophobicity": 4.5,  "charge": 0,  "volume": 166.7},
    "L": {"hydrophobicity": 3.8,  "charge": 0,  "volume": 166.7},
    "K": {"hydrophobicity": -3.9, "charge": 1,  "volume": 168.6},
    "M": {"hydrophobicity": 1.9,  "charge": 0,  "volume": 162.9},
    "F": {"hydrophobicity": 2.8,  "charge": 0,  "volume": 189.9},
    "P": {"hydrophobicity": -1.6, "charge": 0,  "volume": 112.7},
    "S": {"hydrophobicity": -0.8, "charge": 0,  "volume": 89.0},
    "T": {"hydrophobicity": -0.7, "charge": 0,  "volume": 116.1},
    "W": {"hydrophobicity": -0.9, "charge": 0,  "volume": 227.8},
    "Y": {"hydrophobicity": -1.3, "charge": 0,  "volume": 193.6},
    "V": {"hydrophobicity": 4.2,  "charge": 0,  "volume": 140.0},
}


# ============================================================
# PARSE PROTEIN VARIANT
# ============================================================

def parse_protein_variant(value):

    if pd.isna(value):
        return None, None, None

    text = str(value).strip()

    # Three-letter amino-acid notation
    aa_map = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
    }

    pattern = (
        r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|"
        r"Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)"
        r"(\d+)"
        r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|"
        r"Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)"
    )

    match = re.search(pattern, text)

    if match:

        ref = aa_map[match.group(1)]
        position = int(match.group(2))
        alt = aa_map[match.group(3)]

        return ref, position, alt

    # One-letter notation
    match = re.search(
        r"([A-Z])(\d+)([A-Z])",
        text
    )

    if match:

        ref = match.group(1)
        position = int(match.group(2))
        alt = match.group(3)

        if (
            ref in AA_PROPERTIES
            and alt in AA_PROPERTIES
        ):
            return ref, position, alt

    return None, None, None


# ============================================================
# LOAD OVERLAP DATASET
# ============================================================

print("=" * 70)
print("41. BUILD INDEPENDENT FUNCTIONAL ML DATASET")
print("=" * 70)

if not OVERLAP_FILE.exists():
    raise FileNotFoundError(
        f"Missing file:\n{OVERLAP_FILE}"
    )

df = pd.read_csv(OVERLAP_FILE)

print()
print("Input dataset:")
print("Rows:", len(df))
print("Columns:", len(df.columns))

print()
print("Input columns:")

for column in df.columns:
    print(" -", column)


# ============================================================
# VERIFY EXPECTED COLUMNS
# ============================================================

required_columns = [
    "functional_cDNA",
    "functional_protein",
    "normalized_protein",
    "protein_position",
    "functional_class",
    "phenotype",
    "inheritance",
    "deepgene_match",
    "deepgene_match_count",
    "deepgene_variant_ids",
]

missing_columns = [
    column
    for column in required_columns
    if column not in df.columns
]

if missing_columns:

    raise ValueError(
        "Missing expected columns:\n"
        + "\n".join(
            f" - {column}"
            for column in missing_columns
        )
    )


# ============================================================
# KEEP ONLY SUCCESSFULLY MATCHED VARIANTS
# ============================================================

# ============================================================
# KEEP SUCCESSFULLY MATCHED VARIANTS
# ============================================================
#
# Step 38-39 stores the actual matched DeepGene IDs in
# deepgene_variant_ids. Therefore we use that field rather
# than assuming deepgene_match contains the text "YES".
#

df["deepgene_variant_ids"] = (
    df["deepgene_variant_ids"]
    .fillna("")
    .astype(str)
    .str.strip()
)

df["deepgene_match_count"] = pd.to_numeric(
    df["deepgene_match_count"],
    errors="coerce"
).fillna(0)

matched = df[
    (
        df["deepgene_match_count"] > 0
    )
    &
    (
        df["deepgene_variant_ids"] != ""
    )
    &
    (
        df["deepgene_variant_ids"].str.lower() != "nan"
    )
].copy()

print()
print("Matched variants:", len(matched))

if len(matched) == 0:
    raise ValueError(
        "No matched variants detected. "
        "Check deepgene_match_count and deepgene_variant_ids."
    )

print()
print("Matched variants:", len(matched))


# ============================================================
# EXTRACT DEEPGENE VARIANT ID
# ============================================================

# deepgene_variant_ids can contain one or more IDs.
#
# For this dataset, every functional variant should correspond
# to the matched DeepGene variant(s).

matched["variant_id"] = (
    matched["deepgene_variant_ids"]
    .astype(str)
    .str.split("|")
    .str[0]
    .str.strip()
)


# ============================================================
# PARSE NORMALIZED PROTEIN VARIANT
# ============================================================

parsed = matched[
    [
        "variant_id",
        "functional_protein",
        "normalized_protein",
        "protein_position",
        "functional_class",
        "phenotype",
        "inheritance",
    ]
].copy()

parsed.columns = [
    "variant_id",
    "functional_protein",
    "protein_variant",
    "protein_position_original",
    "functional_effect",
    "phenotype",
    "inheritance",
]


parsed[
    [
        "reference_aa",
        "parsed_position",
        "alternate_aa",
    ]
] = parsed["protein_variant"].apply(
    lambda x: pd.Series(
        parse_protein_variant(x)
    )
)


# ============================================================
# USE PARSED POSITION
# ============================================================

parsed["protein_position"] = pd.to_numeric(
    parsed["parsed_position"],
    errors="coerce"
)


# ============================================================
# REMOVE UNPARSED VARIANTS
# ============================================================

before = len(parsed)

parsed = parsed.dropna(
    subset=[
        "reference_aa",
        "protein_position",
        "alternate_aa",
    ]
).copy()

after = len(parsed)

print()
print("Protein parsing:")
print("Before:", before)
print("Successfully parsed:", after)
print("Failed:", before - after)


# ============================================================
# REMOVE DUPLICATE DEEPGENE VARIANTS
# ============================================================

parsed = parsed.drop_duplicates(
    subset=["variant_id"]
).reset_index(drop=True)


# ============================================================
# AMINO ACID FEATURES
# ============================================================

parsed["same_amino_acid"] = (
    parsed["reference_aa"]
    == parsed["alternate_aa"]
).astype(int)


def get_property(
    amino_acid,
    property_name
):

    if amino_acid not in AA_PROPERTIES:
        return np.nan

    return AA_PROPERTIES[
        amino_acid
    ][property_name]


# ------------------------------------------------------------
# Hydrophobicity
# ------------------------------------------------------------

parsed["reference_hydrophobicity"] = (
    parsed["reference_aa"]
    .apply(
        lambda x:
        get_property(
            x,
            "hydrophobicity"
        )
    )
)

parsed["alternate_hydrophobicity"] = (
    parsed["alternate_aa"]
    .apply(
        lambda x:
        get_property(
            x,
            "hydrophobicity"
        )
    )
)

parsed["delta_hydrophobicity"] = (
    parsed["alternate_hydrophobicity"]
    - parsed["reference_hydrophobicity"]
)

parsed["absolute_delta_hydrophobicity"] = (
    parsed["delta_hydrophobicity"]
    .abs()
)


# ------------------------------------------------------------
# Charge
# ------------------------------------------------------------

parsed["reference_charge"] = (
    parsed["reference_aa"]
    .apply(
        lambda x:
        get_property(
            x,
            "charge"
        )
    )
)

parsed["alternate_charge"] = (
    parsed["alternate_aa"]
    .apply(
        lambda x:
        get_property(
            x,
            "charge"
        )
    )
)

parsed["delta_charge"] = (
    parsed["alternate_charge"]
    - parsed["reference_charge"]
)

parsed["absolute_delta_charge"] = (
    parsed["delta_charge"]
    .abs()
)


# ------------------------------------------------------------
# Molecular volume
# ------------------------------------------------------------

parsed["reference_volume"] = (
    parsed["reference_aa"]
    .apply(
        lambda x:
        get_property(
            x,
            "volume"
        )
    )
)

parsed["alternate_volume"] = (
    parsed["alternate_aa"]
    .apply(
        lambda x:
        get_property(
            x,
            "volume"
        )
    )
)

parsed["delta_volume"] = (
    parsed["alternate_volume"]
    - parsed["reference_volume"]
)

parsed["absolute_delta_volume"] = (
    parsed["delta_volume"]
    .abs()
)


# ============================================================
# AMINO ACID CHANGE TYPE
# ============================================================

def classify_change(
    ref,
    alt
):

    if ref == alt:
        return "same"

    ref_charge = AA_PROPERTIES[
        ref
    ]["charge"]

    alt_charge = AA_PROPERTIES[
        alt
    ]["charge"]

    if ref_charge != alt_charge:
        return "charge_change"

    ref_hydro = AA_PROPERTIES[
        ref
    ]["hydrophobicity"]

    alt_hydro = AA_PROPERTIES[
        alt
    ]["hydrophobicity"]

    if (
        (ref_hydro >= 0)
        !=
        (alt_hydro >= 0)
    ):
        return "hydrophobicity_class_change"

    return "physicochemical_change"


parsed["amino_acid_change_type"] = (
    parsed.apply(
        lambda row:
        classify_change(
            row["reference_aa"],
            row["alternate_aa"],
        ),
        axis=1,
    )
)


# ============================================================
# TARGET
# ============================================================

# IMPORTANT:
#
# functional_effect comes from the independent experimental
# functional dataset.
#
# We DO NOT convert it into pathogenic/benign.

parsed["functional_target"] = (
    parsed["functional_effect"]
    .astype(str)
    .str.strip()
)


# ============================================================
# FEATURE COLUMNS
# ============================================================

feature_columns = [

    "protein_position",

    "same_amino_acid",

    "reference_hydrophobicity",
    "alternate_hydrophobicity",
    "delta_hydrophobicity",
    "absolute_delta_hydrophobicity",

    "reference_charge",
    "alternate_charge",
    "delta_charge",
    "absolute_delta_charge",

    "reference_volume",
    "alternate_volume",
    "delta_volume",
    "absolute_delta_volume",

]


# ============================================================
# FINAL DATASET
# ============================================================

final_columns = [

    "variant_id",

    "functional_protein",
    "protein_variant",

    "reference_aa",
    "alternate_aa",

    "protein_position",

    "amino_acid_change_type",

    "functional_effect",
    "functional_target",

    "phenotype",
    "inheritance",

] + feature_columns


final_df = parsed[
    final_columns
].copy()


# ============================================================
# CHECK FEATURE MISSINGNESS
# ============================================================

print()
print("Missing feature values:")

print(
    final_df[
        feature_columns
    ].isna().sum()
)


# ============================================================
# TARGET DISTRIBUTION
# ============================================================

print()
print("=" * 70)
print("FUNCTIONAL TARGET DISTRIBUTION")
print("=" * 70)

print(
    final_df[
        "functional_target"
    ].value_counts(
        dropna=False
    )
)


# ============================================================
# DATASET SUMMARY
# ============================================================

print()
print("=" * 70)
print("FINAL DATASET")
print("=" * 70)

print(
    "Rows:",
    len(final_df)
)

print(
    "Unique variants:",
    final_df[
        "variant_id"
    ].nunique()
)

print(
    "Features:",
    len(feature_columns)
)


# ============================================================
# SAVE DATASET
# ============================================================

OUTPUT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

final_df.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# REPORT
# ============================================================

report = []

report.append(
    "DEEPGENE STEP 41"
)

report.append(
    "INDEPENDENT FUNCTIONAL ML DATASET"
)

report.append("=" * 60)

report.append(
    f"Input overlap rows: {len(df)}"
)

report.append(
    f"Matched rows: {len(matched)}"
)

report.append(
    f"Parsed rows: {len(final_df)}"
)

report.append(
    f"Unique DeepGene variants: "
    f"{final_df['variant_id'].nunique()}"
)

report.append("")

report.append(
    "TARGET"
)

report.append(
    "Independent experimental functional classification"
)

report.append(
    "NOT pathogenic/benign"
)

report.append("")

report.append(
    "FEATURE POLICY"
)

report.append(
    "Variant-derived biological features only."
)

report.append(
    "ClinVar clinical significance excluded."
)

report.append(
    "ClinVar review status excluded."
)

report.append(
    "ClinVar pathogenicity counts excluded."
)

report.append(
    "ClinVar submitter/reliability fields excluded."
)

report.append(
    "ClinVar conflict fields excluded."
)

report.append(
    "Phenotype retained only as source metadata;"
)

report.append(
    "not used as an ML feature."
)

report.append(
    "Inheritance retained only as source metadata;"
)

report.append(
    "not used as an ML feature."
)

report.append("")

report.append(
    "ML FEATURES:"
)

for feature in feature_columns:

    report.append(
        f" - {feature}"
    )

report.append("")

report.append(
    "FUNCTIONAL TARGET DISTRIBUTION:"
)

for category, count in (
    final_df[
        "functional_target"
    ]
    .value_counts(
        dropna=False
    )
    .items()
):

    report.append(
        f" - {category}: {count}"
    )

report.append("")

report.append(
    "STATUS"
)

report.append(
    "Independent experimental dataset: YES"
)

report.append(
    "Biological variant overlap: YES"
)

report.append(
    "Independent target candidate: YES"
)

report.append(
    "ClinVar pathogenicity target: NO"
)

report.append(
    "ML training: NOT STARTED"
)

report.append(
    "Train/test split: NOT STARTED"
)

report.append(
    "Model training: NOT STARTED"
)

report.append(
    "Evaluation: NOT STARTED"
)


REPORT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True
)

REPORT_FILE.write_text(
    "\n".join(report),
    encoding="utf-8"
)


# ============================================================
# FINAL
# ============================================================

print()
print("=" * 70)
print("STEP 41 COMPLETE")
print("=" * 70)

print()
print("Dataset:")
print(OUTPUT_FILE)

print()
print("Report:")
print(REPORT_FILE)

print()
print("STATUS: PASS")

print()
print("ML TRAINING: NOT STARTED")