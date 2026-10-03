# ============================================================
# 35. SCN1A FUNCTIONAL DATASET AUDIT
# ============================================================
#
# Purpose:
#   Inspect the Brunklaus et al. 2020 SCN1A functional
#   supplementary dataset and determine:
#
#       1. Workbook structure
#       2. Number of sheets
#       3. Number of records
#       4. Variant-related columns
#       5. Functional measurement columns
#       6. Missingness
#       7. Possible SCN1A variant representations
#       8. Preliminary overlap with DeepGene V1
#
# IMPORTANT:
#
#   This script does NOT:
#       - train ML
#       - create labels
#       - modify DeepGene V1
#       - convert functional measurements into pathogenicity labels
#       - merge datasets
#
#   Variant overlap calculated here is PRELIMINARY.
#   Final overlap requires careful variant normalization.
#
# ============================================================

from pathlib import Path
import pandas as pd
import re


# ============================================================
# 1. PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[1]

RAW_DIR = BASE_DIR / "raw"
PROCESSED_DIR = BASE_DIR / "processed"
ANALYSIS_DIR = BASE_DIR / "analysis_results"

V1_FILE = (
    PROCESSED_DIR /
    "deepgene_v1_final.csv"
)

FUNCTIONAL_FILE = (
    RAW_DIR /
    "SCN1A_functional_Brunklaus_2020.xlsx"
)

REPORT_FILE = (
    ANALYSIS_DIR /
    "35_functional_dataset_audit.txt"
)


# ============================================================
# 2. START
# ============================================================

print("=" * 70)
print("DEEPGENE — SCN1A FUNCTIONAL DATASET AUDIT")
print("=" * 70)


# ============================================================
# 3. CHECK FUNCTIONAL DATASET
# ============================================================

print("\n" + "=" * 70)
print("1. FUNCTIONAL DATASET")
print("=" * 70)

print("\nExpected file:")
print(FUNCTIONAL_FILE)

if not FUNCTIONAL_FILE.exists():

    print("\nFUNCTIONAL DATASET NOT FOUND.")
    print("--------------------------------")

    print(
        "\nPlease download the Brunklaus et al. 2020 "
        "supplementary Excel spreadsheet manually."
    )

    print(
        "\nPlace it exactly at:"
    )

    print(FUNCTIONAL_FILE)

    print(
        "\nThe script will NOT disable SSL verification "
        "or use an unsafe certificate workaround."
    )

    raise FileNotFoundError(
        "\nMissing functional dataset."
    )

print("\nFunctional dataset found.")

file_size_kb = (
    FUNCTIONAL_FILE.stat().st_size / 1024
)

print(
    f"File size: {file_size_kb:.2f} KB"
)


# ============================================================
# 4. LOAD WORKBOOK
# ============================================================

print("\n" + "=" * 70)
print("2. WORKBOOK STRUCTURE")
print("=" * 70)

try:

    workbook = pd.ExcelFile(
        FUNCTIONAL_FILE
    )

except Exception as exc:

    print("\nCould not open Excel workbook.")
    print(str(exc))
    raise


print(
    f"\nNumber of sheets: "
    f"{len(workbook.sheet_names)}"
)

for sheet in workbook.sheet_names:
    print(f"- {sheet}")


# ============================================================
# 5. LOAD ALL SHEETS
# ============================================================

sheet_data = {}

for sheet in workbook.sheet_names:

    print("\n" + "=" * 70)
    print(f"SHEET: {sheet}")
    print("=" * 70)

    try:

        data = pd.read_excel(
            FUNCTIONAL_FILE,
            sheet_name=sheet
        )

        sheet_data[sheet] = data

        print(f"Rows: {len(data)}")
        print(f"Columns: {len(data.columns)}")

        print("\nColumns:")

        for column in data.columns:

            print(
                f"  - {column}"
            )

    except Exception as exc:

        print(
            f"Could not read sheet: {exc}"
        )


# ============================================================
# 6. DATASET PREVIEW
# ============================================================

print("\n" + "=" * 70)
print("3. DATASET PREVIEW")
print("=" * 70)

for sheet, data in sheet_data.items():

    print(f"\nSHEET: {sheet}")
    print("-" * 60)

    if len(data) == 0:

        print("Empty sheet.")
        continue

    print(
        data.head(5).to_string(
            index=False
        )
    )


# ============================================================
# 7. VARIANT-RELATED COLUMN DETECTION
# ============================================================

print("\n" + "=" * 70)
print("4. VARIANT-RELATED COLUMNS")
print("=" * 70)


variant_keywords = [
    "variant",
    "mutation",
    "hgvs",
    "protein",
    "amino",
    "change",
    "substitution",
    "nucleotide",
    "c.",
    "p.",
]


variant_columns = {}

for sheet, data in sheet_data.items():

    candidates = []

    for column in data.columns:

        column_text = (
            str(column)
            .lower()
        )

        if any(
            keyword in column_text
            for keyword in variant_keywords
        ):

            candidates.append(
                str(column)
            )

    variant_columns[sheet] = candidates

    print(f"\n{sheet}")

    if candidates:

        for column in candidates:
            print(
                f"  - {column}"
            )

    else:

        print(
            "  No obvious variant-related columns detected."
        )


# ============================================================
# 8. FUNCTIONAL / MEASUREMENT COLUMN DETECTION
# ============================================================

print("\n" + "=" * 70)
print("5. FUNCTIONAL / MEASUREMENT COLUMNS")
print("=" * 70)


functional_keywords = [
    "current",
    "function",
    "functional",
    "density",
    "activation",
    "inactivation",
    "voltage",
    "electrophysi",
    "gating",
    "peak",
    "normalized",
    "residual",
    "effect",
    "assay",
]


functional_columns = {}

for sheet, data in sheet_data.items():

    candidates = []

    for column in data.columns:

        column_text = (
            str(column)
            .lower()
        )

        if any(
            keyword in column_text
            for keyword in functional_keywords
        ):

            candidates.append(
                str(column)
            )

    functional_columns[sheet] = candidates

    print(f"\n{sheet}")

    if candidates:

        for column in candidates:

            print(
                f"  - {column}"
            )

    else:

        print(
            "  No obvious functional columns detected."
        )


# ============================================================
# 9. MISSINGNESS AUDIT
# ============================================================

print("\n" + "=" * 70)
print("6. MISSINGNESS AUDIT")
print("=" * 70)


for sheet, data in sheet_data.items():

    print(f"\nSHEET: {sheet}")

    total_cells = (
        data.shape[0] *
        data.shape[1]
    )

    missing_cells = (
        int(data.isna().sum().sum())
    )

    print(
        f"Total cells: {total_cells}"
    )

    print(
        f"Missing cells: {missing_cells}"
    )

    if total_cells > 0:

        print(
            f"Missing percentage: "
            f"{missing_cells / total_cells * 100:.2f}%"
        )

    print("\nColumns with missing values:")

    missing_by_column = (
        data.isna()
        .sum()
    )

    missing_by_column = (
        missing_by_column[
            missing_by_column > 0
        ]
        .sort_values(
            ascending=False
        )
    )

    if len(missing_by_column) == 0:

        print(
            "  None"
        )

    else:

        for column, count in (
            missing_by_column.items()
        ):

            print(
                f"  - {column}: {count}"
            )


# ============================================================
# 10. LOAD DEEPGENE V1
# ============================================================

print("\n" + "=" * 70)
print("7. DEEPGENE V1")
print("=" * 70)


if not V1_FILE.exists():

    raise FileNotFoundError(
        f"DeepGene V1 not found: {V1_FILE}"
    )


v1 = pd.read_csv(
    V1_FILE,
    low_memory=False
)


print(
    f"Rows: {len(v1)}"
)

print(
    f"Unique variants: "
    f"{v1['variant_id'].nunique()}"
)

print(
    f"Unique HGVS names: "
    f"{v1['hgvs_name'].nunique()}"
)


# ============================================================
# 11. VARIANT TEXT NORMALIZATION
# ============================================================

def normalize_variant_text(value):

    if pd.isna(value):

        return ""

    text = str(value).strip()

    # Remove whitespace
    text = re.sub(
        r"\s+",
        "",
        text
    )

    # Uppercase for comparison
    text = text.upper()

    return text


# ============================================================
# 12. V1 HGVS SET
# ============================================================

v1_hgvs = set(
    v1["hgvs_name"]
    .dropna()
    .map(
        normalize_variant_text
    )
)

v1_hgvs = {
    value
    for value in v1_hgvs
    if value
}


print("\nNormalized V1 HGVS representations:")
print(len(v1_hgvs))


# ============================================================
# 13. PRELIMINARY DIRECT MATCHING
# ============================================================

print("\n" + "=" * 70)
print("8. PRELIMINARY VARIANT OVERLAP")
print("=" * 70)


overlap_records = []

functional_variant_values = []


for sheet, data in sheet_data.items():

    for column in variant_columns.get(
        sheet,
        []
    ):

        for value in data[column]:

            normalized = (
                normalize_variant_text(
                    value
                )
            )

            if not normalized:
                continue

            functional_variant_values.append(
                (
                    sheet,
                    column,
                    value,
                    normalized
                )
            )


for (
    sheet,
    column,
    original,
    normalized
) in functional_variant_values:

    if normalized in v1_hgvs:

        overlap_records.append(
            (
                sheet,
                column,
                original,
                normalized
            )
        )


print(
    f"Functional variant values inspected: "
    f"{len(functional_variant_values)}"
)

print(
    f"Direct normalized matches: "
    f"{len(overlap_records)}"
)


# ============================================================
# 14. UNIQUE MATCHES
# ============================================================

unique_overlap_values = {
    item[3]
    for item in overlap_records
}


print(
    f"Unique directly matched representations: "
    f"{len(unique_overlap_values)}"
)


# ============================================================
# 15. SHOW MATCH EXAMPLES
# ============================================================

print("\nDirect match examples:")

if overlap_records:

    for item in overlap_records[:30]:

        sheet, column, original, normalized = item

        print(
            f"- Sheet={sheet} | "
            f"Column={column} | "
            f"Original={original} | "
            f"Normalized={normalized}"
        )

else:

    print(
        "No direct normalized matches found."
    )


# ============================================================
# 16. DUPLICATE FUNCTIONAL VARIANT VALUES
# ============================================================

print("\n" + "=" * 70)
print("9. FUNCTIONAL VARIANT DUPLICATION")
print("=" * 70)


normalized_values = [
    item[3]
    for item in functional_variant_values
]


value_counts = (
    pd.Series(
        normalized_values
    )
    .value_counts()
)


print(
    f"Unique normalized functional variant values: "
    f"{len(value_counts)}"
)

print(
    f"Repeated normalized variant values: "
    f"{(value_counts > 1).sum()}"
)


print("\nMost repeated values:")

for value, count in (
    value_counts.head(20).items()
):

    print(
        f"- {value}: {count}"
    )


# ============================================================
# 17. FUNCTIONAL DATASET SIZE
# ============================================================

print("\n" + "=" * 70)
print("10. FUNCTIONAL DATASET SIZE")
print("=" * 70)


total_rows = sum(
    len(data)
    for data in sheet_data.values()
)


print(
    f"Total workbook rows across sheets: "
    f"{total_rows}"
)


# ============================================================
# 18. TARGET INTERPRETATION
# ============================================================

print("\n" + "=" * 70)
print("11. TARGET INTERPRETATION")
print("=" * 70)


print(
    """
The functional dataset must NOT yet be converted into
Pathogenic / Benign labels.

The purpose of this step is to determine what was actually
measured experimentally.

Potential measurements may include:

- sodium current
- current density
- activation
- inactivation
- voltage dependence
- gating properties
- residual function
- functional effect categories

The biological meaning of each measurement must be understood
before deciding whether it can become a DeepGene prediction target.
"""
)


# ============================================================
# 19. IMPORTANT OVERLAP WARNING
# ============================================================

print("\n" + "=" * 70)
print("12. OVERLAP WARNING")
print("=" * 70)


print(
    """
The overlap calculated in this script is preliminary.

A direct text match does not prove biological identity.

A biological variant may have different:

- transcript versions
- HGVS descriptions
- protein descriptions
- genomic representations
- reference assemblies

Therefore, direct matches are useful for discovery but are NOT
yet sufficient for constructing an ML dataset.
"""
)


# ============================================================
# 20. FINAL STATUS
# ============================================================

print("\n" + "=" * 70)
print("13. STEP 35 STATUS")
print("=" * 70)


print(
    "Functional dataset located: YES"
)

print(
    "Workbook successfully inspected: YES"
)

print(
    "Functional measurements identified: "
    f"{sum(len(x) for x in functional_columns.values())}"
)

print(
    "Preliminary direct overlap calculated: YES"
)

print(
    "Validated biological variant overlap: NOT YET"
)

print(
    "Independent ML target defined: NOT YET"
)

print(
    "\n🚨 ML HAS NOT STARTED."
)


# ============================================================
# 21. NEXT STEP
# ============================================================

print("\n" + "=" * 70)
print("14. NEXT STEP")
print("=" * 70)


print(
    """
The next step is to inspect the actual functional measurements
and determine:

1. What biological quantity was measured?
2. How many variants have usable measurements?
3. Are measurements comparable across experiments?
4. How many variants overlap DeepGene V1 after proper normalization?
5. Can the measurements form a defensible target?
"""
)


# ============================================================
# 22. END
# ============================================================

print("\n" + "=" * 70)
print("STEP 35 COMPLETE")
print("=" * 70)