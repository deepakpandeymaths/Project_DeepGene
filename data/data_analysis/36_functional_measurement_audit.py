# ============================================================
# 36. FUNCTIONAL MEASUREMENT AUDIT
# ============================================================
#
# Purpose:
#   Inspect the actual experimental functional measurements
#   in the Brunklaus 2020 SCN1A dataset.
#
# This step does NOT:
#   - create pathogenic/benign labels
#   - train an ML model
#   - modify DeepGene V1
#   - assume that "Function" is the final target
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
    "36_functional_measurement_audit.txt"
)


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

def section(title):
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def clean_text(value):
    if pd.isna(value):
        return ""

    return str(value).strip()


def numeric_like_count(series):
    """
    Count values that can reasonably be interpreted as numeric.

    Handles values such as:
        40
        40.0
        -15.2
        "40"
        "<5"
        ">10"

    Values containing obvious ranges or text are not counted
    as numeric-like.
    """

    count = 0

    for value in series:

        if pd.isna(value):
            continue

        text = str(value).strip()

        if text == "":
            continue

        # Direct numeric conversion
        try:
            float(text)
            count += 1
            continue
        except ValueError:
            pass

        # Simple numeric values with units or symbols
        if re.fullmatch(
            r"[<>~=]?\s*-?\d+(?:\.\d+)?\s*(?:mV|ms|%)?",
            text,
            flags=re.IGNORECASE
        ):
            count += 1

    return count


def nonempty_count(series):
    return int(series.notna().sum())


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

os.makedirs(RESULTS_DIR, exist_ok=True)

section("DEEPGENE — SCN1A FUNCTIONAL MEASUREMENT AUDIT")

print("Functional dataset:")
print(FUNCTIONAL_FILE)

if not os.path.exists(FUNCTIONAL_FILE):
    print()
    print("ERROR: Functional dataset not found.")
    raise SystemExit(1)

print()
print("Functional dataset found.")


# ------------------------------------------------------------
# LOAD WORKBOOK
# ------------------------------------------------------------

section("1. LOAD FUNCTIONAL WORKBOOK")

try:
    workbook = pd.ExcelFile(FUNCTIONAL_FILE)
except Exception as exc:
    print("ERROR: Could not open workbook.")
    print(exc)
    raise SystemExit(1)

print("Sheets:")
for sheet in workbook.sheet_names:
    print(f"- {sheet}")


# ------------------------------------------------------------
# LOAD MAIN SHEET
# ------------------------------------------------------------

section("2. LOAD MAIN FUNCTIONAL SHEET")

MAIN_SHEET = "Sheet1"

if MAIN_SHEET not in workbook.sheet_names:
    MAIN_SHEET = workbook.sheet_names[0]

df = pd.read_excel(
    FUNCTIONAL_FILE,
    sheet_name=MAIN_SHEET
)

# Remove completely empty rows
df = df.dropna(
    how="all"
).reset_index(drop=True)

print(f"Sheet used: {MAIN_SHEET}")
print(f"Rows after removing empty rows: {len(df)}")
print(f"Columns: {len(df.columns)}")


# ------------------------------------------------------------
# BASIC VARIANT INFORMATION
# ------------------------------------------------------------

section("3. VARIANT INFORMATION")

variant_columns = [
    "cDNA",
    "Variant",
    "Function",
    "Phenotype",
    "Inheritance"
]

for column in variant_columns:

    if column not in df.columns:
        print(f"{column}: NOT FOUND")
        continue

    count = nonempty_count(df[column])

    print(f"{column}: {count} non-empty values")


if "Variant" in df.columns:

    variants = (
        df["Variant"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    variants = variants[variants != ""]

    print()
    print(f"Unique protein-level variant descriptions: {variants.nunique()}")
    print(f"Total non-empty variant descriptions: {len(variants)}")


# ------------------------------------------------------------
# FUNCTION CLASSIFICATION
# ------------------------------------------------------------

section("4. FUNCTION CLASSIFICATION")

if "Function" in df.columns:

    function_values = (
        df["Function"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    function_values = function_values[
        function_values != ""
    ]

    print(f"Non-empty Function values: {len(function_values)}")
    print(f"Unique Function values: {function_values.nunique()}")

    print()
    print("Function categories:")

    counts = function_values.value_counts()

    for value, count in counts.items():
        print(f"- {value}: {count}")

else:
    print("Function column not found.")


# ------------------------------------------------------------
# CANDIDATE EXPERIMENTAL MEASUREMENTS
# ------------------------------------------------------------

section("5. EXPERIMENTAL MEASUREMENT COLUMNS")

# These are the columns that appear to represent actual
# electrophysiological measurements based on their names.

measurement_keywords = [
    "Whole-cell current",
    "Peak current Density",
    "Change in Peak current Density",
    "Persistent Na current",
    "Activation",
    "Shift of Activation",
    "Fast Inactivation",
    "Shift of Fast Inactivation",
    "Recovery from Fast Inactivation",
    "Onset of Slow Inactivation",
    "Slow inactivation",
    "Recovery from Slow inactivation"
]

measurement_columns = []

for column in df.columns:

    column_text = str(column).strip().lower()

    for keyword in measurement_keywords:

        if keyword.lower() in column_text:

            measurement_columns.append(column)
            break


# Remove duplicates while preserving order
measurement_columns = list(
    dict.fromkeys(measurement_columns)
)

print(f"Candidate experimental measurement columns: {len(measurement_columns)}")

for column in measurement_columns:
    print(f"- {column}")


# ------------------------------------------------------------
# MEASUREMENT COVERAGE
# ------------------------------------------------------------

section("6. MEASUREMENT COVERAGE")

measurement_summary = []

for column in measurement_columns:

    nonempty = nonempty_count(df[column])
    numeric_like = numeric_like_count(df[column])

    missing = len(df) - nonempty

    measurement_summary.append(
        {
            "column": column,
            "non_empty": nonempty,
            "missing": missing,
            "numeric_like": numeric_like,
            "coverage_percent": round(
                100 * nonempty / len(df),
                2
            )
        }
    )

    print()
    print(f"Measurement: {column}")
    print(f"  Non-empty: {nonempty}")
    print(f"  Missing: {missing}")
    print(f"  Numeric-like: {numeric_like}")
    print(
        f"  Coverage: "
        f"{100 * nonempty / len(df):.2f}%"
    )


# ------------------------------------------------------------
# MEASUREMENT VALUES PREVIEW
# ------------------------------------------------------------

section("7. MEASUREMENT VALUE PREVIEW")

for column in measurement_columns:

    values = (
        df[column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    values = values[
        values != ""
    ]

    print()
    print(f"--- {column} ---")

    if len(values) == 0:
        print("No values.")
        continue

    unique_values = values.unique()

    print(f"Unique non-empty values: {len(unique_values)}")

    for value in unique_values[:10]:
        print(f"  {value}")

    if len(unique_values) > 10:
        print("  ...")


# ------------------------------------------------------------
# ACMG COLUMNS — DO NOT USE AS EXPERIMENTAL TARGET
# ------------------------------------------------------------

section("8. ACMG INTERPRETATION COLUMNS")

acmg_columns = [
    column
    for column in df.columns
    if "ACMG" in str(column)
]

for column in acmg_columns:

    print()
    print(f"- {column}")
    print(f"  Non-empty: {nonempty_count(df[column])}")

print()
print(
    "ACMG columns are interpretation/evidence fields and "
    "will NOT be treated as independent experimental measurements."
)


# ------------------------------------------------------------
# INDEPENDENCE CHECK
# ------------------------------------------------------------

section("9. POTENTIAL TARGET CATEGORIES")

print()
print("The following categories were detected:")

print()
print("A. QUALITATIVE FUNCTION")
print("   Example: Loss / Mixed / etc.")
print("   Status: REQUIRES PROVENANCE AUDIT")

print()
print("B. CONTINUOUS ELECTROPHYSIOLOGICAL MEASUREMENTS")
print("   Examples:")
print("   - Peak current density")
print("   - Persistent current")
print("   - Activation shift")
print("   - Inactivation shift")
print("   - Recovery measurements")
print("   Status: POTENTIALLY SUITABLE")

print()
print("C. ACMG CLASSIFICATION")
print("   Example: Class 4 / Class 5")
print("   Status: NOT AN INDEPENDENT TARGET")

print()
print("D. IN-SILICO PREDICTIONS")
print("   Examples: SIFT, CADD, REVEL, PolyPhen")
print("   Status: INPUT/EVIDENCE FEATURES, NOT TARGET")


# ------------------------------------------------------------
# DEEPGENE V1 CHECK
# ------------------------------------------------------------

section("10. DEEPGENE V1 CHECK")

if not os.path.exists(V1_FILE):

    print("DeepGene V1 file not found.")
    print(V1_FILE)

else:

    v1 = pd.read_csv(V1_FILE)

    print(f"DeepGene V1 rows: {len(v1)}")
    print(
        f"DeepGene V1 unique variants: "
        f"{v1['variant_id'].nunique()}"
    )

    print()
    print("DeepGene V1 is unchanged by this step.")


# ------------------------------------------------------------
# TARGET DECISION
# ------------------------------------------------------------

section("11. STEP 36 INTERPRETATION")

print()
print("Current conclusion:")

print()
print(
    "The dataset contains genuine experimental "
    "electrophysiological measurements."
)

print()
print(
    "However, no measurement is automatically accepted "
    "as the final DeepGene ML target."
)

print()
print(
    "The next decision requires determining whether one "
    "or more measurements are comparable across variants "
    "and sufficiently independent from the ClinVar-derived "
    "features in DeepGene V1."
)

print()
print(
    "No pathogenic/benign labels were created."
)

print(
    "No ML model was trained."
)

print(
    "DeepGene V1 was not modified."
)


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
        "DEEPGENE — SCN1A FUNCTIONAL MEASUREMENT AUDIT\n"
    )
    report.write("=" * 70 + "\n\n")

    report.write(
        f"Workbook: {FUNCTIONAL_FILE}\n"
    )

    report.write(
        f"Sheet: {MAIN_SHEET}\n"
    )

    report.write(
        f"Rows: {len(df)}\n"
    )

    report.write(
        f"Columns: {len(df.columns)}\n\n"
    )

    report.write(
        "FUNCTION CATEGORIES\n"
    )
    report.write("-" * 70 + "\n")

    if "Function" in df.columns:

        function_values = (
            df["Function"]
            .dropna()
            .astype(str)
            .str.strip()
        )

        function_values = function_values[
            function_values != ""
        ]

        for value, count in function_values.value_counts().items():
            report.write(
                f"{value}: {count}\n"
            )

    report.write("\n")
    report.write(
        "EXPERIMENTAL MEASUREMENTS\n"
    )
    report.write("-" * 70 + "\n")

    for item in measurement_summary:

        report.write(
            f"\n{item['column']}\n"
        )

        report.write(
            f"  non_empty: {item['non_empty']}\n"
        )

        report.write(
            f"  missing: {item['missing']}\n"
        )

        report.write(
            f"  numeric_like: {item['numeric_like']}\n"
        )

        report.write(
            f"  coverage_percent: "
            f"{item['coverage_percent']}\n"
        )

    report.write("\n")
    report.write(
        "TARGET STATUS\n"
    )
    report.write("-" * 70 + "\n")

    report.write(
        "Independent ML target: NOT YET DEFINED\n"
    )

    report.write(
        "Pathogenic/Benign labels: NOT CREATED\n"
    )

    report.write(
        "ML training: NOT STARTED\n"
    )

    report.write(
        "DeepGene V1 modified: NO\n"
    )


# ------------------------------------------------------------
# FINAL STATUS
# ------------------------------------------------------------

section("13. STEP 36 STATUS")

print()
print("Functional measurements inspected: YES")
print(
    f"Candidate measurement columns identified: "
    f"{len(measurement_columns)}"
)
print("Measurement coverage calculated: YES")
print("ACMG fields separated from measurements: YES")
print("Independent ML target defined: NOT YET")

print()
print("🚨 ML HAS NOT STARTED.")

print()
print("Report:")
print(OUTPUT_FILE)

print()
print("STEP 36 COMPLETE")