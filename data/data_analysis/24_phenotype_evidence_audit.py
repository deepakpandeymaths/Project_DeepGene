# ============================================================
# 24. PHENOTYPE EVIDENCE AUDIT
# ============================================================
#
# Purpose:
# Audit phenotype information stored in the final V1 dataset
# and compare it with the underlying ClinVar records.
#
# This step checks:
#   - phenotype availability
#   - empty / placeholder phenotype values
#   - phenotype diversity
#   - phenotype consistency per variant
#   - consistency between SQLite and final V1 dataset
#
# This is an AUDIT step.
# It does not modify the V1 dataset.
# ============================================================

import os
import sqlite3
import pandas as pd


# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

BASE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

DB_PATH = os.path.join(
    BASE_DIR,
    "database",
    "scn1a.db"
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
    "24_phenotype_evidence_audit.txt"
)


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 24")
print("PHENOTYPE EVIDENCE AUDIT")
print("=" * 70)


# ------------------------------------------------------------
# LOAD FINAL DATASET
# ------------------------------------------------------------

print("\nLoading final V1 dataset...")

final_df = pd.read_csv(INPUT_FILE)

print(f"Rows loaded: {len(final_df)}")
print(f"Columns loaded: {len(final_df.columns)}")


# ------------------------------------------------------------
# CONNECT DATABASE
# ------------------------------------------------------------

print("\nConnecting to SQLite database...")

conn = sqlite3.connect(DB_PATH)


# ------------------------------------------------------------
# LOAD PHENOTYPE INFORMATION
# ------------------------------------------------------------

clinvar = pd.read_sql_query(
    """
    SELECT
        variant_id,
        phenotypes
    FROM clinvar_records
    """,
    conn
)

print(
    f"ClinVar phenotype records loaded: "
    f"{len(clinvar)}"
)


# ------------------------------------------------------------
# NORMALIZE PHENOTYPE VALUES
# ------------------------------------------------------------

clinvar["phenotypes_clean"] = (
    clinvar["phenotypes"]
    .fillna("")
    .astype(str)
    .str.strip()
)

final_df["phenotypes_clean"] = (
    final_df["phenotypes"]
    .fillna("")
    .astype(str)
    .str.strip()
)


# ------------------------------------------------------------
# 1. FINAL DATASET PHENOTYPE COVERAGE
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. FINAL DATASET PHENOTYPE COVERAGE")
print("=" * 70)

phenotype_available_counts = (
    final_df["phenotype_available"]
    .value_counts(dropna=False)
    .sort_index()
)

print(
    phenotype_available_counts.to_string()
)

available_count = (
    final_df["phenotype_available"] == 1
).sum()

unavailable_count = (
    final_df["phenotype_available"] == 0
).sum()

print(
    f"\nPhenotype available: "
    f"{available_count}"
)

print(
    f"Phenotype unavailable: "
    f"{unavailable_count}"
)


# ------------------------------------------------------------
# 2. PHENOTYPE FIELD QUALITY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("2. PHENOTYPE FIELD QUALITY")
print("=" * 70)

empty_phenotype_count = (
    final_df["phenotypes_clean"] == ""
).sum()

placeholder_values = {
    "",
    "-",
    "na",
    "n/a",
    "none",
    "not provided"
}

placeholder_count = (
    final_df["phenotypes_clean"]
    .str.lower()
    .isin(placeholder_values)
    .sum()
)

non_empty_count = (
    len(final_df) - empty_phenotype_count
)

print(
    f"Empty phenotype values: "
    f"{empty_phenotype_count}"
)

print(
    f"Placeholder/missing-style phenotype values: "
    f"{placeholder_count}"
)

print(
    f"Non-empty phenotype values: "
    f"{non_empty_count}"
)


# ------------------------------------------------------------
# 3. UNIQUE PHENOTYPE STRINGS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. PHENOTYPE STRING DIVERSITY")
print("=" * 70)

unique_phenotypes = (
    final_df["phenotypes_clean"]
    .nunique()
)

print(
    f"Unique phenotype strings in final dataset: "
    f"{unique_phenotypes}"
)

print("\nMost common phenotype strings:")

phenotype_frequency = (
    final_df["phenotypes_clean"]
    .value_counts()
    .head(20)
)

print(
    phenotype_frequency.to_string()
)


# ------------------------------------------------------------
# 4. PHENOTYPE STRING LENGTH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. PHENOTYPE STRING LENGTH")
print("=" * 70)

phenotype_lengths = (
    final_df["phenotypes_clean"]
    .str.len()
)

print(
    f"Minimum length: "
    f"{phenotype_lengths.min()}"
)

print(
    f"Maximum length: "
    f"{phenotype_lengths.max()}"
)

print(
    f"Mean length: "
    f"{phenotype_lengths.mean():.2f}"
)

print(
    f"Median length: "
    f"{phenotype_lengths.median():.2f}"
)


# ------------------------------------------------------------
# 5. SQLITE PHENOTYPE COVERAGE
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. SQLITE PHENOTYPE COVERAGE")
print("=" * 70)

sqlite_empty_count = (
    clinvar["phenotypes_clean"] == ""
).sum()

sqlite_nonempty_count = (
    len(clinvar) - sqlite_empty_count
)

print(
    f"Empty phenotype records in SQLite: "
    f"{sqlite_empty_count}"
)

print(
    f"Non-empty phenotype records in SQLite: "
    f"{sqlite_nonempty_count}"
)


# ------------------------------------------------------------
# 6. PHENOTYPE CONSISTENCY PER VARIANT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("6. PHENOTYPE CONSISTENCY PER VARIANT")
print("=" * 70)

phenotype_per_variant = (
    clinvar
    .groupby("variant_id")["phenotypes_clean"]
    .nunique()
)

multiple_phenotype_values = (
    phenotype_per_variant > 1
).sum()

single_phenotype_value = (
    phenotype_per_variant == 1
).sum()

print(
    f"Variants with exactly one phenotype value: "
    f"{single_phenotype_value}"
)

print(
    f"Variants with multiple distinct phenotype values: "
    f"{multiple_phenotype_values}"
)


# ------------------------------------------------------------
# 7. CREATE ONE SQLITE PHENOTYPE VALUE PER VARIANT
# ------------------------------------------------------------
#
# Step 22 showed that each variant currently has one distinct
# phenotype value. We therefore create one source value per
# variant for direct comparison with the final dataset.
# ------------------------------------------------------------

sqlite_variant_phenotypes = (
    clinvar[
        [
            "variant_id",
            "phenotypes_clean"
        ]
    ]
    .drop_duplicates(
        subset=[
            "variant_id",
            "phenotypes_clean"
        ]
    )
)

sqlite_variant_phenotypes = (
    sqlite_variant_phenotypes
    .drop_duplicates(
        subset=["variant_id"]
    )
    .rename(
        columns={
            "phenotypes_clean":
            "phenotypes_sqlite"
        }
    )
)


# ------------------------------------------------------------
# 8. FINAL DATASET VS SQLITE
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("7. FINAL DATASET VS SQLITE")
print("=" * 70)

comparison = final_df[
    [
        "variant_id",
        "phenotypes_clean",
        "phenotype_available"
    ]
].copy()

comparison = comparison.rename(
    columns={
        "phenotypes_clean":
        "phenotypes_final"
    }
)

comparison = comparison.merge(
    sqlite_variant_phenotypes,
    on="variant_id",
    how="left"
)

comparison["phenotypes_sqlite"] = (
    comparison["phenotypes_sqlite"]
    .fillna("")
    .astype(str)
    .str.strip()
)

phenotype_mismatches = (
    comparison["phenotypes_final"]
    != comparison["phenotypes_sqlite"]
).sum()

print(
    f"Phenotype value mismatches: "
    f"{phenotype_mismatches}"
)


# ------------------------------------------------------------
# 9. PHENOTYPE AVAILABLE FLAG CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("8. PHENOTYPE AVAILABILITY FLAG")
print("=" * 70)

expected_phenotype_available = (
    comparison["phenotypes_sqlite"] != ""
).astype(int)

actual_phenotype_available = (
    comparison["phenotype_available"]
    .astype(int)
)

availability_mismatches = (
    expected_phenotype_available
    != actual_phenotype_available
).sum()

print(
    f"phenotype_available mismatches: "
    f"{availability_mismatches}"
)


# ------------------------------------------------------------
# 10. SAMPLE PHENOTYPE VALUES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("9. SAMPLE PHENOTYPE VALUES")
print("=" * 70)

sample_values = (
    final_df["phenotypes_clean"]
    .drop_duplicates()
    .head(10)
)

for index, value in enumerate(
    sample_values,
    start=1
):
    print(
        f"{index}. {value}"
    )


# ------------------------------------------------------------
# 11. OVERALL RESULT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("10. OVERALL PHENOTYPE AUDIT")
print("=" * 70)

total_mismatches = (
    phenotype_mismatches
    + availability_mismatches
)

if total_mismatches == 0:

    audit_status = "PASS"

    print("STATUS: PASS")

    print(
        "\nPhenotype values and phenotype availability "
        "flags are consistent with the SQLite source."
    )

else:

    audit_status = "REVIEW REQUIRED"

    print("STATUS: REVIEW REQUIRED")

    print(
        "\nOne or more phenotype consistency checks "
        "failed."
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
        "DEEPGENE STEP 24\n"
        "PHENOTYPE EVIDENCE AUDIT\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Final dataset rows: {len(final_df)}\n"
    )

    f.write(
        f"Phenotype available: "
        f"{available_count}\n"
    )

    f.write(
        f"Phenotype unavailable: "
        f"{unavailable_count}\n"
    )

    f.write(
        f"Empty phenotype values: "
        f"{empty_phenotype_count}\n"
    )

    f.write(
        f"Placeholder/missing-style phenotype values: "
        f"{placeholder_count}\n"
    )

    f.write(
        f"Unique phenotype strings: "
        f"{unique_phenotypes}\n"
    )

    f.write(
        f"Variants with multiple distinct "
        f"phenotype values: "
        f"{multiple_phenotype_values}\n"
    )

    f.write(
        f"Phenotype value mismatches: "
        f"{phenotype_mismatches}\n"
    )

    f.write(
        f"Phenotype availability mismatches: "
        f"{availability_mismatches}\n"
    )

    f.write(
        f"\nSTATUS: {audit_status}\n"
    )


# ------------------------------------------------------------
# CLOSE DATABASE
# ------------------------------------------------------------

conn.close()


# ------------------------------------------------------------
# FINISH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("REPORT SAVED")
print("=" * 70)

print(OUTPUT_FILE)

print("\nStep 24 complete.")