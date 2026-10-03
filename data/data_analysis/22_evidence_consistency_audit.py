# ============================================================
# 22. EVIDENCE CONSISTENCY AUDIT
# ============================================================
#
# Purpose:
# Audit the relationship between:
#   - unique SCN1A variants
#   - ClinVar evidence records
#   - clinical significance categories
#   - review status
#   - genomic representations
#   - GRCh37 / GRCh38 representations
#
# This is an AUDIT step.
# It does not create ML features.
# It does not change the V1 dataset.
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

OUTPUT_FILE = os.path.join(
    RESULTS_DIR,
    "22_evidence_consistency_audit.txt"
)


# ------------------------------------------------------------
# CONNECT TO DATABASE
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 22")
print("EVIDENCE CONSISTENCY AUDIT")
print("=" * 70)

print("\nDATABASE:")
print(DB_PATH)

conn = sqlite3.connect(DB_PATH)


# ------------------------------------------------------------
# BASIC DATABASE COUNTS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. BASIC DATABASE COUNTS")
print("=" * 70)

genes_count = pd.read_sql_query(
    "SELECT COUNT(*) AS count FROM genes",
    conn
).iloc[0]["count"]

variants_count = pd.read_sql_query(
    "SELECT COUNT(*) AS count FROM variants",
    conn
).iloc[0]["count"]

clinvar_count = pd.read_sql_query(
    "SELECT COUNT(*) AS count FROM clinvar_records",
    conn
).iloc[0]["count"]

representation_count = pd.read_sql_query(
    "SELECT COUNT(*) AS count FROM variant_representations",
    conn
).iloc[0]["count"]

print(f"Genes: {genes_count}")
print(f"Unique variants: {variants_count}")
print(f"ClinVar records: {clinvar_count}")
print(f"Genomic representations: {representation_count}")


# ------------------------------------------------------------
# 2. CLINVAR RECORDS PER VARIANT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("2. CLINVAR RECORDS PER VARIANT")
print("=" * 70)

records_per_variant = pd.read_sql_query(
    """
    SELECT
        variant_id,
        COUNT(*) AS clinvar_record_count
    FROM clinvar_records
    GROUP BY variant_id
    """,
    conn
)

print(
    f"Variants with ClinVar records: "
    f"{len(records_per_variant)}"
)

print(
    f"Maximum ClinVar records for one variant: "
    f"{records_per_variant['clinvar_record_count'].max()}"
)

print(
    f"Average ClinVar records per variant: "
    f"{records_per_variant['clinvar_record_count'].mean():.2f}"
)

print("\nDistribution:")

print(
    records_per_variant["clinvar_record_count"]
    .value_counts()
    .sort_index()
    .to_string()
)


# ------------------------------------------------------------
# 3. MULTIPLE CLINVAR RECORDS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. VARIANTS WITH MULTIPLE CLINVAR RECORDS")
print("=" * 70)

multiple_records = records_per_variant[
    records_per_variant["clinvar_record_count"] > 1
]

print(
    f"Variants with more than one ClinVar record: "
    f"{len(multiple_records)}"
)

print(
    f"Variants with exactly one ClinVar record: "
    f"{(records_per_variant['clinvar_record_count'] == 1).sum()}"
)


# ------------------------------------------------------------
# 4. DISTINCT CLINICAL SIGNIFICANCE CATEGORIES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. CLINICAL SIGNIFICANCE CONSISTENCY")
print("=" * 70)

clinical_categories = pd.read_sql_query(
    """
    SELECT
        variant_id,
        COUNT(DISTINCT clinical_significance)
            AS distinct_clinical_significance_count
    FROM clinvar_records
    GROUP BY variant_id
    """,
    conn
)

multiple_clinical = clinical_categories[
    clinical_categories["distinct_clinical_significance_count"] > 1
]

print(
    "Variants with more than one distinct "
    "clinical significance category: "
    f"{len(multiple_clinical)}"
)

print(
    "Variants with exactly one clinical "
    "significance category: "
    f"{(clinical_categories['distinct_clinical_significance_count'] == 1).sum()}"
)


# ------------------------------------------------------------
# 5. SHOW CLINICAL SIGNIFICANCE COMBINATIONS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. CLINICAL SIGNIFICANCE COMBINATIONS")
print("=" * 70)

clinical_combinations = pd.read_sql_query(
    """
    SELECT
        variant_id,
        GROUP_CONCAT(
            DISTINCT clinical_significance
        ) AS clinical_significance_categories
    FROM clinvar_records
    GROUP BY variant_id
    """,
    conn
)

clinical_combinations["clinical_significance_categories"] = (
    clinical_combinations["clinical_significance_categories"]
    .fillna("")
)

multi_category_combinations = clinical_combinations[
    clinical_combinations[
        "clinical_significance_categories"
    ].str.contains(",", regex=False)
]

print(
    f"Variants with multiple distinct clinical "
    f"significance categories: "
    f"{len(multi_category_combinations)}"
)

if len(multi_category_combinations) > 0:

    print("\nMost common combinations:")

    combination_counts = (
        multi_category_combinations[
            "clinical_significance_categories"
        ]
        .value_counts()
        .head(20)
    )

    print(combination_counts.to_string())


# ------------------------------------------------------------
# 6. REVIEW STATUS CONSISTENCY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("6. REVIEW STATUS CONSISTENCY")
print("=" * 70)

review_status = pd.read_sql_query(
    """
    SELECT
        variant_id,
        COUNT(DISTINCT review_status)
            AS distinct_review_status_count
    FROM clinvar_records
    GROUP BY variant_id
    """,
    conn
)

multiple_review_status = review_status[
    review_status["distinct_review_status_count"] > 1
]

print(
    "Variants with more than one distinct review status: "
    f"{len(multiple_review_status)}"
)

print(
    "Variants with exactly one review status: "
    f"{(review_status['distinct_review_status_count'] == 1).sum()}"
)


# ------------------------------------------------------------
# 7. GENOMIC REPRESENTATIONS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("7. GENOMIC REPRESENTATION AUDIT")
print("=" * 70)

representations_per_variant = pd.read_sql_query(
    """
    SELECT
        variant_id,
        COUNT(*) AS representation_count
    FROM variant_representations
    GROUP BY variant_id
    """,
    conn
)

print(
    f"Variants with genomic representations: "
    f"{len(representations_per_variant)}"
)

print(
    f"Maximum representations for one variant: "
    f"{representations_per_variant['representation_count'].max()}"
)

print(
    f"Average representations per variant: "
    f"{representations_per_variant['representation_count'].mean():.2f}"
)

print("\nRepresentation count distribution:")

print(
    representations_per_variant["representation_count"]
    .value_counts()
    .sort_index()
    .to_string()
)


# ------------------------------------------------------------
# 8. MULTIPLE ASSEMBLIES
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("8. ASSEMBLY REPRESENTATION")
print("=" * 70)

assembly_counts = pd.read_sql_query(
    """
    SELECT
        assembly,
        COUNT(*) AS representation_count
    FROM variant_representations
    GROUP BY assembly
    ORDER BY representation_count DESC
    """,
    conn
)

print(assembly_counts.to_string(index=False))


# ------------------------------------------------------------
# 9. VARIANTS WITH BOTH GRCh37 AND GRCh38
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("9. GRCh37 / GRCh38 COVERAGE")
print("=" * 70)

assembly_presence = pd.read_sql_query(
    """
    SELECT
        variant_id,
        MAX(
            CASE
                WHEN assembly = 'GRCh37'
                THEN 1 ELSE 0
            END
        ) AS grch37_present,

        MAX(
            CASE
                WHEN assembly = 'GRCh38'
                THEN 1 ELSE 0
            END
        ) AS grch38_present

    FROM variant_representations
    GROUP BY variant_id
    """,
    conn
)

both_assemblies = assembly_presence[
    (assembly_presence["grch37_present"] == 1) &
    (assembly_presence["grch38_present"] == 1)
]

only_grch37 = assembly_presence[
    (assembly_presence["grch37_present"] == 1) &
    (assembly_presence["grch38_present"] == 0)
]

only_grch38 = assembly_presence[
    (assembly_presence["grch37_present"] == 0) &
    (assembly_presence["grch38_present"] == 1)
]

neither = assembly_presence[
    (assembly_presence["grch37_present"] == 0) &
    (assembly_presence["grch38_present"] == 0)
]

print(f"Both GRCh37 and GRCh38: {len(both_assemblies)}")
print(f"Only GRCh37: {len(only_grch37)}")
print(f"Only GRCh38: {len(only_grch38)}")
print(f"Neither: {len(neither)}")


# ------------------------------------------------------------
# 10. DATABASE DUPLICATE CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("10. DUPLICATE CHECK")
print("=" * 70)

duplicate_variants = pd.read_sql_query(
    """
    SELECT
        variant_id,
        COUNT(*) AS count
    FROM variants
    GROUP BY variant_id
    HAVING COUNT(*) > 1
    """,
    conn
)

duplicate_clinvar_ids = pd.read_sql_query(
    """
    SELECT
        clinvar_record_id,
        COUNT(*) AS count
    FROM clinvar_records
    GROUP BY clinvar_record_id
    HAVING COUNT(*) > 1
    """,
    conn
)

duplicate_representations = pd.read_sql_query(
    """
    SELECT
        representation_id,
        COUNT(*) AS count
    FROM variant_representations
    GROUP BY representation_id
    HAVING COUNT(*) > 1
    """,
    conn
)

print(
    f"Duplicate variant IDs: "
    f"{len(duplicate_variants)}"
)

print(
    f"Duplicate ClinVar record IDs: "
    f"{len(duplicate_clinvar_ids)}"
)

print(
    f"Duplicate representation IDs: "
    f"{len(duplicate_representations)}"
)


# ------------------------------------------------------------
# 11. CONSISTENCY WITH FINAL V1 DATASET
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("11. CONSISTENCY WITH FINAL V1 DATASET")
print("=" * 70)

final_dataset_path = os.path.join(
    PROCESSED_DIR,
    "deepgene_v1_final.csv"
)

final_df = pd.read_csv(final_dataset_path)

print(
    f"Final V1 dataset rows: "
    f"{len(final_df)}"
)

print(
    f"Final V1 dataset columns: "
    f"{len(final_df.columns)}"
)

print(
    f"Unique variant IDs in final dataset: "
    f"{final_df['variant_id'].nunique()}"
)

database_variant_ids = set(
    pd.read_sql_query(
        "SELECT variant_id FROM variants",
        conn
    )["variant_id"]
)

final_variant_ids = set(
    final_df["variant_id"]
)

missing_from_final = database_variant_ids - final_variant_ids
extra_in_final = final_variant_ids - database_variant_ids

print(
    f"Database variants missing from final dataset: "
    f"{len(missing_from_final)}"
)

print(
    f"Final dataset variants not found in database: "
    f"{len(extra_in_final)}"
)


# ------------------------------------------------------------
# 12. AUDIT CONCLUSION
# ------------------------------------------------------------

audit_pass = True

if len(duplicate_variants) > 0:
    audit_pass = False

if len(duplicate_clinvar_ids) > 0:
    audit_pass = False

if len(duplicate_representations) > 0:
    audit_pass = False

if len(missing_from_final) > 0:
    audit_pass = False

if len(extra_in_final) > 0:
    audit_pass = False


print("\n" + "=" * 70)
print("12. AUDIT CONCLUSION")
print("=" * 70)

if audit_pass:
    print("STATUS: PASS")
    print()
    print(
        "No duplicate primary identifiers were detected, "
        "and the final V1 dataset is consistent with "
        "the database variant set."
    )
else:
    print("STATUS: REVIEW REQUIRED")
    print()
    print(
        "One or more structural consistency checks "
        "require investigation."
    )


# ------------------------------------------------------------
# SAVE REPORT
# ------------------------------------------------------------

with open(OUTPUT_FILE, "w", encoding="utf-8") as f:

    f.write("DEEPGENE STEP 22\n")
    f.write("EVIDENCE CONSISTENCY AUDIT\n")
    f.write("=" * 70 + "\n\n")

    f.write(f"Genes: {genes_count}\n")
    f.write(f"Unique variants: {variants_count}\n")
    f.write(f"ClinVar records: {clinvar_count}\n")
    f.write(
        f"Genomic representations: "
        f"{representation_count}\n"
    )

    f.write("\n")
    f.write(
        f"Variants with multiple ClinVar records: "
        f"{len(multiple_records)}\n"
    )

    f.write(
        f"Variants with multiple clinical significance "
        f"categories: {len(multiple_clinical)}\n"
    )

    f.write(
        f"Variants with multiple review statuses: "
        f"{len(multiple_review_status)}\n"
    )

    f.write(
        f"Variants with both GRCh37 and GRCh38: "
        f"{len(both_assemblies)}\n"
    )

    f.write(
        f"Variants with only GRCh37: "
        f"{len(only_grch37)}\n"
    )

    f.write(
        f"Variants with only GRCh38: "
        f"{len(only_grch38)}\n"
    )

    f.write(
        f"Variants with neither assembly: "
        f"{len(neither)}\n"
    )

    f.write("\n")
    f.write(
        f"Duplicate variant IDs: "
        f"{len(duplicate_variants)}\n"
    )

    f.write(
        f"Duplicate ClinVar record IDs: "
        f"{len(duplicate_clinvar_ids)}\n"
    )

    f.write(
        f"Duplicate representation IDs: "
        f"{len(duplicate_representations)}\n"
    )

    f.write("\n")
    f.write(
        f"Final V1 dataset rows: "
        f"{len(final_df)}\n"
    )

    f.write(
        f"Final V1 unique variants: "
        f"{final_df['variant_id'].nunique()}\n"
    )

    f.write(
        f"Database variants missing from final dataset: "
        f"{len(missing_from_final)}\n"
    )

    f.write(
        f"Final dataset variants not found in database: "
        f"{len(extra_in_final)}\n"
    )

    f.write("\n")
    f.write(
        "AUDIT STATUS: "
        + ("PASS" if audit_pass else "REVIEW REQUIRED")
        + "\n"
    )


# ------------------------------------------------------------
# CLOSE DATABASE
# ------------------------------------------------------------

conn.close()

print("\n" + "=" * 70)
print("REPORT SAVED")
print("=" * 70)

print(OUTPUT_FILE)

print("\nStep 22 complete.")