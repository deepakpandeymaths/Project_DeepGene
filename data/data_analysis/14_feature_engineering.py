# ============================================================
# 14. FEATURE ENGINEERING
# Project: DeepGene
# Focus: SCN1A
#
# Purpose:
# Build a clean, transparent evidence dataset from the SQLite
# database without inventing or discarding variant information.
# ============================================================

import sqlite3
import pandas as pd
from pathlib import Path


# ============================================================
# 1. PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data" / "database" / "scn1a.db"

OUTPUT_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_evidence_v1.csv"
)


# ============================================================
# 2. CONNECT TO DATABASE
# ============================================================

print("=" * 60)
print("DEEPGENE FEATURE ENGINEERING")
print("=" * 60)

print(f"\nDatabase: {DB_PATH}")

conn = sqlite3.connect(DB_PATH)


# ============================================================
# 3. LOAD VARIANTS
# ============================================================

variants = pd.read_sql_query(
    """
    SELECT
        variant_id,
        gene_id,
        allele_id,
        variation_id,
        hgvs_name
    FROM variants
    """,
    conn
)

print(f"\nVariants loaded: {len(variants)}")


# ============================================================
# 4. LOAD CLINVAR EVIDENCE
# ============================================================

clinvar = pd.read_sql_query(
    """
    SELECT
        variant_id,
        clinical_significance,
        review_status,
        number_submitters,
        phenotypes
    FROM clinvar_records
    """,
    conn
)

print(f"Raw ClinVar rows loaded: {len(clinvar)}")


# ============================================================
# 5. REMOVE DUPLICATE EVIDENCE
#
# The same ClinVar evidence can appear multiple times because
# the same variant has multiple genomic representations
# (e.g. GRCh37 and GRCh38).
#
# We therefore deduplicate evidence using the actual evidence
# fields rather than simply dropping variants.
# ============================================================

evidence_columns = [
    "variant_id",
    "clinical_significance",
    "review_status",
    "number_submitters",
    "phenotypes"
]

clinvar = clinvar.drop_duplicates(
    subset=evidence_columns
).reset_index(drop=True)

print(
    f"Unique ClinVar evidence rows after deduplication: "
    f"{len(clinvar)}"
)


# ============================================================
# 6. NORMALIZE CLINICAL SIGNIFICANCE
# ============================================================

clinvar["clinical_significance"] = (
    clinvar["clinical_significance"]
    .fillna("")
    .astype(str)
    .str.strip()
)

clinvar["review_status"] = (
    clinvar["review_status"]
    .fillna("")
    .astype(str)
    .str.strip()
)

clinvar["phenotypes"] = (
    clinvar["phenotypes"]
    .fillna("")
    .astype(str)
    .str.strip()
)


# ============================================================
# 7. CREATE CLASSIFICATION FLAGS
# ============================================================

clinvar["is_pathogenic"] = (
    clinvar["clinical_significance"]
    == "Pathogenic"
)

clinvar["is_likely_pathogenic"] = (
    clinvar["clinical_significance"]
    == "Likely pathogenic"
)

clinvar["is_pathogenic_likely_pathogenic"] = (
    clinvar["clinical_significance"]
    == "Pathogenic/Likely pathogenic"
)

clinvar["is_benign"] = (
    clinvar["clinical_significance"]
    == "Benign"
)

clinvar["is_likely_benign"] = (
    clinvar["clinical_significance"]
    == "Likely benign"
)

clinvar["is_benign_likely_benign"] = (
    clinvar["clinical_significance"]
    == "Benign/Likely benign"
)

clinvar["is_vus"] = (
    clinvar["clinical_significance"]
    == "Uncertain significance"
)

clinvar["is_conflicting"] = (
    clinvar["clinical_significance"]
    == "Conflicting classifications of pathogenicity"
)


# ============================================================
# 8. AGGREGATE CLINICAL EVIDENCE PER VARIANT
# ============================================================

clinical_features = (
    clinvar
    .groupby("variant_id")
    .agg(
        pathogenic_count=("is_pathogenic", "sum"),
        likely_pathogenic_count=("is_likely_pathogenic", "sum"),
        pathogenic_likely_pathogenic_count=(
            "is_pathogenic_likely_pathogenic",
            "sum"
        ),

        benign_count=("is_benign", "sum"),
        likely_benign_count=("is_likely_benign", "sum"),
        benign_likely_benign_count=(
            "is_benign_likely_benign",
            "sum"
        ),

        vus_count=("is_vus", "sum"),
        conflicting_count=("is_conflicting", "sum"),

        # Use the largest submitter count among the
        # deduplicated evidence assertions.
        submitter_count=("number_submitters", "max"),
    )
    .reset_index()
)


# ============================================================
# 9. PRIMARY CLINVAR CLASSIFICATION
#
# Keep the classification and review status visible in the
# evidence dataset for transparency.
#
# Current database representation has one deduplicated
# evidence row per variant, so these are directly associated
# with each variant.
# ============================================================

primary_evidence = (
    clinvar[
        [
            "variant_id",
            "clinical_significance",
            "review_status",
            "phenotypes"
        ]
    ]
    .drop_duplicates("variant_id")
)


# ============================================================
# 10. RELIABILITY FEATURES
# ============================================================

reliability_features = (
    clinvar
    .groupby("variant_id")
    .agg(
        multiple_submitters=(
            "number_submitters",
            lambda x: int((x > 1).any())
        ),

        expert_panel_review=(
            "review_status",
            lambda x: int(
                x.astype(str)
                .str.contains(
                    "expert panel",
                    case=False,
                    na=False
                )
                .any()
            )
        ),

        conflict_flag=(
            "clinical_significance",
            lambda x: int(
                x.astype(str)
                .str.contains(
                    "conflicting classifications",
                    case=False,
                    na=False
                )
                .any()
            )
        ),

        phenotype_available=(
            "phenotypes",
            lambda x: int(
                x.astype(str)
                .str.strip()
                .ne("")
                .any()
            )
        )
    )
    .reset_index()
)


# ============================================================
# 11. LOAD GENOMIC REPRESENTATIONS
# ============================================================

representations = pd.read_sql_query(
    """
    SELECT
        variant_id,
        assembly,
        chromosome,
        start,
        stop,
        dbsnp_id
    FROM variant_representations
    """,
    conn
)

print(
    f"Genomic representations loaded: "
    f"{len(representations)}"
)


# ============================================================
# 12. GENOMIC / QC FEATURES
# ============================================================

genomic_features = (
    representations
    .groupby("variant_id")
    .agg(
        grch37_present=(
            "assembly",
            lambda x: int(
                x.astype(str)
                .str.upper()
                .eq("GRCH37")
                .any()
            )
        ),

        grch38_present=(
            "assembly",
            lambda x: int(
                x.astype(str)
                .str.upper()
                .eq("GRCH38")
                .any()
            )
        ),

        dbsnp_present=(
            "dbsnp_id",
            lambda x: int(
                x.astype(str)
                .str.strip()
                .replace(
                    {
                        "": pd.NA,
                        "-1": pd.NA,
                        "nan": pd.NA,
                        "None": pd.NA
                    }
                )
                .notna()
                .any()
            )
        ),

        genomic_coordinates_present=(
            "start",
            lambda x: int(x.notna().any())
        )
    )
    .reset_index()
)


# ============================================================
# 13. MERGE EVERYTHING
# ============================================================

evidence = variants.merge(
    primary_evidence,
    on="variant_id",
    how="left"
)

evidence = evidence.merge(
    clinical_features,
    on="variant_id",
    how="left"
)

evidence = evidence.merge(
    reliability_features,
    on="variant_id",
    how="left"
)

evidence = evidence.merge(
    genomic_features,
    on="variant_id",
    how="left"
)


# ============================================================
# 14. HGVS QUALITY FLAG
# ============================================================

evidence["hgvs_present"] = (
    evidence["hgvs_name"]
    .fillna("")
    .astype(str)
    .str.strip()
    .ne("")
    .astype(int)
)


# ============================================================
# 15. FILL NUMERIC FEATURE MISSING VALUES
# ============================================================

numeric_columns = [
    "pathogenic_count",
    "likely_pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "benign_count",
    "likely_benign_count",
    "benign_likely_benign_count",
    "vus_count",
    "conflicting_count",
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
    "phenotype_available",
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present"
]

for column in numeric_columns:
    evidence[column] = evidence[column].fillna(0)


# ============================================================
# 16. IMPORTANT:
# total_clinvar_records IS NOT INCLUDED
#
# After evidence deduplication, the current database contains
# one unique evidence row per variant.
#
# Therefore total_clinvar_records would be constant (=1)
# and carries no useful information for downstream modelling.
# ============================================================


# ============================================================
# 17. ORDER COLUMNS
# ============================================================

column_order = [
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
    "hgvs_name",

    # Primary evidence
    "clinical_significance",
    "review_status",
    "phenotypes",

    # Clinical evidence
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",

    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",

    "vus_count",
    "conflicting_count",

    # Reliability
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
    "phenotype_available",

    # Genomic / QC
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present"
]

evidence = evidence[column_order]


# ============================================================
# 18. SAVE
# ============================================================

evidence.to_csv(
    OUTPUT_PATH,
    index=False
)


# ============================================================
# 19. VALIDATION
# ============================================================

print("\n" + "=" * 60)
print("VALIDATION")
print("=" * 60)

print(f"\nUnique variants: {evidence['variant_id'].nunique()}")
print(f"Rows in output: {len(evidence)}")
print(f"Columns in output: {len(evidence.columns)}")

print("\nClinical significance distribution:")

print(
    evidence["clinical_significance"]
    .value_counts(dropna=False)
)


print("\nReliability feature distribution:")

for column in [
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag"
]:
    print(f"\n{column}:")
    print(evidence[column].value_counts().sort_index())


print("\nPhenotype availability:")
print(
    evidence["phenotype_available"]
    .value_counts()
    .sort_index()
)


print("\nGenomic/QC features:")

for column in [
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present"
]:
    print(
        f"{column}: "
        f"{evidence[column].sum()} present / "
        f"{len(evidence)} total"
    )


print("\nFirst 5 rows:")
print(
    evidence.head(5).to_string(index=False)
)


print("\n" + "=" * 60)
print("FEATURE ENGINEERING COMPLETE")
print("=" * 60)

print(f"\nOutput saved to:")
print(OUTPUT_PATH)


conn.close()