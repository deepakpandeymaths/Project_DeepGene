# ============================================================
# 31. COMPREHENSIVE V1 DATA ANALYSIS AUDIT
# ============================================================
#
# Purpose:
#   Perform a final report-only audit of the DeepGene V1
#   data-analysis pipeline.
#
# This script:
#   - does NOT modify datasets
#   - checks dataset integrity
#   - checks feature completeness
#   - checks redundancy / constant columns
#   - compares derived datasets
#   - audits phenotype representations
#   - checks SQLite consistency
#   - identifies limitations
#   - identifies the canonical V1 analytical dataset
#   - identifies the natural boundary before ML
#
# ============================================================

from pathlib import Path
import sqlite3
import pandas as pd


# ============================================================
# 1. PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[1]

PROCESSED_DIR = BASE_DIR / "processed"
ANALYSIS_DIR = BASE_DIR / "analysis_results"
DATABASE_DIR = BASE_DIR / "database"

FINAL_DATASET = PROCESSED_DIR / "deepgene_v1_final.csv"
EVIDENCE_DATASET = PROCESSED_DIR / "deepgene_evidence_v1.csv"
VARIANT_PROFILE_DATASET = PROCESSED_DIR / "deepgene_variant_profiles_v1.csv"
PHENOTYPE_PROFILE_DATASET = PROCESSED_DIR / "deepgene_phenotype_profiles_v1.csv"
PHENOTYPE_NORMALIZED_DATASET = (
    PROCESSED_DIR / "deepgene_phenotype_normalized_v1.csv"
)

DATABASE = DATABASE_DIR / "scn1a.db"

REPORT_FILE = ANALYSIS_DIR / "31_comprehensive_v1_audit.txt"


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def section(title):
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def load_dataset(path):
    if not path.exists():
        raise FileNotFoundError(f"Dataset not found: {path}")

    return pd.read_csv(path, low_memory=False)


def dataset_summary(name, df):
    print(f"\n{name}")
    print("-" * len(name))

    print(f"Rows: {len(df)}")
    print(f"Columns: {len(df.columns)}")
    print(f"Unique variant IDs: {df['variant_id'].nunique() if 'variant_id' in df.columns else 'N/A'}")
    print(f"Duplicate rows: {df.duplicated().sum()}")

    if "variant_id" in df.columns:
        print(f"Duplicate variant IDs: {df['variant_id'].duplicated().sum()}")

    missing = int(df.isna().sum().sum())
    print(f"Total missing cells: {missing}")


def column_profile(df):
    rows = []

    for col in df.columns:
        rows.append(
            {
                "column": col,
                "dtype": str(df[col].dtype),
                "unique_values": df[col].nunique(dropna=False),
                "missing_values": int(df[col].isna().sum()),
                "constant": df[col].nunique(dropna=False) <= 1,
            }
        )

    return pd.DataFrame(rows)


def compare_variant_ids(name_a, df_a, name_b, df_b):
    if "variant_id" not in df_a.columns or "variant_id" not in df_b.columns:
        print(f"{name_a} vs {name_b}: variant_id unavailable")
        return

    ids_a = set(df_a["variant_id"])
    ids_b = set(df_b["variant_id"])

    missing_from_b = ids_a - ids_b
    extra_in_b = ids_b - ids_a

    print(f"\n{name_a} vs {name_b}")
    print(f"  Variants in {name_a}: {len(ids_a)}")
    print(f"  Variants in {name_b}: {len(ids_b)}")
    print(f"  Missing from {name_b}: {len(missing_from_b)}")
    print(f"  Extra in {name_b}: {len(extra_in_b)}")

    return len(missing_from_b), len(extra_in_b)


# ============================================================
# 3. LOAD DATASETS
# ============================================================

print("=" * 70)
print("DEEPGENE V1 — COMPREHENSIVE DATA ANALYSIS AUDIT")
print("=" * 70)

print("\nLoading datasets...")

final_df = load_dataset(FINAL_DATASET)
evidence_df = load_dataset(EVIDENCE_DATASET)
variant_profiles_df = load_dataset(VARIANT_PROFILE_DATASET)
phenotype_profiles_df = load_dataset(PHENOTYPE_PROFILE_DATASET)
phenotype_normalized_df = load_dataset(PHENOTYPE_NORMALIZED_DATASET)

print("All datasets loaded successfully.")


# ============================================================
# 4. BASIC DATASET INVENTORY
# ============================================================

section("1. DATASET INVENTORY")

datasets = {
    "Canonical V1 dataset": final_df,
    "Evidence dataset": evidence_df,
    "Variant evidence profiles": variant_profiles_df,
    "Phenotype profiles": phenotype_profiles_df,
    "Normalized phenotype dataset": phenotype_normalized_df,
}

for name, df in datasets.items():
    dataset_summary(name, df)


# ============================================================
# 5. CANONICAL V1 INTEGRITY
# ============================================================

section("2. CANONICAL V1 DATASET INTEGRITY")

print(f"Rows: {len(final_df)}")
print(f"Columns: {len(final_df.columns)}")

if "variant_id" in final_df.columns:
    print(f"Unique variants: {final_df['variant_id'].nunique()}")
    print(
        f"Duplicate variant IDs: "
        f"{final_df['variant_id'].duplicated().sum()}"
    )

print(f"Duplicate rows: {final_df.duplicated().sum()}")
print(f"Total missing cells: {int(final_df.isna().sum().sum())}")


# ============================================================
# 6. FINAL V1 COLUMN INVENTORY
# ============================================================

section("3. CANONICAL V1 COLUMN INVENTORY")

for i, col in enumerate(final_df.columns, start=1):
    print(f"{i:02d}. {col}")


# ============================================================
# 7. FEATURE GROUP CLASSIFICATION
# ============================================================

section("4. FEATURE GROUP CLASSIFICATION")

clinical_features = [
    "clinical_significance",
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
]

reliability_features = [
    "review_status",
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
]

phenotype_features = [
    "phenotypes",
    "phenotype_available",
]

identity_features = [
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
    "hgvs_name",
]

representation_features = [
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",
]

feature_groups = {
    "Clinical evidence": clinical_features,
    "Reliability / review context": reliability_features,
    "Phenotype": phenotype_features,
    "Variant identity": identity_features,
    "Representation / QC": representation_features,
}

for group_name, columns in feature_groups.items():
    present = [c for c in columns if c in final_df.columns]
    missing = [c for c in columns if c not in final_df.columns]

    print(f"\n{group_name}")
    print(f"  Present: {len(present)}")
    print(f"  Missing: {len(missing)}")

    if present:
        print("  Columns:")
        for c in present:
            print(f"    - {c}")

    if missing:
        print("  Missing columns:")
        for c in missing:
            print(f"    - {c}")


# ============================================================
# 8. CONSTANT / LOW-VARIATION FEATURES
# ============================================================

section("5. CONSTANT AND LOW-VARIATION FEATURES")

profile = column_profile(final_df)

constant_columns = profile.loc[
    profile["constant"], "column"
].tolist()

print(f"Constant columns: {len(constant_columns)}")

if constant_columns:
    for col in constant_columns:
        print(f"  - {col}")
else:
    print("  None")

print("\nColumns with <= 2 unique values:")

low_variation = profile.loc[
    profile["unique_values"] <= 2,
    ["column", "unique_values", "missing_values"]
]

if len(low_variation) == 0:
    print("  None")
else:
    for _, row in low_variation.iterrows():
        print(
            f"  - {row['column']}: "
            f"{row['unique_values']} unique values, "
            f"{row['missing_values']} missing"
        )


# ============================================================
# 9. IDENTITY VS ANALYTICAL FEATURES
# ============================================================

section("6. IDENTIFIER VS ANALYTICAL FEATURES")

print("Identifier / representation fields are retained for traceability.")
print("They should not automatically be treated as biological predictors.")

print("\nIdentifier fields:")

for col in identity_features:
    if col in final_df.columns:
        print(f"  - {col}")

print("\nRepresentation / QC fields:")

for col in representation_features:
    if col in final_df.columns:
        print(f"  - {col}")

print("\nInterpretation:")
print(
    "variant_id, gene_id, allele_id, and variation_id primarily identify "
    "records or biological entities."
)

print(
    "HGVS and genomic representation fields describe how the variant is "
    "represented and should not automatically be interpreted as biological "
    "pathogenicity predictors."
)


# ============================================================
# 10. CLINICAL EVIDENCE DISTRIBUTION
# ============================================================

section("7. CLINICAL EVIDENCE DISTRIBUTION")

if "clinical_significance" in final_df.columns:
    counts = final_df["clinical_significance"].value_counts(dropna=False)

    for value, count in counts.items():
        print(f"{value}: {count}")

print("\nEvidence count ranges:")

count_columns = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count",
    "submitter_count",
]

for col in count_columns:
    if col in final_df.columns:
        print(
            f"{col}: "
            f"min={final_df[col].min()}, "
            f"max={final_df[col].max()}, "
            f"mean={final_df[col].mean():.3f}"
        )


# ============================================================
# 11. RELIABILITY / REVIEW CONTEXT
# ============================================================

section("8. RELIABILITY AND REVIEW CONTEXT")

for col in [
    "review_status",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",
]:
    if col in final_df.columns:
        print(f"\n{col}")
        print(final_df[col].value_counts(dropna=False).to_string())


# ============================================================
# 12. PHENOTYPE BRANCH AUDIT
# ============================================================

section("9. PHENOTYPE BRANCH")

print("Original phenotype dataset:")
print(f"  Rows: {len(phenotype_profiles_df)}")
print(f"  Columns: {len(phenotype_profiles_df.columns)}")

print("\nNormalized phenotype dataset:")
print(f"  Rows: {len(phenotype_normalized_df)}")
print(f"  Columns: {len(phenotype_normalized_df.columns)}")

if "phenotype_category" in phenotype_profiles_df.columns:
    print("\nOriginal phenotype categories:")
    print(
        phenotype_profiles_df["phenotype_category"]
        .value_counts(dropna=False)
        .to_string()
    )

if "normalization_changed" in phenotype_normalized_df.columns:
    print("\nNormalization status:")
    print(
        phenotype_normalized_df["normalization_changed"]
        .value_counts(dropna=False)
        .to_string()
    )

if "normalized_phenotype_terms" in phenotype_normalized_df.columns:
    normalized_term_counts = (
        phenotype_normalized_df["normalized_phenotype_terms"]
        .fillna("")
        .apply(
            lambda x: 0
            if not str(x).strip()
            else len([t for t in str(x).split("|") if t.strip()])
        )
    )

    print(
        f"\nNormalized phenotype terms across all variants: "
        f"{normalized_term_counts.sum()}"
    )

    print(
        f"Variants with >=1 normalized phenotype term: "
        f"{(normalized_term_counts > 0).sum()}"
    )

    print(
        f"Variants with zero normalized phenotype terms: "
        f"{(normalized_term_counts == 0).sum()}"
    )


# ============================================================
# 13. PHENOTYPE DATASET RELATIONSHIP
# ============================================================

section("10. PHENOTYPE DATASET RELATIONSHIP")

compare_variant_ids(
    "Canonical V1",
    final_df,
    "Normalized phenotype dataset",
    phenotype_normalized_df,
)

compare_variant_ids(
    "Canonical V1",
    final_df,
    "Phenotype profiles",
    phenotype_profiles_df,
)

print(
    "\nThe normalized phenotype dataset is treated as a derived "
    "representation rather than automatically replacing the canonical "
    "26-column V1 dataset."
)


# ============================================================
# 14. CROSS-DATASET VARIANT CONSISTENCY
# ============================================================

section("11. CROSS-DATASET VARIANT CONSISTENCY")

for name, df in datasets.items():
    if "variant_id" not in df.columns:
        print(f"{name}: variant_id missing")
        continue

    ids = set(df["variant_id"])

    canonical_ids = set(final_df["variant_id"])

    missing = canonical_ids - ids
    extra = ids - canonical_ids

    print(f"\n{name}")
    print(f"  Unique variants: {len(ids)}")
    print(f"  Missing vs canonical V1: {len(missing)}")
    print(f"  Extra vs canonical V1: {len(extra)}")


# ============================================================
# 15. COLUMN OVERLAP / REDUNDANCY
# ============================================================

section("12. COLUMN OVERLAP ACROSS DATASETS")

dataset_columns = {
    name: set(df.columns)
    for name, df in datasets.items()
}

for name, columns in dataset_columns.items():
    print(f"\n{name}: {len(columns)} columns")

    overlap = columns & set(final_df.columns)

    print(
        f"  Columns overlapping canonical V1: "
        f"{len(overlap)}"
    )


# ============================================================
# 16. DUPLICATE INFORMATION CHECK
# ============================================================

section("13. EXACT DUPLICATE COLUMN CONTENT")

duplicate_column_pairs = []

columns = list(final_df.columns)

for i in range(len(columns)):
    for j in range(i + 1, len(columns)):

        col_a = columns[i]
        col_b = columns[j]

        try:
            if final_df[col_a].equals(final_df[col_b]):
                duplicate_column_pairs.append((col_a, col_b))
        except Exception:
            pass

print(
    f"Exact duplicate column pairs: "
    f"{len(duplicate_column_pairs)}"
)

if duplicate_column_pairs:
    for a, b in duplicate_column_pairs:
        print(f"  - {a} == {b}")
else:
    print("  None detected.")


# ============================================================
# 17. PHENOTYPE NORMALIZATION IMPACT
# ============================================================

section("14. PHENOTYPE NORMALIZATION IMPACT")

if (
    "normalization_changed" in phenotype_normalized_df.columns
    and "variant_id" in phenotype_normalized_df.columns
):

    changed = phenotype_normalized_df[
        phenotype_normalized_df["normalization_changed"] == True
    ]

    unchanged = phenotype_normalized_df[
        phenotype_normalized_df["normalization_changed"] == False
    ]

    print(f"Changed variants: {len(changed)}")
    print(f"Unchanged variants: {len(unchanged)}")

    print(
        f"Percentage changed: "
        f"{len(changed) / len(phenotype_normalized_df) * 100:.2f}%"
    )

    if "original_phenotypes" in phenotype_normalized_df.columns:
        original_semicolon = (
            phenotype_normalized_df["original_phenotypes"]
            .fillna("")
            .astype(str)
            .str.contains(";", regex=False)
        )

        print(
            f"Original variants containing semicolon: "
            f"{original_semicolon.sum()}"
        )

        print(
            "Normalization change count matches original semicolon "
            "representation:",
            changed.shape[0] == original_semicolon.sum(),
        )


# ============================================================
# 18. SQLITE CROSS-CHECK
# ============================================================

section("15. SQLITE DATABASE CROSS-CHECK")

if not DATABASE.exists():
    print("SQLite database not found.")
else:

    conn = sqlite3.connect(DATABASE)

    try:
        queries = {
            "Genes": "SELECT COUNT(*) FROM genes",
            "Unique variants": "SELECT COUNT(*) FROM variants",
            "ClinVar records": "SELECT COUNT(*) FROM clinvar_records",
            "Genomic representations": (
                "SELECT COUNT(*) FROM variant_representations"
            ),
        }

        db_counts = {}

        for name, query in queries.items():
            value = conn.execute(query).fetchone()[0]
            db_counts[name] = value
            print(f"{name}: {value}")

        print("\nDatabase vs canonical V1:")

        print(
            "  Unique variants:",
            db_counts["Unique variants"],
            "vs",
            len(final_df)
        )

        print(
            "  Variant count match:",
            db_counts["Unique variants"] == len(final_df)
        )

        print(
            "  ClinVar records:",
            db_counts["ClinVar records"]
        )

        print(
            "  Genomic representations:",
            db_counts["Genomic representations"]
        )

    finally:
        conn.close()


# ============================================================
# 19. INFORMATION AVAILABILITY SUMMARY
# ============================================================

section("16. INFORMATION AVAILABLE IN V1")

available_information = {
    "SCN1A variant identity": True,
    "HGVS representation": "hgvs_name" in final_df.columns,
    "GRCh37 representation": "grch37_present" in final_df.columns,
    "GRCh38 representation": "grch38_present" in final_df.columns,
    "dbSNP representation": "dbsnp_present" in final_df.columns,
    "Clinical significance": "clinical_significance" in final_df.columns,
    "Pathogenic evidence counts": "pathogenic_count" in final_df.columns,
    "Likely pathogenic evidence counts": "likely_pathogenic_count" in final_df.columns,
    "Benign evidence counts": "benign_count" in final_df.columns,
    "Likely benign evidence counts": "likely_benign_count" in final_df.columns,
    "VUS evidence counts": "vus_count" in final_df.columns,
    "Conflict evidence": "conflicting_count" in final_df.columns,
    "Submitter count": "submitter_count" in final_df.columns,
    "Multiple submitter context": "multiple_submitters" in final_df.columns,
    "Expert panel context": "expert_panel_review" in final_df.columns,
    "Phenotype information": "phenotypes" in final_df.columns,
    "Phenotype availability flag": "phenotype_available" in final_df.columns,
}

for item, available in available_information.items():
    print(f"{item}: {available}")


# ============================================================
# 20. INFORMATION NOT PRESENT
# ============================================================

section("17. IMPORTANT INFORMATION NOT PRESENT IN V1")

missing_external_information = [
    "Independent population frequency information",
    "Independent functional assay evidence",
    "Independent clinical literature evidence",
    "Independent curated gene-disease evidence beyond ClinVar",
    "Independent phenotype ontology mapping",
    "Independent pathogenicity benchmark labels",
    "Independent external validation dataset",
]

for item in missing_external_information:
    print(f"- {item}")


# ============================================================
# 21. SCIENTIFIC LIMITATIONS
# ============================================================

section("18. SCIENTIFIC LIMITATIONS OF V1")

limitations = [
    (
        "ClinVar-only evidence",
        "V1 is based on ClinVar and therefore does not independently "
        "represent population, functional, literature, or other external evidence."
    ),
    (
        "Clinical significance dependency",
        "Clinical significance is itself ClinVar-derived and should not "
        "be treated as an independent ML benchmark label when ClinVar "
        "evidence is used as an input."
    ),
    (
        "Phenotype incompleteness",
        "Phenotype information is present in the source records but is "
        "not uniformly informative across all variants."
    ),
    (
        "Representation complexity",
        "A biological variant may have multiple genomic representations "
        "and ClinVar records."
    ),
    (
        "No independent validation",
        "V1 does not contain an independent external benchmark dataset."
    ),
    (
        "No clinical decision capability",
        "V1 is an evidence organization and prioritization framework, "
        "not a clinical diagnostic system."
    ),
]

for title, description in limitations:
    print(f"\n{title}")
    print(f"  {description}")


# ============================================================
# 22. REVIEW OF DATA-ANALYSIS COMPLETENESS
# ============================================================

section("19. DATA-ANALYSIS COMPLETENESS REVIEW")

checks = []

checks.append(
    (
        "Canonical V1 exists",
        FINAL_DATASET.exists()
    )
)

checks.append(
    (
        "Canonical V1 has unique variants",
        final_df["variant_id"].nunique() == len(final_df)
        if "variant_id" in final_df.columns
        else False
    )
)

checks.append(
    (
        "Canonical V1 has no duplicate rows",
        final_df.duplicated().sum() == 0
    )
)

checks.append(
    (
        "Evidence dataset aligns with canonical V1",
        set(final_df["variant_id"]) == set(evidence_df["variant_id"])
    )
)

checks.append(
    (
        "Variant profiles align with canonical V1",
        set(final_df["variant_id"])
        == set(variant_profiles_df["variant_id"])
    )
)

checks.append(
    (
        "Phenotype profiles align with canonical V1",
        set(final_df["variant_id"])
        == set(phenotype_profiles_df["variant_id"])
    )
)

checks.append(
    (
        "Normalized phenotype dataset aligns with canonical V1",
        set(final_df["variant_id"])
        == set(phenotype_normalized_df["variant_id"])
    )
)

checks.append(
    (
        "Canonical V1 has no missing cells",
        final_df.isna().sum().sum() == 0
    )
)

for check_name, passed in checks:
    print(f"{'PASS' if passed else 'FAIL'}: {check_name}")


# ============================================================
# 23. NATURAL ML BOUNDARY
# ============================================================

section("20. NATURAL BOUNDARY BEFORE ML")

print(
    "The current V1 data-analysis pipeline has reached a point where "
    "the available ClinVar evidence can be represented, validated, "
    "audited, and organized into variant-level analytical profiles."
)

print("\nBefore supervised ML, the following are still required:")

ml_requirements = [
    "Define an independent prediction target.",
    "Obtain labels that are independent of the input evidence used by the model.",
    "Define train/validation/test separation without variant leakage.",
    "Decide which feature groups are biologically appropriate predictors.",
    "Exclude identifiers and non-predictive representation fields where appropriate.",
    "Define evaluation metrics before training.",
    "Establish a reproducible baseline.",
    "Perform model validation and error analysis.",
]

for item in ml_requirements:
    print(f"- {item}")


# ============================================================
# 24. FINAL ARCHITECTURE DECISION
# ============================================================

section("21. V1 DATA ARCHITECTURE")

print("CANONICAL ANALYTICAL DATASET")
print("--------------------------------")
print("data/processed/deepgene_v1_final.csv")

print("\nSUPPORTING DERIVED DATASETS")
print("--------------------------------")
print("1. deepgene_evidence_v1.csv")
print("2. deepgene_variant_profiles_v1.csv")
print("3. deepgene_phenotype_profiles_v1.csv")
print("4. deepgene_phenotype_normalized_v1.csv")

print("\nDATABASE SOURCE")
print("--------------------------------")
print("data/database/scn1a.db")

print("\nINTERPRETATION")
print("--------------------------------")
print(
    "The canonical V1 dataset remains the 26-column validated evidence "
    "dataset. The normalized phenotype dataset is retained as a separate "
    "derived representation so that the original phenotype representation "
    "is not lost or silently replaced."
)


# ============================================================
# 25. FINAL CONCLUSION
# ============================================================

section("22. FINAL CONCLUSION")

print(
    "DeepGene V1 currently provides a validated SCN1A ClinVar-based "
    "variant evidence representation with clinical significance, "
    "evidence counts, reliability/review context, phenotype information, "
    "variant identity, and genomic representation/QC fields."
)

print(
    "\nThe data-analysis pipeline should not be extended with arbitrary "
    "additional features simply to increase dataset size."
)

print(
    "\nThe next major methodological boundary is the definition of an "
    "independent prediction target and validation strategy."
)

print(
    "\nNo machine-learning model is created by this audit."
)

print(
    "\nThe audit itself is report-only and does not modify any dataset."
)


# ============================================================
# 26. SAVE REPORT
# ============================================================

ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

# Re-run output capture is intentionally avoided.
# The report is generated by redirecting stdout from this script
# externally if desired.

print()
print("=" * 70)
print("AUDIT COMPLETE")
print("=" * 70)
print(f"Report target: {REPORT_FILE}")
print("=" * 70)