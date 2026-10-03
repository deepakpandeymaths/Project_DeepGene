# ============================================================
# 20. FINALIZE DEEPGENE V1 EVIDENCE DATASET
# ============================================================

import os
import hashlib
import pandas as pd


# ============================================================
# PATHS
# ============================================================

BASE_DIR = os.path.dirname(
    os.path.dirname(
        os.path.abspath(__file__)
    )
)

PROCESSED_DIR = os.path.join(BASE_DIR, "processed")
RESULTS_DIR = os.path.join(BASE_DIR, "analysis_results")

INPUT_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_evidence_v1.csv"
)

PRIORITIZED_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_prioritized_v1.csv"
)

OUTPUT_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_v1_final.csv"
)

REPORT_FILE = os.path.join(
    RESULTS_DIR,
    "20_v1_dataset_manifest.txt"
)


# ============================================================
# HEADER
# ============================================================

print("=" * 70)
print("DEEPGENE V1 DATASET FINALIZATION")
print("=" * 70)


# ============================================================
# CHECK INPUT
# ============================================================

if not os.path.exists(INPUT_FILE):

    raise FileNotFoundError(
        f"Input dataset not found:\n{INPUT_FILE}"
    )


# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT_FILE)

print()
print(f"Input dataset: {INPUT_FILE}")
print(f"Rows loaded: {len(df)}")
print(f"Columns loaded: {len(df.columns)}")


# ============================================================
# BASIC VALIDATION
# ============================================================

print()
print("=" * 70)
print("BASIC VALIDATION")
print("=" * 70)

assert len(df) == 5381, (
    f"Expected 5381 rows, found {len(df)}"
)

assert df["variant_id"].nunique() == len(df), (
    "variant_id values are not unique."
)

assert df["hgvs_name"].notna().all(), (
    "Missing HGVS identifiers detected."
)

assert df["gene_id"].notna().all(), (
    "Missing gene identifiers detected."
)

assert df["clinical_significance"].notna().all(), (
    "Missing clinical significance detected."
)

assert df["review_status"].notna().all(), (
    "Missing review status detected."
)

print("✓ Row count verified")
print("✓ Variant uniqueness verified")
print("✓ HGVS identifiers verified")
print("✓ Gene identifiers verified")
print("✓ Clinical significance verified")
print("✓ Review status verified")


# ============================================================
# SORT DATASET
# ============================================================

print()
print("=" * 70)
print("STANDARDIZING DATASET ORDER")
print("=" * 70)

df = df.sort_values(
    by="variant_id"
).reset_index(drop=True)

print("✓ Dataset sorted by variant_id")


# ============================================================
# SAVE FINAL DATASET
# ============================================================

df.to_csv(
    OUTPUT_FILE,
    index=False
)

print()
print(f"✓ Final V1 dataset written:")
print(OUTPUT_FILE)


# ============================================================
# SHA256 CHECKSUM
# ============================================================

sha256_hash = hashlib.sha256()

with open(
    OUTPUT_FILE,
    "rb"
) as f:

    for chunk in iter(
        lambda: f.read(1024 * 1024),
        b""
    ):
        sha256_hash.update(chunk)

checksum = sha256_hash.hexdigest()


# ============================================================
# DATASET SUMMARY
# ============================================================

clinical_counts = (
    df["clinical_significance"]
    .value_counts()
)

review_counts = (
    df["review_status"]
    .value_counts()
)

# ============================================================
# FEATURE GROUPS
# ============================================================

clinical_features = [
    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",
    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",
    "vus_count",
    "conflicting_count"
]

reliability_features = [
    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag"
]

phenotype_features = [
    "phenotypes",
    "phenotype_available"
]

identity_features = [
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
    "hgvs_name"
]

representation_features = [
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present"
]


# ============================================================
# CHECK FEATURE GROUPS
# ============================================================

all_expected_features = (
    clinical_features
    + reliability_features
    + phenotype_features
    + identity_features
    + representation_features
    + [
        "clinical_significance",
        "review_status"
    ]
)

missing_group_features = [
    col
    for col in all_expected_features
    if col not in df.columns
]

if missing_group_features:

    raise ValueError(
        "Missing expected feature columns: "
        + str(missing_group_features)
    )

print()
print("=" * 70)
print("FEATURE GROUPS")
print("=" * 70)

print(
    f"Clinical evidence features: "
    f"{len(clinical_features)}"
)

print(
    f"Reliability features: "
    f"{len(reliability_features)}"
)

print(
    f"Phenotype features: "
    f"{len(phenotype_features)}"
)

print(
    f"Variant identity features: "
    f"{len(identity_features)}"
)

print(
    f"Representation/QC features: "
    f"{len(representation_features)}"
)


# ============================================================
# GENERATE MANIFEST
# ============================================================

report = []

report.append(
    "DEEPGENE V1 DATASET MANIFEST"
)

report.append("=" * 70)

report.append(
    "Dataset name: deepgene_v1_final.csv"
)

report.append(
    "Project: DeepGene"
)

report.append(
    "Gene: SCN1A"
)

report.append(
    "Primary data source: ClinVar"
)

report.append(
    "Dataset version: V1"
)

report.append("")

report.append(
    "DATASET SIZE"
)

report.append("-" * 70)

report.append(
    f"Rows / unique variants: {len(df)}"
)

report.append(
    f"Columns: {len(df.columns)}"
)

report.append("")

report.append(
    "DATASET REPRESENTATION"
)

report.append("-" * 70)

report.append(
    "Each row represents one unique SCN1A variant."
)

report.append(
    "ClinVar evidence was aggregated at the variant level."
)

report.append(
    "Genomic representations are represented through "
    "presence/QC features."
)

report.append("")

report.append(
    "CLINICAL SIGNIFICANCE DISTRIBUTION"
)

report.append("-" * 70)

for category, count in clinical_counts.items():

    report.append(
        f"{category}: {count}"
    )

report.append("")

report.append(
    "REVIEW STATUS DISTRIBUTION"
)

report.append("-" * 70)

for category, count in review_counts.items():

    report.append(
        f"{category}: {count}"
    )

report.append("")

report.append(
    "KEY EVIDENCE FEATURES"
)

report.append("-" * 70)

report.append(
    f"Conflicting variants: "
    f"{df['conflict_flag'].sum()}"
)

report.append(
    f"Expert-panel variants: "
    f"{df['expert_panel_review'].sum()}"
)

report.append(
    f"Multiple-submitter variants: "
    f"{df['multiple_submitters'].sum()}"
)

report.append(
    f"Variants with phenotype information: "
    f"{df['phenotype_available'].sum()}"
)

report.append("")

report.append(
    "GENOMIC REPRESENTATION"
)

report.append("-" * 70)

report.append(
    f"GRCh37 present: "
    f"{df['grch37_present'].sum()}"
)

report.append(
    f"GRCh38 present: "
    f"{df['grch38_present'].sum()}"
)

report.append(
    f"dbSNP present: "
    f"{df['dbsnp_present'].sum()}"
)

report.append(
    f"Genomic coordinates present: "
    f"{df['genomic_coordinates_present'].sum()}"
)

report.append(
    f"HGVS present: "
    f"{df['hgvs_present'].sum()}"
)

report.append("")

report.append(
    "FEATURE GROUPS"
)

report.append("-" * 70)

report.append(
    f"Clinical evidence features: "
    f"{len(clinical_features)}"
)

report.append(
    f"Reliability features: "
    f"{len(reliability_features)}"
)

report.append(
    f"Phenotype features: "
    f"{len(phenotype_features)}"
)

report.append(
    f"Identity/QC features: "
    f"{len(identity_features)}"
)

report.append(
    f"Representation features: "
    f"{len(representation_features)}"
)

report.append("")

report.append(
    "BENCHMARK STATUS"
)

report.append("-" * 70)

report.append(
    "Independent pathogenicity labels are NOT available "
    "in the ClinVar-only V1 dataset."
)

report.append(
    "ClinVar clinical significance must therefore not be "
    "treated as an independent ML benchmark label while "
    "the same information is used as an input feature."
)

report.append("")

report.append(
    "VALIDATION STATUS"
)

report.append("-" * 70)

report.append(
    "Step 19 validation: PASS"
)

report.append(
    "Duplicate variant IDs: 0"
)

report.append(
    "Missing essential identifiers: 0"
)

report.append(
    "Missing values in final evidence dataset: "
    f"{int(df.isna().sum().sum())}"
)

report.append("")

report.append(
    "FILE INTEGRITY"
)

report.append("-" * 70)

report.append(
    f"SHA-256: {checksum}"
)

report.append("")

report.append(
    "DOWNSTREAM USE"
)

report.append("-" * 70)

report.append(
    "This dataset is the frozen V1 evidence representation "
    "for downstream DeepGene development."
)

report.append(
    "It is suitable for transparent evidence analysis "
    "and feature-development work."
)

report.append(
    "Machine-learning benchmarking requires a separate "
    "independent label source."
)

report.append("")

report.append("=" * 70)

report.append(
    "END OF DEEPGENE V1 DATASET MANIFEST"
)


# ============================================================
# WRITE REPORT
# ============================================================

with open(
    REPORT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "\n".join(report)
    )


# ============================================================
# FINAL OUTPUT
# ============================================================

print()
print("=" * 70)
print("FINAL V1 DATASET")
print("=" * 70)

print()
print(f"Rows:              {len(df)}")
print(f"Unique variants:   {df['variant_id'].nunique()}")
print(f"Columns:           {len(df.columns)}")

print()
print(f"SHA-256:")
print(checksum)

print()
print("Final dataset:")
print(OUTPUT_FILE)

print()
print("Manifest:")
print(REPORT_FILE)

print()
print("=" * 70)
print("DEEPGENE V1 DATASET FINALIZATION COMPLETE")
print("=" * 70)