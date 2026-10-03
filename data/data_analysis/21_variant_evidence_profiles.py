# ============================================================
# 21. GENERATE DEEPGENE VARIANT EVIDENCE PROFILES
# ============================================================

import os
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
    "deepgene_v1_final.csv"
)

OUTPUT_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_variant_profiles_v1.csv"
)

REPORT_FILE = os.path.join(
    RESULTS_DIR,
    "21_variant_evidence_profiles.txt"
)


# ============================================================
# LOAD DATA
# ============================================================

print("=" * 70)
print("DEEPGENE VARIANT EVIDENCE PROFILES")
print("=" * 70)

print()
print(f"Input: {INPUT_FILE}")

if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(
        f"Input dataset not found:\n{INPUT_FILE}"
    )

df = pd.read_csv(INPUT_FILE)

print()
print(f"Variants loaded: {len(df)}")
print(f"Columns loaded: {len(df.columns)}")


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def has_pathogenic_evidence(row):
    return (
        row["pathogenic_count"] > 0
        or row["pathogenic_likely_pathogenic_count"] > 0
    )


def has_likely_pathogenic_evidence(row):
    return row["likely_pathogenic_count"] > 0


def has_benign_evidence(row):
    return (
        row["benign_count"] > 0
        or row["benign_likely_benign_count"] > 0
    )


def has_likely_benign_evidence(row):
    return row["likely_benign_count"] > 0


def has_vus_evidence(row):
    return row["vus_count"] > 0


def has_conflicting_evidence(row):
    return row["conflicting_count"] > 0


# ============================================================
# CREATE EVIDENCE FLAGS
# ============================================================

print()
print("=" * 70)
print("CREATING EVIDENCE FLAGS")
print("=" * 70)

df["has_pathogenic_evidence"] = df.apply(
    has_pathogenic_evidence,
    axis=1
).astype(int)

df["has_likely_pathogenic_evidence"] = df.apply(
    has_likely_pathogenic_evidence,
    axis=1
).astype(int)

df["has_benign_evidence"] = df.apply(
    has_benign_evidence,
    axis=1
).astype(int)

df["has_likely_benign_evidence"] = df.apply(
    has_likely_benign_evidence,
    axis=1
).astype(int)

df["has_vus_evidence"] = df.apply(
    has_vus_evidence,
    axis=1
).astype(int)

df["has_conflicting_evidence"] = df.apply(
    has_conflicting_evidence,
    axis=1
).astype(int)

print("✓ Pathogenic evidence flag")
print("✓ Likely pathogenic evidence flag")
print("✓ Benign evidence flag")
print("✓ Likely benign evidence flag")
print("✓ VUS evidence flag")
print("✓ Conflicting evidence flag")


# ============================================================
# EVIDENCE PROFILE
# ============================================================

print()
print("=" * 70)
print("BUILDING EVIDENCE PROFILE")
print("=" * 70)


def build_evidence_profile(row):

    categories = []

    if row["has_pathogenic_evidence"] == 1:
        categories.append("Pathogenic")

    if row["has_likely_pathogenic_evidence"] == 1:
        categories.append("Likely pathogenic")

    if row["has_benign_evidence"] == 1:
        categories.append("Benign")

    if row["has_likely_benign_evidence"] == 1:
        categories.append("Likely benign")

    if row["has_vus_evidence"] == 1:
        categories.append("Uncertain significance")

    if row["has_conflicting_evidence"] == 1:
        categories.append("Conflicting")

    if len(categories) == 0:
        return "No classified evidence category detected"

    return "; ".join(categories)


df["evidence_profile"] = df.apply(
    build_evidence_profile,
    axis=1
)


# ============================================================
# RELIABILITY PROFILE
# ============================================================

def build_reliability_profile(row):

    features = []

    if row["expert_panel_review"] == 1:
        features.append("Expert panel")

    if row["multiple_submitters"] == 1:
        features.append("Multiple submitters")
    else:
        features.append("Single submitter")

    if row["conflict_flag"] == 1:
        features.append("Conflicting classifications")

    if len(features) == 0:
        return "No additional reliability context"

    return "; ".join(features)


df["reliability_profile"] = df.apply(
    build_reliability_profile,
    axis=1
)


# ============================================================
# REPRESENTATION PROFILE
# ============================================================

def build_representation_profile(row):

    features = []

    if row["hgvs_present"] == 1:
        features.append("HGVS")

    if row["genomic_coordinates_present"] == 1:
        features.append("Genomic coordinates")

    if row["grch37_present"] == 1:
        features.append("GRCh37")

    if row["grch38_present"] == 1:
        features.append("GRCh38")

    if row["dbsnp_present"] == 1:
        features.append("dbSNP")

    return "; ".join(features)


df["representation_profile"] = df.apply(
    build_representation_profile,
    axis=1
)


# ============================================================
# PHENOTYPE PROFILE
# ============================================================

def build_phenotype_profile(row):

    if row["phenotype_available"] == 1:

        phenotype_text = str(row["phenotypes"])

        if phenotype_text.strip() == "":
            return "Phenotype information available"

        return (
            "Phenotype information available: "
            + phenotype_text
        )

    return "No phenotype information available"


df["phenotype_profile"] = df.apply(
    build_phenotype_profile,
    axis=1
)


# ============================================================
# REVIEW CONTEXT
# ============================================================

def build_review_context(row):

    status = str(row["review_status"])

    submitters = int(row["submitter_count"])

    if row["conflict_flag"] == 1:

        return (
            "Conflicting ClinVar classifications; "
            f"{submitters} submitter(s); "
            f"review status: {status}"
        )

    if row["expert_panel_review"] == 1:

        return (
            "Expert-panel review; "
            f"{submitters} submitter(s); "
            f"review status: {status}"
        )

    if row["multiple_submitters"] == 1:

        return (
            f"{submitters} submitter(s); "
            f"review status: {status}"
        )

    return (
        "Single submitter; "
        f"review status: {status}"
    )


df["review_context"] = df.apply(
    build_review_context,
    axis=1
)


# ============================================================
# EVIDENCE SUMMARY
# ============================================================

def build_evidence_summary(row):

    summary = []

    summary.append(
        f"ClinVar classification: "
        f"{row['clinical_significance']}"
    )

    summary.append(
        f"Evidence profile: "
        f"{row['evidence_profile']}"
    )

    summary.append(
        f"Reliability context: "
        f"{row['reliability_profile']}"
    )

    summary.append(
        f"Phenotype: "
        f"{'available' if row['phenotype_available'] == 1 else 'not available'}"
    )

    summary.append(
        f"Representation: "
        f"{row['representation_profile']}"
    )

    return " | ".join(summary)


df["evidence_summary"] = df.apply(
    build_evidence_summary,
    axis=1
)


# ============================================================
# SELECT OUTPUT COLUMNS
# ============================================================

profile_columns = [
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
    "hgvs_name",

    "clinical_significance",
    "review_status",

    "pathogenic_count",
    "pathogenic_likely_pathogenic_count",
    "likely_pathogenic_count",

    "benign_count",
    "benign_likely_benign_count",
    "likely_benign_count",

    "vus_count",
    "conflicting_count",

    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",

    "phenotype_available",
    "phenotypes",

    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present",

    "has_pathogenic_evidence",
    "has_likely_pathogenic_evidence",
    "has_benign_evidence",
    "has_likely_benign_evidence",
    "has_vus_evidence",
    "has_conflicting_evidence",

    "evidence_profile",
    "reliability_profile",
    "representation_profile",
    "phenotype_profile",
    "review_context",
    "evidence_summary"
]

profiles = df[profile_columns].copy()


# ============================================================
# VALIDATION
# ============================================================

print()
print("=" * 70)
print("PROFILE VALIDATION")
print("=" * 70)

assert len(profiles) == 5381

assert (
    profiles["variant_id"].nunique()
    == 5381
)

assert (
    profiles["hgvs_name"].notna().all()
)

assert (
    profiles["evidence_profile"].notna().all()
)

assert (
    profiles["reliability_profile"].notna().all()
)

assert (
    profiles["representation_profile"].notna().all()
)

assert (
    profiles["phenotype_profile"].notna().all()
)

assert (
    profiles["review_context"].notna().all()
)

assert (
    profiles["evidence_summary"].notna().all()
)

print("✓ 5,381 profiles generated")
print("✓ Variant IDs remain unique")
print("✓ Evidence profiles present")
print("✓ Reliability profiles present")
print("✓ Representation profiles present")
print("✓ Phenotype profiles present")
print("✓ Review contexts present")
print("✓ Evidence summaries present")


# ============================================================
# PROFILE DISTRIBUTIONS
# ============================================================

print()
print("=" * 70)
print("EVIDENCE PROFILE DISTRIBUTION")
print("=" * 70)

print(
    profiles["evidence_profile"]
    .value_counts()
    .to_string()
)


print()
print("=" * 70)
print("RELIABILITY PROFILE DISTRIBUTION")
print("=" * 70)

print(
    profiles["reliability_profile"]
    .value_counts()
    .to_string()
)


# ============================================================
# SAVE DATASET
# ============================================================

profiles.to_csv(
    OUTPUT_FILE,
    index=False
)

print()
print("=" * 70)
print("OUTPUT")
print("=" * 70)

print()
print("Variant profiles:")
print(OUTPUT_FILE)


# ============================================================
# GENERATE REPORT
# ============================================================

report = []

report.append(
    "DEEPGENE VARIANT EVIDENCE PROFILES"
)

report.append("=" * 70)

report.append(
    f"Variants profiled: {len(profiles)}"
)

report.append(
    f"Profile columns: {len(profiles.columns)}"
)

report.append("")

report.append(
    "PURPOSE"
)

report.append("-" * 70)

report.append(
    "This dataset converts the validated DeepGene V1 "
    "evidence representation into an explainable "
    "variant-level profile."
)

report.append(
    "It does not assign an arbitrary pathogenicity score."
)

report.append(
    "It preserves the underlying ClinVar evidence "
    "and explicitly separates evidence, reliability, "
    "phenotype, representation, and review context."
)

report.append("")

report.append(
    "EVIDENCE PROFILE DISTRIBUTION"
)

report.append("-" * 70)

for category, count in (
    profiles["evidence_profile"]
    .value_counts()
    .items()
):

    report.append(
        f"{category}: {count}"
    )

report.append("")

report.append(
    "RELIABILITY PROFILE DISTRIBUTION"
)

report.append("-" * 70)

for category, count in (
    profiles["reliability_profile"]
    .value_counts()
    .items()
):

    report.append(
        f"{category}: {count}"
    )

report.append("")

report.append(
    "IMPORTANT METHODOLOGICAL NOTE"
)

report.append("-" * 70)

report.append(
    "The evidence profile is descriptive."
)

report.append(
    "It should not be interpreted as an independent "
    "pathogenicity prediction."
)

report.append(
    "ClinVar clinical significance remains part of "
    "the V1 evidence representation and therefore "
    "cannot simultaneously serve as an independent "
    "benchmark label for machine-learning evaluation."
)

report.append("")

report.append("=" * 70)

report.append(
    "END OF REPORT"
)


with open(
    REPORT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "\n".join(report)
    )


print()
print("Report:")
print(REPORT_FILE)

print()
print("=" * 70)
print("VARIANT EVIDENCE PROFILE GENERATION COMPLETE")
print("=" * 70)