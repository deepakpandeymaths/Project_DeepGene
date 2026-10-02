# ============================================================
# 16. TRANSPARENT BASELINE
# Project: DeepGene
# Focus: SCN1A
#
# Purpose:
# Create a transparent, rule-based interpretation of the
# evidence already present in deepgene_evidence_v1.csv.
#
# IMPORTANT:
# This is NOT a clinical diagnostic model.
# It is NOT an ML model.
# It does NOT assign arbitrary numerical pathogenicity scores.
# ============================================================

import pandas as pd
from pathlib import Path


# ============================================================
# 1. PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_evidence_v1.csv"
)

OUTPUT_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_baseline_v1.csv"
)

RESULTS_DIR = (
    PROJECT_ROOT
    / "data"
    / "analysis_results"
)

RESULTS_PATH = (
    RESULTS_DIR
    / "16_transparent_baseline.txt"
)


# ============================================================
# 2. LOAD DATA
# ============================================================

print("=" * 60)
print("DEEPGENE TRANSPARENT BASELINE")
print("=" * 60)

print(f"\nInput: {INPUT_PATH}")

df = pd.read_csv(INPUT_PATH)

print(f"Variants loaded: {len(df)}")


# ============================================================
# 3. DEFINE EVIDENCE CATEGORY
#
# These categories describe the existing ClinVar evidence.
# They are NOT newly inferred biological labels.
# ============================================================

def classify_evidence(row):

    clinical_significance = str(
        row["clinical_significance"]
    ).strip()

    # --------------------------------------------------------
    # Conflicting evidence
    # --------------------------------------------------------

    if (
        clinical_significance
        == "Conflicting classifications of pathogenicity"
    ):
        return "Conflicting evidence"

    # --------------------------------------------------------
    # Pathogenic evidence
    # --------------------------------------------------------

    if clinical_significance == "Pathogenic":
        return "Pathogenic evidence"

    if clinical_significance == "Likely pathogenic":
        return "Likely pathogenic evidence"

    if clinical_significance == "Pathogenic/Likely pathogenic":
        return "Pathogenic/Likely pathogenic evidence"

    # --------------------------------------------------------
    # Uncertain evidence
    # --------------------------------------------------------

    if clinical_significance == "Uncertain significance":
        return "Uncertain evidence"

    # --------------------------------------------------------
    # Benign evidence
    # --------------------------------------------------------

    if clinical_significance == "Benign":
        return "Benign evidence"

    if clinical_significance == "Likely benign":
        return "Likely benign evidence"

    if clinical_significance == "Benign/Likely benign":
        return "Benign/Likely benign evidence"

    # --------------------------------------------------------
    # Other classifications
    # --------------------------------------------------------

    return "Other / insufficiently classified"


df["evidence_category"] = df.apply(
    classify_evidence,
    axis=1
)


# ============================================================
# 4. DESCRIBE EVIDENCE RELIABILITY
#
# This is descriptive context, not a numerical score.
# ============================================================

def describe_reliability(row):

    if row["expert_panel_review"] == 1:
        return "Expert panel reviewed"

    if row["multiple_submitters"] == 1:
        return "Multiple submitters"

    return "Single submitter / no multiple-submitter evidence"


df["reliability_context"] = df.apply(
    describe_reliability,
    axis=1
)


# ============================================================
# 5. CREATE CONFLICT / REVIEW FLAG
#
# Conflicting evidence should remain visibly identifiable
# rather than being hidden inside another category.
# ============================================================

df["requires_evidence_review"] = (
    df["conflict_flag"] == 1
).astype(int)


# ============================================================
# 6. CREATE HUMAN-READABLE EXPLANATION
# ============================================================

def create_explanation(row):

    category = row["evidence_category"]
    reliability = row["reliability_context"]

    explanation = (
        f"{category}; "
        f"reliability context: {reliability}"
    )

    if row["conflict_flag"] == 1:
        explanation += (
            "; conflicting classification flag is present"
        )

    if row["phenotype_available"] == 1:
        explanation += (
            "; phenotype information is available"
        )

    return explanation


df["evidence_explanation"] = df.apply(
    create_explanation,
    axis=1
)


# ============================================================
# 7. CREATE TRANSPARENT REVIEW GROUP
#
# This is intended for evidence review, not clinical
# diagnosis.
# ============================================================

def review_group(row):

    if row["conflict_flag"] == 1:
        return "Review conflicting evidence"

    if row["evidence_category"] in [
        "Pathogenic evidence",
        "Likely pathogenic evidence",
        "Pathogenic/Likely pathogenic evidence"
    ]:
        return "Pathogenicity evidence present"

    if row["evidence_category"] == "Uncertain evidence":
        return "Uncertain evidence"

    if row["evidence_category"] in [
        "Benign evidence",
        "Likely benign evidence",
        "Benign/Likely benign evidence"
    ]:
        return "Benign evidence"

    return "Other evidence"


df["review_group"] = df.apply(
    review_group,
    axis=1
)


# ============================================================
# 8. SELECT OUTPUT COLUMNS
# ============================================================

output_columns = [
    "variant_id",
    "hgvs_name",

    "clinical_significance",
    "review_status",
    "phenotypes",

    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",

    "phenotype_available",

    "evidence_category",
    "reliability_context",
    "requires_evidence_review",
    "review_group",
    "evidence_explanation"
]

baseline = df[output_columns].copy()


# ============================================================
# 9. SAVE BASELINE DATASET
# ============================================================

baseline.to_csv(
    OUTPUT_PATH,
    index=False
)


# ============================================================
# 10. CREATE ANALYSIS REPORT
# ============================================================

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True
)

results = []

results.append("=" * 60)
results.append("DEEPGENE TRANSPARENT BASELINE")
results.append("=" * 60)

results.append("")
results.append(
    f"Variants analysed: {len(baseline)}"
)

results.append("")
results.append("=" * 60)
results.append("EVIDENCE CATEGORY DISTRIBUTION")
results.append("=" * 60)

category_counts = (
    baseline["evidence_category"]
    .value_counts()
)

for category, count in category_counts.items():

    percentage = (
        count / len(baseline) * 100
    )

    results.append(
        f"{category}: {count} "
        f"({percentage:.2f}%)"
    )


results.append("")
results.append("=" * 60)
results.append("RELIABILITY CONTEXT")
results.append("=" * 60)

reliability_counts = (
    baseline["reliability_context"]
    .value_counts()
)

for category, count in reliability_counts.items():

    percentage = (
        count / len(baseline) * 100
    )

    results.append(
        f"{category}: {count} "
        f"({percentage:.2f}%)"
    )


results.append("")
results.append("=" * 60)
results.append("EVIDENCE REVIEW FLAGS")
results.append("=" * 60)

review_count = (
    baseline["requires_evidence_review"]
    .sum()
)

results.append(
    f"Variants with conflicting evidence: "
    f"{review_count}"
)

results.append(
    f"Variants without conflicting evidence: "
    f"{len(baseline) - review_count}"
)


results.append("")
results.append("=" * 60)
results.append("REVIEW GROUP DISTRIBUTION")
results.append("=" * 60)

review_groups = (
    baseline["review_group"]
    .value_counts()
)

for group, count in review_groups.items():

    percentage = (
        count / len(baseline) * 100
    )

    results.append(
        f"{group}: {count} "
        f"({percentage:.2f}%)"
    )


results.append("")
results.append("=" * 60)
results.append("BASELINE DESIGN")
results.append("=" * 60)

results.append(
    "No arbitrary numerical pathogenicity score was assigned."
)

results.append(
    "Evidence categories directly describe the existing "
    "clinical significance classification."
)

results.append(
    "Reliability context is reported separately from "
    "clinical significance."
)

results.append(
    "Conflicting evidence remains explicitly identifiable."
)

results.append(
    "Phenotype information is preserved for future "
    "phenotype-based analysis."
)

results.append(
    "This baseline is an evidence-prioritization framework, "
    "not a clinical diagnostic system."
)


# ============================================================
# 11. SAVE REPORT
# ============================================================

with open(
    RESULTS_PATH,
    "w",
    encoding="utf-8"
) as file:

    file.write(
        "\n".join(results)
    )


# ============================================================
# 12. PRINT VALIDATION
# ============================================================

print("\n" + "=" * 60)
print("EVIDENCE CATEGORY DISTRIBUTION")
print("=" * 60)

print(category_counts)


print("\n" + "=" * 60)
print("RELIABILITY CONTEXT")
print("=" * 60)

print(reliability_counts)


print("\n" + "=" * 60)
print("REVIEW FLAGS")
print("=" * 60)

print(
    f"Conflicting evidence: {review_count}"
)

print(
    f"No conflicting evidence: "
    f"{len(baseline) - review_count}"
)


print("\n" + "=" * 60)
print("OUTPUT")
print("=" * 60)

print(f"\nBaseline dataset:")
print(OUTPUT_PATH)

print(f"\nAnalysis report:")
print(RESULTS_PATH)


print("\n" + "=" * 60)
print("TRANSPARENT BASELINE COMPLETE")
print("=" * 60)