# ============================================================
# 18. EVIDENCE PRIORITIZATION
# Project: DeepGene
# Focus: SCN1A
#
# Purpose:
# Create a transparent evidence-review priority using the
# evidence already available in DeepGene V1.
#
# IMPORTANT:
# This is NOT a clinical pathogenicity score.
# This is NOT a probability.
# This is NOT an ML prediction.
#
# The output identifies variants that may deserve earlier
# evidence review based on documented ClinVar evidence and
# reliability/context signals.
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
    / "deepgene_prioritized_v1.csv"
)

RESULTS_DIR = (
    PROJECT_ROOT
    / "data"
    / "analysis_results"
)

RESULTS_PATH = (
    RESULTS_DIR
    / "18_evidence_prioritization.txt"
)


# ============================================================
# 2. LOAD DATA
# ============================================================

print("=" * 60)
print("DEEPGENE EVIDENCE PRIORITIZATION")
print("=" * 60)

print(f"\nInput: {INPUT_PATH}")

df = pd.read_csv(INPUT_PATH)

print(f"Variants loaded: {len(df)}")


# ============================================================
# 3. DEFINE EVIDENCE PRIORITY
#
# IMPORTANT:
# This is a REVIEW PRIORITY CATEGORY.
#
# It is deliberately not called:
# - pathogenicity score
# - risk score
# - probability
#
# The categories are based on explicit evidence states.
# ============================================================

def determine_priority(row):

    clinical = str(
        row["clinical_significance"]
    ).strip()

    conflict = row["conflict_flag"]

    expert_panel = row["expert_panel_review"]

    multiple_submitters = row["multiple_submitters"]

    # --------------------------------------------------------
    # 1. Conflicting evidence
    #
    # Conflicting variants deserve explicit review because
    # their evidence does not agree.
    # --------------------------------------------------------

    if conflict == 1:
        return "HIGH REVIEW PRIORITY"

    # --------------------------------------------------------
    # 2. Expert-panel reviewed evidence
    #
    # Keep this explicitly identifiable.
    # --------------------------------------------------------

    if expert_panel == 1:
        return "HIGH REVIEW PRIORITY"

    # --------------------------------------------------------
    # 3. Pathogenicity evidence with multiple submitters
    # --------------------------------------------------------

    if clinical in [
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic/Likely pathogenic"
    ] and multiple_submitters == 1:

        return "HIGH REVIEW PRIORITY"

    # --------------------------------------------------------
    # 4. Pathogenicity evidence
    # --------------------------------------------------------

    if clinical in [
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic/Likely pathogenic"
    ]:

        return "MEDIUM REVIEW PRIORITY"

    # --------------------------------------------------------
    # 5. Uncertain evidence
    # --------------------------------------------------------

    if clinical == "Uncertain significance":

        return "MEDIUM REVIEW PRIORITY"

    # --------------------------------------------------------
    # 6. Benign evidence
    #
    # These remain in the dataset but are lower priority for
    # pathogenicity-focused evidence review.
    # --------------------------------------------------------

    if clinical in [
        "Benign",
        "Likely benign",
        "Benign/Likely benign"
    ]:

        return "LOW REVIEW PRIORITY"

    # --------------------------------------------------------
    # 7. Other classifications
    # --------------------------------------------------------

    return "REVIEW CLASSIFICATION"


df["review_priority"] = df.apply(
    determine_priority,
    axis=1
)


# ============================================================
# 4. CREATE PRIORITY REASON
# ============================================================

def determine_reason(row):

    clinical = str(
        row["clinical_significance"]
    ).strip()

    reasons = []

    if row["conflict_flag"] == 1:
        reasons.append(
            "conflicting classifications"
        )

    if row["expert_panel_review"] == 1:
        reasons.append(
            "expert panel review"
        )

    if (
        row["multiple_submitters"] == 1
        and row["submitter_count"] > 1
    ):
        reasons.append(
            f"{int(row['submitter_count'])} submitters"
        )

    if clinical in [
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic/Likely pathogenic"
    ]:
        reasons.append(
            f"ClinVar classification: {clinical}"
        )

    elif clinical == "Uncertain significance":
        reasons.append(
            "ClinVar classification: uncertain significance"
        )

    elif clinical in [
        "Benign",
        "Likely benign",
        "Benign/Likely benign"
    ]:
        reasons.append(
            f"ClinVar classification: {clinical}"
        )

    if row["phenotype_available"] == 1:
        reasons.append(
            "phenotype information available"
        )

    if len(reasons) == 0:
        return "No specific prioritization signal identified"

    return "; ".join(reasons)


df["priority_reason"] = df.apply(
    determine_reason,
    axis=1
)


# ============================================================
# 5. CREATE EXPLICIT REVIEW FLAGS
# ============================================================

df["conflicting_evidence_review"] = (
    df["conflict_flag"] == 1
).astype(int)

df["expert_panel_evidence"] = (
    df["expert_panel_review"] == 1
).astype(int)

df["multiple_submitter_evidence"] = (
    df["multiple_submitters"] == 1
).astype(int)


# ============================================================
# 6. SORT FOR REVIEW
#
# Order:
#   1. High review priority
#   2. Medium review priority
#   3. Low review priority
#
# Within each group, variants with more submitters appear
# first as an evidence-context ordering.
#
# This is an ordering for review, NOT a pathogenicity ranking.
# ============================================================

priority_order = {
    "HIGH REVIEW PRIORITY": 1,
    "MEDIUM REVIEW PRIORITY": 2,
    "LOW REVIEW PRIORITY": 3,
    "REVIEW CLASSIFICATION": 4
}

df["priority_order"] = (
    df["review_priority"]
    .map(priority_order)
)


df = df.sort_values(
    by=[
        "priority_order",
        "submitter_count",
        "variant_id"
    ],
    ascending=[
        True,
        False,
        True
    ]
).reset_index(drop=True)


# ============================================================
# 7. SELECT OUTPUT COLUMNS
# ============================================================

output_columns = [
    "variant_id",
    "hgvs_name",

    "clinical_significance",
    "review_status",

    "submitter_count",
    "multiple_submitters",
    "expert_panel_review",
    "conflict_flag",

    "phenotype_available",

    "review_priority",
    "priority_reason",

    "conflicting_evidence_review",
    "expert_panel_evidence",
    "multiple_submitter_evidence"
]

prioritized = df[output_columns].copy()


# ============================================================
# 8. SAVE PRIORITIZED DATASET
# ============================================================

prioritized.to_csv(
    OUTPUT_PATH,
    index=False
)


# ============================================================
# 9. CREATE REPORT
# ============================================================

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True
)

results = []

results.append("=" * 60)
results.append("DEEPGENE EVIDENCE PRIORITIZATION")
results.append("=" * 60)

results.append("")
results.append(
    f"Variants analysed: {len(prioritized)}"
)

results.append("")
results.append("=" * 60)
results.append("REVIEW PRIORITY DISTRIBUTION")
results.append("=" * 60)

priority_counts = (
    prioritized["review_priority"]
    .value_counts()
)

for priority, count in priority_counts.items():

    percentage = (
        count / len(prioritized) * 100
    )

    results.append(
        f"{priority}: {count} "
        f"({percentage:.2f}%)"
    )


results.append("")
results.append("=" * 60)
results.append("CONFLICTING EVIDENCE")
results.append("=" * 60)

conflicting_count = (
    prioritized["conflicting_evidence_review"]
    .sum()
)

results.append(
    f"Variants with conflicting evidence: "
    f"{conflicting_count}"
)


results.append("")
results.append("=" * 60)
results.append("EXPERT PANEL EVIDENCE")
results.append("=" * 60)

expert_count = (
    prioritized["expert_panel_evidence"]
    .sum()
)

results.append(
    f"Variants with expert panel review: "
    f"{expert_count}"
)


results.append("")
results.append("=" * 60)
results.append("MULTIPLE SUBMITTER EVIDENCE")
results.append("=" * 60)

multiple_count = (
    prioritized["multiple_submitter_evidence"]
    .sum()
)

results.append(
    f"Variants with multiple submitters: "
    f"{multiple_count}"
)


results.append("")
results.append("=" * 60)
results.append("DESIGN NOTES")
results.append("=" * 60)

results.append(
    "This output is a review-priority framework."
)

results.append(
    "It is not a pathogenicity probability."
)

results.append(
    "It is not a clinical diagnostic classification."
)

results.append(
    "No arbitrary numerical pathogenicity weights were used."
)

results.append(
    "Conflicting classifications are explicitly prioritized "
    "for review."
)

results.append(
    "Expert-panel evidence is explicitly retained."
)

results.append(
    "Multiple-submitter evidence is retained as reliability "
    "context."
)

results.append(
    "Benign evidence remains in the dataset and is not "
    "discarded."
)


# ============================================================
# 10. SAVE REPORT
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
# 11. PRINT VALIDATION
# ============================================================

print("\n" + "=" * 60)
print("REVIEW PRIORITY DISTRIBUTION")
print("=" * 60)

print(priority_counts)


print("\n" + "=" * 60)
print("KEY REVIEW FLAGS")
print("=" * 60)

print(
    f"Conflicting evidence: {conflicting_count}"
)

print(
    f"Expert panel review: {expert_count}"
)

print(
    f"Multiple submitters: {multiple_count}"
)


print("\n" + "=" * 60)
print("TOP 20 REVIEW ITEMS")
print("=" * 60)

print(
    prioritized[
        [
            "variant_id",
            "hgvs_name",
            "clinical_significance",
            "submitter_count",
            "review_priority",
            "priority_reason"
        ]
    ]
    .head(20)
    .to_string(index=False)
)


print("\n" + "=" * 60)
print("OUTPUT")
print("=" * 60)

print("\nPrioritized dataset:")
print(OUTPUT_PATH)

print("\nAnalysis report:")
print(RESULTS_PATH)


print("\n" + "=" * 60)
print("EVIDENCE PRIORITIZATION COMPLETE")
print("=" * 60)