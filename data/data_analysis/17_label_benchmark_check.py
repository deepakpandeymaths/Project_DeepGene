# ============================================================
# 17. LABEL / BENCHMARK FEASIBILITY CHECK
# Project: DeepGene
# Focus: SCN1A
#
# Purpose:
# Determine whether the current ClinVar-only dataset contains
# enough independent information to construct benchmark labels
# without leakage.
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


# ============================================================
# 2. LOAD DATA
# ============================================================

print("=" * 60)
print("DEEPGENE LABEL / BENCHMARK FEASIBILITY CHECK")
print("=" * 60)

print(f"\nInput: {INPUT_PATH}")

df = pd.read_csv(INPUT_PATH)

print(f"Variants loaded: {len(df)}")


# ============================================================
# 3. AVAILABLE CANDIDATE LABEL INFORMATION
# ============================================================

print("\n" + "=" * 60)
print("AVAILABLE CLINVAR CLASSIFICATION INFORMATION")
print("=" * 60)

print(
    df["clinical_significance"]
    .value_counts(dropna=False)
)


# ============================================================
# 4. REVIEW STATUS
# ============================================================

print("\n" + "=" * 60)
print("REVIEW STATUS")
print("=" * 60)

print(
    df["review_status"]
    .value_counts(dropna=False)
)


# ============================================================
# 5. CHECK WHETHER AN INDEPENDENT LABEL SOURCE EXISTS
# ============================================================

print("\n" + "=" * 60)
print("INDEPENDENT LABEL SOURCE CHECK")
print("=" * 60)

print(
    "Current V1 data source: ClinVar"
)

print(
    "Independent external pathogenicity labels: NOT PRESENT"
)

print(
    "Therefore, ClinVar clinical significance cannot be "
    "treated as an independent ML benchmark label while "
    "the same information is also used as an input feature."
)


# ============================================================
# 6. CHECK PHENOTYPE INFORMATION
# ============================================================

print("\n" + "=" * 60)
print("PHENOTYPE INFORMATION")
print("=" * 60)

phenotype_available = (
    df["phenotype_available"]
    .sum()
)

print(
    f"Variants with phenotype information: "
    f"{phenotype_available}"
)

print(
    f"Variants without phenotype information: "
    f"{len(df) - phenotype_available}"
)


# ============================================================
# 7. CHECK REVIEW-BASED SUBSETS
# ============================================================

print("\n" + "=" * 60)
print("REVIEW-BASED SUBSETS")
print("=" * 60)

expert_panel = (
    df["expert_panel_review"]
    .sum()
)

multiple_submitters = (
    df["multiple_submitters"]
    .sum()
)

conflicting = (
    df["conflict_flag"]
    .sum()
)

print(
    f"Expert panel reviewed: {expert_panel}"
)

print(
    f"Multiple submitters: {multiple_submitters}"
)

print(
    f"Conflicting classifications: {conflicting}"
)


# ============================================================
# 8. BENCHMARK FEASIBILITY CONCLUSION
# ============================================================

print("\n" + "=" * 60)
print("BENCHMARK CONCLUSION")
print("=" * 60)

print(
    "STATUS: INDEPENDENT LABELS NOT AVAILABLE IN V1"
)

print(
    "\nReason:"
)

print(
    "The current dataset is ClinVar-only."
)

print(
    "ClinVar clinical significance is already part of "
    "the evidence representation."
)

print(
    "Using the same classification as both input evidence "
    "and target label would introduce target leakage."
)

print(
    "\nTherefore:"
)

print(
    "1. Preserve the current evidence dataset."
)

print(
    "2. Do not create artificial independent labels."
)

print(
    "3. Keep ML benchmarking separate until an appropriate "
    "independent label source is available."
)

print(
    "\nFor V1, DeepGene can continue with transparent "
    "evidence representation and rule-based analysis."
)

print("\n" + "=" * 60)
print("LABEL / BENCHMARK CHECK COMPLETE")
print("=" * 60)