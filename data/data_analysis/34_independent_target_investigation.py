# ============================================================
# 34. DEEPGENE INDEPENDENT TARGET INVESTIGATION
# ============================================================
#
# Purpose:
#   Document and compare candidate independent targets for the
#   future DeepGene ML phase.
#
# IMPORTANT:
#   This script does NOT:
#       - train ML
#       - create labels
#       - download datasets
#       - modify V1
#
#   It is a methodological decision record.
#
# ============================================================

from pathlib import Path
import pandas as pd


# ============================================================
# 1. PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[1]

PROCESSED_DIR = BASE_DIR / "processed"
ANALYSIS_DIR = BASE_DIR / "analysis_results"

V1_FILE = PROCESSED_DIR / "deepgene_v1_final.csv"

REPORT_FILE = (
    ANALYSIS_DIR /
    "34_independent_target_investigation.txt"
)


# ============================================================
# 2. LOAD CURRENT V1 DATASET
# ============================================================

print("=" * 70)
print("DEEPGENE — INDEPENDENT TARGET INVESTIGATION")
print("=" * 70)

if not V1_FILE.exists():
    raise FileNotFoundError(
        f"Canonical V1 dataset not found: {V1_FILE}"
    )

df = pd.read_csv(
    V1_FILE,
    low_memory=False
)

print("\nCanonical V1 loaded.")
print(f"Rows: {len(df)}")
print(f"Columns: {len(df.columns)}")
print(
    f"Unique variants: "
    f"{df['variant_id'].nunique()}"
)


# ============================================================
# 3. PURPOSE
# ============================================================

print("\n" + "=" * 70)
print("1. PURPOSE OF THIS STEP")
print("=" * 70)

print(
    """
The purpose of Step 34 is to identify a scientifically defensible
prediction target for a future DeepGene machine-learning phase.

The current V1 dataset is ClinVar-based. Therefore, ClinVar-derived
clinical significance should not automatically be used as an
independent ML target while the same ClinVar-derived evidence is
used as model input.
"""
)


# ============================================================
# 4. CANDIDATE TARGET TYPES
# ============================================================

print("\n" + "=" * 70)
print("2. CANDIDATE TARGET TYPES")
print("=" * 70)

candidate_targets = [
    {
        "target": "Pathogenicity classification",
        "example": "Pathogenic vs non-pathogenic",
        "current_status": "NOT READY",
        "main_issue": (
            "Requires an independent target because ClinVar "
            "clinical significance is already part of V1."
        ),
    },
    {
        "target": "Disease association",
        "example": "Disease-associated vs not disease-associated",
        "current_status": "NOT READY",
        "main_issue": (
            "Requires an independently curated disease-association "
            "target with clear provenance."
        ),
    },
    {
        "target": "Functional effect",
        "example": "Normal / partial loss / loss of function",
        "current_status": "PROMISING — INVESTIGATE",
        "main_issue": (
            "Experimental functional datasets exist, but their "
            "size, measurements, variant overlap, and target "
            "definition must be audited."
        ),
    },
    {
        "target": "Continuous functional phenotype",
        "example": "Residual sodium current or other electrophysiological measurement",
        "current_status": "PROMISING — INVESTIGATE",
        "main_issue": (
            "Measurements may be biologically meaningful but "
            "different studies may use different assays and scales."
        ),
    },
    {
        "target": "Evidence-based prioritization",
        "example": "Review priority based on available evidence",
        "current_status": "POSSIBLE WITHOUT ML",
        "main_issue": (
            "Useful as a transparent evidence system, but it is "
            "not an independent biological prediction target."
        ),
    },
]

for item in candidate_targets:

    print(f"\nTarget: {item['target']}")
    print(f"Example: {item['example']}")
    print(f"Status: {item['current_status']}")
    print(f"Issue: {item['main_issue']}")


# ============================================================
# 5. CANDIDATE EXTERNAL SOURCES
# ============================================================

print("\n" + "=" * 70)
print("3. CANDIDATE EXTERNAL SOURCES")
print("=" * 70)

sources = [
    {
        "source": "ClinGen Epilepsy Sodium Channel VCEP",
        "type": "Expert-curated variant interpretation",
        "potential_use": (
            "Reference / validation / evidence provenance"
        ),
        "independent_target": "NOT AUTOMATICALLY",
        "reason": (
            "ClinGen classifications integrate multiple evidence "
            "types and may overlap conceptually or directly with "
            "ClinVar-derived evidence."
        ),
    },
    {
        "source": "Brunklaus et al. 2020 functional dataset",
        "type": "Experimental electrophysiology",
        "potential_use": (
            "Potential functional target"
        ),
        "independent_target": "POTENTIALLY",
        "reason": (
            "Contains experimentally measured SCN1A functional "
            "information rather than simply ClinVar classifications. "
            "Requires variant-level overlap and measurement audit."
        ),
    },
    {
        "source": "ClinGen SCN1A functional-control data",
        "type": "Voltage-clamp functional evidence",
        "potential_use": (
            "Potential functional benchmark / target"
        ),
        "independent_target": "POTENTIALLY",
        "reason": (
            "VCEP documentation reports 23 benign-control and "
            "63 pathogenic/likely-pathogenic control variants with "
            "electrophysiological data. Underlying provenance must "
            "be examined before use."
        ),
    },
    {
        "source": "Recent SCN1A functional studies",
        "type": "Experimental electrophysiology",
        "potential_use": (
            "Supplementary functional evidence"
        ),
        "independent_target": "POTENTIALLY",
        "reason": (
            "Recent studies provide quantitative functional data, "
            "but individual studies may contain relatively few variants."
        ),
    },
    {
        "source": "Historical SCN1A mutation databases",
        "type": "Literature/database compilation",
        "potential_use": (
            "Context / secondary validation"
        ),
        "independent_target": "LIMITED",
        "reason": (
            "These databases combine published classifications, "
            "phenotypes, and functional information rather than "
            "representing a single independent experimental assay."
        ),
    },
]


for item in sources:

    print(f"\nSource: {item['source']}")
    print(f"Type: {item['type']}")
    print(f"Potential use: {item['potential_use']}")
    print(f"Independent target status: {item['independent_target']}")
    print(f"Reason: {item['reason']}")


# ============================================================
# 6. FUNCTIONAL TARGET REQUIREMENTS
# ============================================================

print("\n" + "=" * 70)
print("4. REQUIREMENTS FOR A FUNCTIONAL ML TARGET")
print("=" * 70)

requirements = [
    "Variant identity must be reliably matched to DeepGene V1.",
    "The experimental measurement must be clearly defined.",
    "The assay system must be documented.",
    "The measurement scale must be interpretable.",
    "Wild-type normalization must be understood where applicable.",
    "Different functional assays must not be silently merged.",
    "Repeated measurements must be handled explicitly.",
    "Target labels must not simply reproduce ClinVar classification.",
    "The number of independently measured variants must be sufficient.",
    "Train/test splitting must prevent information leakage.",
]

for item in requirements:
    print(f"- {item}")


# ============================================================
# 7. FUNCTIONAL TARGET FORMS
# ============================================================

print("\n" + "=" * 70)
print("5. POSSIBLE FUNCTIONAL TARGET FORMS")
print("=" * 70)

functional_forms = [
    (
        "Continuous",
        "Predict a measured functional quantity, such as normalized "
        "whole-cell sodium current."
    ),
    (
        "Binary",
        "Predict broadly normal-function vs abnormal-function."
    ),
    (
        "Multiclass",
        "Predict categories such as normal, partial loss, "
        "or severe loss of function."
    ),
    (
        "Phenotype-linked",
        "Predict functional class together with an associated "
        "clinical phenotype."
    ),
]

for name, description in functional_forms:

    print(f"\n{name}")
    print(f"  {description}")


# ============================================================
# 8. WHY WE SHOULD NOT CREATE A LABEL YET
# ============================================================

print("\n" + "=" * 70)
print("6. WHY NO LABEL IS CREATED IN STEP 34")
print("=" * 70)

print(
    """
A functional label cannot safely be created simply by converting
ClinVar categories into:

    Pathogenic = loss of function
    Benign = normal function

SCN1A biology is more complicated than that. Functional effects
can include loss-of-function, partial loss-of-function, gain-of-
function, mixed effects, and phenotype-dependent consequences.

Therefore, the experimental measurements must be inspected before
choosing a target representation.
"""
)


# ============================================================
# 9. CURRENT BEST INVESTIGATION PATH
# ============================================================

print("\n" + "=" * 70)
print("7. CURRENT INVESTIGATION PATH")
print("=" * 70)

investigation_path = [
    "Audit the Brunklaus 2020 functional supplementary dataset.",
    "Audit the ClinGen SCN1A functional-control data.",
    "Normalize variant representations.",
    "Measure exact and carefully validated variant overlap with V1.",
    "Inspect the functional measurement types.",
    "Determine whether measurements can be harmonized.",
    "Determine the number of usable independently measured variants.",
    "Only then decide whether a functional target is suitable."
]

for number, item in enumerate(
    investigation_path,
    start=1
):
    print(f"{number}. {item}")


# ============================================================
# 10. CURRENT DECISION
# ============================================================

print("\n" + "=" * 70)
print("8. CURRENT DECISION")
print("=" * 70)

print(
    """
The functional-effect route is currently the most promising
candidate for investigation because it can provide a biological
measurement that is conceptually different from ClinVar clinical
significance.

However, it is NOT yet accepted as the DeepGene ML target.

The decision remains pending until the underlying functional
datasets have been inspected and variant overlap has been measured.
"""
)


# ============================================================
# 11. WHAT IS NOT ACCEPTED AS TARGET
# ============================================================

print("\n" + "=" * 70)
print("9. TARGETS NOT ACCEPTED AT THIS STAGE")
print("=" * 70)

not_accepted = [
    "ClinVar clinical_significance copied directly as an ML target.",
    "Artificial pathogenic/benign labels created from ClinVar evidence.",
    "Review priority converted into an ML target.",
    "Expert-panel status converted into an ML target.",
    "Phenotype strings converted into arbitrary binary labels.",
]

for item in not_accepted:
    print(f"- {item}")


# ============================================================
# 12. ML BOUNDARY
# ============================================================

print("\n" + "=" * 70)
print("10. ML BOUNDARY")
print("=" * 70)

print(
    """
🚨 ML HAS NOT STARTED.

Step 34 is still target-definition and data-source investigation.

ML will begin only after:

1. An independent target is selected.
2. The target data are acquired.
3. Variant matching is validated.
4. Leakage is audited.
5. Train/validation/test strategy is defined.
6. Evaluation metrics are defined.
"""
)


# ============================================================
# 13. FINAL STATUS
# ============================================================

print("\n" + "=" * 70)
print("FINAL STATUS")
print("=" * 70)

print("PHASE 1 — DATA ANALYSIS: COMPLETE")
print("PHASE 2 — PREDICTION PROBLEM DEFINITION: IN PROGRESS")
print("STEP 34 — INDEPENDENT TARGET INVESTIGATION: COMPLETE")
print("PHASE 3 — MACHINE LEARNING: NOT STARTED")

print(
    "\nNext candidate investigation:"
)

print(
    "Step 35 — Functional dataset acquisition and audit"
)

print(
    "\nNo dataset was modified.")
print("No labels were created.")
print("No ML model was trained.")


# ============================================================
# 14. END
# ============================================================

print("\n" + "=" * 70)
print("STEP 34 COMPLETE")
print("=" * 70)