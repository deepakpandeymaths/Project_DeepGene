# ============================================================
# 40. FUNCTIONAL TARGET DEFINITION + LEAKAGE AUDIT
# ============================================================
#
# DeepGene — SCN1A
#
# PURPOSE
# -------
# This is the final scientific gate before ML.
#
# It determines:
#   1. Which functional measurements are actually usable.
#   2. Whether the matched variants have sufficient target data.
#   3. Whether the qualitative Function field is suitable.
#   4. Which DeepGene features must be excluded because they
#      are derived from ClinVar interpretation/evidence.
#   5. Whether a defensible independent functional target exists.
#
# IMPORTANT:
#   - NO ML model is trained.
#   - NO pathogenic/benign labels are created.
#   - DeepGene V1 is NOT modified.
#
# ============================================================

import os
import re
import pandas as pd
import numpy as np


# ============================================================
# 1. PATHS
# ============================================================

PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)

FUNCTIONAL_FILE = os.path.join(
    PROJECT_ROOT,
    "data",
    "raw",
    "SCN1A_functional_Brunklaus_2020.xlsx"
)

V1_FILE = os.path.join(
    PROJECT_ROOT,
    "data",
    "processed",
    "deepgene_v1_final.csv"
)

OVERLAP_FILE = os.path.join(
    PROJECT_ROOT,
    "data",
    "processed",
    "deepgene_functional_overlap_v1.csv"
)

RESULTS_DIR = os.path.join(
    PROJECT_ROOT,
    "data",
    "analysis_results"
)

PROCESSED_DIR = os.path.join(
    PROJECT_ROOT,
    "data",
    "processed"
)

REPORT_FILE = os.path.join(
    RESULTS_DIR,
    "40_functional_target_definition.txt"
)

TARGET_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_functional_target_v1.csv"
)


# ============================================================
# 2. HELPERS
# ============================================================

def section(title):

    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def clean(value):

    if pd.isna(value):
        return ""

    return str(value).strip()


def numeric_value(value):

    if pd.isna(value):
        return np.nan

    text = str(value).strip()

    if text == "":
        return np.nan

    # Direct numeric
    try:
        return float(text)
    except ValueError:
        pass

    # Numeric with common symbols/units
    match = re.fullmatch(
        r"[<>~=]?\s*(-?\d+(?:\.\d+)?)\s*(?:mV|ms|%)?",
        text,
        flags=re.IGNORECASE
    )

    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return np.nan

    return np.nan


# ============================================================
# 3. START
# ============================================================

os.makedirs(
    RESULTS_DIR,
    exist_ok=True
)

os.makedirs(
    PROCESSED_DIR,
    exist_ok=True
)

section(
    "DEEPGENE — STEP 40 FUNCTIONAL TARGET DEFINITION "
    "+ LEAKAGE AUDIT"
)


# ============================================================
# 4. LOAD FUNCTIONAL DATA
# ============================================================

section("1. LOAD FUNCTIONAL DATASET")

if not os.path.exists(FUNCTIONAL_FILE):

    print("ERROR: Functional dataset not found.")
    print(FUNCTIONAL_FILE)

    raise SystemExit(1)


functional = pd.read_excel(
    FUNCTIONAL_FILE,
    sheet_name="Sheet1"
)

functional = functional.dropna(
    how="all"
).reset_index(
    drop=True
)

print(
    f"Functional workbook rows: "
    f"{len(functional)}"
)


# ============================================================
# 5. LOAD OVERLAP DATASET
# ============================================================

section("2. LOAD VALIDATED OVERLAP DATASET")

if not os.path.exists(OVERLAP_FILE):

    print("ERROR: Step 38–39 overlap file not found.")
    print(OVERLAP_FILE)

    raise SystemExit(1)


overlap = pd.read_csv(
    OVERLAP_FILE
)

print(
    f"Overlap rows: "
    f"{len(overlap)}"
)

matched = overlap[
    overlap["deepgene_match"] == True
].copy()

print(
    f"Matched DeepGene variants: "
    f"{len(matched)}"
)


# ============================================================
# 6. LOAD DEEPGENE V1
# ============================================================

section("3. LOAD DEEPGENE V1")

if not os.path.exists(V1_FILE):

    print("ERROR: DeepGene V1 not found.")
    print(V1_FILE)

    raise SystemExit(1)


v1 = pd.read_csv(
    V1_FILE
)

print(
    f"DeepGene V1 rows: "
    f"{len(v1)}"
)

print(
    f"DeepGene V1 unique variants: "
    f"{v1['variant_id'].nunique()}"
)


# ============================================================
# 7. FUNCTIONAL CLASS TARGET AUDIT
# ============================================================

section("4. QUALITATIVE FUNCTION TARGET")

if "functional_class" not in matched.columns:

    print(
        "ERROR: functional_class column missing."
    )

    raise SystemExit(1)


function_counts = (
    matched["functional_class"]
    .fillna("")
    .astype(str)
    .str.strip()
    .value_counts()
)

print(
    "Functional classification among matched variants:"
)

for category, count in function_counts.items():

    print(
        f"- {category}: {count}"
    )


print()
print(
    "Interpretation:"
)

print(
    "The Function field is experimental functional "
    "classification information."
)

print(
    "It is NOT automatically equivalent to "
    "Pathogenic/Benign."
)

print(
    "Therefore it will be retained as a candidate "
    "experimental target, but not converted into "
    "pathogenicity labels."
)


# ============================================================
# 8. CANDIDATE CONTINUOUS TARGETS
# ============================================================

section("5. CANDIDATE CONTINUOUS FUNCTIONAL TARGETS")


candidate_measurements = [
    "Peak current Density            (% of WT)",
    "Change in Peak current Density         (% of WT)",
    "Persistent Na current (% of peak current)",
    "Persistent Na current Ratio vs WT",
    "WT Activation (mV)",
    "Variant Activation (mV)",
    "Shift of Activation curve (mV)",
    "WT Fast Inactivation (mV)",
    "Variant Fast Inactivation (mV)",
    "Shift of Fast Inactivation (mV)",
    "Recovery from Fast Inactivation",
    "Onset of Slow Inactivation",
    "WT Slow inactivation (mV)",
    "Variant Slow inactivation (mV)",
    "Slow inactivation (mV)",
    "Recovery from Slow inactivation"
]


measurement_summary = []


for column in candidate_measurements:

    if column not in functional.columns:

        print(
            f"{column}: NOT FOUND"
        )

        continue


    numeric_series = (
        functional[column]
        .apply(numeric_value)
    )

    usable = int(
        numeric_series.notna().sum()
    )

    total = len(functional)

    coverage = (
        100 * usable / total
        if total > 0
        else 0
    )

    if usable > 0:

        minimum = numeric_series.min()
        maximum = numeric_series.max()
        median = numeric_series.median()

    else:

        minimum = np.nan
        maximum = np.nan
        median = np.nan


    measurement_summary.append(
        {
            "measurement": column,
            "usable_values": usable,
            "missing_values": total - usable,
            "coverage_percent": round(
                coverage,
                2
            ),
            "minimum": minimum,
            "maximum": maximum,
            "median": median
        }
    )

    print()
    print(
        f"Measurement: {column}"
    )

    print(
        f"  Usable numeric values: {usable}"
    )

    print(
        f"  Missing/unusable: {total - usable}"
    )

    print(
        f"  Coverage: {coverage:.2f}%"
    )

    if usable > 0:

        print(
            f"  Range: {minimum} to {maximum}"
        )

        print(
            f"  Median: {median}"
        )


# ============================================================
# 9. MATCHED VARIANT TARGET COVERAGE
# ============================================================

section(
    "6. TARGET COVERAGE AMONG MATCHED VARIANTS"
)

# Recreate functional lookup using protein representation
# from the overlap dataset.

target_columns = [
    item["measurement"]
    for item in measurement_summary
]


# Build functional lookup
functional_lookup = {}

for _, row in functional.iterrows():

    protein = clean(
        row.get("Variant", "")
    )

    if not protein:
        continue

    # Normalize same basic representation used by overlap
    protein_norm = protein.upper()

    functional_lookup[
        protein_norm
    ] = row


print(
    f"Matched variants available for assessment: "
    f"{len(matched)}"
)


# ============================================================
# 10. TARGET COVERAGE USING MATCH TABLE
# ============================================================

coverage_rows = []

for column in target_columns:

    if column not in functional.columns:

        continue

    values = []

    for _, row in matched.iterrows():

        protein = clean(
            row.get(
                "functional_protein",
                ""
            )
        )

        source_row = functional[
            functional["Variant"]
            .astype(str)
            .str.strip()
            == protein
        ]

        if len(source_row) == 0:
            continue

        value = numeric_value(
            source_row.iloc[0][column]
        )

        values.append(value)


    values = pd.Series(values)

    usable = int(
        values.notna().sum()
    )

    coverage = (
        100 * usable / len(matched)
        if len(matched) > 0
        else 0
    )

    coverage_rows.append(
        {
            "measurement": column,
            "matched_variants_with_value": usable,
            "matched_variants": len(matched),
            "coverage_percent": round(
                coverage,
                2
            )
        }
    )

    print(
        f"- {column}: "
        f"{usable}/{len(matched)} "
        f"({coverage:.2f}%)"
    )


# ============================================================
# 11. PRIMARY TARGET CANDIDATE
# ============================================================

section("7. PRIMARY TARGET CANDIDATE")

# We do NOT automatically choose a target based solely
# on coverage. We identify the most complete experimentally
# measured quantities for scientific review.

coverage_df = pd.DataFrame(
    coverage_rows
)

if len(coverage_df) > 0:

    coverage_df = coverage_df.sort_values(
        "matched_variants_with_value",
        ascending=False
    )

    print(
        "Candidate measurements ranked by usable "
        "matched-variant coverage:"
    )

    for _, row in coverage_df.iterrows():

        print(
            f"- {row['measurement']}: "
            f"{row['matched_variants_with_value']}/"
            f"{row['matched_variants']} "
            f"({row['coverage_percent']}%)"
        )


# ============================================================
# 12. LEAKAGE AUDIT — DEEPGENE FEATURES
# ============================================================

section(
    "8. DEEPGENE FEATURE LEAKAGE AUDIT"
)

print(
    "DeepGene V1 columns:"
)

for column in v1.columns:

    print(
        f"- {column}"
    )


# ------------------------------------------------------------
# ClinVar-derived features
# ------------------------------------------------------------

clinvar_derived = [
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
    "conflict_flag"
]


print()
print(
    "ClinVar-derived features requiring exclusion "
    "from a functional-effect model:"
)

for column in clinvar_derived:

    if column in v1.columns:

        print(
            f"- {column}"
        )


# ------------------------------------------------------------
# Identity features
# ------------------------------------------------------------

identity_features = [
    "variant_id",
    "allele_id",
    "gene_id",
    "hgvs_name",
    "variation_id"
]


print()
print(
    "Identifier features:"
)

for column in identity_features:

    if column in v1.columns:

        print(
            f"- {column}"
        )


# ------------------------------------------------------------
# Representation / QC
# ------------------------------------------------------------

representation_features = [
    "grch37_present",
    "grch38_present",
    "dbsnp_present",
    "genomic_coordinates_present",
    "hgvs_present"
]


print()
print(
    "Representation/QC features:"
)

for column in representation_features:

    if column in v1.columns:

        print(
            f"- {column}"
        )


# ------------------------------------------------------------
# Phenotype
# ------------------------------------------------------------

phenotype_features = [
    "phenotypes",
    "phenotype_available"
]


print()
print(
    "Phenotype features:"
)

for column in phenotype_features:

    if column in v1.columns:

        print(
            f"- {column}"
        )


# ============================================================
# 13. PROPOSED SAFE FEATURE SET
# ============================================================

section(
    "9. PROPOSED SAFE FEATURE SET FOR FUNCTIONAL ML"
)

print()
print(
    "For a functional-effect prediction model, "
    "ClinVar-derived interpretation/evidence fields "
    "should not be used as predictors."
)

print()
print(
    "Identifier fields should also not be treated "
    "as biological predictors."
)

print()
print(
    "Representation/QC fields are metadata/context "
    "and should not automatically be biological features."
)

print()
print(
    "Phenotype information requires separate encoding "
    "and should not be included until its representation "
    "and leakage implications are explicitly resolved."
)

print()
print(
    "Therefore, based on the current V1 dataset alone, "
    "there is NOT YET a sufficiently independent "
    "biological feature matrix."
)


# ============================================================
# 14. IMPORTANT CONSEQUENCE
# ============================================================

section(
    "10. CRITICAL CONSEQUENCE"
)

print()
print(
    "The current DeepGene V1 dataset was deliberately "
    "built from ClinVar evidence."
)

print()
print(
    "The Brunklaus dataset provides an independent "
    "experimental target candidate."
)

print()
print(
    "However, the current V1 feature set contains mostly "
    "ClinVar-derived evidence and metadata."
)

print()
print(
    "Therefore we must NOT simply train:"
)

print(
    "ClinVar evidence features"
)

print(
    "        ↓"
)

print(
    "Brunklaus functional measurement"
)

print()
print(
    "without a separate feature-engineering decision."
)


# ============================================================
# 15. TARGET DECISION
# ============================================================

section(
    "11. STEP 40 TARGET DECISION"
)

print()
print(
    "Independent experimental dataset: YES"
)

print(
    "Validated biological overlap: YES"
)

print(
    f"Matched variants: {len(matched)}"
)

print(
    "Experimental functional measurements: YES"
)

print(
    "Pathogenic/Benign target: NO"
)

print(
    "Independent functional target candidate: YES"
)

print(
    "Final target accepted for ML: NOT YET"
)

print(
    "Reason: feature independence and measurement "
    "harmonization still require explicit finalization."
)


# ============================================================
# 16. NO LABEL CREATION
# ============================================================

section(
    "12. LABEL SAFETY CHECK"
)

print()
print(
    "No artificial pathogenic/benign labels created."
)

print(
    "No ClinVar clinical significance copied into "
    "the functional target."
)

print(
    "No expert-panel label created."
)

print(
    "No review-priority label created."
)

print(
    "No ML model trained."
)


# ============================================================
# 17. SAVE TARGET AUDIT DATA
# ============================================================

section(
    "13. SAVE TARGET AUDIT DATA"
)

coverage_df.to_csv(
    TARGET_FILE,
    index=False
)

print(
    "Saved target audit table:"
)

print(
    TARGET_FILE
)


# ============================================================
# 18. SAVE REPORT
# ============================================================

section(
    "14. SAVE REPORT"
)

with open(
    REPORT_FILE,
    "w",
    encoding="utf-8"
) as report:

    report.write(
        "DEEPGENE — STEP 40 FUNCTIONAL TARGET "
        "DEFINITION + LEAKAGE AUDIT\n"
    )

    report.write("=" * 70 + "\n\n")

    report.write(
        f"Functional workbook variants: "
        f"{len(functional)}\n"
    )

    report.write(
        f"Validated matched variants: "
        f"{len(matched)}\n"
    )

    report.write(
        "\nTARGET CANDIDATES\n"
    )

    report.write("-" * 70 + "\n")

    for _, row in coverage_df.iterrows():

        report.write(
            f"{row['measurement']}: "
            f"{row['matched_variants_with_value']}/"
            f"{row['matched_variants']} "
            f"({row['coverage_percent']}%)\n"
        )

    report.write(
        "\nTARGET DECISION\n"
    )

    report.write("-" * 70 + "\n")

    report.write(
        "Independent experimental dataset: YES\n"
    )

    report.write(
        "Validated biological overlap: YES\n"
    )

    report.write(
        "Independent functional target candidate: YES\n"
    )

    report.write(
        "Pathogenic/Benign target: NO\n"
    )

    report.write(
        "Final ML target accepted: NOT YET\n"
    )

    report.write(
        "\nLEAKAGE\n"
    )

    report.write("-" * 70 + "\n")

    report.write(
        "ClinVar-derived evidence fields require exclusion "
        "from a functional-effect prediction model.\n"
    )

    report.write(
        "Identifiers are not biological predictors.\n"
    )

    report.write(
        "Representation/QC fields are metadata/context.\n"
    )

    report.write(
        "Phenotype representation requires separate review.\n"
    )

    report.write(
        "\nML training: NOT STARTED\n"
    )

    report.write(
        "DeepGene V1 modified: NO\n"
    )


# ============================================================
# 19. FINAL STATUS
# ============================================================

section(
    "15. STEP 40 STATUS"
)

print()
print(
    "Independent experimental dataset: YES"
)

print(
    "Validated biological overlap: YES"
)

print(
    f"Matched variants: {len(matched)}"
)

print(
    "Functional measurements available: YES"
)

print(
    "Pathogenic/Benign labels created: NO"
)

print(
    "Feature leakage audit: COMPLETED"
)

print(
    "Final ML target: NOT YET ACCEPTED"
)

print()
print(
    "🚨 ML HAS NOT STARTED."
)

print()
print(
    "STEP 40 COMPLETE"
)