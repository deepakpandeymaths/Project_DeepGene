# ============================================================
# 23. EVIDENCE COUNT CONSISTENCY AUDIT
# ============================================================
#
# Purpose:
# Verify that the evidence/count features in
# deepgene_v1_final.csv agree with the underlying
# ClinVar records stored in SQLite.
#
# This is an AUDIT step.
# It does not modify the V1 dataset.
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

INPUT_FILE = os.path.join(
    PROCESSED_DIR,
    "deepgene_v1_final.csv"
)

OUTPUT_FILE = os.path.join(
    RESULTS_DIR,
    "23_evidence_count_consistency.txt"
)


# ------------------------------------------------------------
# START
# ------------------------------------------------------------

print("=" * 70)
print("DEEPGENE STEP 23")
print("EVIDENCE COUNT CONSISTENCY AUDIT")
print("=" * 70)


# ------------------------------------------------------------
# LOAD FINAL DATASET
# ------------------------------------------------------------

print("\nLoading final V1 dataset...")

final_df = pd.read_csv(INPUT_FILE)

print(f"Rows loaded: {len(final_df)}")
print(f"Columns loaded: {len(final_df.columns)}")


# ------------------------------------------------------------
# CONNECT DATABASE
# ------------------------------------------------------------

print("\nConnecting to SQLite database...")

conn = sqlite3.connect(DB_PATH)


# ------------------------------------------------------------
# LOAD CLINVAR RECORDS
# ------------------------------------------------------------

clinvar = pd.read_sql_query(
    """
    SELECT
        variant_id,
        clinical_significance,
        review_status,
        number_submitters
    FROM clinvar_records
    """,
    conn
)

print(
    f"ClinVar records loaded: {len(clinvar)}"
)


# ------------------------------------------------------------
# DEDUPLICATE EVIDENCE
# ------------------------------------------------------------
#
# The database contains one ClinVar record per variant in
# the current representation. We nevertheless keep this
# deduplication logic explicit so that evidence is not
# accidentally counted twice.
# ------------------------------------------------------------

clinvar = clinvar.drop_duplicates(
    subset=[
        "variant_id",
        "clinical_significance",
        "review_status",
        "number_submitters"
    ]
).reset_index(drop=True)

print(
    f"Unique evidence rows after deduplication: "
    f"{len(clinvar)}"
)


# ------------------------------------------------------------
# NORMALIZE CLINICAL SIGNIFICANCE
# ------------------------------------------------------------

clinvar["clinical_significance"] = (
    clinvar["clinical_significance"]
    .fillna("")
    .astype(str)
    .str.strip()
)


# ------------------------------------------------------------
# CALCULATE EXPECTED EVIDENCE COUNTS
# ------------------------------------------------------------

print("\nCalculating expected evidence counts...")


def count_category(df, category):
    return (
        df["clinical_significance"]
        .eq(category)
        .groupby(df["variant_id"])
        .sum()
    )


pathogenic = count_category(
    clinvar,
    "Pathogenic"
)

likely_pathogenic = count_category(
    clinvar,
    "Likely pathogenic"
)

pathogenic_likely_pathogenic = count_category(
    clinvar,
    "Pathogenic/Likely pathogenic"
)

benign = count_category(
    clinvar,
    "Benign"
)

likely_benign = count_category(
    clinvar,
    "Likely benign"
)

benign_likely_benign = count_category(
    clinvar,
    "Benign/Likely benign"
)

vus = count_category(
    clinvar,
    "Uncertain significance"
)

conflicting = count_category(
    clinvar,
    "Conflicting classifications of pathogenicity"
)


# ------------------------------------------------------------
# BUILD EXPECTED DATAFRAME
# ------------------------------------------------------------

expected = pd.DataFrame({
    "variant_id": final_df["variant_id"]
}).drop_duplicates()

expected = expected.set_index("variant_id")

expected["pathogenic_count_expected"] = (
    pathogenic
)

expected["likely_pathogenic_count_expected"] = (
    likely_pathogenic
)

expected["pathogenic_likely_pathogenic_count_expected"] = (
    pathogenic_likely_pathogenic
)

expected["benign_count_expected"] = (
    benign
)

expected["likely_benign_count_expected"] = (
    likely_benign
)

expected["benign_likely_benign_count_expected"] = (
    benign_likely_benign
)

expected["vus_count_expected"] = (
    vus
)

expected["conflicting_count_expected"] = (
    conflicting
)

expected = expected.fillna(0)

expected = expected.reset_index()


# ------------------------------------------------------------
# MERGE WITH FINAL DATASET
# ------------------------------------------------------------

comparison = final_df[
    [
        "variant_id",
        "pathogenic_count",
        "likely_pathogenic_count",
        "pathogenic_likely_pathogenic_count",
        "benign_count",
        "likely_benign_count",
        "benign_likely_benign_count",
        "vus_count",
        "conflicting_count"
    ]
].merge(
    expected,
    on="variant_id",
    how="left"
)


# ------------------------------------------------------------
# CHECK EACH EVIDENCE COUNT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("1. EVIDENCE COUNT COMPARISON")
print("=" * 70)

count_pairs = [
    (
        "pathogenic_count",
        "pathogenic_count_expected"
    ),
    (
        "likely_pathogenic_count",
        "likely_pathogenic_count_expected"
    ),
    (
        "pathogenic_likely_pathogenic_count",
        "pathogenic_likely_pathogenic_count_expected"
    ),
    (
        "benign_count",
        "benign_count_expected"
    ),
    (
        "likely_benign_count",
        "likely_benign_count_expected"
    ),
    (
        "benign_likely_benign_count",
        "benign_likely_benign_count_expected"
    ),
    (
        "vus_count",
        "vus_count_expected"
    ),
    (
        "conflicting_count",
        "conflicting_count_expected"
    )
]

count_results = {}

for actual_col, expected_col in count_pairs:

    matches = (
        comparison[actual_col]
        == comparison[expected_col]
    )

    mismatch_count = (~matches).sum()

    count_results[actual_col] = mismatch_count

    print(
        f"{actual_col:<45} "
        f"mismatches: {mismatch_count}"
    )


# ------------------------------------------------------------
# SUBMITTER COUNT CHECK
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("2. SUBMITTER COUNT CONSISTENCY")
print("=" * 70)

submitter_expected = (
    clinvar.groupby("variant_id")[
        "number_submitters"
    ]
    .max()
    .rename("submitter_count_expected")
    .reset_index()
)

submitter_comparison = final_df[
    [
        "variant_id",
        "submitter_count"
    ]
].merge(
    submitter_expected,
    on="variant_id",
    how="left"
)

submitter_comparison[
    "submitter_count_expected"
] = submitter_comparison[
    "submitter_count_expected"
].fillna(0)

submitter_mismatches = (
    submitter_comparison["submitter_count"]
    != submitter_comparison["submitter_count_expected"]
).sum()

print(
    f"submitter_count mismatches: "
    f"{submitter_mismatches}"
)


# ------------------------------------------------------------
# MULTIPLE SUBMITTER FLAG
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("3. MULTIPLE-SUBMITTER FLAG")
print("=" * 70)

expected_multiple_submitters = (
    submitter_comparison[
        "submitter_count_expected"
    ] > 1
).astype(int)

actual_multiple_submitters = (
    final_df["multiple_submitters"]
    .astype(int)
)

multiple_submitter_mismatches = (
    actual_multiple_submitters
    != expected_multiple_submitters
).sum()

print(
    f"multiple_submitters mismatches: "
    f"{multiple_submitter_mismatches}"
)


# ------------------------------------------------------------
# CONFLICT FLAG
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("4. CONFLICT FLAG")
print("=" * 70)

expected_conflict_flag = (
    comparison["conflicting_count_expected"] > 0
).astype(int)

actual_conflict_flag = (
    final_df["conflict_flag"]
    .astype(int)
)

conflict_flag_mismatches = (
    actual_conflict_flag
    != expected_conflict_flag
).sum()

print(
    f"conflict_flag mismatches: "
    f"{conflict_flag_mismatches}"
)


# ------------------------------------------------------------
# EXPERT PANEL FLAG
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("5. EXPERT PANEL FLAG")
print("=" * 70)

expert_panel_expected = (
    clinvar.groupby("variant_id")[
        "review_status"
    ]
    .apply(
        lambda x: int(
            x.astype(str)
             .str.contains(
                 "reviewed by expert panel",
                 case=False,
                 na=False
             )
             .any()
        )
    )
    .rename("expert_panel_expected")
    .reset_index()
)

expert_panel_comparison = final_df[
    [
        "variant_id",
        "expert_panel_review"
    ]
].merge(
    expert_panel_expected,
    on="variant_id",
    how="left"
)

expert_panel_comparison[
    "expert_panel_expected"
] = expert_panel_comparison[
    "expert_panel_expected"
].fillna(0)

expert_panel_mismatches = (
    expert_panel_comparison[
        "expert_panel_review"
    ].astype(int)
    !=
    expert_panel_comparison[
        "expert_panel_expected"
    ].astype(int)
).sum()

print(
    f"expert_panel_review mismatches: "
    f"{expert_panel_mismatches}"
)


# ------------------------------------------------------------
# OVERALL RESULT
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("6. OVERALL CONSISTENCY RESULT")
print("=" * 70)

total_count_mismatches = sum(
    count_results.values()
)

total_mismatches = (
    total_count_mismatches
    + submitter_mismatches
    + multiple_submitter_mismatches
    + conflict_flag_mismatches
    + expert_panel_mismatches
)

print(
    f"Total evidence-count mismatches: "
    f"{total_count_mismatches}"
)

print(
    f"Total consistency mismatches: "
    f"{total_mismatches}"
)


if total_mismatches == 0:

    audit_status = "PASS"

    print("\nSTATUS: PASS")

    print(
        "\nAll tested evidence counts and "
        "derived reliability flags are consistent "
        "with the underlying ClinVar records."
    )

else:

    audit_status = "REVIEW REQUIRED"

    print("\nSTATUS: REVIEW REQUIRED")

    print(
        "\nOne or more evidence features do not "
        "match the underlying ClinVar records."
    )


# ------------------------------------------------------------
# SAVE REPORT
# ------------------------------------------------------------

with open(
    OUTPUT_FILE,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "DEEPGENE STEP 23\n"
        "EVIDENCE COUNT CONSISTENCY AUDIT\n"
    )

    f.write("=" * 70 + "\n\n")

    f.write(
        f"Final dataset rows: {len(final_df)}\n"
    )

    f.write(
        f"Final dataset columns: {len(final_df.columns)}\n"
    )

    f.write(
        f"ClinVar evidence rows checked: "
        f"{len(clinvar)}\n\n"
    )

    f.write(
        "EVIDENCE COUNT MISMATCHES\n"
    )

    f.write("-" * 70 + "\n")

    for column, mismatch_count in count_results.items():

        f.write(
            f"{column}: "
            f"{mismatch_count}\n"
        )

    f.write("\n")

    f.write(
        f"submitter_count mismatches: "
        f"{submitter_mismatches}\n"
    )

    f.write(
        f"multiple_submitters mismatches: "
        f"{multiple_submitter_mismatches}\n"
    )

    f.write(
        f"conflict_flag mismatches: "
        f"{conflict_flag_mismatches}\n"
    )

    f.write(
        f"expert_panel_review mismatches: "
        f"{expert_panel_mismatches}\n"
    )

    f.write("\n")

    f.write(
        f"TOTAL MISMATCHES: "
        f"{total_mismatches}\n"
    )

    f.write(
        f"STATUS: {audit_status}\n"
    )


# ------------------------------------------------------------
# CLOSE DATABASE
# ------------------------------------------------------------

conn.close()


# ------------------------------------------------------------
# FINISH
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("REPORT SAVED")
print("=" * 70)

print(OUTPUT_FILE)

print("\nStep 23 complete.")