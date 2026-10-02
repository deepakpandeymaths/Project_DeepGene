import sqlite3
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data" / "database" / "scn1a.db"
RESULTS_DIR = PROJECT_ROOT / "data" / "analysis_results"
OUTPUT_PATH = RESULTS_DIR / "10_evidence_quality.txt"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# Database connection
# --------------------------------------------------

conn = sqlite3.connect(DB_PATH)
cursor = conn.cursor()

results = []


def section(title):
    results.append("")
    results.append("=" * 60)
    results.append(title)
    results.append("=" * 60)


def line(text):
    results.append(str(text))


# --------------------------------------------------
# 1. Overall counts
# --------------------------------------------------

section("1. OVERALL EVIDENCE COUNTS")

cursor.execute("SELECT COUNT(*) FROM variants")
unique_variants = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM clinvar_records")
clinvar_records = cursor.fetchone()[0]

line(f"Unique variants: {unique_variants}")
line(f"ClinVar records: {clinvar_records}")


# --------------------------------------------------
# 2. Review status
# --------------------------------------------------

section("2. REVIEW STATUS")

cursor.execute("""
    SELECT
        review_status,
        COUNT(*) AS record_count,
        COUNT(DISTINCT variant_id) AS variant_count
    FROM clinvar_records
    GROUP BY review_status
    ORDER BY variant_count DESC
""")

for review_status, record_count, variant_count in cursor.fetchall():

    review_status = review_status if review_status else "NULL"

    line(
        f"{review_status}: "
        f"{variant_count} unique variants, "
        f"{record_count} records"
    )


# --------------------------------------------------
# 3. Number of submitters
# --------------------------------------------------

section("3. NUMBER OF SUBMITTERS")

cursor.execute("""
    SELECT
        number_submitters,
        COUNT(*) AS record_count,
        COUNT(DISTINCT variant_id) AS variant_count
    FROM clinvar_records
    GROUP BY number_submitters
    ORDER BY number_submitters
""")

for submitters, record_count, variant_count in cursor.fetchall():

    submitters = submitters if submitters is not None else "NULL"

    line(
        f"{submitters} submitter(s): "
        f"{variant_count} unique variants, "
        f"{record_count} records"
    )


# --------------------------------------------------
# 4. Single vs multiple submitters
# --------------------------------------------------

section("4. SINGLE VS MULTIPLE SUBMITTERS")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE number_submitters = 1
""")

single_submitter = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE number_submitters > 1
""")

multiple_submitters = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE number_submitters IS NULL
""")

unknown_submitters = cursor.fetchone()[0]

line(f"Variants with one submitter: {single_submitter}")
line(f"Variants with multiple submitters: {multiple_submitters}")
line(f"Variants with unknown submitter count: {unknown_submitters}")


# --------------------------------------------------
# 5. Expert panel
# --------------------------------------------------

section("5. EXPERT PANEL EVIDENCE")

cursor.execute("""
    SELECT
        clinical_significance,
        COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE review_status = 'reviewed by expert panel'
    GROUP BY clinical_significance
    ORDER BY COUNT(DISTINCT variant_id) DESC
""")

expert_results = cursor.fetchall()

expert_total = sum(count for _, count in expert_results)

line(f"Total unique variants with expert panel review: {expert_total}")
line("")

for significance, count in expert_results:

    significance = significance if significance else "NULL"

    line(f"{significance}: {count}")


# --------------------------------------------------
# 6. Conflicting classifications
# --------------------------------------------------

section("6. CONFLICTING CLASSIFICATIONS")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE review_status =
          'criteria provided, conflicting classifications'
""")

conflicting_review = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE clinical_significance =
          'Conflicting classifications of pathogenicity'
""")

conflicting_significance = cursor.fetchone()[0]

line(
    "Variants with conflicting review status: "
    f"{conflicting_review}"
)

line(
    "Variants with conflicting clinical significance: "
    f"{conflicting_significance}"
)


# --------------------------------------------------
# 7. Review status × submitters
# --------------------------------------------------

section("7. REVIEW STATUS × NUMBER OF SUBMITTERS")

cursor.execute("""
    SELECT
        review_status,
        number_submitters,
        COUNT(DISTINCT variant_id)
    FROM clinvar_records
    GROUP BY review_status, number_submitters
    ORDER BY review_status, number_submitters
""")

current_status = None

for review_status, submitters, count in cursor.fetchall():

    review_status = review_status if review_status else "NULL"
    submitters = submitters if submitters is not None else "NULL"

    if review_status != current_status:
        results.append("")
        line(f"[{review_status}]")
        current_status = review_status

    line(
        f"  {submitters} submitter(s): "
        f"{count} unique variants"
    )


# --------------------------------------------------
# 8. Evidence-quality categories
# --------------------------------------------------

section("8. BROAD EVIDENCE CATEGORIES")

evidence_categories = {
    "Expert panel reviewed":
        "reviewed by expert panel",

    "Multiple submitters, no conflicts":
        "criteria provided, multiple submitters, no conflicts",

    "Single submitter":
        "criteria provided, single submitter",

    "Conflicting":
        "criteria provided, conflicting classifications",

    "No assertion criteria":
        "no assertion criteria provided",

    "No classification":
        "no classification provided",

    "Unspecified":
        "-",
}

for label, status in evidence_categories.items():

    cursor.execute("""
        SELECT COUNT(DISTINCT variant_id)
        FROM clinvar_records
        WHERE review_status = ?
    """, (status,))

    count = cursor.fetchone()[0]

    line(f"{label}: {count}")


# --------------------------------------------------
# 9. Potential evidence signals
# --------------------------------------------------

section("9. POTENTIAL EVIDENCE SIGNALS FOR DEEPGENE")

line("The current database contains several fields that")
line("could later contribute to a transparent prioritization")
line("system:")
line("")
line("- Clinical significance")
line("- Review status")
line("- Number of submitters")
line("- Expert panel review")
line("- Conflicting classification status")
line("- Phenotype information")
line("")
line("These fields should be treated as evidence signals,")
line("not as an automatic disease-causality score.")
line("")
line("This analysis does not assign weights or rank variants.")


# --------------------------------------------------
# 10. Summary
# --------------------------------------------------

section("10. ANALYSIS SUMMARY")

line(
    "Evidence quality was examined using ClinVar review status "
    "and submitter information."
)
line("")
line(
    "Review status provides information about how ClinVar "
    "classifications were supported, while NumberSubmitters "
    "provides an additional description of the available evidence."
)
line("")
line(
    "These fields can later be incorporated into a transparent "
    "variant-prioritization framework."
)
line("")
line(
    "No pathogenicity score or ranking is produced by this analysis."
)


# --------------------------------------------------
# Save results
# --------------------------------------------------

conn.close()

OUTPUT_PATH.write_text(
    "\n".join(results),
    encoding="utf-8"
)

print("\n".join(results))

print("\n" + "=" * 60)
print("Analysis saved to:")
print(OUTPUT_PATH)
print("=" * 60)