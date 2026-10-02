import sqlite3
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data" / "database" / "scn1a.db"
RESULTS_DIR = PROJECT_ROOT / "data" / "analysis_results"
OUTPUT_PATH = RESULTS_DIR / "09_clinical_significance.txt"

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
# 1. Overall database counts
# --------------------------------------------------

section("1. OVERALL CLINVAR COUNTS")

cursor.execute("SELECT COUNT(*) FROM variants")
unique_variants = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM clinvar_records")
clinvar_records = cursor.fetchone()[0]

line(f"Unique variants: {unique_variants}")
line(f"ClinVar records: {clinvar_records}")


# --------------------------------------------------
# 2. Clinical significance at record level
# --------------------------------------------------

section("2. CLINICAL SIGNIFICANCE — CLINVAR RECORD LEVEL")

cursor.execute("""
    SELECT
        clinical_significance,
        COUNT(*) AS record_count
    FROM clinvar_records
    GROUP BY clinical_significance
    ORDER BY record_count DESC
""")

record_results = cursor.fetchall()

for significance, count in record_results:
    significance = significance if significance else "NULL"
    line(f"{significance}: {count}")


# --------------------------------------------------
# 3. Clinical significance at unique variant level
# --------------------------------------------------

section("3. CLINICAL SIGNIFICANCE — UNIQUE VARIANT LEVEL")

cursor.execute("""
    SELECT
        clinical_significance,
        COUNT(DISTINCT variant_id) AS variant_count
    FROM clinvar_records
    GROUP BY clinical_significance
    ORDER BY variant_count DESC
""")

variant_results = cursor.fetchall()

for significance, count in variant_results:
    significance = significance if significance else "NULL"
    line(f"{significance}: {count}")


# --------------------------------------------------
# 4. Variants with multiple classifications
# --------------------------------------------------

section("4. VARIANTS WITH MULTIPLE CLINVAR CLASSIFICATIONS")

cursor.execute("""
    SELECT COUNT(*)
    FROM (
        SELECT variant_id
        FROM clinvar_records
        GROUP BY variant_id
        HAVING COUNT(DISTINCT clinical_significance) > 1
    )
""")

multiple_classification_variants = cursor.fetchone()[0]

line(
    f"Variants with more than one clinical significance: "
    f"{multiple_classification_variants}"
)


# --------------------------------------------------
# 5. Conflicting classifications
# --------------------------------------------------

section("5. CONFLICTING CLASSIFICATIONS")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE clinical_significance =
          'Conflicting classifications of pathogenicity'
""")

conflicting_variants = cursor.fetchone()[0]

line(
    f"Variants classified as conflicting: "
    f"{conflicting_variants}"
)


# --------------------------------------------------
# 6. Pathogenicity-related categories
# --------------------------------------------------

section("6. PATHOGENICITY-RELATED CATEGORIES")

categories = {
    "Pathogenic": "Pathogenic",
    "Likely pathogenic": "Likely pathogenic",
    "Pathogenic/Likely pathogenic": "Pathogenic/Likely pathogenic",
    "Uncertain significance": "Uncertain significance",
    "Likely benign": "Likely benign",
    "Benign": "Benign",
    "Benign/Likely benign": "Benign/Likely benign",
    "Conflicting classifications": "Conflicting classifications of pathogenicity",
}

for label, value in categories.items():

    cursor.execute("""
        SELECT COUNT(DISTINCT variant_id)
        FROM clinvar_records
        WHERE clinical_significance = ?
    """, (value,))

    count = cursor.fetchone()[0]

    line(f"{label}: {count}")


# --------------------------------------------------
# 7. Variants with no clear classification
# --------------------------------------------------

section("7. OTHER / UNRESOLVED CLASSIFICATIONS")

cursor.execute("""
    SELECT
        clinical_significance,
        COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE clinical_significance NOT IN (
        'Pathogenic',
        'Likely pathogenic',
        'Pathogenic/Likely pathogenic',
        'Uncertain significance',
        'Likely benign',
        'Benign',
        'Benign/Likely benign',
        'Conflicting classifications of pathogenicity'
    )
    GROUP BY clinical_significance
    ORDER BY COUNT(DISTINCT variant_id) DESC
""")

other_results = cursor.fetchall()

for significance, count in other_results:
    significance = significance if significance else "NULL"
    line(f"{significance}: {count}")


# --------------------------------------------------
# 8. Review status vs clinical significance
# --------------------------------------------------

section("8. REVIEW STATUS AND CLINICAL SIGNIFICANCE")

cursor.execute("""
    SELECT
        review_status,
        clinical_significance,
        COUNT(DISTINCT variant_id)
    FROM clinvar_records
    GROUP BY review_status, clinical_significance
    ORDER BY review_status, COUNT(DISTINCT variant_id) DESC
""")

review_results = cursor.fetchall()

current_review_status = None

for review_status, significance, count in review_results:

    review_status = review_status if review_status else "NULL"
    significance = significance if significance else "NULL"

    if review_status != current_review_status:
        results.append("")
        line(f"[{review_status}]")
        current_review_status = review_status

    line(f"  {significance}: {count}")


# --------------------------------------------------
# 9. Variants with expert panel review
# --------------------------------------------------

section("9. EXPERT PANEL REVIEW")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE review_status = 'reviewed by expert panel'
""")

expert_panel_variants = cursor.fetchone()[0]

line(
    f"Unique variants with expert panel review: "
    f"{expert_panel_variants}"
)


# --------------------------------------------------
# 10. Summary
# --------------------------------------------------

section("10. ANALYSIS SUMMARY")

line(f"Total unique SCN1A variants: {unique_variants}")
line(f"Total ClinVar records: {clinvar_records}")
line("")
line(
    "Clinical significance was examined at both the ClinVar-record "
    "level and the unique-variant level."
)
line("")
line(
    "Unique-variant counts are used for downstream interpretation "
    "so that multiple genomic representations of the same variant "
    "do not automatically create duplicate biological variants."
)
line("")
line(
    "Clinical significance is an evidence field from ClinVar. "
    "This analysis does not independently determine whether a "
    "variant is disease-causing."
)


# --------------------------------------------------
# Save
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