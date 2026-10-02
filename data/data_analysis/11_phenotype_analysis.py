import sqlite3
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data" / "database" / "scn1a.db"
RESULTS_DIR = PROJECT_ROOT / "data" / "analysis_results"
OUTPUT_PATH = RESULTS_DIR / "11_phenotype_analysis.txt"

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
# 1. Overall phenotype coverage
# --------------------------------------------------

section("1. PHENOTYPE COVERAGE")

cursor.execute("SELECT COUNT(*) FROM variants")
unique_variants = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM clinvar_records")
total_records = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE phenotypes IS NOT NULL
      AND TRIM(phenotypes) != ''
      AND phenotypes != '-'
""")

variants_with_phenotypes = cursor.fetchone()[0]

variants_without_phenotypes = (
    unique_variants - variants_with_phenotypes
)

line(f"Unique variants: {unique_variants}")
line(f"ClinVar records: {total_records}")
line(f"Variants with phenotype information: {variants_with_phenotypes}")
line(f"Variants without phenotype information: {variants_without_phenotypes}")

if unique_variants > 0:
    coverage = (
        variants_with_phenotypes / unique_variants
    ) * 100

    line(f"Phenotype coverage: {coverage:.2f}%")


# --------------------------------------------------
# 2. Phenotype coverage at record level
# --------------------------------------------------

section("2. PHENOTYPE COVERAGE — RECORD LEVEL")

cursor.execute("""
    SELECT COUNT(*)
    FROM clinvar_records
    WHERE phenotypes IS NOT NULL
      AND TRIM(phenotypes) != ''
      AND phenotypes != '-'
""")

records_with_phenotypes = cursor.fetchone()[0]

records_without_phenotypes = (
    total_records - records_with_phenotypes
)

line(f"Records with phenotype information: {records_with_phenotypes}")
line(f"Records without phenotype information: {records_without_phenotypes}")

if total_records > 0:
    coverage = (
        records_with_phenotypes / total_records
    ) * 100

    line(f"Record-level phenotype coverage: {coverage:.2f}%")


# --------------------------------------------------
# 3. Most frequent phenotype descriptions
# --------------------------------------------------

section("3. MOST FREQUENT PHENOTYPE DESCRIPTIONS")

cursor.execute("""
    SELECT
        phenotypes,
        COUNT(DISTINCT variant_id) AS variant_count
    FROM clinvar_records
    WHERE phenotypes IS NOT NULL
      AND TRIM(phenotypes) != ''
      AND phenotypes != '-'
    GROUP BY phenotypes
    ORDER BY variant_count DESC
    LIMIT 30
""")

phenotype_results = cursor.fetchall()

for phenotype, count in phenotype_results:
    line(f"{count}: {phenotype}")


# --------------------------------------------------
# 4. Variants with multiple phenotype descriptions
# --------------------------------------------------

section("4. VARIANTS WITH MULTIPLE PHENOTYPE DESCRIPTIONS")

cursor.execute("""
    SELECT COUNT(*)
    FROM (
        SELECT variant_id
        FROM clinvar_records
        WHERE phenotypes IS NOT NULL
          AND TRIM(phenotypes) != ''
          AND phenotypes != '-'
        GROUP BY variant_id
        HAVING COUNT(DISTINCT phenotypes) > 1
    )
""")

multiple_phenotypes = cursor.fetchone()[0]

line(
    f"Variants with multiple distinct phenotype descriptions: "
    f"{multiple_phenotypes}"
)


# --------------------------------------------------
# 5. Number of phenotype descriptions per variant
# --------------------------------------------------

section("5. PHENOTYPE COUNT PER VARIANT")

cursor.execute("""
    SELECT
        phenotype_count,
        COUNT(*) AS variant_count
    FROM (
        SELECT
            variant_id,
            COUNT(DISTINCT phenotypes) AS phenotype_count
        FROM clinvar_records
        WHERE phenotypes IS NOT NULL
          AND TRIM(phenotypes) != ''
          AND phenotypes != '-'
        GROUP BY variant_id
    )
    GROUP BY phenotype_count
    ORDER BY phenotype_count
""")

for phenotype_count, variant_count in cursor.fetchall():
    line(
        f"{phenotype_count} phenotype description(s): "
        f"{variant_count} variants"
    )


# --------------------------------------------------
# 6. Phenotypes associated with clinical significance
# --------------------------------------------------

section("6. PHENOTYPE INFORMATION BY CLINICAL SIGNIFICANCE")

cursor.execute("""
    SELECT
        clinical_significance,
        COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE phenotypes IS NOT NULL
      AND TRIM(phenotypes) != ''
      AND phenotypes != '-'
    GROUP BY clinical_significance
    ORDER BY COUNT(DISTINCT variant_id) DESC
""")

for significance, count in cursor.fetchall():

    significance = (
        significance
        if significance
        else "NULL"
    )

    line(
        f"{significance}: "
        f"{count} variants with phenotype information"
    )


# --------------------------------------------------
# 7. Phenotype information for conflicting variants
# --------------------------------------------------

section("7. PHENOTYPE INFORMATION FOR CONFLICTING VARIANTS")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE clinical_significance =
          'Conflicting classifications of pathogenicity'
      AND phenotypes IS NOT NULL
      AND TRIM(phenotypes) != ''
      AND phenotypes != '-'
""")

conflicting_with_phenotype = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE clinical_significance =
          'Conflicting classifications of pathogenicity'
""")

total_conflicting = cursor.fetchone()[0]

line(
    f"Conflicting variants with phenotype information: "
    f"{conflicting_with_phenotype}"
)

line(
    f"Total conflicting variants: "
    f"{total_conflicting}"
)


# --------------------------------------------------
# 8. Example phenotype records
# --------------------------------------------------

section("8. EXAMPLE PHENOTYPE RECORDS")

cursor.execute("""
    SELECT DISTINCT
        v.variant_id,
        v.allele_id,
        v.hgvs_name,
        c.clinical_significance,
        c.phenotypes
    FROM variants v
    JOIN clinvar_records c
        ON v.variant_id = c.variant_id
    WHERE c.phenotypes IS NOT NULL
      AND TRIM(c.phenotypes) != ''
      AND c.phenotypes != '-'
    ORDER BY v.variant_id
    LIMIT 20
""")

for variant_id, allele_id, hgvs_name, significance, phenotype in cursor.fetchall():

    line("")
    line(f"Variant ID: {variant_id}")
    line(f"Allele ID: {allele_id}")
    line(f"HGVS: {hgvs_name}")
    line(f"Clinical significance: {significance}")
    line(f"Phenotype: {phenotype}")


# --------------------------------------------------
# 9. Potential DeepGene use
# --------------------------------------------------

section("9. POTENTIAL USE IN DEEPGENE")

line(
    "Phenotype information can provide clinical context for "
    "SCN1A variants."
)
line("")
line(
    "The current database preserves the phenotype information "
    "reported in ClinVar."
)
line("")
line(
    "Phenotype information should not be treated as independent "
    "proof of pathogenicity."
)
line("")
line(
    "For future prioritization, phenotype matching could become "
    "one evidence component alongside clinical significance, "
    "review status, submitter information, and other available "
    "evidence."
)


# --------------------------------------------------
# 10. Summary
# --------------------------------------------------

section("10. ANALYSIS SUMMARY")

line(f"Unique SCN1A variants: {unique_variants}")
line(f"Variants with phenotype information: {variants_with_phenotypes}")
line(f"Variants without phenotype information: {variants_without_phenotypes}")
line("")
line(
    "This analysis describes phenotype coverage and the phenotype "
    "descriptions present in the current ClinVar dataset."
)
line("")
line(
    "It does not determine whether a phenotype is caused by a "
    "specific variant."
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