import sqlite3
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data" / "database" / "scn1a.db"
RESULTS_DIR = PROJECT_ROOT / "data" / "analysis_results"
OUTPUT_PATH = RESULTS_DIR / "08_variant_statistics.txt"


# --------------------------------------------------
# Setup
# --------------------------------------------------

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

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
# 1. Database overview
# --------------------------------------------------

section("1. DATABASE OVERVIEW")

cursor.execute("SELECT COUNT(*) FROM genes")
line(f"Genes: {cursor.fetchone()[0]}")

cursor.execute("SELECT COUNT(*) FROM variants")
line(f"Unique variants: {cursor.fetchone()[0]}")

cursor.execute("SELECT COUNT(*) FROM variant_representations")
line(f"Genomic representations: {cursor.fetchone()[0]}")

cursor.execute("SELECT COUNT(*) FROM clinvar_records")
line(f"ClinVar records: {cursor.fetchone()[0]}")


# --------------------------------------------------
# 2. Variant type
# --------------------------------------------------

section("2. VARIANT TYPES")

# Type is not currently stored in the database.
# This section is intentionally documented so we don't
# incorrectly derive variant type from HGVS text.

line("Variant type is not currently stored as a separate")
line("database field.")
line("")
line("The current database schema stores:")
line("- HGVS name")
line("- genomic representation")
line("- reference allele")
line("- alternate allele")
line("- ClinVar clinical evidence")
line("")
line("Therefore, variant-type statistics will be added")
line("later if/when the database schema includes this field.")


# --------------------------------------------------
# 3. dbSNP coverage
# --------------------------------------------------

section("3. dbSNP COVERAGE")

cursor.execute("""
    SELECT COUNT(*)
    FROM variant_representations
    WHERE dbsnp_id IS NOT NULL
      AND TRIM(dbsnp_id) != ''
      AND dbsnp_id != '-'
""")

with_dbsnp = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(*)
    FROM variant_representations
""")

total_representations = cursor.fetchone()[0]

without_dbsnp = total_representations - with_dbsnp

line(f"Representations with dbSNP ID: {with_dbsnp}")
line(f"Representations without dbSNP ID: {without_dbsnp}")

if total_representations > 0:
    percentage = (with_dbsnp / total_representations) * 100
    line(f"dbSNP coverage: {percentage:.2f}%")


# --------------------------------------------------
# 4. Genomic coordinates
# --------------------------------------------------

section("4. GENOMIC COORDINATE COVERAGE")

cursor.execute("""
    SELECT COUNT(*)
    FROM variant_representations
    WHERE chromosome IS NOT NULL
      AND TRIM(chromosome) != ''
      AND chromosome != 'na'
      AND start IS NOT NULL
      AND stop IS NOT NULL
""")

with_coordinates = cursor.fetchone()[0]

without_coordinates = total_representations - with_coordinates

line(f"Representations with genomic coordinates: {with_coordinates}")
line(f"Representations without genomic coordinates: {without_coordinates}")

if total_representations > 0:
    percentage = (with_coordinates / total_representations) * 100
    line(f"Coordinate coverage: {percentage:.2f}%")


# --------------------------------------------------
# 5. Assembly distribution
# --------------------------------------------------

section("5. GENOME ASSEMBLY DISTRIBUTION")

cursor.execute("""
    SELECT assembly, COUNT(*)
    FROM variant_representations
    GROUP BY assembly
    ORDER BY COUNT(*) DESC
""")

for assembly, count in cursor.fetchall():
    line(f"{assembly}: {count}")


# --------------------------------------------------
# 6. Variants represented on both assemblies
# --------------------------------------------------

section("6. VARIANT REPRESENTATION ACROSS ASSEMBLIES")

cursor.execute("""
    SELECT
        COUNT(*) AS total_variants,
        SUM(CASE WHEN assembly = 'GRCh37' THEN 1 ELSE 0 END),
        SUM(CASE WHEN assembly = 'GRCh38' THEN 1 ELSE 0 END)
    FROM (
        SELECT variant_id, assembly
        FROM variant_representations
        GROUP BY variant_id, assembly
    )
""")

# Calculate using variant-level queries instead
cursor.execute("""
    SELECT COUNT(*)
    FROM variants v
    WHERE EXISTS (
        SELECT 1
        FROM variant_representations r
        WHERE r.variant_id = v.variant_id
          AND r.assembly = 'GRCh37'
    )
    AND EXISTS (
        SELECT 1
        FROM variant_representations r
        WHERE r.variant_id = v.variant_id
          AND r.assembly = 'GRCh38'
    )
""")

both = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(*)
    FROM variants v
    WHERE EXISTS (
        SELECT 1
        FROM variant_representations r
        WHERE r.variant_id = v.variant_id
          AND r.assembly = 'GRCh37'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM variant_representations r
        WHERE r.variant_id = v.variant_id
          AND r.assembly = 'GRCh38'
    )
""")

grch37_only = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(*)
    FROM variants v
    WHERE EXISTS (
        SELECT 1
        FROM variant_representations r
        WHERE r.variant_id = v.variant_id
          AND r.assembly = 'GRCh38'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM variant_representations r
        WHERE r.variant_id = v.variant_id
          AND r.assembly = 'GRCh37'
    )
""")

grch38_only = cursor.fetchone()[0]

line(f"Variants represented on both GRCh37 and GRCh38: {both}")
line(f"Variants represented only on GRCh37: {grch37_only}")
line(f"Variants represented only on GRCh38: {grch38_only}")


# --------------------------------------------------
# 7. Reference / alternate allele information
# --------------------------------------------------

section("7. ALLELE INFORMATION")

cursor.execute("""
    SELECT COUNT(*)
    FROM variant_representations
    WHERE reference_allele IS NOT NULL
      AND alternate_allele IS NOT NULL
      AND TRIM(reference_allele) != ''
      AND TRIM(alternate_allele) != ''
""")

complete_alleles = cursor.fetchone()[0]

line(f"Representations with reference and alternate alleles: {complete_alleles}")
line(f"Representations without complete allele information: "
     f"{total_representations - complete_alleles}")


# --------------------------------------------------
# 8. Basic consistency checks
# --------------------------------------------------

section("8. DATABASE CONSISTENCY CHECKS")

cursor.execute("""
    SELECT COUNT(*)
    FROM variants
    WHERE gene_id IS NULL
""")

line(f"Variants without gene_id: {cursor.fetchone()[0]}")

cursor.execute("""
    SELECT COUNT(*)
    FROM variant_representations r
    LEFT JOIN variants v
        ON r.variant_id = v.variant_id
    WHERE v.variant_id IS NULL
""")

line(f"Orphan genomic representations: {cursor.fetchone()[0]}")

cursor.execute("""
    SELECT COUNT(*)
    FROM clinvar_records c
    LEFT JOIN variants v
        ON c.variant_id = v.variant_id
    WHERE v.variant_id IS NULL
""")

line(f"Orphan ClinVar records: {cursor.fetchone()[0]}")


# --------------------------------------------------
# 9. Summary
# --------------------------------------------------

section("9. ANALYSIS SUMMARY")

cursor.execute("SELECT COUNT(*) FROM variants")
unique_variants = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM variant_representations")
representations = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM clinvar_records")
clinvar_records = cursor.fetchone()[0]

line(f"The database contains {unique_variants} unique SCN1A variants.")
line(f"These variants have {representations} genomic representations.")
line(f"The database contains {clinvar_records} ClinVar records.")
line("")
line("The genomic representation layer allows the same biological")
line("variant to be represented separately on GRCh37 and GRCh38.")
line("")
line("This analysis describes the structure and completeness of")
line("the current SCN1A ClinVar database. It does not assign")
line("pathogenicity or prioritize variants.")


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
print(f"Analysis saved to:")
print(OUTPUT_PATH)
print("=" * 60)