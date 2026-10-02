import sqlite3
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data" / "database" / "scn1a.db"
RESULTS_DIR = PROJECT_ROOT / "data" / "analysis_results"
OUTPUT_PATH = RESULTS_DIR / "12_variant_representation.txt"

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
# 1. Overall representation counts
# --------------------------------------------------

section("1. OVERALL VARIANT REPRESENTATION")

cursor.execute("SELECT COUNT(*) FROM variants")
unique_variants = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM variant_representations")
total_representations = cursor.fetchone()[0]

line(f"Unique biological variants: {unique_variants}")
line(f"Genomic representations: {total_representations}")


# --------------------------------------------------
# 2. Assembly distribution
# --------------------------------------------------

section("2. ASSEMBLY DISTRIBUTION")

cursor.execute("""
    SELECT
        assembly,
        COUNT(DISTINCT variant_id) AS variant_count,
        COUNT(*) AS representation_count
    FROM variant_representations
    GROUP BY assembly
    ORDER BY representation_count DESC
""")

for assembly, variant_count, representation_count in cursor.fetchall():
    assembly = assembly if assembly else "NULL"

    line(
        f"{assembly}: "
        f"{variant_count} variants, "
        f"{representation_count} representations"
    )


# --------------------------------------------------
# 3. Variants on both assemblies
# --------------------------------------------------

section("3. VARIANTS REPRESENTED ON BOTH ASSEMBLIES")

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

cursor.execute("""
    SELECT COUNT(*)
    FROM variants v
    WHERE NOT EXISTS (
        SELECT 1
        FROM variant_representations r
        WHERE r.variant_id = v.variant_id
          AND r.assembly IN ('GRCh37', 'GRCh38')
    )
""")

neither = cursor.fetchone()[0]

line(f"Both GRCh37 and GRCh38: {both}")
line(f"GRCh37 only: {grch37_only}")
line(f"GRCh38 only: {grch38_only}")
line(f"Neither GRCh37 nor GRCh38: {neither}")


# --------------------------------------------------
# 4. Missing coordinate information
# --------------------------------------------------

section("4. MISSING GENOMIC COORDINATES")

cursor.execute("""
    SELECT COUNT(*)
    FROM variant_representations
    WHERE chromosome IS NULL
       OR TRIM(chromosome) = ''
       OR chromosome = 'na'
       OR start IS NULL
       OR stop IS NULL
""")

missing_coordinates = cursor.fetchone()[0]

line(f"Representations without complete coordinates: {missing_coordinates}")

cursor.execute("""
    SELECT
        assembly,
        COUNT(*)
    FROM variant_representations
    WHERE chromosome IS NULL
       OR TRIM(chromosome) = ''
       OR chromosome = 'na'
       OR start IS NULL
       OR stop IS NULL
    GROUP BY assembly
""")

for assembly, count in cursor.fetchall():
    assembly = assembly if assembly else "NULL"
    line(f"{assembly}: {count}")


# --------------------------------------------------
# 5. dbSNP consistency
# --------------------------------------------------

section("5. dbSNP REPRESENTATION CONSISTENCY")

cursor.execute("""
    SELECT COUNT(*)
    FROM (
        SELECT variant_id
        FROM variant_representations
        WHERE dbsnp_id IS NOT NULL
          AND TRIM(dbsnp_id) != ''
          AND dbsnp_id != '-'
        GROUP BY variant_id
        HAVING COUNT(DISTINCT dbsnp_id) > 1
    )
""")

multiple_dbsnp = cursor.fetchone()[0]

line(
    f"Variants with multiple dbSNP IDs across representations: "
    f"{multiple_dbsnp}"
)

cursor.execute("""
    SELECT COUNT(*)
    FROM (
        SELECT variant_id
        FROM variant_representations
        GROUP BY variant_id
        HAVING COUNT(*) > 1
           AND COUNT(DISTINCT reference_allele) > 1
    )
""")

different_reference = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(*)
    FROM (
        SELECT variant_id
        FROM variant_representations
        GROUP BY variant_id
        HAVING COUNT(*) > 1
           AND COUNT(DISTINCT alternate_allele) > 1
    )
""")

different_alternate = cursor.fetchone()[0]

line(
    f"Variants with different reference alleles across representations: "
    f"{different_reference}"
)

line(
    f"Variants with different alternate alleles across representations: "
    f"{different_alternate}"
)


# --------------------------------------------------
# 6. Coordinate examples across assemblies
# --------------------------------------------------

section("6. EXAMPLES OF MULTI-ASSEMBLY REPRESENTATIONS")

cursor.execute("""
    SELECT
        v.variant_id,
        v.allele_id,
        v.hgvs_name,
        r.assembly,
        r.chromosome,
        r.start,
        r.stop,
        r.reference_allele,
        r.alternate_allele,
        r.dbsnp_id
    FROM variants v
    JOIN variant_representations r
        ON v.variant_id = r.variant_id
    WHERE v.variant_id IN (
        SELECT variant_id
        FROM variant_representations
        GROUP BY variant_id
        HAVING COUNT(DISTINCT assembly) >= 2
    )
    ORDER BY v.variant_id, r.assembly
    LIMIT 20
""")

for row in cursor.fetchall():

    (
        variant_id,
        allele_id,
        hgvs_name,
        assembly,
        chromosome,
        start,
        stop,
        reference,
        alternate,
        dbsnp
    ) = row

    line("")
    line(f"Variant ID: {variant_id}")
    line(f"Allele ID: {allele_id}")
    line(f"HGVS: {hgvs_name}")
    line(f"Assembly: {assembly}")
    line(f"Coordinate: chr{chromosome}:{start}-{stop}")
    line(f"Reference: {reference}")
    line(f"Alternate: {alternate}")
    line(f"dbSNP: {dbsnp}")


# --------------------------------------------------
# 7. Representation model validation
# --------------------------------------------------

section("7. REPRESENTATION MODEL VALIDATION")

cursor.execute("""
    SELECT COUNT(*)
    FROM variants v
    WHERE NOT EXISTS (
        SELECT 1
        FROM variant_representations r
        WHERE r.variant_id = v.variant_id
    )
""")

variants_without_representation = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(*)
    FROM variant_representations r
    LEFT JOIN variants v
        ON r.variant_id = v.variant_id
    WHERE v.variant_id IS NULL
""")

orphan_representations = cursor.fetchone()[0]

line(
    f"Variants without any genomic representation: "
    f"{variants_without_representation}"
)

line(
    f"Orphan genomic representations: "
    f"{orphan_representations}"
)


# --------------------------------------------------
# 8. Potential DeepGene use
# --------------------------------------------------

section("8. POTENTIAL USE IN DEEPGENE")

line(
    "The representation layer separates biological variants "
    "from genome-assembly-specific coordinates."
)
line("")
line(
    "This prevents GRCh37 and GRCh38 representations of the "
    "same variant from being treated as separate biological variants."
)
line("")
line(
    "Assembly-specific coordinates can later support genomic "
    "annotation and external database integration."
)
line("")
line(
    "This analysis does not assign pathogenicity or prioritize variants."
)


# --------------------------------------------------
# 9. Summary
# --------------------------------------------------

section("9. ANALYSIS SUMMARY")

line(f"Unique biological variants: {unique_variants}")
line(f"Genomic representations: {total_representations}")
line(f"Represented on both GRCh37 and GRCh38: {both}")
line(f"GRCh37-only variants: {grch37_only}")
line(f"GRCh38-only variants: {grch38_only}")
line(f"Representations without coordinates: {missing_coordinates}")
line("")
line(
    "The representation layer successfully separates biological "
    "variant identity from assembly-specific genomic coordinates."
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