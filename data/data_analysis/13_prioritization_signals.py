import sqlite3
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data" / "database" / "scn1a.db"
RESULTS_DIR = PROJECT_ROOT / "data" / "analysis_results"
OUTPUT_PATH = RESULTS_DIR / "13_prioritization_signals.txt"

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
# 1. Dataset overview
# --------------------------------------------------

section("1. DATASET OVERVIEW")

cursor.execute("SELECT COUNT(*) FROM variants")
unique_variants = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM clinvar_records")
clinvar_records = cursor.fetchone()[0]

line(f"Unique variants: {unique_variants}")
line(f"ClinVar records: {clinvar_records}")


# --------------------------------------------------
# 2. Clinical significance signals
# --------------------------------------------------

section("2. CLINICAL SIGNIFICANCE SIGNAL")

cursor.execute("""
    SELECT
        clinical_significance,
        COUNT(DISTINCT variant_id) AS variant_count
    FROM clinvar_records
    GROUP BY clinical_significance
    ORDER BY variant_count DESC
""")

for significance, count in cursor.fetchall():
    line(f"{significance}: {count}")


# --------------------------------------------------
# 3. Review status signals
# --------------------------------------------------

section("3. REVIEW STATUS SIGNAL")

cursor.execute("""
    SELECT
        review_status,
        COUNT(DISTINCT variant_id) AS variant_count
    FROM clinvar_records
    GROUP BY review_status
    ORDER BY variant_count DESC
""")

for review_status, count in cursor.fetchall():
    line(f"{review_status}: {count}")


# --------------------------------------------------
# 4. Submitter evidence
# --------------------------------------------------

section("4. SUBMITTER COUNT SIGNAL")

cursor.execute("""
    SELECT
        number_submitters,
        COUNT(DISTINCT variant_id) AS variant_count
    FROM clinvar_records
    GROUP BY number_submitters
    ORDER BY number_submitters
""")

for submitters, count in cursor.fetchall():
    line(f"{submitters} submitter(s): {count} variants")


# --------------------------------------------------
# 5. Single vs multiple submitters
# --------------------------------------------------

section("5. SINGLE VS MULTIPLE SUBMITTERS")

cursor.execute("""
    SELECT
        CASE
            WHEN number_submitters = 1 THEN 'one submitter'
            WHEN number_submitters > 1 THEN 'multiple submitters'
            ELSE 'unknown'
        END AS category,
        COUNT(DISTINCT variant_id)
    FROM clinvar_records
    GROUP BY category
""")

for category, count in cursor.fetchall():
    line(f"{category}: {count}")


# --------------------------------------------------
# 6. Expert panel signal
# --------------------------------------------------

section("6. EXPERT PANEL SIGNAL")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE LOWER(review_status) = 'reviewed by expert panel'
""")

expert_panel_variants = cursor.fetchone()[0]

line(f"Variants reviewed by expert panel: {expert_panel_variants}")


# --------------------------------------------------
# 7. Conflicting evidence signal
# --------------------------------------------------

section("7. CONFLICTING EVIDENCE SIGNAL")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE LOWER(review_status)
          LIKE '%conflicting classifications%'
""")

conflicting_variants = cursor.fetchone()[0]

line(f"Variants with conflicting review status: {conflicting_variants}")


# --------------------------------------------------
# 8. Phenotype signal
# --------------------------------------------------

section("8. PHENOTYPE SIGNAL")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE phenotypes IS NOT NULL
      AND TRIM(phenotypes) != ''
""")

variants_with_phenotype = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM clinvar_records
    WHERE phenotypes IS NULL
       OR TRIM(phenotypes) = ''
""")

variants_without_phenotype = cursor.fetchone()[0]

line(f"Variants with phenotype information: {variants_with_phenotype}")
line(f"Variants without phenotype information: {variants_without_phenotype}")


# --------------------------------------------------
# 9. Genomic representation signal
# --------------------------------------------------

section("9. GENOMIC REPRESENTATION SIGNAL")

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM variant_representations
    WHERE assembly = 'GRCh37'
""")

grch37_variants = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM variant_representations
    WHERE assembly = 'GRCh38'
""")

grch38_variants = cursor.fetchone()[0]

cursor.execute("""
    SELECT COUNT(DISTINCT variant_id)
    FROM variant_representations
    WHERE dbsnp_id IS NOT NULL
      AND TRIM(dbsnp_id) != ''
""")

dbsnp_variants = cursor.fetchone()[0]

line(f"Variants represented on GRCh37: {grch37_variants}")
line(f"Variants represented on GRCh38: {grch38_variants}")
line(f"Variants with dbSNP identifiers: {dbsnp_variants}")


# --------------------------------------------------
# 10. Potential prioritization signals
# --------------------------------------------------

section("10. AVAILABLE PRIORITIZATION SIGNALS")

signals = [
    "Clinical significance",
    "Review status",
    "Number of submitters",
    "Expert panel review",
    "Conflicting classification status",
    "Phenotype information",
    "Genome assembly representation",
    "Genomic coordinates",
    "dbSNP identifier",
    "HGVS variant description"
]

for signal in signals:
    line(f"- {signal}")


# --------------------------------------------------
# 11. Signals NOT currently available
# --------------------------------------------------

section("11. IMPORTANT SIGNALS NOT CURRENTLY AVAILABLE")

missing_signals = [
    "Population allele frequency",
    "Functional assay results",
    "Protein structure information",
    "Conservation scores",
    "Splice prediction scores",
    "Protein impact prediction scores",
    "Gene constraint metrics",
    "Individual patient-level genotype/phenotype data",
    "Independent population databases",
    "Experimental functional evidence beyond ClinVar annotations"
]

for signal in missing_signals:
    line(f"- {signal}")


# --------------------------------------------------
# 12. Evidence interpretation
# --------------------------------------------------

section("12. EVIDENCE INTERPRETATION")

line(
    "The current dataset contains several signals that can describe "
    "the strength and context of existing ClinVar evidence."
)

line("")
line(
    "Clinical significance provides the existing ClinVar classification "
    "associated with a variant."
)

line("")
line(
    "Review status and submitter count provide information about how "
    "the classification was supported and reviewed."
)

line("")
line(
    "Expert panel review and conflicting classifications are additional "
    "evidence-quality signals."
)

line("")
line(
    "Phenotype information provides disease or phenotype context, "
    "but should not be treated as independent proof of pathogenicity."
)

line("")
line(
    "Genomic representation and dbSNP identifiers primarily support "
    "variant identity and integration rather than directly measuring "
    "pathogenicity."
)


# --------------------------------------------------
# 13. Future DeepGene model
# --------------------------------------------------

section("13. FUTURE DEEPGENE MODEL CONSIDERATIONS")

line(
    "The current dataset is sufficient to build an initial evidence-aware "
    "variant prioritization framework."
)

line("")
line(
    "However, a robust predictive model should eventually combine "
    "ClinVar evidence with additional independent biological and "
    "population-level evidence."
)

line("")
line(
    "The current analysis does not assign scores, weights, rankings, "
    "or pathogenicity predictions."
)


# --------------------------------------------------
# 14. Summary
# --------------------------------------------------

section("14. ANALYSIS SUMMARY")

line(f"Unique variants available: {unique_variants}")
line(f"ClinVar records available: {clinvar_records}")
line(f"Variants with expert panel review: {expert_panel_variants}")
line(f"Variants with conflicting review status: {conflicting_variants}")
line(f"Variants with phenotype information: {variants_with_phenotype}")
line(f"Variants represented on GRCh37: {grch37_variants}")
line(f"Variants represented on GRCh38: {grch38_variants}")
line("")
line(
    "The current ClinVar database contains multiple useful evidence "
    "signals for future DeepGene prioritization."
)
line(
    "Additional biological and population-level datasets will be "
    "needed for a more comprehensive model."
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