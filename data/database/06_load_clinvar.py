import sqlite3
import pandas as pd

csv_file = "data/raw/SCN1A_clinvar.tsv"
database_file = "data/database/scn1a.db"

print("Loading ClinVar data...")

df = pd.read_csv(
    csv_file,
    sep="\t",
    low_memory=False
)

print(f"ClinVar records found: {len(df)}")

connection = sqlite3.connect(database_file)
cursor = connection.cursor()

# Get SCN1A gene
cursor.execute("""
SELECT gene_id
FROM genes
WHERE gene_symbol = 'SCN1A';
""")

gene_id = cursor.fetchone()[0]

# Create one variant per AlleleID
variant_ids = {}

for allele_id, group in df.groupby("#AlleleID", sort=False):

    first = group.iloc[0]

    cursor.execute("""
    INSERT INTO variants (
        gene_id,
        allele_id,
        variation_id,
        hgvs_name
    )
    VALUES (?, ?, ?, ?);
    """, (
        gene_id,
        int(allele_id),
        int(first["VariationID"]),
        first["Name"]
    ))

    variant_id = cursor.lastrowid
    variant_ids[allele_id] = variant_id

    # Add every assembly-specific representation
    for _, row in group.iterrows():

        cursor.execute("""
        INSERT INTO variant_representations (
            variant_id,
            assembly,
            chromosome,
            start,
            stop,
            reference_allele,
            alternate_allele,
            dbsnp_id
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?);
        """, (
            variant_id,
            row["Assembly"],
            row["Chromosome"],
            row["Start"],
            row["Stop"],
            row["ReferenceAllele"],
            row["AlternateAllele"],
            row["RS# (dbSNP)"]
        ))

        # Clinical evidence
        cursor.execute("""
        INSERT INTO clinvar_records (
            variant_id,
            clinical_significance,
            review_status,
            number_submitters,
            phenotypes
        )
        VALUES (?, ?, ?, ?, ?);
        """, (
            variant_id,
            row["ClinicalSignificance"],
            row["ReviewStatus"],
            row["NumberSubmitters"],
            row["PhenotypeList"]
        ))

connection.commit()

# Verification
cursor.execute("SELECT COUNT(*) FROM variants;")
variant_count = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM variant_representations;")
representation_count = cursor.fetchone()[0]

cursor.execute("SELECT COUNT(*) FROM clinvar_records;")
record_count = cursor.fetchone()[0]

connection.close()

print("\nImport complete!")
print(f"Unique variants inserted: {variant_count}")
print(f"Genomic representations inserted: {representation_count}")
print(f"ClinVar records inserted: {record_count}")