import sqlite3

database_file = "data/database/scn1a.db"

connection = sqlite3.connect(database_file)
cursor = connection.cursor()

# Remove old tables
cursor.execute("DROP TABLE IF EXISTS clinvar_records")
cursor.execute("DROP TABLE IF EXISTS variant_representations")
cursor.execute("DROP TABLE IF EXISTS variants")
cursor.execute("DROP TABLE IF EXISTS genes")

# 1. Genes
cursor.execute("""
CREATE TABLE genes (
    gene_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT NOT NULL UNIQUE,
    gene_name TEXT
);
""")

# 2. Biological ClinVar variants
cursor.execute("""
CREATE TABLE variants (
    variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    allele_id INTEGER UNIQUE,
    variation_id INTEGER,
    hgvs_name TEXT,
    
    FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
);
""")

# 3. Assembly-specific genomic representations
cursor.execute("""
CREATE TABLE variant_representations (
    representation_id INTEGER PRIMARY KEY AUTOINCREMENT,
    variant_id INTEGER NOT NULL,
    assembly TEXT,
    chromosome TEXT,
    start INTEGER,
    stop INTEGER,
    reference_allele TEXT,
    alternate_allele TEXT,
    dbsnp_id TEXT,
    
    FOREIGN KEY (variant_id) REFERENCES variants(variant_id)
);
""")

# 4. ClinVar clinical evidence
cursor.execute("""
CREATE TABLE clinvar_records (
    clinvar_record_id INTEGER PRIMARY KEY AUTOINCREMENT,
    variant_id INTEGER NOT NULL,
    clinical_significance TEXT,
    review_status TEXT,
    number_submitters INTEGER,
    phenotypes TEXT,
    
    FOREIGN KEY (variant_id) REFERENCES variants(variant_id)
);
""")

# Insert SCN1A
cursor.execute("""
INSERT INTO genes (gene_symbol, gene_name)
VALUES (?, ?);
""", (
    "SCN1A",
    "Sodium voltage-gated channel alpha subunit 1"
))

connection.commit()

print("Database created successfully!")

print("\nTables:")
cursor.execute("""
SELECT name
FROM sqlite_master
WHERE type = 'table'
ORDER BY name;
""")

for table in cursor.fetchall():
    print(" -", table[0])

connection.close()