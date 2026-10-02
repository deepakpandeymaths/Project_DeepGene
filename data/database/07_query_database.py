import sqlite3

database_file = "data/database/scn1a.db"

connection = sqlite3.connect(database_file)
cursor = connection.cursor()

print("\n=== DATABASE SUMMARY ===")

cursor.execute("SELECT COUNT(*) FROM genes;")
print("Genes:", cursor.fetchone()[0])

cursor.execute("SELECT COUNT(*) FROM variants;")
print("Unique variants:", cursor.fetchone()[0])

cursor.execute("SELECT COUNT(*) FROM variant_representations;")
print("Genomic representations:", cursor.fetchone()[0])

cursor.execute("SELECT COUNT(*) FROM clinvar_records;")
print("ClinVar records:", cursor.fetchone()[0])


print("\n=== CLINICAL SIGNIFICANCE ===")

cursor.execute("""
SELECT clinical_significance, COUNT(*)
FROM clinvar_records
GROUP BY clinical_significance
ORDER BY COUNT(*) DESC;
""")

for significance, count in cursor.fetchall():
    print(f"{significance}: {count}")


print("\n=== ASSEMBLY ===")

cursor.execute("""
SELECT assembly, COUNT(*)
FROM variant_representations
GROUP BY assembly
ORDER BY COUNT(*) DESC;
""")

for assembly, count in cursor.fetchall():
    print(f"{assembly}: {count}")


print("\n=== EXAMPLE VARIANTS ===")

cursor.execute("""
SELECT DISTINCT
    v.variant_id,
    v.allele_id,
    v.variation_id,
    v.hgvs_name,
    r.assembly,
    r.chromosome,
    r.start,
    r.stop,
    c.clinical_significance
FROM variants v
JOIN variant_representations r
    ON v.variant_id = r.variant_id
JOIN clinvar_records c
    ON v.variant_id = c.variant_id
LIMIT 10;
""")

for row in cursor.fetchall():
    print(row)


connection.close()