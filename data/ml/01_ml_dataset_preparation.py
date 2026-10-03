import pandas as pd
from pathlib import Path

# ==========================================================
# Paths
# ==========================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = (
    PROJECT_ROOT /
    "data" /
    "processed" /
    "deepgene_v1_final.csv"
)

OUTPUT_DIR = PROJECT_ROOT / "data" / "processed"
REPORT_DIR = PROJECT_ROOT / "data" / "ml_results"

OUTPUT_DIR.mkdir(exist_ok=True)
REPORT_DIR.mkdir(exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "deepgene_ml_ready_v1.csv"
REPORT_FILE = REPORT_DIR / "01_ml_dataset_preparation.txt"

# ==========================================================
# Load
# ==========================================================

df = pd.read_csv(INPUT_FILE)

report = []

def line(x):
    report.append(str(x))

line("=" * 70)
line("PROJECT DEEPGENE")
line("STEP 01 — ML DATASET PREPARATION")
line("=" * 70)

line(f"Input rows: {len(df)}")
line(f"Input columns: {len(df.columns)}")
line("")

# ==========================================================
# Remove identifiers
# ==========================================================

remove_cols = [
    "variant_id",
    "gene_id",
    "allele_id",
    "variation_id",
    "hgvs_name",
    "phenotypes"
]

existing = [c for c in remove_cols if c in df.columns]

df_ml = df.drop(columns=existing)

line("Removed columns:")
for c in existing:
    line(f"  - {c}")

line("")
line(f"Remaining features: {len(df_ml.columns)}")

# ==========================================================
# Missing values
# ==========================================================

missing = df_ml.isnull().sum().sum()

line("")
line(f"Total missing values: {missing}")

# ==========================================================
# Save
# ==========================================================

df_ml.to_csv(OUTPUT_FILE, index=False)

with open(REPORT_FILE, "w") as f:
    for r in report:
        f.write(r + "\n")

# ==========================================================
# Console output
# ==========================================================

for r in report:
    print(r)

print("")
print("ML-ready dataset:")
print(OUTPUT_FILE)

print("")
print("Report:")
print(REPORT_FILE)

print("")
print("STATUS: PASS")
print("=" * 70)