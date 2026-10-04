"""Normalize the public SCN1A supplementary spreadsheet to a training CSV."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/SCN1A_Brunklaus_2020_supplementary.xlsx"
OUT = ROOT / "data/processed/deepgene_scn1a_brunklaus_2020_functional.csv"
SOURCE_URL = "https://discovery.ucl.ac.uk/id/eprint/10091164/10/Schorge_SCN1A%20variants%20from%20bench%20to%20bedside-improved%20clinical%20prediction%20from%20functional%20characterization_Spreadsheet.xlsx"
AA = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PATTERN = re.compile(r"^(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2})$")

source = pd.read_excel(RAW, sheet_name="Sheet1")
rows = []
for _, row in source[source["Function"].isin(["Loss", "Gain"])].iterrows():
    match = PATTERN.match(str(row["Variant"]).strip())
    if not match or match.group("ref") not in AA or match.group("alt") not in AA:
        continue
    rows.append({"id": f"brunklaus2020_{len(rows) + 1:03d}", "gene": "SCN1A", "y": "LOF" if row["Function"] == "Loss" else "GOF", "pheno": str(row.get("Phenotype", "")), "aa1": AA[match.group("ref")], "aa2": AA[match.group("alt")], "pos": int(match.group("pos")), "variant_key": str(row["Variant"]).strip(), "functional_label": str(row["Function"]), "label_status": "published_experimental_function", "source_name": "Brunklaus et al. 2020 supplementary functional characterization", "source_url": SOURCE_URL})
pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"Wrote {len(rows)} labelled rows to {OUT}")
