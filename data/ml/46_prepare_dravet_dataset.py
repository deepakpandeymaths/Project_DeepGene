"""Build leakage-controlled Dravet training features from curated labels."""

from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
LABELS = ROOT / "data" / "raw" / "dravet_clinical_labels.csv"
VARIANTS = ROOT / "data" / "processed" / "deepgene_v1_final.csv"
OUT = ROOT / "data" / "processed" / "deepgene_dravet_ml_ready_v1.csv"

HYDRO = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
CHARGE = {"D":-1,"E":-1,"H":1,"K":1,"R":1}
VOLUME = {"A":88.6,"C":108.5,"D":111.1,"E":138.4,"F":189.9,"G":60.1,"H":153.2,"I":166.7,"K":168.6,"L":166.7,"M":162.9,"N":114.1,"P":112.7,"Q":143.8,"R":173.4,"S":89.0,"T":116.1,"V":140.0,"W":227.8,"Y":193.6}

def normalize(value):
    value = str(value).replace(" ", "")
    m = re.search(r"(?:SCN1A:)?(?:p\.)?([A-Za-z]{1,3})([0-9]+)([A-Za-z]{1,3})$", value)
    if not m:
        return None
    names = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
    a, p, b = m.groups()
    a, b = names.get(a.title(), a.upper()), names.get(b.title(), b.upper())
    return f"SCN1A:p.{a}{int(p)}{b}" if a in HYDRO and b in HYDRO else None

def main():
    labels = pd.read_csv(LABELS)
    required = {"variant_key", "dravet_label", "label_source", "cohort_id", "reviewer_status"}
    missing = required - set(labels.columns)
    if missing:
        raise ValueError(f"Missing label columns: {sorted(missing)}")
    labels["variant_key"] = labels["variant_key"].map(normalize)
    labels = labels[(labels["variant_key"].notna()) & (labels["dravet_label"].isin([0, 1])) & (labels["reviewer_status"].str.lower().str.contains("review", na=False))].copy()
    if labels.empty:
        print("BLOCKED: no reviewed Dravet labels are available. Fill data/raw/dravet_clinical_labels.csv first.")
        return
    if labels["variant_key"].duplicated().any():
        raise ValueError("Duplicate variant keys require adjudication before modeling.")
    variants = pd.read_csv(VARIANTS)
    rows = []
    for _, row in variants.iterrows():
        m = re.search(r"\(p\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})\)", str(row["hgvs_name"]))
        if not m:
            continue
        names = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
        ref, pos, alt = names.get(m.group(1)), int(m.group(2)), names.get(m.group(3))
        if not ref or not alt or ref not in HYDRO or alt not in HYDRO: continue
        rows.append({"variant_key": f"SCN1A:p.{ref}{pos}{alt}", "protein_position": pos, "reference_hydrophobicity": HYDRO[ref], "alternate_hydrophobicity": HYDRO[alt], "delta_hydrophobicity": HYDRO[alt]-HYDRO[ref], "absolute_delta_hydrophobicity": abs(HYDRO[alt]-HYDRO[ref]), "reference_charge": CHARGE.get(ref,0), "alternate_charge": CHARGE.get(alt,0), "delta_charge": CHARGE.get(alt,0)-CHARGE.get(ref,0), "absolute_delta_charge": abs(CHARGE.get(alt,0)-CHARGE.get(ref,0)), "reference_volume": VOLUME[ref], "alternate_volume": VOLUME[alt], "delta_volume": VOLUME[alt]-VOLUME[ref], "absolute_delta_volume": abs(VOLUME[alt]-VOLUME[ref])})
    features = pd.DataFrame(rows).merge(labels[["variant_key", "dravet_label", "cohort_id", "label_source"]], on="variant_key", how="inner")
    if features.empty:
        print("BLOCKED: no curated labels matched the DeepGene SCN1A variant table.")
        return
    OUT.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(OUT, index=False)
    print(f"Prepared {len(features)} labelled variants: {OUT}")

if __name__ == "__main__":
    main()
