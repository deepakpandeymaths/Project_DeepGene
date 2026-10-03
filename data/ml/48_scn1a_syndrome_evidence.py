"""Report documented SCN1A syndrome evidence without inventing probabilities.

This is an evidence lookup, not a Dravet predictor. It deliberately reports
what the local ClinVar-derived dataset says for a matching variant and leaves
clinical expression probabilities unavailable.
"""

from pathlib import Path
import argparse
import json
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "processed" / "deepgene_v1_final.csv"
THREE = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
ONE_TO_THREE = {value: key for key, value in THREE.items()}

def normalize(value):
    value = str(value).replace(" ", "")
    match = re.search(r"(?:p\.)?([A-Za-z]{1,3})([0-9]+)([A-Za-z]{1,3})$", value)
    if not match:
        raise ValueError("Use a protein variant such as p.Arg1648His or R1648H")
    ref, pos, alt = match.groups()
    ref = THREE.get(ref.title(), ref.upper())
    alt = THREE.get(alt.title(), alt.upper())
    return f"p.{ref}{int(pos)}{alt}"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", required=True)
    args = parser.parse_args()
    query = normalize(args.variant)
    query_match = f"p.{ONE_TO_THREE[query[2]]}{query[3:-1]}{ONE_TO_THREE[query[-1]]}"
    df = pd.read_csv(INPUT)
    hgvs = df["hgvs_name"].astype(str)
    matches = df[hgvs.str.contains(re.escape(query), case=False, regex=True, na=False) | hgvs.str.contains(re.escape(query_match), case=False, regex=True, na=False)].copy()
    if matches.empty:
        output = {"variant": query, "found_in_local_dataset": False, "documented_syndrome_mentions": [], "dravet_expression_probability": None, "interpretation": "No matching local ClinVar-derived record. No probability was estimated."}
        print(json.dumps(output, indent=2))
        return
    phenotype_text = "|".join(matches["phenotypes"].fillna("").astype(str).tolist())
    terms = sorted({term.strip() for term in phenotype_text.split("|") if term.strip() and term.strip().lower() not in {"not provided", "nan"}})
    dravet_terms = [term for term in terms if "dravet" in term.lower() or "severe myoclonic epilepsy in infancy" in term.lower()]
    output = {"variant": query, "found_in_local_dataset": True, "matching_records": len(matches), "clinical_significance_values": sorted(matches["clinical_significance"].dropna().astype(str).unique().tolist()), "documented_syndrome_mentions": terms, "dravet_term_mentioned": bool(dravet_terms), "dravet_related_terms": dravet_terms, "dravet_expression_probability": None, "interpretation": "These are documented phenotype/evidence mentions only. They are not a patient-specific diagnosis, penetrance estimate, or probability of developing Dravet syndrome."}
    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
