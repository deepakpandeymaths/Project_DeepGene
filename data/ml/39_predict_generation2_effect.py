"""Predict LOF/GOF with the Generation-2 cross-channel candidate."""

from pathlib import Path
import argparse
import json
import pickle
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "data" / "model" / "deepgene_generation2_gbm.pkl"
HYDRO = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
CHARGE = {"D":-1,"E":-1,"H":1,"K":1,"R":1}
VOLUME = {"A":88.6,"C":108.5,"D":111.1,"E":138.4,"F":189.9,"G":60.1,"H":153.2,"I":166.7,"K":168.6,"L":166.7,"M":162.9,"N":114.1,"P":112.7,"Q":143.8,"R":173.4,"S":89.0,"T":116.1,"V":140.0,"W":227.8,"Y":193.6}
THREE = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}


def parse(value):
    import re
    m = re.fullmatch(r"(?:p\.)?([A-Za-z]{1,3})([0-9]+)([A-Za-z]{1,3})", value.replace(" ", ""))
    if not m:
        raise ValueError("Use p.Arg1234Gly or R1234G")
    ref, pos, alt = m.groups()
    ref, alt = THREE.get(ref.title(), ref.upper()), THREE.get(alt.title(), alt.upper())
    if ref not in HYDRO or alt not in HYDRO or ref == alt:
        raise ValueError("Invalid or identical amino-acid residues")
    return ref, int(pos), alt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--variant", required=True)
    args = parser.parse_args()
    with ARTIFACT.open("rb") as handle:
        bundle = pickle.load(handle)
    ref, pos, alt = parse(args.variant)
    row = {"protein_position": pos, "same_amino_acid": int(ref == alt)}
    for name, table in [("hydro", HYDRO), ("charge", CHARGE), ("volume", VOLUME)]:
        row[f"reference_{name}"] = table.get(ref, 0)
        row[f"alternate_{name}"] = table.get(alt, 0)
        row[f"delta_{name}"] = row[f"alternate_{name}"] - row[f"reference_{name}"]
        row[f"absolute_delta_{name}"] = abs(row[f"delta_{name}"])
    for gene in bundle["genes"]:
        row[f"gene_{gene}"] = int(args.gene == gene)
    X = pd.DataFrame([row]).reindex(columns=bundle["features"], fill_value=0)
    probs = bundle["model"].predict_proba(X)[0]
    classes = list(bundle["model"].classes_)
    output = {"gene": args.gene, "variant": f"p.{ref}{pos}{alt}", "predicted_functional_effect": classes[int(probs.argmax())], "probabilities": {c: float(p) for c, p in zip(classes, probs)}, "model_status": bundle["status"], "deployment_allowed": False, "interpretation": "Generation-2 cross-channel research estimate; not a clinical diagnosis or pathogenicity prediction."}
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
