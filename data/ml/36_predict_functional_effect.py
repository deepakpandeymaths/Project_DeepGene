"""Predict SCN1A functional effect for one missense substitution.

Example:
    python data/ml/36_predict_functional_effect.py --variant p.Arg1234Gly
"""

from pathlib import Path
import argparse
import json
import re
import pickle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "data" / "model" / "deepgene_functional_exploratory.pkl"

HYDRO = {"A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4,
         "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5,
         "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7,
         "V": 4.2, "W": -0.9, "Y": -1.3}
CHARGE = {"D": -1, "E": -1, "H": 1, "K": 1, "R": 1}
VOLUME = {"A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9,
          "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7,
          "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4,
          "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6}
THREE_TO_ONE = {"Ala":"A", "Arg":"R", "Asn":"N", "Asp":"D", "Cys":"C", "Gln":"Q",
                "Glu":"E", "Gly":"G", "His":"H", "Ile":"I", "Leu":"L", "Lys":"K",
                "Met":"M", "Phe":"F", "Pro":"P", "Ser":"S", "Thr":"T", "Trp":"W",
                "Tyr":"Y", "Val":"V"}


def parse_variant(raw):
    value = raw.strip().replace(" ", "").replace("p.", "")
    match = re.fullmatch(r"([A-Za-z]{1,3})([0-9]+)([A-Za-z]{1,3})", value)
    if not match:
        raise ValueError("Use a protein missense variant such as p.Arg1234Gly or R1234G.")
    ref, position, alt = match.groups()
    ref = THREE_TO_ONE.get(ref.title(), ref.upper())
    alt = THREE_TO_ONE.get(alt.title(), alt.upper())
    if ref not in HYDRO or alt not in HYDRO or ref == alt:
        raise ValueError("Reference and alternate residues must be different valid amino acids.")
    return ref, int(position), alt


def feature_row(ref, position, alt):
    rh, ah = HYDRO[ref], HYDRO[alt]
    rc, ac = CHARGE.get(ref, 0), CHARGE.get(alt, 0)
    rv, av = VOLUME[ref], VOLUME[alt]
    return pd.DataFrame([{
        "protein_position": position, "same_amino_acid": int(ref == alt),
        "reference_hydrophobicity": rh, "alternate_hydrophobicity": ah,
        "delta_hydrophobicity": ah - rh, "absolute_delta_hydrophobicity": abs(ah - rh),
        "reference_charge": rc, "alternate_charge": ac, "delta_charge": ac - rc,
        "absolute_delta_charge": abs(ac - rc), "reference_volume": rv,
        "alternate_volume": av, "delta_volume": av - rv, "absolute_delta_volume": abs(av - rv),
    }])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", required=True)
    args = parser.parse_args()
    if not ARTIFACT.exists():
        raise FileNotFoundError(f"Model not found. Run: python data/ml/35_build_exploratory_predictor.py")
    with ARTIFACT.open("rb") as handle:
        bundle = pickle.load(handle)
    ref, position, alt = parse_variant(args.variant)
    X = feature_row(ref, position, alt)[bundle["features"]]
    probabilities = bundle["model"].predict_proba(X)[0]
    classes = list(bundle["model"].classes_)
    probability_map = {label: float(probabilities[i]) for i, label in enumerate(classes)}
    entropy = float(-(probabilities * np.log(np.clip(probabilities, 1e-12, 1))).sum() / np.log(len(classes)))
    top_class = classes[int(np.argmax(probabilities))]
    top_probability = float(np.max(probabilities))
    output = {
        "variant": f"p.{ref}{position}{alt}", "predicted_functional_effect": top_class,
        "probabilities": probability_map, "top_probability": top_probability,
        "normalized_entropy": entropy,
        "abstain": bool(top_probability < 0.60 or entropy > 0.80),
        "deployment_allowed": False,
        "model_status": bundle["status"], "selected_model": bundle["selected_model"],
        "interpretation": "Research-only functional-effect estimate; do not interpret as pathogenicity, disease risk, diagnosis, or treatment advice.",
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
