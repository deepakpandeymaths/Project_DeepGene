"""Create the leakage-audited Generation-2 sodium-channel dataset.

The SCION table is a published, openly available functional-effect dataset.
This step preserves source provenance, normalizes labels, creates an exact
protein-level identity key, and removes exact overlap with DeepGene V1 from
the Gen-2 training matrix.
"""

from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data" / "raw" / "scion_clean_tbl.csv"
V1 = ROOT / "data" / "processed" / "deepgene_ml_ready_v1.csv"
OUT = ROOT / "data" / "processed" / "deepgene_gen2_scion_harmonized_v1.csv"
ML_OUT = ROOT / "data" / "processed" / "deepgene_gen2_ml_ready_v1.csv"
REPORT = ROOT / "data" / "analysis_results" / "43_generation2_data_expansion.txt"
MANIFEST = ROOT / "data" / "raw" / "generation2_source_manifest.json"


def key(gene, ref, pos, alt):
    return f"{gene}:p.{ref}{int(pos)}{alt}"


def main():
    source = pd.read_csv(SOURCE)
    v1 = pd.read_csv(V1)
    required = {"id", "gene", "y", "aa1", "aa2", "pos"}
    if not required.issubset(source.columns):
        raise ValueError(f"SCION source is missing: {sorted(required - set(source.columns))}")

    source = source.copy()
    source["variant_key"] = [key(*row) for row in source[["gene", "aa1", "pos", "aa2"]].itertuples(index=False, name=None)]
    source["functional_label"] = source["y"].map({"LOF": "LOF", "GOF": "GOF"})
    source["label_status"] = source["functional_label"].map({"LOF": "primary", "GOF": "primary"}).fillna("excluded")
    v1_keys = {key("SCN1A", r, p, a) for r, p, a in v1[["reference_aa", "protein_position", "alternate_aa"]].itertuples(index=False, name=None)}
    source["overlap_with_deepgene_v1"] = source["variant_key"].isin(v1_keys)
    source["source_name"] = "SCION / Bosselmann et al. 2023"
    source["source_url"] = "https://github.com/christianbosselmann/SCION"

    harmonized = source.drop_duplicates("variant_key").reset_index(drop=True)
    train = harmonized[(harmonized["label_status"] == "primary") & (~harmonized["overlap_with_deepgene_v1"])].copy()

    # Minimal reproducible richer representation available directly from the
    # curated table: gene/task identity, position, residue properties, and
    # substitution deltas. More expensive conservation/structure features can
    # be added later without changing the identity or label policy.
    hydro = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
    charge = {"D":-1,"E":-1,"H":1,"K":1,"R":1}
    volume = {"A":88.6,"C":108.5,"D":111.1,"E":138.4,"F":189.9,"G":60.1,"H":153.2,"I":166.7,"K":168.6,"L":166.7,"M":162.9,"N":114.1,"P":112.7,"Q":143.8,"R":173.4,"S":89.0,"T":116.1,"V":140.0,"W":227.8,"Y":193.6}
    for aa, table in [("hydro", hydro), ("charge", charge), ("volume", volume)]:
        train[f"reference_{aa}"] = train["aa1"].map(table).fillna(0)
        train[f"alternate_{aa}"] = train["aa2"].map(table).fillna(0)
        train[f"delta_{aa}"] = train[f"alternate_{aa}"] - train[f"reference_{aa}"]
        train[f"absolute_delta_{aa}"] = train[f"delta_{aa}"].abs()
    train["protein_position"] = train["pos"].astype(int)
    train["same_amino_acid"] = (train["aa1"] == train["aa2"]).astype(int)
    gene_dummies = pd.get_dummies(train["gene"], prefix="gene", dtype=int)
    ml = pd.concat([train[["variant_key", "gene", "aa1", "aa2", "protein_position", "functional_label"]], train.filter(regex="^(reference_|alternate_|delta_|absolute_delta_|same_amino_acid$)"), gene_dummies], axis=1)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    harmonized.to_csv(OUT, index=False)
    ml.to_csv(ML_OUT, index=False)
    manifest = {"source": "SCION", "url": "https://github.com/christianbosselmann/SCION", "source_file": str(SOURCE), "source_rows": len(source), "source_genes": sorted(source.gene.unique().tolist()), "exact_v1_overlap": int(source.overlap_with_deepgene_v1.sum()), "gen2_training_rows_after_overlap_exclusion": len(train), "label_policy": {"LOF": "primary", "GOF": "primary", "mixed_or_unclear": "excluded"}, "independence_note": "SCION and DeepGene V1 are related literature-derived functional sources; exact variant overlaps are excluded, but source-level dependence cannot be ruled out."}
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    lines = ["DEEPGENE GENERATION-2 DATA EXPANSION", "", f"Source rows: {len(source)}", f"Unique variant keys: {len(harmonized)}", f"Genes: {', '.join(sorted(harmonized.gene.unique()))}", f"LOF: {(harmonized.functional_label == 'LOF').sum()}", f"GOF: {(harmonized.functional_label == 'GOF').sum()}", f"Exact overlaps with DeepGene V1: {int(source.overlap_with_deepgene_v1.sum())}", f"Gen-2 training rows after exact-overlap exclusion: {len(train)}", "", "IMPORTANT: Exact overlap was excluded. The source is still literature-related to V1, so this is expanded cross-channel training data, not a fully independent external test set.", "", f"Harmonized table: {OUT}", f"ML matrix: {ML_OUT}", f"Manifest: {MANIFEST}", "", "STATUS: PASS"]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
