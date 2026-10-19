"""Build transparent, research-only V2 variant prioritization."""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
ML_RESULTS = ROOT / "data" / "ml_results"
ANALYSIS = ROOT / "data" / "analysis_results"


AA_ONE_TO_THREE = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}
AA_THREE_TO_THREE = {value.lower(): value for value in AA_ONE_TO_THREE.values()}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def number(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def canonical_amino_acid(token: str) -> str:
    token = clean(token)
    if len(token) == 1 and token.upper() in AA_ONE_TO_THREE:
        return AA_ONE_TO_THREE[token.upper()]
    if token.lower() in AA_THREE_TO_THREE:
        return AA_THREE_TO_THREE[token.lower()]
    if token in {"Ter", "*"}:
        return "Ter"
    return ""


def canonical_variant_key(value: str) -> str:
    value = clean(value)
    match = re.match(
        r"^([^:]+):p\.([A-Z][a-z]{2}|[A-Z])(\d+)(Ter|\*|[A-Z][a-z]{2}|[A-Z])$",
        value,
    )
    if not match:
        return value
    reference = canonical_amino_acid(match.group(2))
    alternate = canonical_amino_acid(match.group(4))
    if not reference or not alternate:
        return value
    return f"{match.group(1)}:p.{reference}{match.group(3)}{alternate}"


def main() -> None:
    predictions = read_csv(ML_RESULTS / "v2_uncertainty_predictions.csv")
    master = read_csv(PROCESSED / "deepgene_v2_master_variants.csv")
    completeness = read_csv(PROCESSED / "deepgene_v2_evidence_completeness.csv")

    clinical_by_key: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in master:
        key = clean(row.get("protein_key"))
        if key:
            clinical_by_key[key].append(row)
    completeness_by_id = {clean(row.get("variant_id")): row for row in completeness}

    selected = [row for row in predictions if clean(row.get("model")) == "gradient_boosting"]
    output_fields = [
        "rank", "source_row_id", "variant_key", "split", "functional_source_ids", "clinical_variant_ids",
        "observed_target", "predicted_target", "gof_probability", "confidence", "uncertainty_level",
        "evidence_completeness_fraction", "evidence_completeness_category", "functional_signal_component",
        "uncertainty_component", "evidence_gap_component", "conflict_component", "research_priority_score",
        "priority_band", "population_evidence_status", "population_max_af", "supporting_evidence",
        "missing_evidence", "interpretation_note",
    ]
    prioritized: list[dict[str, object]] = []
    for row in selected:
        key = canonical_variant_key(row.get("variant_key", ""))
        clinical_rows = clinical_by_key.get(key, [])
        clinical_ids = sorted({clean(item.get("variant_id")) for item in clinical_rows if clean(item.get("variant_id"))})
        completeness_rows = [completeness_by_id[item] for item in clinical_ids if item in completeness_by_id]
        completeness_fraction = max((number(item.get("evidence_completeness_fraction")) for item in completeness_rows), default=0.0)
        completeness_category = "|".join(sorted({clean(item.get("evidence_completeness_category")) for item in completeness_rows})) or "not_joined"
        conflict = int(any(clean(item.get("conflict_flag")).lower() in {"1", "true", "yes"} for item in clinical_rows))
        probability = number(row.get("gof_probability"))
        confidence = number(row.get("confidence"))
        entropy = number(row.get("entropy_bits"))
        functional_signal = confidence
        uncertainty_component = max(0.0, min(1.0, entropy))
        evidence_gap = max(0.0, min(1.0, 1.0 - completeness_fraction))
        conflict_component = float(conflict)
        score = 0.45 * functional_signal + 0.20 * uncertainty_component + 0.25 * evidence_gap + 0.10 * conflict_component
        if score >= 0.75:
            band = "high_research_priority"
        elif score >= 0.55:
            band = "moderate_research_priority"
        else:
            band = "lower_research_priority"
        functional_sources = sorted({clean(item.get("functional_source_ids")) for item in clinical_rows if clean(item.get("functional_source_ids"))})
        supporting = ["explicit functional assay label", "calibrated functional-effect probability"]
        if clinical_rows:
            supporting.append("ClinVar representation and evidence context")
        if conflict:
            supporting.append("clinical conflict flag")
        missing = []
        if not clinical_rows:
            missing.append("clinical variant join")
        else:
            for completeness_row in completeness_rows:
                missing.extend(clean(completeness_row.get("missing_evidence_categories")).split("|"))
            missing = sorted({item for item in missing if item})
        population_statuses = sorted({clean(item.get("population_evidence_status")) for item in completeness_rows if clean(item.get("population_evidence_status"))})
        population_status = "|".join(population_statuses) or "not_joined"
        population_values = [number(item.get("population_max_af"), -1.0) for item in completeness_rows]
        population_values = [value for value in population_values if 0.0 <= value <= 1.0]
        prioritized.append({
            "rank": 0,
            "source_row_id": clean(row.get("source_row_id")),
            "variant_key": key,
            "split": clean(row.get("split")),
            "functional_source_ids": "|".join(functional_sources),
            "clinical_variant_ids": "|".join(clinical_ids),
            "observed_target": clean(row.get("observed_target")),
            "predicted_target": clean(row.get("predicted_target")),
            "gof_probability": clean(row.get("gof_probability")),
            "confidence": clean(row.get("confidence")),
            "uncertainty_level": clean(row.get("uncertainty_level")),
            "evidence_completeness_fraction": f"{completeness_fraction:.6f}",
            "evidence_completeness_category": completeness_category,
            "functional_signal_component": f"{functional_signal:.6f}",
            "uncertainty_component": f"{uncertainty_component:.6f}",
            "evidence_gap_component": f"{evidence_gap:.6f}",
            "conflict_component": f"{conflict_component:.6f}",
            "research_priority_score": f"{score:.6f}",
            "priority_band": band,
            "population_evidence_status": population_status,
            "population_max_af": f"{max(population_values):.12g}" if population_values else "",
            "supporting_evidence": "|".join(supporting),
            "missing_evidence": "|".join(missing),
            "interpretation_note": "Research prioritization only; not a clinical classification, diagnosis, risk estimate, or treatment recommendation.",
        })

    prioritized.sort(key=lambda item: (-number(item["research_priority_score"]), clean(item["variant_key"])))
    for rank, row in enumerate(prioritized, start=1):
        row["rank"] = rank
    write_csv(PROCESSED / "deepgene_v2_research_prioritization.csv", prioritized, output_fields)

    explanation_rows = []
    for row in prioritized[:25]:
        explanation_rows.append({
            "rank": row["rank"],
            "variant_key": row["variant_key"],
            "priority_band": row["priority_band"],
            "score": row["research_priority_score"],
            "explanation": f"Ranked using calibrated functional signal ({row['functional_signal_component']}), model uncertainty ({row['uncertainty_component']}), evidence gap ({row['evidence_gap_component']}), and clinical conflict component ({row['conflict_component']}). Missing evidence remains unobserved, not negative.",
        })
    (ML_RESULTS / "v2_prioritization_explanations.json").write_text(json.dumps({
        "created": date.today().isoformat(),
        "model_reference": "gradient_boosting calibrated functional baseline",
        "rows_ranked": len(prioritized),
        "top_ranked_explanations": explanation_rows,
    }, indent=2), encoding="utf-8")

    report = {
        "created": date.today().isoformat(),
        "status": "RESEARCH_ONLY_PRIORITIZATION",
        "rows_ranked": len(prioritized),
        "clinical_joined_rows": sum(1 for row in prioritized if clean(row["clinical_variant_ids"])),
        "population_available_rows": sum(1 for row in prioritized if "available" in clean(row["population_evidence_status"])),
        "priority_band_counts": {band: sum(1 for row in prioritized if row["priority_band"] == band) for band in sorted({row["priority_band"] for row in prioritized})},
        "split_counts": {split: sum(1 for row in prioritized if row["split"] == split) for split in sorted({row["split"] for row in prioritized})},
        "formula": "0.45 functional signal + 0.20 uncertainty + 0.25 evidence gap + 0.10 conflict component",
        "join_policy": "Canonicalize one-letter and three-letter protein HGVS tokens before exact protein-key matching; do not use position-only or amino-acid-only matches.",
        "interpretation": "This ranks variants for research review. It does not estimate pathogenicity, disease risk, diagnosis, treatment response, or clinical action.",
        "limitations": [
            "Only evaluated functional rows with calibrated predictions are ranked.",
            "Clinical joins use protein-level keys and may contain multiple nucleotide representations.",
            "Missing external evidence is not interpreted as benign or absent disease.",
            "The functional model is source-dependent and research-only.",
        ],
    }
    (ANALYSIS / "v2_prioritization_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "rows_ranked": len(prioritized),
        "priority_band_counts": report["priority_band_counts"],
        "outputs": [
            "data/processed/deepgene_v2_research_prioritization.csv",
            "data/ml_results/v2_prioritization_explanations.json",
            "data/analysis_results/v2_prioritization_report.json",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
