"""Lift unmatched GRCh37 representations to GRCh38 for targeted gnomAD queries."""

from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import date
from pathlib import Path

from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"
CHAIN = RAW / "hg19ToHg38.over.chain.gz"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def main() -> None:
    lift = LiftOver(str(CHAIN))
    annotations = read_csv(PROCESSED / "deepgene_v2_population_annotations.csv")
    targets = []
    seen = set()
    audit = Counter()
    for row in annotations:
        if clean(row.get("join_status")) != "unmatched" or clean(row.get("assembly")) != "GRCh37":
            continue
        chromosome = clean(row.get("chromosome"))
        position = clean(row.get("position"))
        reference = clean(row.get("reference_allele")).upper()
        alternate = clean(row.get("alternate_allele")).upper()
        try:
            lifted = lift.convert_coordinate(f"chr{chromosome}", int(position) - 1)
        except (TypeError, ValueError):
            lifted = None
        if not lifted:
            audit["no_liftover_result"] += 1
            continue
        if len(lifted) != 1:
            audit["ambiguous_liftover"] += 1
            continue
        target_chromosome, target_position_zero, target_strand, _ = lifted[0]
        if target_strand != "+":
            audit["negative_strand_excluded"] += 1
            continue
        target = f"{target_chromosome.removeprefix('chr')}-{target_position_zero + 1}-{reference}-{alternate}"
        if target in seen:
            audit["duplicate_lifted_target"] += 1
            continue
        seen.add(target)
        targets.append({
            "source_variant_id": clean(row.get("variant_id")),
            "source_variation_id": clean(row.get("variation_id")),
            "source_assembly": "GRCh37",
            "source_chromosome": chromosome,
            "source_position": position,
            "reference_allele": reference,
            "alternate_allele": alternate,
            "lifted_variant_id": target,
            "target_assembly": "GRCh38",
            "chain_file": str(CHAIN.relative_to(ROOT)),
        })
        audit["lifted_unique_targets"] += 1
    PROCESSED.mkdir(parents=True, exist_ok=True)
    (PROCESSED / "deepgene_v2_grch37_liftover_targets.json").write_text(json.dumps(targets, indent=2), encoding="utf-8")
    report = {
        "created": date.today().isoformat(),
        "status": "LIFTOVER_COMPLETE",
        "chain_file": str(CHAIN.relative_to(ROOT)),
        "input_unmatched_grch37_rows": sum(clean(row.get("join_status")) == "unmatched" and clean(row.get("assembly")) == "GRCh37" for row in annotations),
        "audit": dict(audit),
        "target_output": str((PROCESSED / "deepgene_v2_grch37_liftover_targets.json").relative_to(ROOT)),
        "safety_rule": "Lifted coordinates are candidates for API matching only; reference and alternate alleles must still agree.",
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "v2_grch37_liftover_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
