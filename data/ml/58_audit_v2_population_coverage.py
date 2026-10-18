"""Audit safe coverage limits of the current gnomAD population export."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"


def read(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def ids(value: object) -> set[str]:
    result = set()
    for item in re.split(r"[|;,\s]+", clean(value)):
        item = item.lower().strip()
        if item and item not in {"na", "-", "."}:
            result.add(item.removeprefix("rs"))
    return result


def main() -> None:
    clinvar = read(RAW / "SCN1A_clinvar.tsv", "\t")
    gnomad = read(RAW / "gnomad_scn1a_population.csv", ",")
    annotations = read(PROCESSED / "deepgene_v2_population_annotations.csv", ",")
    completeness = read(PROCESSED / "deepgene_v2_evidence_completeness.csv", ",")
    gnomad_by_rsid: dict[str, list[dict[str, str]]] = {}
    for row in gnomad:
        for rsid in ids(row.get("rsIDs")):
            gnomad_by_rsid.setdefault(rsid, []).append(row)
    clinvar_by_allele = {}
    for row in clinvar:
        clinvar_by_allele.setdefault(clean(row.get("#AlleleID")), row)

    reasons = Counter()
    unmatched_by_assembly = Counter()
    shared_rsid = 0
    allele_mismatch = 0
    for row in annotations:
        if clean(row.get("join_status")) != "unmatched":
            continue
        unmatched_by_assembly[clean(row.get("assembly")) or "missing_assembly"] += 1
        source = clinvar_by_allele.get(clean(row.get("variant_id")), {})
        candidates = []
        for rsid in ids(source.get("RS# (dbSNP)")):
            candidates.extend(gnomad_by_rsid.get(rsid, []))
        if not candidates:
            reasons["no_shared_rsid"] += 1
            continue
        shared_rsid += 1
        allele = (clean(row.get("reference_allele")).upper(), clean(row.get("alternate_allele")).upper())
        if not any((clean(item.get("Reference")).upper(), clean(item.get("Alternate")).upper()) == allele for item in candidates):
            reasons["shared_rsid_allele_mismatch"] += 1
            allele_mismatch += 1
        else:
            reasons["shared_rsid_requires_review"] += 1

    audit = {
        "created": date.today().isoformat(),
        "status": "API_ENRICHED_PARTIAL_SOURCE_SCOPE_LIMITED",
        "source": "gnomAD v4.1.2 local SCN1A export",
        "source_rows": len(gnomad),
        "source_chromosomes": sorted({clean(row.get("Chromosome")) for row in gnomad}),
        "source_transcripts": len({clean(row.get("Transcript")) for row in gnomad if clean(row.get("Transcript"))}),
        "representation_rows": len(annotations),
        "matched_representation_rows": sum(clean(row.get("join_status")) == "matched" for row in annotations),
        "representation_coverage": sum(clean(row.get("join_status")) == "matched" for row in annotations) / len(annotations) if annotations else 0.0,
        "variant_level_rows": len(completeness),
        "variant_level_matched_rows": sum(clean(row.get("population_evidence_status")) == "available" for row in completeness),
        "variant_level_coverage": sum(clean(row.get("population_evidence_status")) == "available" for row in completeness) / len(completeness) if completeness else 0.0,
        "unmatched_rows": sum(unmatched_by_assembly.values()),
        "unmatched_by_assembly": dict(unmatched_by_assembly),
        "unmatched_reason_counts": dict(reasons),
        "unmatched_with_shared_rsid": shared_rsid,
        "shared_rsid_allele_mismatches": allele_mismatch,
        "safe_additional_matches_from_current_export": 0,
        "reason": "The local export and targeted API/liftover enrichment are partial; remaining shared rsIDs with allele disagreement are excluded from matching.",
        "next_action": "Retry failed API batches with backoff and evaluate additional transcript or regional gnomAD coverage.",
        "safety_rule": "Do not match by rsID alone when reference and alternate alleles disagree; unmatched does not mean rare or benign.",
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "v2_population_unmatched_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
