"""Merge allele-consistent gnomAD matches found after GRCh37 liftover."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def number(value: object) -> float | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if 0 <= value <= 1 else None


def main() -> None:
    path = PROCESSED / "deepgene_v2_population_annotations.csv"
    matches = json.loads((PROCESSED / "deepgene_v2_liftover_gnomad_matches.json").read_text(encoding="utf-8"))
    by_source = {clean(item.get("source_variant_id")): item for item in matches}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0]) if rows else []
    merged = 0
    for row in rows:
        if clean(row.get("assembly")) != "GRCh37" or clean(row.get("join_status")) == "matched":
            continue
        item = by_source.get(clean(row.get("variant_id")))
        if not item:
            continue
        api = item.get("api_variant") or {}
        datasets = [api.get("exome") or {}, api.get("genome") or {}]
        frequencies = [number(dataset.get("af")) for dataset in datasets]
        frequencies = [value for value in frequencies if value is not None]
        if not frequencies:
            continue
        ac = sum(int(dataset.get("ac")) for dataset in datasets if str(dataset.get("ac", "")).isdigit())
        an = sum(int(dataset.get("an")) for dataset in datasets if str(dataset.get("an", "")).isdigit())
        row.update({
            "allele_frequency": f"{max(frequencies):.12g}",
            "max_population_af": f"{max(frequencies):.12g}",
            "allele_count": str(ac),
            "allele_number": str(an),
            "join_status": "matched",
            "match_method": "gnomad_api_liftover",
            "source_version": "gnomAD_v4_api_liftover",
            "missing_value_note": "Frequency retrieved after validated GRCh37-to-GRCh38 liftover; annotation only, not a training label.",
        })
        merged += 1
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"api_liftover_matches_available": len(matches), "merged_rows": merged, "output": str(path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
