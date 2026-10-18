"""Merge allele-consistent gnomAD API results into the annotation table."""

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
    annotation_path = PROCESSED / "deepgene_v2_population_annotations.csv"
    matches = json.loads((PROCESSED / "deepgene_v2_gnomad_api_matches.json").read_text(encoding="utf-8"))
    by_id = {clean(item.get("requested_variant_id")): item for item in matches}
    with annotation_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0]) if rows else []
    merged = 0
    for row in rows:
        item = by_id.get("-".join([clean(row.get("chromosome")), clean(row.get("position")), clean(row.get("reference_allele")), clean(row.get("alternate_allele"))]))
        if not item or item.get("api_status") != "matched":
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
            "match_method": "gnomad_api_exact",
            "source_version": "gnomAD_v4_api",
            "missing_value_note": "Frequency retrieved from the official gnomAD GraphQL API; annotation only, not a training label.",
        })
        merged += 1
    with annotation_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"api_matches_available": len(matches), "merged_rows": merged, "output": str(annotation_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
