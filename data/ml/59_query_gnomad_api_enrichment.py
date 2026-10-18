"""Query gnomAD v4 for unmatched GRCh38 SCN1A representations.

The script uses the public GraphQL API in small batches, preserves raw API
responses, and accepts a result only when the returned genomic alleles agree
with the requested variant. It never treats an API miss as a negative label.
"""

from __future__ import annotations

import csv
import json
import time
from datetime import date
from pathlib import Path
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"
API_URL = "https://gnomad.broadinstitute.org/api"
BATCH_SIZE = 20


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def variant_id(row: dict[str, str]) -> str:
    return "-".join([
        clean(row.get("chromosome")),
        clean(row.get("position")),
        clean(row.get("reference_allele")),
        clean(row.get("alternate_allele")),
    ])


def query_batch(ids: list[str]) -> dict:
    fields = "variantId reference_genome chrom pos ref alt rsids exome { ac an af } genome { ac an af }"
    selections = []
    for index, value in enumerate(ids):
        safe = value.replace('\\', '\\\\').replace('"', '\\"')
        selections.append(f'v{index}: variant(variantId: "{safe}", dataset: gnomad_r4) {{ {fields} }}')
    payload = json.dumps({"query": "query { " + " ".join(selections) + " }"}).encode()
    request = Request(API_URL, data=payload, headers={"Content-Type": "application/json", "User-Agent": "DeepGene-V2/1.0"})
    with urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def main() -> None:
    annotations = read_csv(PROCESSED / "deepgene_v2_population_annotations.csv")
    targets = []
    seen = set()
    for row in annotations:
        if clean(row.get("join_status")) == "unmatched" and clean(row.get("assembly")) == "GRCh38":
            requested = variant_id(row)
            if requested not in seen and all(requested.split("-")):
                seen.add(requested)
                targets.append({"requested_variant_id": requested, **row})

    raw_path = RAW / "gnomad_scn1a_api_enrichment.jsonl"
    results = []
    errors = []
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    with raw_path.open("w", encoding="utf-8") as raw_handle:
        for offset in range(0, len(targets), BATCH_SIZE):
            batch = targets[offset:offset + BATCH_SIZE]
            ids = [item["requested_variant_id"] for item in batch]
            try:
                response = query_batch(ids)
                raw_handle.write(json.dumps({"requested": ids, "response": response}) + "\n")
                data = response.get("data", {}) if isinstance(response, dict) else {}
                for index, item in enumerate(batch):
                    returned = data.get(f"v{index}")
                    if not returned:
                        results.append({"requested_variant_id": item["requested_variant_id"], "api_status": "not_found", "api_variant": None})
                        continue
                    requested_alleles = (clean(item.get("reference_allele")).upper(), clean(item.get("alternate_allele")).upper())
                    returned_alleles = (clean(returned.get("ref")).upper(), clean(returned.get("alt")).upper())
                    status = "matched" if requested_alleles == returned_alleles else "allele_conflict"
                    results.append({"requested_variant_id": item["requested_variant_id"], "api_status": status, "api_variant": returned})
            except Exception as exc:
                errors.append({"batch_start": offset, "requested": ids, "error": str(exc)})
                raw_handle.write(json.dumps({"requested": ids, "error": str(exc)}) + "\n")
            if offset + BATCH_SIZE < len(targets):
                time.sleep(0.15)

    matched = [item for item in results if item["api_status"] == "matched"]
    report = {
        "created": date.today().isoformat(),
        "status": "API_ENRICHMENT_COMPLETE" if not errors else "API_ENRICHMENT_PARTIAL",
        "api": API_URL,
        "dataset": "gnomad_r4",
        "target_rows": len(targets),
        "result_rows": len(results),
        "allele_consistent_matches": len(matched),
        "not_found": sum(item["api_status"] == "not_found" for item in results),
        "allele_conflicts": sum(item["api_status"] == "allele_conflict" for item in results),
        "request_errors": len(errors),
        "raw_response_path": str(raw_path.relative_to(ROOT)),
        "next_action": "Merge only allele-consistent API matches into the population annotation table and rerun coverage audit.",
        "training_label_policy": "API frequency annotations are not training labels.",
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "v2_gnomad_api_enrichment_report.json").write_text(json.dumps({**report, "errors": errors}, indent=2), encoding="utf-8")
    (PROCESSED / "deepgene_v2_gnomad_api_matches.json").write_text(json.dumps(matched, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
