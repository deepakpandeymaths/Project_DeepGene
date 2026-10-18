"""Retry failed gnomAD API batches by splitting them into smaller requests."""

from __future__ import annotations

import json
import time
from pathlib import Path
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
API_URL = "https://gnomad.broadinstitute.org/api"


def query(ids: list[str]) -> dict:
    fields = "variantId reference_genome chrom pos ref alt rsids exome { ac an af } genome { ac an af }"
    selections = []
    for index, value in enumerate(ids):
        safe = value.replace('\\', '\\\\').replace('"', '\\"')
        selections.append(f'v{index}: variant(variantId: "{safe}", dataset: gnomad_r4) {{ {fields} }}')
    body = json.dumps({"query": "query { " + " ".join(selections) + " }"}).encode()
    request = Request(API_URL, data=body, headers={"Content-Type": "application/json", "User-Agent": "DeepGene-V2/1.0"})
    with urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def retry(ids: list[str], results: list[dict], errors: list[dict]) -> None:
    try:
        response = query(ids)
        data = response.get("data", {})
        for index, requested in enumerate(ids):
            returned = data.get(f"v{index}")
            if returned:
                requested_parts = requested.split("-")
                requested_alleles = (requested_parts[2].upper(), requested_parts[3].upper())
                returned_alleles = (str(returned.get("ref", "")).upper(), str(returned.get("alt", "")).upper())
                if requested_alleles == returned_alleles:
                    results.append({"requested_variant_id": requested, "api_status": "matched", "api_variant": returned})
                else:
                    errors.append({"requested_variant_id": requested, "status": "allele_conflict"})
            else:
                errors.append({"requested_variant_id": requested, "status": "not_found"})
        return
    except Exception as exc:
        if len(ids) == 1:
            errors.append({"requested_variant_id": ids[0], "status": "request_error", "error": str(exc)})
            return
        middle = max(1, len(ids) // 2)
        retry(ids[:middle], results, errors)
        retry(ids[middle:], results, errors)


def main() -> None:
    failed_batches = []
    with (RAW / "gnomad_scn1a_api_enrichment.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            item = json.loads(line)
            if item.get("error") and item.get("requested"):
                failed_batches.append(item["requested"])
    targets = [variant for batch in failed_batches for variant in batch]
    results: list[dict] = []
    errors: list[dict] = []
    for index, batch in enumerate(failed_batches):
        retry(batch, results, errors)
        if index + 1 < len(failed_batches):
            time.sleep(0.15)
    (PROCESSED / "deepgene_v2_gnomad_api_retry_matches.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (PROCESSED / "deepgene_v2_gnomad_api_retry_report.json").write_text(json.dumps({"failed_batches": len(failed_batches), "retry_targets": len(targets), "matches": len(results), "errors": len(errors), "error_details": errors}, indent=2), encoding="utf-8")
    print(json.dumps({"failed_batches": len(failed_batches), "retry_targets": len(targets), "retry_matches": len(results), "retry_errors": len(errors)}, indent=2))


if __name__ == "__main__":
    main()
