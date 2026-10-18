"""Query gnomAD v4 for validated GRCh37-to-GRCh38 lifted targets."""

from __future__ import annotations

import json
import time
from pathlib import Path
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"
API_URL = "https://gnomad.broadinstitute.org/api"
BATCH_SIZE = 20


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


def main() -> None:
    targets = json.loads((PROCESSED / "deepgene_v2_grch37_liftover_targets.json").read_text(encoding="utf-8"))
    results = []
    errors = []
    for offset in range(0, len(targets), BATCH_SIZE):
        batch = targets[offset:offset + BATCH_SIZE]
        ids = [item["lifted_variant_id"] for item in batch]
        try:
            response = query(ids)
            data = response.get("data", {})
            for index, item in enumerate(batch):
                returned = data.get(f"v{index}")
                if not returned:
                    continue
                alleles = (item["reference_allele"].upper(), item["alternate_allele"].upper())
                returned_alleles = (str(returned.get("ref", "")).upper(), str(returned.get("alt", "")).upper())
                if alleles == returned_alleles:
                    results.append({**item, "api_status": "matched", "api_variant": returned})
        except Exception as exc:
            errors.append({"batch_start": offset, "requested": ids, "error": str(exc)})
        if offset + BATCH_SIZE < len(targets):
            time.sleep(0.15)
    (PROCESSED / "deepgene_v2_liftover_gnomad_matches.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    report = {"status": "LIFTOVER_API_PARTIAL" if errors else "LIFTOVER_API_COMPLETE", "target_rows": len(targets), "allele_consistent_matches": len(results), "request_errors": len(errors), "errors": errors, "dataset": "gnomad_r4", "api": API_URL}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "v2_liftover_gnomad_api_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key != "errors"}, indent=2))


if __name__ == "__main__":
    main()
