"""Resumable, rate-limited retry runner for lifted gnomAD API requests.

The runner is intentionally conservative: it retries only request errors from
the recorded liftover audit, writes one checkpoint record per attempt, and
never treats a missing or allele-conflicting API response as a population
negative. Use ``--dry-run`` to inspect the pending workload without network
requests.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"
API_URL = "https://gnomad.broadinstitute.org/api"
AUDIT = RESULTS / "v2_liftover_gnomad_api_report.json"
TARGETS = PROCESSED / "deepgene_v2_grch37_liftover_targets.json"
CHECKPOINT = PROCESSED / "v2_liftover_api_checkpoint.jsonl"
MATCHES = PROCESSED / "deepgene_v2_liftover_api_retry_matches.json"
REPORT = RESULTS / "v2_liftover_api_retry_report.json"


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def failed_ids(audit: dict) -> set[str]:
    return {value for error in audit.get("errors", []) for value in error.get("requested", [])}


def load_latest_checkpoint(path: Path = CHECKPOINT) -> dict[str, dict]:
    latest = {}
    if not path.exists():
        return latest
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            try:
                item = json.loads(line)
            except json.JSONDecodeError:
                continue
            key = item.get("lifted_variant_id")
            if key:
                latest[key] = item
    return latest


def pending_targets(targets: list[dict], failed: set[str], latest: dict[str, dict], max_attempts: int) -> list[dict]:
    pending = []
    for item in targets:
        key = item.get("lifted_variant_id")
        if key not in failed:
            continue
        previous = latest.get(key, {})
        if previous.get("status") in {"matched", "not_found", "allele_conflict"}:
            continue
        if int(previous.get("attempt", 0)) >= max_attempts:
            continue
        pending.append(item)
    return sorted(pending, key=lambda item: item.get("lifted_variant_id", ""))


def query_one(item: dict, timeout: int) -> dict:
    value = item["lifted_variant_id"]
    safe = value.replace("\\", "\\\\").replace('"', '\\"')
    fields = "variantId reference_genome chrom pos ref alt rsids exome { ac an af } genome { ac an af }"
    body = json.dumps({"query": f'query {{ v: variant(variantId: "{safe}", dataset: gnomad_r4) {{ {fields} }} }}'}).encode()
    request = Request(API_URL, data=body, headers={"Content-Type": "application/json", "User-Agent": "DeepGene-V2/1.0"})
    try:
        with urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
        returned = (payload.get("data") or {}).get("v")
        if not returned:
            return {"status": "not_found"}
        expected = (item["reference_allele"].upper(), item["alternate_allele"].upper())
        actual = (str(returned.get("ref", "")).upper(), str(returned.get("alt", "")).upper())
        if expected != actual:
            return {"status": "allele_conflict", "api_variant": returned}
        return {"status": "matched", "api_variant": returned}
    except Exception as exc:
        return {"status": "request_error", "error": f"{type(exc).__name__}: {exc}"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-requests", type=int, default=0, help="Maximum new requests this run; 0 means all pending targets.")
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--sleep", type=float, default=0.5, help="Seconds between requests.")
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    audit = load_json(AUDIT)
    targets = load_json(TARGETS)
    latest = load_latest_checkpoint()
    pending = pending_targets(targets, failed_ids(audit), latest, args.max_attempts)
    selected = pending[:args.max_requests] if args.max_requests else pending
    plan = {"pending_targets": len(pending), "selected_targets": len(selected), "checkpoint": str(CHECKPOINT.relative_to(ROOT)), "dry_run": args.dry_run}
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return

    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    with CHECKPOINT.open("a", encoding="utf-8") as handle:
        for index, item in enumerate(selected):
            key = item["lifted_variant_id"]
            previous_attempt = int(latest.get(key, {}).get("attempt", 0))
            result = query_one(item, args.timeout)
            record = {
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "lifted_variant_id": key,
                "attempt": previous_attempt + 1,
                **result,
            }
            handle.write(json.dumps(record) + "\n")
            handle.flush()
            latest[key] = record
            if index + 1 < len(selected):
                time.sleep(max(0.0, args.sleep))

    matched = []
    for item in targets:
        record = latest.get(item.get("lifted_variant_id"), {})
        if record.get("status") == "matched":
            matched.append({**item, "api_status": "matched", "api_variant": record.get("api_variant")})
    MATCHES.write_text(json.dumps(matched, indent=2), encoding="utf-8")
    statuses = Counter(item.get("status", "unknown") for item in latest.values())
    report = {**plan, "status_counts": dict(statuses), "matches": len(matched), "next_action": "Run data/ml/64_merge_liftover_gnomad_matches.py only after reviewing the checkpoint and audit."}
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
