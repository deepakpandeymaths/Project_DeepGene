"""Create deterministic V2 grouped splits and cross-source overlap audits."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"
ML_RESULTS = ROOT / "data" / "ml_results"


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


def split_for_group(group_key: str) -> str:
    digest = hashlib.sha256(group_key.encode("utf-8")).hexdigest()
    bucket = int(digest[:8], 16) % 100
    if bucket < 80:
        return "train"
    if bucket < 90:
        return "validation"
    return "test"


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    ML_RESULTS.mkdir(parents=True, exist_ok=True)
    clinical = read_csv(PROCESSED / "deepgene_v2_clinical_evidence.csv")
    master = read_csv(PROCESSED / "deepgene_v2_master_variants.csv")
    functional = read_csv(PROCESSED / "deepgene_v2_functional_evidence.csv")

    master_by_id = {clean(row.get("variant_id")): row for row in master}
    split_rows: list[dict[str, object]] = []
    group_products: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))

    for row in clinical:
        variant_id = clean(row.get("variant_id"))
        master_row = master_by_id.get(variant_id, {})
        protein_key = clean(master_row.get("protein_key"))
        group_key = protein_key or f"CLINVAR_VARIANT:{variant_id}"
        split = split_for_group(group_key)
        group_products[group_key]["clinical"].add(variant_id)
        split_rows.append({
            "product": "clinical_evidence",
            "row_id": variant_id,
            "source_id": clean(row.get("clinical_source_id")),
            "group_key": group_key,
            "split": split,
            "target_policy": "clinical labels are not used as features for same-target prediction",
        })

    for row in functional:
        source_id = clean(row.get("source_id"))
        row_id = clean(row.get("source_row_id"))
        group_key = clean(row.get("variant_key")) or f"{source_id}:{row_id}"
        split = split_for_group(group_key)
        group_products[group_key]["functional"].add(row_id)
        split_rows.append({
            "product": "functional_evidence",
            "row_id": row_id,
            "source_id": source_id,
            "group_key": group_key,
            "split": split,
            "target_policy": "functional labels are targets and are excluded from feature matrices",
        })

    fields = ["product", "row_id", "source_id", "group_key", "split", "target_policy"]
    write_csv(ML_RESULTS / "v2_grouped_split_ids.csv", split_rows, fields)

    cross_product_groups = {
        key: {product: sorted(ids) for product, ids in products.items()}
        for key, products in group_products.items()
        if len(products) > 1
    }
    split_group_sets = defaultdict(set)
    for row in split_rows:
        split_group_sets[row["split"]].add(row["group_key"])
    cross_split_groups = set()
    for group_key in set(row["group_key"] for row in split_rows):
        assigned = {row["split"] for row in split_rows if row["group_key"] == group_key}
        if len(assigned) > 1:
            cross_split_groups.add(group_key)

    source_groups = defaultdict(set)
    for row in functional:
        source_groups[clean(row.get("source_id"))].add(clean(row.get("variant_key")))
    source_ids = sorted(source_groups)
    source_overlap_counts = {}
    for index, left in enumerate(source_ids):
        for right in source_ids[index + 1:]:
            source_overlap_counts[f"{left}__{right}"] = len(source_groups[left] & source_groups[right])

    overlap_audit = {
        "audit_date": date.today().isoformat(),
        "status": "PASS_WITH_REPORTED_OVERLAPS" if not cross_split_groups else "FAIL",
        "source_row_counts": {
            "clinical_evidence": len(clinical),
            "functional_evidence": len(functional),
        },
        "source_pair_overlap_counts": source_overlap_counts,
        "clinical_functional_shared_group_count": sum(1 for products in cross_product_groups.values() if "clinical" in products and "functional" in products),
        "clinical_functional_shared_groups": cross_product_groups,
        "groups_assigned_to_multiple_splits": sorted(cross_split_groups),
        "interpretation": "Shared variants are retained for audit but assigned to one deterministic split across products. They must not be treated as independent observations.",
    }
    (RESULTS / "v2_cross_source_overlap_audit.json").write_text(json.dumps(overlap_audit, indent=2), encoding="utf-8")

    split_counts = defaultdict(lambda: defaultdict(int))
    for row in split_rows:
        split_counts[row["product"]][row["split"]] += 1
    split_policy = {
        "created": date.today().isoformat(),
        "algorithm": "SHA-256 of normalized biological group key; first 32 bits modulo 100",
        "assignment": {"train": "0-79", "validation": "80-89", "test": "90-99"},
        "grouping": {
            "clinical_evidence": "protein_key when available, otherwise ClinVar variant_id",
            "functional_evidence": "variant_key when available, otherwise source_id and source_row_id",
        },
        "guarantees": [
            "All rows for one group_key receive one split.",
            "The split is deterministic and reproducible without a mutable random seed.",
            "Clinical and functional products remain separate even when they share a group_key.",
            "No clinical or functional target value is used to assign a split.",
        ],
        "row_counts": {product: dict(counts) for product, counts in split_counts.items()},
        "group_count_by_split": {split: len(groups) for split, groups in split_group_sets.items()},
        "next_validation_requirement": "Future model scripts must consume these split IDs or document a stronger source-aware external split.",
    }
    (RESULTS / "v2_split_policy.json").write_text(json.dumps(split_policy, indent=2), encoding="utf-8")

    print(json.dumps({
        "split_rows": len(split_rows),
        "clinical_rows": len(clinical),
        "functional_rows": len(functional),
        "shared_group_count": overlap_audit["clinical_functional_shared_group_count"],
        "groups_in_multiple_splits": len(cross_split_groups),
        "outputs": [
            "data/ml_results/v2_grouped_split_ids.csv",
            "data/analysis_results/v2_cross_source_overlap_audit.json",
            "data/analysis_results/v2_split_policy.json",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
