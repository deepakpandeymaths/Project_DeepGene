"""Stage conservative HPO mappings from the existing phenotype descriptions.

Only exact, high-confidence seed mappings are emitted as HPO IDs. Disease
names, clinical labels, and unmapped terms remain pending until an official,
versioned HPO mapping resource is added.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"

# Conservative exact-term seed. This is not a replacement for an official
# HPO release and is kept separate from model features until externally checked.
EXACT_HPO_SEED = {
    "seizure": ("HP:0001250", "Seizure"),
}


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


def normalize_term(value: str) -> str:
    return re.sub(r"\s+", " ", clean(value)).strip().lower()


def main() -> None:
    source = PROCESSED / "deepgene_phenotype_normalized_v1.csv"
    rows = read_csv(source)
    annotations = []
    status_counts = Counter()
    unique_variants = set()
    mapped_variants = set()
    for row in rows:
        variant_id = clean(row.get("variant_id"))
        terms = [term for term in clean(row.get("meaningful_phenotype_terms")).split("|") if clean(term)]
        if not terms:
            terms = [term for term in clean(row.get("phenotypes")).split("|") if clean(term)]
        for term in terms:
            normalized = normalize_term(term)
            mapping = EXACT_HPO_SEED.get(normalized)
            status = "exact_seed_candidate" if mapping else "unmapped_pending_official_hpo"
            hpo_id = mapping[0] if mapping else ""
            hpo_label = mapping[1] if mapping else ""
            annotations.append({
                "variant_id": variant_id,
                "gene": "SCN1A",
                "source_phenotype_term": term,
                "normalized_term": normalized,
                "hpo_id": hpo_id,
                "hpo_label": hpo_label,
                "mapping_status": status,
                "mapping_source": "DeepGene exact-term seed; requires official HPO release verification" if mapping else "none",
                "hpo_release": "pending_official_release",
                "safe_for_training": 0,
                "note": "Candidate mapping is not model-ready until verified against a versioned HPO resource." if mapping else "No HPO ID inferred from this term.",
            })
            status_counts[status] += 1
            if variant_id:
                unique_variants.add(variant_id)
                if mapping:
                    mapped_variants.add(variant_id)

    fields = ["variant_id", "gene", "source_phenotype_term", "normalized_term", "hpo_id", "hpo_label", "mapping_status", "mapping_source", "hpo_release", "safe_for_training", "note"]
    write_csv(PROCESSED / "deepgene_v2_hpo_annotations.csv", annotations, fields)
    audit = {
        "created": date.today().isoformat(),
        "status": "STAGED_PENDING_OFFICIAL_HPO_VERIFICATION",
        "input_file": str(source.relative_to(ROOT)),
        "input_rows": len(rows),
        "annotation_rows": len(annotations),
        "unique_variants_seen": len(unique_variants),
        "unique_variants_with_seed_mapping": len(mapped_variants),
        "mapping_status_counts": dict(status_counts),
        "seed_mappings": {term: {"hpo_id": value[0], "label": value[1]} for term, value in EXACT_HPO_SEED.items()},
        "rules": [
            "Disease names and clinical classifications are not automatically treated as HPO phenotypes.",
            "Unmapped terms remain visible and are not assigned guessed ontology IDs.",
            "Seed mappings are not safe for training until verified against a versioned official HPO resource.",
            "HPO normalization does not create clinical labels or functional LOF/GOF labels.",
        ],
        "next_gate": "Add an official HPO release or validated mapping file, then audit mapped and unmatched terms before updating evidence completeness.",
    }
    (RESULTS / "v2_hpo_join_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps({
        "annotation_rows": len(annotations),
        "unique_variants_seen": len(unique_variants),
        "unique_variants_with_seed_mapping": len(mapped_variants),
        "mapping_status_counts": dict(status_counts),
        "outputs": [
            "data/processed/deepgene_v2_hpo_annotations.csv",
            "data/analysis_results/v2_hpo_join_audit.json",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
