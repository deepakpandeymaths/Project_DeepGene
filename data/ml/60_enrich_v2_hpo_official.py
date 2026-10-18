"""Enrich staged phenotype terms with exact matches from an official HPO OBO release.

Only exact normalized HPO names and exact normalized synonyms are accepted.
Ambiguous terms, disease names, and approximate matches remain unresolved.
The resulting HPO IDs are annotations only and are not model targets.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"
DEFAULT_ONTOLOGY = ROOT / "data" / "raw" / "hp.obo"


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def normalize(value: str) -> str:
    return re.sub(r"\s+", " ", clean(value)).lower().strip()


def parse_obo(path: Path) -> tuple[dict[str, list[tuple[str, str]]], str]:
    candidates: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)
    release = "unknown"
    current: dict[str, object] = {}

    def emit() -> None:
        if current.get("obsolete") or not current.get("id") or not current.get("name"):
            return
        hpo_id = str(current["id"])
        names = [str(current["name"])] + [str(item) for item in current.get("synonyms", [])]
        for label in names:
            key = normalize(label)
            if key:
                candidates[key].append((hpo_id, str(current["name"])))

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("data-version:"):
                release = line.split(":", 1)[1].strip()
            if line == "[Term]":
                emit()
                current = {"synonyms": []}
                continue
            if not current or not line or line.startswith("!"):
                continue
            if line.startswith("id: HP:"):
                current["id"] = line.split(":", 1)[1].strip()
            elif line.startswith("name:"):
                current["name"] = line.split(":", 1)[1].strip()
            elif line.startswith("synonym:"):
                match = re.search(r'"([^"\\]*(?:\\.[^"\\]*)*)"', line)
                if match:
                    current["synonyms"].append(match.group(1).replace('\\"', '"'))
            elif line.startswith("is_obsolete: true"):
                current["obsolete"] = True
        emit()

    return dict(candidates), release


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "variant_id", "gene", "source_phenotype_term", "normalized_term",
        "hpo_id", "hpo_label", "mapping_status", "mapping_source",
        "hpo_release", "safe_for_training", "note",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ontology", type=Path, default=DEFAULT_ONTOLOGY)
    args = parser.parse_args()
    if not args.ontology.exists():
        raise SystemExit(f"HPO ontology not found: {args.ontology}")

    source = PROCESSED / "deepgene_phenotype_normalized_v1.csv"
    terms = read_rows(source)
    ontology, release = parse_obo(args.ontology)
    output = []
    status_counts = Counter()
    variant_ids = set()
    mapped_variants = set()
    seen_terms = set()

    for row in terms:
        variant_id = clean(row.get("variant_id"))
        raw_terms = [x for x in clean(row.get("meaningful_phenotype_terms")).split("|") if clean(x)]
        if not raw_terms:
            raw_terms = [x for x in clean(row.get("phenotypes")).split("|") if clean(x)]
        for term in raw_terms:
            normalized = normalize(term)
            matches = ontology.get(normalized, [])
            unique_matches = sorted(set(matches))
            if len(unique_matches) == 1:
                hpo_id, hpo_label = unique_matches[0]
                status = "official_exact_match"
                source_name = "Official HPO name/synonym exact match"
                note = "Annotation only; not a clinical or functional training label."
                mapped_variants.add(variant_id)
            elif len(unique_matches) > 1:
                hpo_id = hpo_label = ""
                status = "ambiguous_official_match"
                source_name = "Official HPO exact match is ambiguous"
                note = "No ID assigned because the normalized term maps to multiple HPO terms."
            else:
                hpo_id = hpo_label = ""
                status = "unmapped_pending_curated_review"
                source_name = "none"
                note = "No exact HPO name or synonym match; no approximate mapping was attempted."
            output.append({
                "variant_id": variant_id, "gene": "SCN1A",
                "source_phenotype_term": term, "normalized_term": normalized,
                "hpo_id": hpo_id, "hpo_label": hpo_label,
                "mapping_status": status, "mapping_source": source_name,
                "hpo_release": release, "safe_for_training": 0, "note": note,
            })
            status_counts[status] += 1
            variant_ids.add(variant_id)
            seen_terms.add(normalized)

    output_path = PROCESSED / "deepgene_v2_hpo_annotations.csv"
    write_rows(output_path, output)
    audit = {
        "created": date.today().isoformat(),
        "status": "OFFICIAL_HPO_EXACT_MATCH_ANNOTATION_ONLY",
        "ontology_file": str(args.ontology.relative_to(ROOT)) if args.ontology.is_relative_to(ROOT) else str(args.ontology),
        "hpo_release": release,
        "input_file": str(source.relative_to(ROOT)),
        "input_rows": len(terms),
        "annotation_rows": len(output),
        "unique_terms_seen": len(seen_terms),
        "unique_variants_seen": len(variant_ids),
        "unique_variants_with_official_exact_match": len(mapped_variants),
        "mapping_status_counts": dict(status_counts),
        "rules": [
            "Only exact normalized HPO names and synonyms are accepted.",
            "Approximate, disease-name, and ambiguous mappings remain unresolved.",
            "HPO annotations are not LOF/GOF labels, pathogenicity labels, or diagnoses.",
            "safe_for_training remains 0 until a separate feature policy and validation study are approved.",
        ],
    }
    (RESULTS / "v2_hpo_join_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
