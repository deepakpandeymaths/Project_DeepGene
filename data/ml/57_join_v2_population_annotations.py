"""Join an optional versioned gnomAD export to SCN1A ClinVar representations.

The script never derives population frequencies from clinical classifications.
When the expected gnomAD source is absent, it writes a blocked audit and an
explicit source-missing table instead of fabricated values.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"
METADATA = ROOT / "data" / "metadata"
GNOMAD_SOURCE = RAW / "gnomad_scn1a_population.tsv"
CLINVAR_SOURCE = RAW / "SCN1A_clinvar.tsv"


def read_delimited(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def first(row: dict[str, str], names: list[str]) -> str:
    for name in names:
        if clean(row.get(name)) not in {"", "na", "nan", ".", "-"}:
            return clean(row.get(name))
    return ""


def key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (first(row, ["assembly", "Assembly"]), first(row, ["chromosome", "Chromosome"]), first(row, ["position", "pos", "Start", "PositionVCF"]), first(row, ["reference_allele", "ref", "ReferenceAlleleVCF", "ReferenceAllele"]), first(row, ["alternate_allele", "alt", "AlternateAlleleVCF", "AlternateAllele"]))


def main() -> None:
    output_fields = [
        "variant_id", "assembly", "chromosome", "position", "reference_allele", "alternate_allele",
        "allele_frequency", "max_population_af", "allele_count", "allele_number", "homozygote_count",
        "join_status", "source_id", "source_version", "missing_value_note",
    ]
    clinvar = read_delimited(CLINVAR_SOURCE)
    representations = {}
    for row in clinvar:
        variant_id = clean(row.get("#AlleleID"))
        representation = {
            "variant_id": variant_id,
            "assembly": clean(row.get("Assembly")),
            "chromosome": clean(row.get("Chromosome")),
            "position": clean(row.get("PositionVCF")) or clean(row.get("Start")),
            "reference_allele": clean(row.get("ReferenceAlleleVCF")) or clean(row.get("ReferenceAllele")),
            "alternate_allele": clean(row.get("AlternateAlleleVCF")) or clean(row.get("AlternateAllele")),
        }
        representations.setdefault((variant_id, key(representation)), representation)

    output = []
    audit = {
        "created": date.today().isoformat(),
        "source_id": "SRC-GNOMAD",
        "expected_source": str(GNOMAD_SOURCE.relative_to(ROOT)),
        "clinical_representation_rows": len(representations),
        "join_key": "assembly + chromosome + position + reference_allele + alternate_allele",
        "frequency_is_annotation": True,
        "used_as_training_label": False,
    }
    if not GNOMAD_SOURCE.exists():
        for representation in representations.values():
            output.append({**representation, "join_status": "source_missing", "source_id": "SRC-GNOMAD", "source_version": "not_available", "missing_value_note": "No local gnomAD export is present; missing does not mean rare, benign, or absent."})
        audit.update({
            "status": "STATUS: BLOCKED",
            "reason": "The expected gnomAD source file is not present in data/raw.",
            "matched_rows": 0,
            "unmatched_rows": len(output),
            "next_action": "Provide a versioned gnomAD TSV with the documented genomic join columns, then rerun this script.",
        })
    else:
        gnomad = read_delimited(GNOMAD_SOURCE)
        source_version = clean(gnomad[0].get("source_version")) if gnomad else "version_not_declared"
        by_key = {}
        duplicate_keys = 0
        for row in gnomad:
            k = key(row)
            if not all(k):
                continue
            if k in by_key:
                duplicate_keys += 1
            by_key[k] = row
        matched = 0
        for representation in representations.values():
            match = by_key.get(key(representation))
            if match:
                matched += 1
                output.append({**representation, "allele_frequency": first(match, ["allele_frequency", "AF", "AF_nfe"]), "max_population_af": first(match, ["max_population_af", "popmax", "AF_popmax"]), "allele_count": first(match, ["allele_count", "AC"]), "allele_number": first(match, ["allele_number", "AN"]), "homozygote_count": first(match, ["homozygote_count", "nhomalt"]), "join_status": "matched", "source_id": "SRC-GNOMAD", "source_version": source_version or "version_not_declared", "missing_value_note": "Blank frequency fields remain unavailable and are not interpreted as zero."})
            else:
                output.append({**representation, "join_status": "unmatched", "source_id": "SRC-GNOMAD", "source_version": source_version or "version_not_declared", "missing_value_note": "No gnomAD row matched this genomic representation."})
        audit.update({
            "status": "PASS_WITH_COVERAGE_AUDIT",
            "source_rows": len(gnomad),
            "source_duplicate_join_keys": duplicate_keys,
            "matched_rows": matched,
            "unmatched_rows": len(output) - matched,
            "coverage_fraction": matched / len(output) if output else 0.0,
            "next_action": "Review coverage, assembly, and allele semantics before updating V2 evidence completeness.",
        })

    write_csv(PROCESSED / "deepgene_v2_population_annotations.csv", output, output_fields)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "v2_population_join_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    manifest = {
        "source_id": "SRC-GNOMAD",
        "path": str(GNOMAD_SOURCE.relative_to(ROOT)),
        "status": audit["status"],
        "join_key": audit["join_key"],
        "required_columns": ["assembly", "chromosome", "position", "reference_allele", "alternate_allele", "allele_frequency", "source_version"],
        "sha256": hashlib.sha256(GNOMAD_SOURCE.read_bytes()).hexdigest() if GNOMAD_SOURCE.exists() else None,
        "note": "Do not use this annotation as a clinical or functional label.",
    }
    METADATA.mkdir(parents=True, exist_ok=True)
    (METADATA / "v2_population_source_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({
        "status": audit["status"],
        "representation_rows": len(output),
        "matched_rows": audit.get("matched_rows", 0),
        "unmatched_rows": audit.get("unmatched_rows", len(output)),
        "outputs": [
            "data/processed/deepgene_v2_population_annotations.csv",
            "data/analysis_results/v2_population_join_audit.json",
            "data/metadata/v2_population_source_manifest.json",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
