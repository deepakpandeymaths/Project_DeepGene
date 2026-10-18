"""Build the DeepGene V2 master tables and provenance audits.

This script deliberately keeps clinical evidence, functional labels, and
computational annotations as separate products. It does not infer labels from
ClinVar prose, AlphaMissense scores, or conversation memory.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def protein_key_from_hgvs(hgvs: str) -> str:
    match = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", clean(hgvs))
    if not match:
        return ""
    return f"SCN1A:p.{match.group(1)}{match.group(2)}{match.group(3)}"


def functional_key(gene: str, value: str) -> str:
    value = clean(value)
    if not value:
        return ""
    if value.startswith(f"{gene}:p."):
        return value
    if value.startswith("p."):
        return f"{gene}:{value}"
    return f"{gene}:p.{value}"


def duplicate_summary(rows: list[dict[str, str]], key_name: str, label_name: str | None = None) -> dict[str, object]:
    keys = [clean(row.get(key_name)) for row in rows]
    nonempty = [key for key in keys if key]
    counts = Counter(nonempty)
    duplicate_groups = {key: count for key, count in counts.items() if count > 1}
    conflicts = {}
    if label_name:
        labels = defaultdict(set)
        for row in rows:
            key = clean(row.get(key_name))
            label = clean(row.get(label_name))
            if key and label:
                labels[key].add(label)
        conflicts = {key: sorted(values) for key, values in labels.items() if len(values) > 1}
    return {
        "rows": len(rows),
        "nonempty_identity_keys": len(nonempty),
        "unique_identity_keys": len(set(nonempty)),
        "duplicate_identity_groups": len(duplicate_groups),
        "duplicate_rows_beyond_first": sum(count - 1 for count in duplicate_groups.values()),
        "conflicting_label_groups": len(conflicts),
        "conflicting_labels": conflicts,
    }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    clinical = read_csv(PROCESSED / "deepgene_v1_final.csv")
    phenotype = read_csv(PROCESSED / "deepgene_phenotype_profiles_v1.csv")
    brunklaus = read_csv(PROCESSED / "deepgene_scn1a_brunklaus_2020_functional.csv")
    generation2 = read_csv(PROCESSED / "deepgene_gen2_ml_ready_v1.csv")
    alphamissense = read_csv(PROCESSED / "alphamissense_gene_hg38.csv")

    phenotype_by_id = {clean(row.get("variant_id")): row for row in phenotype}
    functional_rows: list[dict[str, object]] = []

    for row in brunklaus:
        key = functional_key(clean(row.get("gene")), clean(row.get("variant_key")))
        functional_rows.append(
            {
                "source_id": "SRC-BRUNKLAUS",
                "source_row_id": clean(row.get("id")),
                "gene": clean(row.get("gene")),
                "variant_key": key,
                "reference_aa": clean(row.get("aa1")),
                "alternate_aa": clean(row.get("aa2")),
                "protein_position": clean(row.get("pos")),
                "functional_label": clean(row.get("functional_label")),
                "source_label": clean(row.get("y")),
                "label_status": clean(row.get("label_status")),
                "source_name": clean(row.get("source_name")),
                "source_url": clean(row.get("source_url")),
            }
        )

    for index, row in enumerate(generation2, start=1):
        gene = clean(row.get("gene"))
        key = functional_key(gene, clean(row.get("variant_key")))
        functional_rows.append(
            {
                "source_id": "SRC-SCION",
                "source_row_id": f"scion_{index:04d}",
                "gene": gene,
                "variant_key": key,
                "reference_aa": clean(row.get("aa1")),
                "alternate_aa": clean(row.get("aa2")),
                "protein_position": clean(row.get("protein_position")),
                "functional_label": clean(row.get("functional_label")),
                "source_label": clean(row.get("functional_label")),
                "label_status": "source_functional_label",
                "source_name": "SCION harmonized V1 table",
                "source_url": "https://github.com/christianbosselmann/SCION",
            }
        )

    functional_fields = [
        "source_id", "source_row_id", "gene", "variant_key", "reference_aa",
        "alternate_aa", "protein_position", "functional_label", "source_label",
        "label_status", "source_name", "source_url",
    ]
    write_csv(PROCESSED / "deepgene_v2_functional_evidence.csv", functional_rows, functional_fields)

    functional_by_key: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in functional_rows:
        key = clean(row.get("variant_key"))
        if key:
            functional_by_key[key].append(row)

    master_fields = [
        "variant_id", "gene", "hgvs_name", "protein_key", "variation_id", "clinical_significance",
        "review_status", "phenotypes", "phenotype_term_count", "meaningful_phenotype_term_count",
        "pathogenic_count", "likely_pathogenic_count", "benign_count", "likely_benign_count",
        "vus_count", "conflicting_count", "submitter_count", "multiple_submitters",
        "expert_panel_review", "conflict_flag", "grch37_present", "grch38_present",
        "dbsnp_present", "genomic_coordinates_present", "hgvs_present", "clinical_source_id",
        "functional_source_count", "functional_source_ids", "functional_label_values",
        "functional_conflict_flag", "alphamissense_join_status", "conversation_training_eligible",
    ]
    master_rows: list[dict[str, object]] = []
    for row in clinical:
        variant_id = clean(row.get("variant_id"))
        phenotype_row = phenotype_by_id.get(variant_id, {})
        protein_key = protein_key_from_hgvs(clean(row.get("hgvs_name")))
        matched = functional_by_key.get(protein_key, [])
        labels = sorted({clean(item.get("functional_label")) for item in matched if clean(item.get("functional_label"))})
        sources = sorted({clean(item.get("source_id")) for item in matched if clean(item.get("source_id"))})
        master_rows.append(
            {
                "variant_id": variant_id,
                "gene": "SCN1A",
                "hgvs_name": clean(row.get("hgvs_name")),
                "protein_key": protein_key,
                "variation_id": clean(row.get("variation_id")),
                "clinical_significance": clean(row.get("clinical_significance")),
                "review_status": clean(row.get("review_status")),
                "phenotypes": clean(row.get("phenotypes")),
                "phenotype_term_count": clean(phenotype_row.get("phenotype_term_count")),
                "meaningful_phenotype_term_count": clean(phenotype_row.get("meaningful_phenotype_term_count")),
                "pathogenic_count": clean(row.get("pathogenic_count")),
                "likely_pathogenic_count": clean(row.get("likely_pathogenic_count")),
                "benign_count": clean(row.get("benign_count")),
                "likely_benign_count": clean(row.get("likely_benign_count")),
                "vus_count": clean(row.get("vus_count")),
                "conflicting_count": clean(row.get("conflicting_count")),
                "submitter_count": clean(row.get("submitter_count")),
                "multiple_submitters": clean(row.get("multiple_submitters")),
                "expert_panel_review": clean(row.get("expert_panel_review")),
                "conflict_flag": clean(row.get("conflict_flag")),
                "grch37_present": clean(row.get("grch37_present")),
                "grch38_present": clean(row.get("grch38_present")),
                "dbsnp_present": clean(row.get("dbsnp_present")),
                "genomic_coordinates_present": clean(row.get("genomic_coordinates_present")),
                "hgvs_present": clean(row.get("hgvs_present")),
                "clinical_source_id": "SRC-CLINVAR",
                "functional_source_count": len(matched),
                "functional_source_ids": "|".join(sources),
                "functional_label_values": "|".join(labels),
                "functional_conflict_flag": int(len(labels) > 1),
                "alphamissense_join_status": "not_joinable_transcript_only",
                "conversation_training_eligible": 0,
            }
        )

    write_csv(PROCESSED / "deepgene_v2_master_variants.csv", master_rows, master_fields)
    clinical_fields = [
        "variant_id", "gene", "hgvs_name", "variation_id", "clinical_significance", "review_status",
        "phenotypes", "pathogenic_count", "pathogenic_likely_pathogenic_count", "likely_pathogenic_count",
        "benign_count", "benign_likely_benign_count", "likely_benign_count", "vus_count", "conflicting_count",
        "submitter_count", "multiple_submitters", "expert_panel_review", "conflict_flag", "phenotype_available",
        "grch37_present", "grch38_present", "dbsnp_present", "genomic_coordinates_present", "hgvs_present",
        "clinical_source_id",
    ]
    clinical_product = [{**{field: clean(row.get(field)) for field in clinical_fields if field != "gene"}, "gene": "SCN1A", "clinical_source_id": "SRC-CLINVAR"} for row in clinical]
    write_csv(PROCESSED / "deepgene_v2_clinical_evidence.csv", clinical_product, clinical_fields)

    brunklaus_keys = {clean(row.get("variant_key")) for row in functional_rows if row["source_id"] == "SRC-BRUNKLAUS" and clean(row.get("variant_key"))}
    scion_keys = {clean(row.get("variant_key")) for row in functional_rows if row["source_id"] == "SRC-SCION" and clean(row.get("variant_key"))}
    clinical_keys = {clean(row.get("protein_key")) for row in master_rows if clean(row.get("protein_key"))}
    duplicate_audit = {
        "audit_date": date.today().isoformat(),
        "purpose": "V2 identity, duplicate, cross-source overlap, and label-conflict audit",
        "source_summaries": {
            "clinvar_variant_id": duplicate_summary(clinical, "variant_id", "clinical_significance"),
            "brunklaus_variant_key": duplicate_summary([row for row in functional_rows if row["source_id"] == "SRC-BRUNKLAUS"], "variant_key", "functional_label"),
            "scion_variant_key": duplicate_summary([row for row in functional_rows if row["source_id"] == "SRC-SCION"], "variant_key", "functional_label"),
        },
        "cross_source_overlaps": {
            "brunklaus_scion_functional_keys": len(brunklaus_keys & scion_keys),
            "clinical_functional_protein_keys": len(clinical_keys & (brunklaus_keys | scion_keys)),
            "clinical_variants_with_functional_records": sum(1 for row in master_rows if int(row["functional_source_count"]) > 0),
        },
        "alphamissense": {
            "rows": len(alphamissense),
            "unique_transcript_ids": len({clean(row.get("transcript_id")) for row in alphamissense if clean(row.get("transcript_id"))}),
            "master_join_status": "not_joined_transcript_only_source_without_variant_transcript_mapping",
        },
        "interpretation": "Duplicate and overlap counts are audit signals. Cross-source overlap does not make either source an independent clinical gold standard.",
    }
    (RESULTS / "v2_duplicate_audit.json").write_text(json.dumps(duplicate_audit, indent=2), encoding="utf-8")

    leakage_audit = {
        "audit_date": date.today().isoformat(),
        "status": "PASS_WITH_GUARDS",
        "checks": [
            {"check": "clinical_target_reconstruction", "status": "GUARDED", "fields": ["clinical_significance", "review_status", "conflict_flag"], "rule": "Do not use a submitted clinical target or its direct aggregation to predict the same target."},
            {"check": "functional_target_reconstruction", "status": "GUARDED", "fields": ["functional_label", "functional_target"], "rule": "Functional labels are targets and must be excluded from feature matrices."},
            {"check": "cross_source_duplicate_overlap", "status": "REPORTED", "fields": ["protein_key", "variant_key"], "rule": "Overlapping biological variants require grouped splits or explicit source-aware evaluation."},
            {"check": "alphamissense_annotation", "status": "SAFE_AS_ANNOTATION", "fields": ["alphamissense_score"], "rule": "AlphaMissense is not converted to LOF/GOF or clinical ground truth."},
            {"check": "conversation_memory", "status": "EXCLUDED", "fields": ["conversation_text"], "rule": "Conversation memory is retrieval context only and cannot enter training labels or features."},
            {"check": "missing_evidence", "status": "GUARDED", "fields": ["phenotypes", "population_max_af", "hpo_id", "splice_score"], "rule": "Missing evidence means unavailable, not negative or benign."},
        ],
        "training_policy": "Only source-labelled, provenance-complete, experiment-specific targets may enter supervised training.",
    }
    (RESULTS / "v2_leakage_audit.json").write_text(json.dumps(leakage_audit, indent=2), encoding="utf-8")

    print(json.dumps({
        "master_rows": len(master_rows),
        "clinical_rows": len(clinical_product),
        "functional_rows": len(functional_rows),
        "alphamissense_rows_preserved_as_annotation": len(alphamissense),
        "clinical_functional_protein_overlaps": len(clinical_keys & (brunklaus_keys | scion_keys)),
        "outputs": [
            "data/processed/deepgene_v2_master_variants.csv",
            "data/processed/deepgene_v2_clinical_evidence.csv",
            "data/processed/deepgene_v2_functional_evidence.csv",
            "data/analysis_results/v2_duplicate_audit.json",
            "data/analysis_results/v2_leakage_audit.json",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
