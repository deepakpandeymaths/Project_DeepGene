"""Build a V2 evidence-completeness profile without inventing missing evidence."""

from __future__ import annotations

import csv
import json
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "data" / "analysis_results"


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


def yes(value: object) -> int:
    return int(clean(value).lower() in {"1", "true", "yes", "y"})


def main() -> None:
    master = read_csv(PROCESSED / "deepgene_v2_master_variants.csv")
    alpha = read_csv(PROCESSED / "alphamissense_gene_hg38.csv")
    population = read_csv(PROCESSED / "deepgene_v2_population_annotations.csv")
    population_by_variant: dict[str, list[dict[str, str]]] = {}
    for item in population:
        population_by_variant.setdefault(clean(item.get("variant_id")), []).append(item)
        variation_id = clean(item.get("variation_id"))
        if variation_id:
            population_by_variant.setdefault(variation_id, []).append(item)
    evidence_fields = [
        "clinical", "phenotype", "functional", "alphamissense_variant_join",
        "hpo_normalized", "population_frequency", "conservation", "splice",
        "structural_domain",
    ]
    fields = [
        "variant_id", "gene", "protein_key", "clinical_significance", "review_status",
        "clinical_evidence_present", "phenotype_evidence_present", "meaningful_phenotype_present",
        "functional_evidence_present", "alphamissense_variant_join_present", "hpo_normalized_present",
        "population_frequency_present", "conservation_present", "splice_evidence_present",
        "structural_domain_present", "evidence_present_count", "evidence_category_count",
        "evidence_completeness_fraction", "evidence_completeness_category", "missing_evidence_categories",
        "population_frequency", "population_max_af", "population_evidence_status",
        "population_match_methods", "population_source_version", "missing_means_unavailable", "training_note",
    ]
    rows: list[dict[str, object]] = []
    category_counts = Counter()
    for source in master:
        clinical_present = int(bool(clean(source.get("clinical_significance")) or clean(source.get("review_status"))))
        phenotype_present = int(bool(clean(source.get("phenotypes"))))
        meaningful_present = int(clean(source.get("meaningful_phenotype_term_count")) not in {"", "0"})
        functional_present = int(clean(source.get("functional_source_count")) not in {"", "0"})
        population_rows = population_by_variant.get(clean(source.get("variant_id")), [])
        population_rows += population_by_variant.get(clean(source.get("variation_id")), [])
        population_rows = list({id(item): item for item in population_rows}.values())
        matched_population = [item for item in population_rows if clean(item.get("join_status")) == "matched"]
        ambiguous_population = any(clean(item.get("join_status")) == "ambiguous" for item in population_rows)
        population_frequency_values = []
        population_max_values = []
        for item in matched_population:
            for field, target in (("allele_frequency", population_frequency_values), ("max_population_af", population_max_values)):
                try:
                    value = float(clean(item.get(field)))
                except (TypeError, ValueError):
                    continue
                if 0 <= value <= 1:
                    target.append(value)
        population_present = int(bool(population_frequency_values or population_max_values))
        if ambiguous_population:
            population_status = "ambiguous"
        elif population_present:
            population_status = "available"
        else:
            population_status = "missing"
        population_methods = sorted({clean(item.get("match_method")) for item in matched_population if clean(item.get("match_method"))})
        population_versions = sorted({clean(item.get("source_version")) for item in population_rows if clean(item.get("source_version"))})
        values = {
            "clinical": clinical_present,
            "phenotype": phenotype_present,
            "functional": functional_present,
            "alphamissense_variant_join": 0,
            "hpo_normalized": 0,
            "population_frequency": population_present,
            "conservation": 0,
            "splice": 0,
            "structural_domain": 0,
        }
        present_count = sum(values.values())
        fraction = present_count / len(evidence_fields)
        if present_count <= 2:
            category = "limited"
        elif present_count <= 4:
            category = "partial"
        else:
            category = "broad"
        category_counts[category] += 1
        missing = [name for name, present in values.items() if not present]
        rows.append({
            "variant_id": clean(source.get("variant_id")),
            "gene": clean(source.get("gene")) or "SCN1A",
            "protein_key": clean(source.get("protein_key")),
            "clinical_significance": clean(source.get("clinical_significance")),
            "review_status": clean(source.get("review_status")),
            "clinical_evidence_present": clinical_present,
            "phenotype_evidence_present": phenotype_present,
            "meaningful_phenotype_present": meaningful_present,
            "functional_evidence_present": functional_present,
            "alphamissense_variant_join_present": 0,
            "hpo_normalized_present": 0,
            "population_frequency_present": population_present,
            "conservation_present": 0,
            "splice_evidence_present": 0,
            "structural_domain_present": 0,
            "evidence_present_count": present_count,
            "evidence_category_count": len(evidence_fields),
            "evidence_completeness_fraction": f"{fraction:.6f}",
            "evidence_completeness_category": category,
            "missing_evidence_categories": "|".join(missing),
            "population_frequency": f"{max(population_frequency_values):.12g}" if population_frequency_values else "",
            "population_max_af": f"{max(population_max_values):.12g}" if population_max_values else "",
            "population_evidence_status": population_status,
            "population_match_methods": "|".join(population_methods),
            "population_source_version": "|".join(population_versions),
            "missing_means_unavailable": 1,
            "training_note": "Population frequency is annotation evidence only; missing evidence is unavailable, not negative or benign; use explicit missingness indicators.",
        })

    write_csv(PROCESSED / "deepgene_v2_evidence_completeness.csv", rows, fields)

    plan = {
        "created": date.today().isoformat(),
        "status": "BASELINE_COMPLETE_EXTERNAL_ENRICHMENT_PENDING",
        "variant_rows_profiled": len(rows),
        "evidence_categories": {
            "clinical": {"status": "available", "source": "ClinVar V1 evidence layer", "safe_for_training": "conditional"},
            "phenotype": {"status": "available_as_raw_text", "source": "ClinVar phenotype descriptions", "safe_for_training": "conditional"},
            "functional": {"status": "available_for_source_specific_cohorts", "source": "SCION and Brunklaus supplementary data", "safe_for_training": "yes_for_experiment_target_only"},
            "alphamissense_variant_join": {"status": "pending_mapping", "source": "AlphaMissense gene/transcript table", "safe_for_training": "annotation_only_until_joined_and_validated"},
            "hpo_normalized": {"status": "pending", "source": "Human Phenotype Ontology", "safe_for_training": "conditional"},
            "population_frequency": {"status": "available_with_coverage_audit", "source": "gnomAD v4.1.2", "safe_for_training": "annotation_only_until_target_is_defined"},
            "conservation": {"status": "pending", "source": "versioned conservation resource", "safe_for_training": "conditional"},
            "splice": {"status": "pending", "source": "versioned splice-impact resource", "safe_for_training": "conditional"},
            "structural_domain": {"status": "pending", "source": "SCN1A domain and structure resources", "safe_for_training": "conditional"},
        },
        "current_category_counts": dict(category_counts),
        "rules": [
            "Evidence completeness is descriptive and is not a pathogenicity score.",
            "Missing evidence is not evidence of benignity, safety, or absence of disease.",
            "Computational annotations remain annotations until a separate validated target is defined.",
            "New sources must include release/version, accession or URL, join key, coverage report, and missing-value semantics.",
        ],
        "next_gate": "Add validated HPO mappings and preserve population-frequency missingness before any model-feature experiment.",
        "alphamissense_rows_preserved": len(alpha),
        "alphamissense_join_note": "The current AlphaMissense table contains transcript-level rows but no validated mapping into the V2 variant master; all variant-level join indicators remain zero.",
        "population_join_note": "gnomAD frequency is joined as annotation evidence using exact coordinates and allele-aware rsID fallbacks; it is never converted into a training label.",
    }
    (RESULTS / "v2_evidence_completeness_report.json").write_text(json.dumps(plan, indent=2), encoding="utf-8")
    (RESULTS / "v2_evidence_enrichment_plan.json").write_text(json.dumps({
        "created": date.today().isoformat(),
        "priority_order": [
            "HPO phenotype normalization with versioned mappings",
            "gnomAD population frequency with assembly and allele matching",
            "AlphaMissense variant-level mapping with transcript validation",
            "conservation annotations",
            "splice-impact annotations",
            "SCN1A domain and structural annotations",
        ],
        "acceptance_gates": [
            "No source is merged without a documented join key.",
            "Coverage and unmatched-row counts are reported.",
            "Source release and checksum are recorded.",
            "Missing values remain explicit and are not converted to negative evidence.",
            "The enriched source is not used as a label without independent validation.",
        ],
    }, indent=2), encoding="utf-8")

    print(json.dumps({
        "profile_rows": len(rows),
        "category_counts": dict(category_counts),
        "alphamissense_rows_preserved_as_annotation": len(alpha),
        "outputs": [
            "data/processed/deepgene_v2_evidence_completeness.csv",
            "data/analysis_results/v2_evidence_completeness_report.json",
            "data/analysis_results/v2_evidence_enrichment_plan.json",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
