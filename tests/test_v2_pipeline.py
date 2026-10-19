import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class V2PipelineContractTests(unittest.TestCase):
    def read_csv(self, relative):
        with (ROOT / relative).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))

    def test_v2_master_and_completeness_have_same_variant_count(self):
        master = self.read_csv("data/processed/deepgene_v2_master_variants.csv")
        completeness = self.read_csv("data/processed/deepgene_v2_evidence_completeness.csv")
        self.assertEqual(len(master), 5381)
        self.assertEqual(len(master), len(completeness))
        self.assertEqual({row["variant_id"] for row in master}, {row["variant_id"] for row in completeness})

    def test_functional_targets_are_explicit_and_separate(self):
        rows = self.read_csv("data/processed/deepgene_v2_functional_evidence.csv")
        self.assertEqual(len(rows), 375)
        self.assertTrue({row["functional_label"] for row in rows} <= {"LOF", "GOF", "Loss", "Gain"})
        self.assertTrue(all(row["source_id"] in {"SRC-BRUNKLAUS", "SRC-SCION"} for row in rows))

    def test_grouped_splits_do_not_cross_boundaries(self):
        rows = self.read_csv("data/ml_results/v2_grouped_split_ids.csv")
        groups = {}
        for row in rows:
            groups.setdefault(row["group_key"], set()).add(row["split"])
        self.assertTrue(all(len(splits) == 1 for splits in groups.values()))

    def test_leakage_audit_excludes_conversation_and_annotations_as_targets(self):
        report = json.loads((ROOT / "data/analysis_results/v2_leakage_audit.json").read_text(encoding="utf-8"))
        checks = {item["check"]: item for item in report["checks"]}
        self.assertEqual(checks["conversation_memory"]["status"], "EXCLUDED")
        self.assertEqual(checks["alphamissense_annotation"]["status"], "SAFE_AS_ANNOTATION")

    def test_priority_join_normalizes_scientific_protein_hgvs(self):
        rows = self.read_csv("data/processed/deepgene_v2_research_prioritization.csv")
        scn1a = [row for row in rows if row["variant_key"].startswith("SCN1A:")]
        self.assertEqual(len(scn1a), 14)
        self.assertTrue(all(row["clinical_variant_ids"] for row in scn1a))
        report = json.loads((ROOT / "data/analysis_results/v2_prioritization_report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["clinical_joined_rows"], 14)


if __name__ == "__main__":
    unittest.main()
