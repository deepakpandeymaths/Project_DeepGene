import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class V2HPOContractTests(unittest.TestCase):
    def test_hpo_audit_is_official_annotation_only(self):
        audit = json.loads((ROOT / "data/analysis_results/v2_hpo_join_audit.json").read_text(encoding="utf-8"))
        self.assertEqual(audit["status"], "OFFICIAL_HPO_EXACT_MATCH_ANNOTATION_ONLY")
        self.assertGreater(audit["unique_variants_with_official_exact_match"], 0)
        self.assertTrue(any("safe_for_training" in rule for rule in audit["rules"]))

    def test_hpo_output_preserves_unmatched_terms_and_disables_training(self):
        path = ROOT / "data/processed/deepgene_v2_hpo_annotations.csv"
        with path.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertTrue(rows)
        self.assertTrue(any(row["mapping_status"] == "official_exact_match" for row in rows))
        self.assertTrue(any(row["mapping_status"] == "unmapped_pending_curated_review" for row in rows))
        self.assertEqual({row["safe_for_training"] for row in rows}, {"0"})


if __name__ == "__main__":
    unittest.main()
