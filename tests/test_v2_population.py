import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class V2PopulationJoinTests(unittest.TestCase):
    def test_population_join_uses_audited_fallbacks(self):
        audit = json.loads(
            (ROOT / "data/analysis_results/v2_population_join_audit.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(audit["status"], "PASS_WITH_COVERAGE_AUDIT")
        self.assertEqual(audit["exact_coordinate_matches"], 1407)
        self.assertEqual(audit["rsid_allele_matches"], 1127)
        self.assertEqual(audit["ambiguous_rsid_rows"], 0)
        self.assertFalse(audit["used_as_training_label"])

    def test_population_output_records_match_method(self):
        with (ROOT / "data/processed/deepgene_v2_population_annotations.csv").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 10696)
        self.assertEqual(
            {row["match_method"] for row in rows},
            {"exact_coordinate", "rsid_allele", "none"},
        )
        self.assertTrue(
            all(
                row["match_method"] != "rsid_allele"
                or row["join_status"] == "matched"
                for row in rows
            )
        )


if __name__ == "__main__":
    unittest.main()
