import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class V2PopulationScopeTests(unittest.TestCase):
    def test_current_export_scope_is_explicit(self):
        report = json.loads(
            (ROOT / "data/analysis_results/v2_population_unmatched_audit.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(report["status"], "API_ENRICHED_PARTIAL_SOURCE_SCOPE_LIMITED")
        self.assertEqual(report["source_chromosomes"], ["2"])
        self.assertEqual(report["safe_additional_matches_from_current_export"], 0)
        self.assertEqual(report["shared_rsid_allele_mismatches"], 600)
        self.assertEqual(report["variant_level_matched_rows"], 1493)
        self.assertAlmostEqual(report["variant_level_coverage"], 0.2774577216, places=6)


if __name__ == "__main__":
    unittest.main()
