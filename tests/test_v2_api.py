import unittest

import app.server as server


class V2ApiContractTests(unittest.TestCase):
    def test_summary_is_research_only_and_has_ranked_rows(self):
        payload = server.v2_summary()
        self.assertEqual(payload["status"], "research_only")
        self.assertGreater(payload["prioritized_rows"], 0)

    def test_prioritization_limit_and_query(self):
        rows = server.v2_prioritization(2)
        self.assertLessEqual(len(rows), 2)
        self.assertTrue(rows)
        key = rows[0]["variant_key"]
        filtered = server.v2_prioritization(20, key)
        self.assertTrue(any(item["variant_key"] == key for item in filtered))

    def test_coverage_contract_is_research_only(self):
        payload = server.v2_coverage()
        self.assertEqual(payload["status"], "research_only")
        self.assertGreater(payload["population"]["representation_total"], 0)
        self.assertGreater(payload["hpo"]["term_rows"], 0)
        self.assertTrue(payload["hpo"]["annotation_only"])

    def test_calibration_contract_is_research_only(self):
        payload = server.v2_calibration()
        self.assertTrue(payload["research_only"])
        self.assertGreater(payload["test_rows"], 0)
        self.assertTrue(payload["models"])
        self.assertIn("expected_calibration_error", payload["models"][0])


if __name__ == "__main__":
    unittest.main()
