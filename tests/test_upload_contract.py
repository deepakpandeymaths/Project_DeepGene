import gzip
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

import app.server as server


class UploadContractTests(unittest.TestCase):
    def test_parser_rejects_invalid_gzip(self):
        with self.assertRaisesRegex(ValueError, "not a valid gzip"):
            server.parse_table("dataset.tsv.gz", b"not gzip")

    def test_parser_rejects_unsupported_extension(self):
        with self.assertRaisesRegex(ValueError, "Only CSV"):
            server.parse_table("dataset.xlsx", b"a,b\n1,2\n")

    def test_parser_reads_gzip_tsv_and_skips_comments(self):
        payload = gzip.compress(b"# metadata\ngene\taa1\taa2\tpos\ty\nSCN1A\tR\tQ\t1\tLOF\n")
        frame = server.parse_table("dataset.tsv.gz", payload)
        self.assertEqual(list(frame.columns), ["gene", "aa1", "aa2", "pos", "y"])
        self.assertEqual(len(frame), 1)

    def test_functional_upload_reports_missing_columns(self):
        frame = pd.DataFrame([{"gene": "SCN1A", "aa1": "R", "aa2": "Q", "pos": 1}])
        with self.assertRaisesRegex(ValueError, "Functional file missing columns"):
            server.retrain_functional(frame)

    def test_functional_upload_rejects_unlabelled_rows(self):
        frame = pd.DataFrame([{"gene": "SCN1A", "aa1": "R", "aa2": "Q", "pos": 1, "y": "unknown"}])
        with self.assertRaisesRegex(ValueError, "No valid LOF/GOF"):
            server.retrain_functional(frame)

    def test_alphamissense_scores_must_be_in_range(self):
        frame = pd.DataFrame([{"transcript_id": "TX1", "mean_am_pathogenicity": 1.2}])
        with self.assertRaisesRegex(ValueError, "between 0 and 1"):
            server.ingest_alphamissense_gene(frame)

    def test_local_ai_fallback_without_api_key(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            memory = Path(temp_dir) / "conversations.jsonl"
            with patch.object(server, "CONVERSATIONS", memory), patch.dict(os.environ, {}, clear=True):
                result = server.ai_answer("How many variants are in the database?", session_id="test")
        self.assertEqual(result["provider"], "local")
        self.assertIn("variants", result["answer"].lower())


if __name__ == "__main__":
    unittest.main()
