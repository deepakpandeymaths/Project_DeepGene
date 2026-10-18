import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "data/ml/66_retry_liftover_api_checkpointed.py"
SPEC = importlib.util.spec_from_file_location("retry_checkpointed", SCRIPT)
retry = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(retry)


class GnomadRetryPolicyTests(unittest.TestCase):
    def test_pending_targets_excludes_terminal_and_maxed_attempts(self):
        targets = [{"lifted_variant_id": "a"}, {"lifted_variant_id": "b"}, {"lifted_variant_id": "c"}]
        latest = {
            "a": {"status": "matched", "attempt": 1},
            "b": {"status": "request_error", "attempt": 3},
            "c": {"status": "request_error", "attempt": 1},
        }
        result = retry.pending_targets(targets, {"a", "b", "c"}, latest, 3)
        self.assertEqual([item["lifted_variant_id"] for item in result], ["c"])

    def test_checkpoint_loader_keeps_latest_attempt(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "checkpoint.jsonl"
            path.write_text(
                json.dumps({"lifted_variant_id": "a", "attempt": 1, "status": "request_error"}) + "\n"
                + json.dumps({"lifted_variant_id": "a", "attempt": 2, "status": "matched"}) + "\n",
                encoding="utf-8",
            )
            latest = retry.load_latest_checkpoint(path)
        self.assertEqual(latest["a"]["attempt"], 2)
        self.assertEqual(latest["a"]["status"], "matched")


if __name__ == "__main__":
    unittest.main()

