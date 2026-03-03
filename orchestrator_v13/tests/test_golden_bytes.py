import unittest
from pathlib import Path

from orchestrator_v13 import run_orchestrator_v13
from orchestrator_v13.io import load_fixture
from orchestrator_v13.util import canonical_json_dumps


FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures"
GOLDEN_ROOT = Path(__file__).resolve().parent / "goldens"


class TestGoldenBytes(unittest.TestCase):
    def test_golden_bytes(self):
        fixture_files = sorted(FIXTURE_DIR.glob("T*.json"))
        self.assertGreaterEqual(len(fixture_files), 8)

        for fixture_path in fixture_files:
            with self.subTest(fixture=fixture_path.name):
                config, raw_log, _expected = load_fixture(fixture_path)
                out = run_orchestrator_v13(raw_log, config)
                target_dir = GOLDEN_ROOT / fixture_path.stem
                self.assertTrue(target_dir.exists(), f"missing golden dir: {target_dir}")

                expected = {
                    "frontier_state.json": canonical_json_dumps(out.frontier_state),
                    "partition_plan.json": canonical_json_dumps(out.partition_plan),
                    "repair_queue.json": canonical_json_dumps(out.repair_queue),
                    "v13_obligations.json": canonical_json_dumps(out.v13_obligations),
                }
                for fname, blob in expected.items():
                    got = (target_dir / fname).read_bytes()
                    self.assertEqual(
                        got,
                        blob,
                        msg=f"golden mismatch for {fixture_path.stem}/{fname}",
                    )


if __name__ == "__main__":
    unittest.main()
