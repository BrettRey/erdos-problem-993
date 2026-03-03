import tempfile
import unittest
from pathlib import Path

from orchestrator_v13 import emit_artifacts, run_orchestrator_v13, verify_replay
from orchestrator_v13.io import load_fixture
from orchestrator_v13.util import canonical_json_dumps


FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures"


class TestFixtureMatrix(unittest.TestCase):
    def test_t0_to_t7_matrix(self):
        fixture_files = sorted(FIXTURE_DIR.glob("T*.json"))
        self.assertGreaterEqual(len(fixture_files), 8)

        for path in fixture_files:
            with self.subTest(fixture=path.name):
                config, raw_log, expected = load_fixture(path)
                out = run_orchestrator_v13(raw_log, config)

                self.assertEqual(
                    out.partition_plan.selected_partition_id.value,
                    expected["selected_partition_id"],
                )
                self.assertEqual(
                    out.v13_obligations.global_closure.status.value,
                    expected["closure_status"],
                )

                with tempfile.TemporaryDirectory() as td:
                    out_dir = Path(td)
                    emit_artifacts(out, out_dir)
                    errs = verify_replay(raw_log, config, out_dir)
                    self.assertEqual(errs, [])

                    # Byte-stability check: second run must match exactly.
                    out2 = run_orchestrator_v13(raw_log, config)
                    self.assertEqual(
                        canonical_json_dumps(out.frontier_state),
                        canonical_json_dumps(out2.frontier_state),
                    )
                    self.assertEqual(
                        canonical_json_dumps(out.partition_plan),
                        canonical_json_dumps(out2.partition_plan),
                    )
                    self.assertEqual(
                        canonical_json_dumps(out.repair_queue),
                        canonical_json_dumps(out2.repair_queue),
                    )
                    self.assertEqual(
                        canonical_json_dumps(out.v13_obligations),
                        canonical_json_dumps(out2.v13_obligations),
                    )


if __name__ == "__main__":
    unittest.main()
