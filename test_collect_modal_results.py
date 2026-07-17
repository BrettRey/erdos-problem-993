"""Focused offline tests for the Modal result collector."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.collect_modal_results import (
    audit_partition_keys,
    build_raw_export,
    load_raw_export,
    require_exact_partition_keys,
    sorted_raw_rows,
    summarize_unimodality,
)


def unimodality_row(count: int) -> dict[str, object]:
    return {"count": count, "counterexample": None}


class TestModalPartitionAudit(unittest.TestCase):
    def test_equal_length_wrong_key_set_is_not_complete(self) -> None:
        rows = {
            "0/2": unimodality_row(3),
            "9/2": unimodality_row(99),
        }
        audit = audit_partition_keys(rows, workers=2)
        self.assertFalse(audit["exact_match"])
        self.assertEqual(audit["missing_keys"], ["1/2"])
        self.assertEqual(audit["unexpected_keys"], ["9/2"])

        summary = summarize_unimodality(rows, n=2, workers=2, expected=4)
        self.assertFalse(summary["complete"])
        self.assertEqual(summary["completed"], 1)
        self.assertEqual(summary["total_trees_so_far"], 3)
        with self.assertRaises(ValueError):
            require_exact_partition_keys(rows, workers=2)

    def test_raw_rows_are_sorted_by_numeric_residue(self) -> None:
        rows = {
            "10/12": unimodality_row(10),
            "2/12": unimodality_row(2),
            "0/12": unimodality_row(0),
        }
        self.assertEqual(
            [row["key"] for row in sorted_raw_rows(rows)],
            ["0/12", "2/12", "10/12"],
        )


class TestModalRawExport(unittest.TestCase):
    def test_round_trip_and_hash_validation(self) -> None:
        rows = {
            "0/2": unimodality_row(3),
            "1/2": unimodality_row(4),
        }
        payload = build_raw_export(
            rows,
            kind="unimodality",
            n=2,
            workers=2,
            expected=7,
            name="test-modal-dict",
        )
        self.assertTrue(payload["partition_key_audit"]["exact_match"])

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "raw.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            loaded_rows, loaded_payload = load_raw_export(path)
            self.assertEqual(loaded_rows, rows)
            self.assertEqual(loaded_payload["rows_sha256"], payload["rows_sha256"])

            tampered = json.loads(path.read_text(encoding="utf-8"))
            tampered["rows"][0]["value"]["count"] = 999
            path.write_text(json.dumps(tampered), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "content_sha256 mismatch"):
                load_raw_export(path)


if __name__ == "__main__":
    unittest.main()
