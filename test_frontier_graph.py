"""Tests for the project-local proof-frontier graph."""

from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

from scripts.frontier_graph import GraphValidationError, validate_graph


GRAPH_PATH = Path(__file__).parent / "proof_graph" / "erdos993_frontier.json"


def load_live_graph() -> dict:
    return json.loads(GRAPH_PATH.read_text(encoding="utf-8"))


class TestFrontierGraph(unittest.TestCase):
    def test_live_graph_validates(self) -> None:
        warnings = validate_graph(load_live_graph())
        self.assertTrue(any("status=open" in warning for warning in warnings))

    def test_dependency_cycle_is_rejected(self) -> None:
        graph = copy.deepcopy(load_live_graph())
        graph["edges"].append(
            {
                "from": "depth3_margin_nonnegative",
                "to": "depth3_rung_lc",
                "type": "depends_on",
            }
        )
        with self.assertRaisesRegex(GraphValidationError, "dependency cycle"):
            validate_graph(graph)

    def test_grounded_status_requires_evidence(self) -> None:
        graph = copy.deepcopy(load_live_graph())
        claim = next(item for item in graph["claims"] if item["id"] == "extendable_pascal_reserve")
        claim["evidence"] = []
        with self.assertRaisesRegex(GraphValidationError, "requires grounded evidence"):
            validate_graph(graph)

    def test_incomplete_bridge_is_rejected(self) -> None:
        graph = copy.deepcopy(load_live_graph())
        del graph["bridges"][0]["kill_test"]
        with self.assertRaisesRegex(GraphValidationError, "kill_test"):
            validate_graph(graph)


if __name__ == "__main__":
    unittest.main()
