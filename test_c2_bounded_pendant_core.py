"""Regression tests for the bounded-pendant-core C2 certificate."""

from __future__ import annotations

import unittest

from indpoly import independence_poly
from verify_c2_bounded_pendant_core_20260808 import (
    connector_adjacency,
    connector_formula,
    integer_partitions,
    rooted_arm_state,
    verify,
)


class TestC2BoundedPendantCore(unittest.TestCase):
    def test_integer_partitions(self) -> None:
        self.assertEqual(
            list(integer_partitions(5)),
            [
                (1, 1, 1, 1, 1),
                (1, 1, 1, 2),
                (1, 1, 3),
                (1, 2, 2),
                (1, 4),
                (2, 3),
                (5,),
            ],
        )

    def test_formula_against_tree_dp(self) -> None:
        left = (1, 2, 4)
        right = (2, 3)
        left_state = (left, *rooted_arm_state(left))
        right_state = (right, *rooted_arm_state(right))
        for internal in range(7):
            adjacency = connector_adjacency(left, right, internal)
            self.assertEqual(
                independence_poly(len(adjacency), adjacency),
                connector_formula(left_state, right_state, internal),
            )

    def test_small_complete_certificate(self) -> None:
        report = verify(8, 5, 3)
        self.assertEqual(report["maximum_total_pendant_vertices"], 8)
        self.assertEqual(report["unordered_arm_vector_pairs"], 223)
        self.assertEqual(
            len(report["canonical_enumeration_sha256"]),
            64,
        )


if __name__ == "__main__":
    unittest.main()
