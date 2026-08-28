"""Tests for the Vatter-style residual-tail probe."""

from __future__ import annotations

import os
import sys
import unittest

_ROOT = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.join(_ROOT, "scripts")
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from indpoly import independence_poly
from probe_vatter_residual_tail import (
    check_all_leaf_orders,
    first_strict_descent,
    leaf_peeling_order,
    residual_polynomials,
    tail_ascent,
)
from trees import trees


class TestVatterResidualTail(unittest.TestCase):
    def test_first_descent_and_tail_ascent(self) -> None:
        self.assertEqual(first_strict_descent([1, 5, 4, 4, 2]), 2)
        self.assertIsNone(first_strict_descent([1, 2, 2, 3]))
        self.assertEqual(tail_ascent([1, 6, 5, 7, 2], 2), (3, 5, 7))
        self.assertIsNone(tail_ascent([1, 6, 5, 5, 2], 2))

    def test_least_vertex_coefficient_identity(self) -> None:
        # Path 0--1--2--3, with a nontrivial leaf-peeling order.
        adj = [[1], [0, 2], [1, 3], [2]]
        order = leaf_peeling_order(adj, max)
        parent = independence_poly(len(adj), adj)
        residuals = residual_polynomials(adj, order)
        for k in range(1, len(parent)):
            total = sum(
                poly[k - 1] if k - 1 < len(poly) else 0
                for _, _, poly in residuals
            )
            self.assertEqual(parent[k], total)

    def test_all_leaf_orders_through_n8(self) -> None:
        for n in range(1, 9):
            for _, adj in trees(n):
                failure, _ = check_all_leaf_orders(adj)
                self.assertIsNone(failure)


if __name__ == "__main__":
    unittest.main()
