"""Regression tests for the once-subdivided double-star theorem."""

import unittest

from indpoly import independence_poly, is_log_concave
from verify_once_subdivided_double_star_lc_20260808 import (
    index_three_gap,
    once_subdivided_adj,
    once_subdivided_factorization,
    once_subdivided_formula,
    symmetric_scaled_gap,
)


class TestOnceSubdividedDoubleStar(unittest.TestCase):
    def test_closed_formula_matches_tree_dp(self) -> None:
        for p in range(7):
            for q in range(7):
                adj = once_subdivided_adj(p, q)
                self.assertEqual(
                    independence_poly(len(adj), adj),
                    once_subdivided_formula(p, q),
                )

    def test_proof_factorization(self) -> None:
        for p in range(31):
            for q in range(31):
                self.assertEqual(
                    once_subdivided_formula(p, q),
                    once_subdivided_factorization(p, q),
                )

    def test_exact_parameter_box(self) -> None:
        for p in range(81):
            for q in range(81):
                self.assertTrue(is_log_concave(once_subdivided_formula(p, q)))

    def test_symmetric_gap_identity(self) -> None:
        for p in range(101):
            for q in range(101):
                self.assertEqual(
                    144 * index_three_gap(p, q),
                    symmetric_scaled_gap(p, q),
                )


if __name__ == "__main__":
    unittest.main()
