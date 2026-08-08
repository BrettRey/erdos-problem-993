"""Regression tests for the exact double-star lemma."""

import unittest

from indpoly import independence_poly, is_log_concave
from verify_double_star_log_concavity_20260808 import (
    double_star_adj,
    double_star_factorization,
    double_star_formula,
    exceptional_gap,
)


class TestDoubleStarLogConcavity(unittest.TestCase):
    def test_closed_formula_matches_tree_dp(self) -> None:
        for p in range(7):
            for q in range(7):
                adj = double_star_adj(p, q)
                self.assertEqual(
                    independence_poly(len(adj), adj),
                    double_star_formula(p, q),
                )

    def test_proof_factorization(self) -> None:
        for p in range(31):
            for q in range(31):
                self.assertEqual(
                    double_star_formula(p, q),
                    double_star_factorization(p, q),
                )

    def test_exact_parameter_box(self) -> None:
        for p in range(81):
            for q in range(81):
                self.assertTrue(is_log_concave(double_star_formula(p, q)))

    def test_endpoint_gap_identities(self) -> None:
        for q in range(201):
            self.assertEqual(
                12 * exceptional_gap(q, 0),
                q * (q - 1) * (q * q - 3 * q + 8),
            )
            self.assertEqual(
                12 * exceptional_gap(q, q),
                q * (q + 1) * (q * q + q + 4),
            )


if __name__ == "__main__":
    unittest.main()
