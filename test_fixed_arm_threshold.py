"""Tests for the effective fixed-arm unimodality theorem."""

from itertools import combinations_with_replacement
import unittest

from indpoly import independence_poly, is_unimodal
from scripts.fixed_arm_unimodality_threshold import (
    certified_threshold,
    early_difference,
    fixed_arm_polynomial,
    fixed_arm_terms,
    multi_arm_adjacency,
    path_independence_poly,
)


class TestFixedArmThreshold(unittest.TestCase):
    def test_path_polynomials(self) -> None:
        self.assertEqual(path_independence_poly(0), [1])
        self.assertEqual(path_independence_poly(1), [1, 1])
        self.assertEqual(path_independence_poly(5), [1, 5, 6, 1])

    def test_formula_matches_tree_dp(self) -> None:
        for arms in ((1,), (2, 3), (5, 5, 4, 2)):
            for s in range(5):
                adjacency = multi_arm_adjacency(s, arms)
                self.assertEqual(
                    fixed_arm_polynomial(s, arms),
                    independence_poly(len(adjacency), adjacency),
                )

    def test_certified_threshold_replay(self) -> None:
        for arms in ((1,), (2, 3, 6), (5, 5, 4, 2), (8,)):
            q, e = fixed_arm_terms(arms)
            threshold, crude, _ = certified_threshold(e)
            self.assertLessEqual(threshold, crude)
            r = len(e) - 1
            for s in (threshold, threshold + 1, threshold + 19):
                self.assertTrue(
                    all(early_difference(s, k, q, e) >= 0 for k in range(r + 1))
                )

    def test_small_vectors_falsification_sweep(self) -> None:
        for arity in range(1, 4):
            for arms in combinations_with_replacement(range(1, 6), arity):
                q, e = fixed_arm_terms(arms)
                threshold, _, _ = certified_threshold(e)
                r = len(e) - 1
                self.assertTrue(all(
                    early_difference(threshold, k, q, e) >= 0
                    for k in range(r + 1)
                ))
                for s in range(threshold + 3):
                    self.assertTrue(is_unimodal(fixed_arm_polynomial(s, arms)))


if __name__ == "__main__":
    unittest.main()
