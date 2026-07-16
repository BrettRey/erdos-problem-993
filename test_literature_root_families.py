"""Exact non-root tests for the literature-defined spectral stress corpus."""

import unittest

from indpoly import independence_poly, is_unimodal
from scripts.stress_literature_root_families import (
    balanced_tree_from_word,
    circle_family_poly,
    circle_recurrence_residual,
    compressed_state_from_word,
    consecutive_family_poly,
    consecutive_recurrence_residual,
    jerrum_patel_word,
    lc_failure_indices,
    pattern_repeat,
    poly_add,
    poly_trim,
    s2_tree,
    t1_tree,
    validate_tree,
)
from scripts.analyze_prufer_corpus import make_li_tree


class TestJerrumPatelConstruction(unittest.TestCase):
    def test_exact_word_order_and_compressed_dp(self) -> None:
        for subdivisions_k, height in ((0, 4), (1, 3), (2, 2)):
            word = jerrum_patel_word(subdivisions_k, height)
            adj = balanced_tree_from_word(word)
            expected_order = (
                (2 * subdivisions_k + 1) * 2 ** (height + 1)
                - (4 * subdivisions_k + 1)
            )
            self.assertEqual(len(adj), expected_order)
            self.assertEqual(validate_tree(adj)["max_degree"], 3)
            excluded, included = compressed_state_from_word(word)
            self.assertEqual(
                independence_poly(len(adj), adj),
                poly_add(excluded, included),
            )


class TestBautistaFamilies(unittest.TestCase):
    def test_circle_formula_and_recurrence(self) -> None:
        for k, offset in ((4, 0), (5, 1), (6, 2)):
            adj = make_li_tree(k, k + offset)
            poly = circle_family_poly(k, offset)
            self.assertEqual(independence_poly(len(adj), adj), poly)
            self.assertEqual(len(poly) - 1, 2 * k + offset + 6)
            self.assertEqual(poly_trim(circle_recurrence_residual(k, offset)), [0])

    def test_consecutive_threshold_cases(self) -> None:
        specs = (
            (2, 4, 2, 9, [92]),
            (3, 4, 2, 16, [162, 163]),
            (7, 5, 2, 13, [161, 162, 163]),
        )
        for ell, arms, width, depth, expected_breaks in specs:
            with self.subTest(ell=ell, arms=arms, width=width, depth=depth):
                base_adj, base_root = t1_tree(ell)
                pendant_adj, pendant_root = s2_tree(arms)
                adj = pattern_repeat(
                    base_adj,
                    base_root,
                    pendant_adj,
                    pendant_root,
                    width,
                    depth,
                )
                poly = consecutive_family_poly(ell, arms, width, depth)
                self.assertEqual(independence_poly(len(adj), adj), poly)
                self.assertEqual(lc_failure_indices(poly), expected_breaks)
                self.assertTrue(is_unimodal(poly))
                self.assertEqual(
                    poly_trim(
                        consecutive_recurrence_residual(
                            ell, arms, width, depth
                        )
                    ),
                    [0],
                )


if __name__ == "__main__":
    unittest.main()
