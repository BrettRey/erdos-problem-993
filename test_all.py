"""Tests for the Erdős Problem #993 search code."""

import contextlib
from fractions import Fraction
import importlib.util
import io
import math
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

from graph6 import parse_graph6
from erdos993_hunt import independence_poly_tree
from indpoly import (
    independence_poly,
    is_log_concave,
    is_unimodal,
    log_concavity_ratio,
    near_miss_ratio,
)
from nm_optimizer import _selection_score, _update_archive, run_optimizer
from trees import trees

_ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS_DIR = os.path.join(_ROOT_DIR, "scripts")
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

_LC_BREAKER_SPEC = importlib.util.spec_from_file_location(
    "lc_breaker_optimizer",
    os.path.join(_SCRIPTS_DIR, "lc_breaker_optimizer.py"),
)
if _LC_BREAKER_SPEC is None or _LC_BREAKER_SPEC.loader is None:
    raise RuntimeError("Unable to load scripts/lc_breaker_optimizer.py")
lc_breaker_optimizer = importlib.util.module_from_spec(_LC_BREAKER_SPEC)
_LC_BREAKER_SPEC.loader.exec_module(lc_breaker_optimizer)

_ANALYZE_CORPUS_SPEC = importlib.util.spec_from_file_location(
    "analyze_prufer_corpus",
    os.path.join(_SCRIPTS_DIR, "analyze_prufer_corpus.py"),
)
if _ANALYZE_CORPUS_SPEC is None or _ANALYZE_CORPUS_SPEC.loader is None:
    raise RuntimeError("Unable to load scripts/analyze_prufer_corpus.py")
analyze_prufer_corpus = importlib.util.module_from_spec(_ANALYZE_CORPUS_SPEC)
sys.modules[_ANALYZE_CORPUS_SPEC.name] = analyze_prufer_corpus
_ANALYZE_CORPUS_SPEC.loader.exec_module(analyze_prufer_corpus)

_ANALYZE_SIGNED_SPEC = importlib.util.spec_from_file_location(
    "analyze_signed_conditionals",
    os.path.join(_SCRIPTS_DIR, "analyze_signed_conditionals.py"),
)
if _ANALYZE_SIGNED_SPEC is None or _ANALYZE_SIGNED_SPEC.loader is None:
    raise RuntimeError("Unable to load scripts/analyze_signed_conditionals.py")
analyze_signed_conditionals = importlib.util.module_from_spec(_ANALYZE_SIGNED_SPEC)
sys.modules[_ANALYZE_SIGNED_SPEC.name] = analyze_signed_conditionals
_ANALYZE_SIGNED_SPEC.loader.exec_module(analyze_signed_conditionals)

from probe_signed_pb_reserve import first_descent, signed_metric  # noqa: E402


class TestGraph6(unittest.TestCase):
    """Test graph6 parsing."""

    def test_single_vertex(self):
        # graph6 for K_1: n=1, no edges
        n, adj = parse_graph6(b"@")
        self.assertEqual(n, 1)
        self.assertEqual(adj, [[]])

    def test_single_edge(self):
        # graph6 for K_2: n=2, one edge
        n, adj = parse_graph6(b"A_")
        self.assertEqual(n, 2)
        self.assertEqual(adj, [[1], [0]])

    def test_path_3(self):
        # P_3: 0-1-2, graph6 = "Bg"
        n, adj = parse_graph6(b"Bg")
        self.assertEqual(n, 3)
        # edges: 0-1, 1-2
        self.assertIn(1, adj[0])
        self.assertIn(0, adj[1])
        self.assertIn(2, adj[1])
        self.assertIn(1, adj[2])

    def test_star_k13(self):
        # K_{1,3}: vertex 0 connected to 1,2,3. graph6 = "Cs"
        n, adj = parse_graph6(b"Cs")
        self.assertEqual(n, 4)
        self.assertEqual(sorted(adj[0]), [1, 2, 3])
        for v in [1, 2, 3]:
            self.assertEqual(adj[v], [0])

    def test_trailing_newline(self):
        n, adj = parse_graph6(b"@\n")
        self.assertEqual(n, 1)


class TestAnalyzeCorpusHarness(unittest.TestCase):
    """Focused checks for adversarial corpus and family harness helpers."""

    def test_graph6_encoder_long_order_roundtrip(self):
        adj = [[] for _ in range(63)]
        encoded = analyze_prufer_corpus.graph6_from_adj(adj)
        n, parsed = parse_graph6(encoded.encode("ascii"))
        self.assertEqual(n, 63)
        self.assertEqual(parsed, adj)

    def test_galvin_generator_order_and_alpha(self):
        adj = analyze_prufer_corpus.make_galvin_tree(4, 3)
        poly = independence_poly(len(adj), adj)
        self.assertEqual(len(adj), 29)
        self.assertEqual(len(poly) - 1, 16)

    def test_bautista_ramos_figure_example(self):
        adj = analyze_prufer_corpus.make_bautista_ramos_tree(2, 5)
        poly = independence_poly(len(adj), adj)
        defects = {defect["k"] for defect in analyze_prufer_corpus.lc_defects(poly)}
        self.assertEqual(len(adj), 70)
        self.assertEqual(len(poly) - 1, 37)
        self.assertEqual(defects, {34, 36})
        self.assertTrue(is_unimodal(poly))

    def test_li_kadrawi_levit_witnesses(self):
        cases = [
            (analyze_prufer_corpus.make_li_tree(4, 4, starred=False), 14),
            (analyze_prufer_corpus.make_li_tree(3, 4, starred=True), 14),
        ]
        for adj, alpha in cases:
            poly = independence_poly(len(adj), adj)
            defects = {defect["k"] for defect in analyze_prufer_corpus.lc_defects(poly)}
            self.assertEqual(len(adj), 26)
            self.assertEqual(len(poly) - 1, alpha)
            self.assertEqual(defects, {13})
            self.assertTrue(is_unimodal(poly))


class TestSignedConditionalReduction(unittest.TestCase):
    """Regression checks for corrected signed conditional reductions."""

    @staticmethod
    def _convolve_fraction(
        left: list[Fraction],
        right: list[Fraction],
    ) -> list[Fraction]:
        out = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
        for i, left_value in enumerate(left):
            for j, right_value in enumerate(right):
                out[i + j] += left_value * right_value
        return out

    @classmethod
    def _exact_grouped_pmf(
        cls,
        blocks: list[tuple[int, Fraction]],
    ) -> list[Fraction]:
        pmf = [Fraction(1)]
        for count, p in blocks:
            block = [
                Fraction(math.comb(count, k)) * p**k * (1 - p) ** (count - k)
                for k in range(count + 1)
            ]
            pmf = cls._convolve_fraction(pmf, block)
        return pmf

    @staticmethod
    def _exact_first_descent(
        seq: list[Fraction],
        *,
        weak: bool = False,
    ) -> int | None:
        for k in range(1, len(seq)):
            if seq[k] < seq[k - 1] or (weak and seq[k] == seq[k - 1]):
                return k
        return None

    @classmethod
    def _exact_signed_profile(
        cls,
        x_blocks: list[tuple[int, Fraction]],
        y_blocks: list[tuple[int, Fraction]],
        *,
        weak: bool = False,
    ) -> tuple[list[Fraction], list[Fraction], list[Fraction], int, int]:
        x_pmf = cls._exact_grouped_pmf(x_blocks)
        y_pmf = cls._exact_grouped_pmf(y_blocks)
        signed = cls._convolve_fraction(x_pmf, list(reversed(y_pmf)))
        descent = cls._exact_first_descent(signed, weak=weak)
        if descent is None:
            raise AssertionError("exact signed profile has no descent")
        support_start = -(len(y_pmf) - 1)
        return x_pmf, y_pmf, signed, descent, support_start + descent

    @staticmethod
    def _exact_ratio(seq: list[Fraction], index: int) -> Fraction:
        if index < 0 or index >= len(seq) - 1:
            return Fraction(0)
        return seq[index + 1] / seq[index]

    @staticmethod
    def _exact_reverse_ratio(seq: list[Fraction], index: int) -> Fraction:
        if index <= 0 or index >= len(seq):
            return Fraction(0)
        return seq[index - 1] / seq[index]

    @staticmethod
    def _exact_expectation(
        weighted_values: list[tuple[int, Fraction]],
        fn,
    ) -> Fraction:
        total_weight = sum((weight for _, weight in weighted_values), Fraction(0))
        return sum(
            (weight * fn(value) for value, weight in weighted_values),
            Fraction(0),
        ) / total_weight

    @classmethod
    def _exact_side_reduction_components(
        cls,
        x_pmf: list[Fraction],
        y_pmf: list[Fraction],
        z: int,
    ) -> dict[str, Fraction]:
        lo = max(0, -z)
        hi = min(len(y_pmf) - 1, len(x_pmf) - 1 - z)
        conditional = [
            (y, x_pmf[z + y] * y_pmf[y])
            for y in range(lo, hi + 1)
        ]
        c_z = sum((weight for _, weight in conditional), Fraction(0))

        x_lower_index = -z - 1
        x_boundary = (
            x_pmf[0] * y_pmf[x_lower_index] / c_z
            if 0 <= x_lower_index < len(y_pmf)
            else Fraction(0)
        )
        y_upper_index = len(y_pmf) - 1
        y_upper_x_index = z + 1 + y_upper_index
        y_boundary = (
            x_pmf[y_upper_x_index] * y_pmf[y_upper_index] / c_z
            if 0 <= y_upper_x_index < len(x_pmf)
            else Fraction(0)
        )
        x_upper_index = len(x_pmf) - 1
        x_upper_y_index = x_upper_index + 1 - z
        x_reciprocal_beta = (
            x_pmf[x_upper_index] * y_pmf[x_upper_y_index] / c_z
            if 0 <= x_upper_y_index < len(y_pmf)
            else Fraction(0)
        )
        y_lower_x_index = z - 1
        y_reciprocal_beta = (
            x_pmf[y_lower_x_index] * y_pmf[0] / c_z
            if 0 <= y_lower_x_index < len(x_pmf)
            else Fraction(0)
        )

        x_forward_mean = cls._exact_expectation(
            conditional,
            lambda y: cls._exact_ratio(x_pmf, z + y),
        )
        x_reciprocal_core = cls._exact_expectation(
            conditional,
            lambda y: cls._exact_reverse_ratio(x_pmf, z + y),
        )
        x_forward_reciprocal_mean = cls._exact_expectation(
            conditional,
            lambda y: (
                cls._exact_ratio(x_pmf, z + y)
                * cls._exact_reverse_ratio(x_pmf, z + y)
            ),
        )
        x_inverse_index_gain = cls._exact_expectation(
            conditional,
            lambda y: Fraction(1, z + y + 1) if z + y + 1 > 0 else Fraction(0),
        )
        x_dispersion_penalty = (
            x_forward_mean * x_reciprocal_core - x_forward_reciprocal_mean
        )
        x_lower_boundary_penalty = x_boundary * x_reciprocal_core
        x_upper_boundary_penalty = (
            x_forward_mean + x_boundary
        ) * x_reciprocal_beta
        x_reduction_bound = (
            x_inverse_index_gain
            - x_dispersion_penalty
            - x_lower_boundary_penalty
            - x_upper_boundary_penalty
        )

        y_backward_mean = cls._exact_expectation(
            conditional,
            lambda y: cls._exact_reverse_ratio(y_pmf, y),
        )
        y_reciprocal_core = cls._exact_expectation(
            conditional,
            lambda y: cls._exact_ratio(y_pmf, y),
        )
        y_backward_reciprocal_mean = cls._exact_expectation(
            conditional,
            lambda y: (
                cls._exact_reverse_ratio(y_pmf, y)
                * cls._exact_ratio(y_pmf, y)
            ),
        )
        y_inverse_index_gain = cls._exact_expectation(
            conditional,
            lambda y: Fraction(1, y + 1),
        )
        y_dispersion_penalty = (
            y_backward_mean * y_reciprocal_core - y_backward_reciprocal_mean
        )
        y_lower_reciprocal_penalty = y_backward_mean * y_reciprocal_beta
        y_upper_boundary_penalty = y_boundary * (
            y_reciprocal_core + y_reciprocal_beta
        )
        y_reduction_bound = (
            y_inverse_index_gain
            - y_dispersion_penalty
            - y_lower_reciprocal_penalty
            - y_upper_boundary_penalty
        )

        return {
            "x_inverse_index_gain": x_inverse_index_gain,
            "x_dispersion_penalty": x_dispersion_penalty,
            "x_lower_boundary_penalty": x_lower_boundary_penalty,
            "x_upper_boundary_penalty": x_upper_boundary_penalty,
            "x_reduction_bound": x_reduction_bound,
            "y_inverse_index_gain": y_inverse_index_gain,
            "y_dispersion_penalty": y_dispersion_penalty,
            "y_lower_reciprocal_penalty": y_lower_reciprocal_penalty,
            "y_upper_boundary_penalty": y_upper_boundary_penalty,
            "y_reduction_bound": y_reduction_bound,
        }

    @staticmethod
    def _exact_effective_drop(
        signed: list[Fraction],
        descent: int,
    ) -> Fraction:
        previous_ratio = signed[descent] / signed[descent - 1]
        pressure = signed[descent + 1] / signed[descent]
        return 1 - pressure / previous_ratio

    @classmethod
    def _exact_signed_scaled_slacks(
        cls,
        x_blocks: list[tuple[int, Fraction]],
        y_blocks: list[tuple[int, Fraction]],
    ) -> tuple[Fraction, Fraction]:
        _, _, signed, descent, _ = cls._exact_signed_profile(
            x_blocks,
            y_blocks,
        )
        pressure = signed[descent + 1] / signed[descent]
        variance = sum(count * p * (1 - p) for count, p in [*x_blocks, *y_blocks])
        return variance * cls._exact_effective_drop(signed, descent), variance * (
            1 - pressure
        )

    def test_fair_binomial_signed_fallback_constants(self):
        best_delta = None
        best_reserve = None
        for total_count in range(4, 41):
            probabilities = [
                Fraction(math.comb(total_count, k), 2**total_count)
                for k in range(total_count + 1)
            ]
            descent = next(
                k
                for k in range(1, total_count + 1)
                if probabilities[k] < probabilities[k - 1]
            )
            previous_ratio = probabilities[descent] / probabilities[descent - 1]
            pressure = (
                probabilities[descent + 1] / probabilities[descent]
                if descent + 1 <= total_count
                else Fraction(0)
            )
            effective_drop = 1 - pressure / previous_ratio
            variance = Fraction(total_count, 4)
            scaled_effective = variance * effective_drop
            scaled_reserve = variance * (1 - pressure)
            best_delta = (
                (total_count, scaled_effective)
                if best_delta is None or scaled_effective < best_delta[1]
                else best_delta
            )
            best_reserve = (
                (total_count, scaled_reserve)
                if best_reserve is None or scaled_reserve < best_reserve[1]
                else best_reserve
            )

        self.assertEqual(best_delta, (4, Fraction(5, 8)))
        self.assertEqual(best_reserve, (4, Fraction(3, 4)))

    def test_first_descent_treats_roundoff_plateau_as_tie(self):
        seq = np.array([0.4, 0.5, np.nextafter(0.5, 0.0), 0.4])
        self.assertEqual(first_descent(seq), 3)

    def test_analyzer_counts_terminal_descent_rows(self):
        metric_counters = {
            "no_descent": 0,
            "terminal_descent": 0,
            "nonpositive_descent_mass": 0,
            "zero_variance": 0,
        }
        self.assertIsNone(
            signed_metric(
                x_blocks=[(1, 0.2)],
                y_blocks=[(0, 0.0)],
                kind="terminal_descent_regression",
                counters=metric_counters,
            )
        )
        self.assertEqual(metric_counters["terminal_descent"], 1)

        counters = {
            "no_descent": 0,
            "terminal_descent": 0,
            "nonpositive_descent_mass": 0,
        }
        analysis = analyze_signed_conditionals.analyze_row(
            {
                "x_blocks": [[1, 0.2]],
                "y_blocks": [[0, 0.0]],
                "x_n": 1,
                "y_n": 0,
                "x_mean": 0.2,
                "y_mean": 0.0,
                "x_variance": 0.16,
                "y_variance": 0.0,
                "variance": 0.16,
            },
            "terminal_descent_regression",
            counters=counters,
        )
        self.assertIsNone(analysis)
        self.assertEqual(counters["terminal_descent"], 1)

    def test_x_reduction_keeps_upper_boundary_beta(self):
        row = signed_metric(
            x_blocks=[(1, 0.5)],
            y_blocks=[(10, 0.5)],
            kind="x_boundary_counterexample",
        )
        self.assertIsNotNone(row)
        analysis = analyze_signed_conditionals.analyze_row(
            row,
            "x_boundary_counterexample",
        )
        self.assertIsNotNone(analysis)

        old_omitted_bound = (
            analysis["x_inverse_index_gain"]
            - analysis["x_dispersion_penalty"]
            - analysis["x_lower_boundary_penalty"]
        )
        self.assertEqual(analysis["first_descent_value"], -3)
        self.assertAlmostEqual(analysis["effective_ratio_drop"], 3.0 / 10.0)
        self.assertAlmostEqual(old_omitted_bound, 4.0 / 11.0)
        self.assertAlmostEqual(
            analysis["x_reciprocal_upper_boundary_beta"],
            42.0 / 55.0,
        )
        self.assertAlmostEqual(analysis["x_reduction_bound"], -1.0 / 55.0)
        self.assertLessEqual(
            analysis["y_reduction_bound"],
            analysis["effective_ratio_drop"] + 1e-12,
        )

    def test_half_heavy_dust_near_misses_are_not_failures(self):
        cases = [
            (
                [(1, 0.25), (6, 0.5)],
                [(4, 0.5)],
                0.23415870316380516,
                0.7852688694525427,
                0.9506802721088434,
            ),
            (
                [(5, 0.004), (7, 0.5)],
                [(1, 0.02), (4, 0.5)],
                0.24306444597986393,
                0.7896792287593093,
                0.789746777744473,
            ),
        ]
        for x_blocks, y_blocks, side_bound, effective_drop, reserve in cases:
            row = signed_metric(
                x_blocks=x_blocks,
                y_blocks=y_blocks,
                kind="half_heavy_dust_near_miss",
            )
            self.assertIsNotNone(row)
            analysis = analyze_signed_conditionals.analyze_row(
                row,
                "half_heavy_dust_near_miss",
            )
            self.assertIsNotNone(analysis)
            self.assertLess(analysis["variance_times_best_side_reduction_bound"], 0.25)
            self.assertAlmostEqual(
                analysis["variance_times_best_side_reduction_bound"],
                side_bound,
            )
            self.assertAlmostEqual(
                analysis["variance_times_effective_ratio_drop"],
                effective_drop,
            )
            self.assertAlmostEqual(analysis["variance_times_reserve"], reserve)
            self.assertGreater(analysis["variance_times_effective_ratio_drop"], 0.75)
            self.assertGreater(analysis["variance_times_reserve"], 0.75)

    def test_half_heavy_side_reduction_disproof_exactly(self):
        x_blocks = [(1, Fraction(1, 4)), (6, Fraction(1, 2))]
        y_blocks = [(4, Fraction(1, 2))]
        x_pmf, y_pmf, signed, descent, z = self._exact_signed_profile(
            x_blocks,
            y_blocks,
        )
        variance = sum(
            count * p * (1 - p) for count, p in [*x_blocks, *y_blocks]
        )
        components = self._exact_side_reduction_components(x_pmf, y_pmf, z)
        x_bound = (
            components["x_inverse_index_gain"]
            - components["x_dispersion_penalty"]
            - components["x_lower_boundary_penalty"]
            - components["x_upper_boundary_penalty"]
        )
        y_bound = (
            components["y_inverse_index_gain"]
            - components["y_dispersion_penalty"]
            - components["y_lower_reciprocal_penalty"]
            - components["y_upper_boundary_penalty"]
        )

        self.assertEqual(z, 2)
        self.assertEqual(variance, Fraction(43, 16))
        self.assertEqual(x_bound, components["x_reduction_bound"])
        self.assertEqual(y_bound, components["y_reduction_bound"])
        self.assertEqual(x_bound, Fraction(399961, 4594590))
        self.assertEqual(y_bound, Fraction(8069, 92610))

        scaled_best = variance * max(x_bound, y_bound)
        self.assertEqual(scaled_best, Fraction(346967, 1481760))
        self.assertLess(scaled_best, Fraction(1, 4))
        self.assertEqual(Fraction(1, 4) - scaled_best, Fraction(23473, 1481760))

        pressure = signed[descent + 1] / signed[descent]
        self.assertEqual(
            variance * self._exact_effective_drop(signed, descent),
            Fraction(19393, 24696),
        )
        self.assertEqual(variance * (1 - pressure), Fraction(559, 588))

    def test_strict_descent_discontinuity_at_binomial_plateau_exactly(self):
        x_blocks = [(5, Fraction(1, 2))]
        y_at_plateau = [(1, Fraction(0))]
        _, _, base_signed, base_descent, base_z = self._exact_signed_profile(
            x_blocks,
            y_at_plateau,
        )
        _, _, base_signed_weak, base_weak_descent, base_weak_z = (
            self._exact_signed_profile(x_blocks, y_at_plateau, weak=True)
        )
        self.assertEqual(base_z, 4)
        self.assertEqual(
            self._exact_effective_drop(base_signed, base_descent),
            Fraction(3, 5),
        )
        self.assertEqual(base_weak_z, 3)
        self.assertEqual(
            self._exact_effective_drop(base_signed_weak, base_weak_descent),
            Fraction(1, 2),
        )

        q = Fraction(1, 1000)
        _, _, signed, descent, z = self._exact_signed_profile(
            x_blocks,
            [(1, q)],
        )
        _, _, signed_weak, weak_descent, weak_z = self._exact_signed_profile(
            x_blocks,
            [(1, q)],
            weak=True,
        )
        perturbed_drop = self._exact_effective_drop(signed, descent)
        self.assertEqual(z, 3)
        self.assertEqual(weak_z, 3)
        self.assertEqual(descent, weak_descent)
        self.assertEqual(signed, signed_weak)
        self.assertEqual(
            perturbed_drop,
            1 - Fraction(10) * (5 - 4 * q) / (10 - 5 * q) ** 2,
        )
        self.assertEqual(perturbed_drop, Fraction(1997601, 3996001))
        self.assertEqual(
            Fraction(1, 2) - perturbed_drop,
            Fraction(799, 7992002),
        )
        self.assertLess(
            abs(perturbed_drop - Fraction(1, 2)),
            abs(perturbed_drop - Fraction(3, 5)),
        )

    def test_half_heavy_dust_point_disproves_eight_tenths_exactly(self):
        scaled_effective, scaled_reserve = self._exact_signed_scaled_slacks(
            x_blocks=[(5, Fraction(1, 250)), (7, Fraction(1, 2))],
            y_blocks=[(1, Fraction(1, 50)), (4, Fraction(1, 2))],
        )
        self.assertGreater(scaled_effective, Fraction(3, 4))
        self.assertGreater(scaled_reserve, Fraction(3, 4))
        self.assertLess(scaled_effective, Fraction(4, 5))
        self.assertLess(scaled_reserve, Fraction(4, 5))


class TestIndependencePoly(unittest.TestCase):
    """Test independence polynomial computation."""

    def test_single_vertex(self):
        # P_1: i_0=1, i_1=1
        poly = independence_poly(1, [[]])
        self.assertEqual(poly, [1, 1])

    def test_path_2(self):
        # P_2 (single edge): i_0=1, i_1=2, i_2=0 → [1, 2]
        # Independent sets: {}, {0}, {1}. By size: i_0=1, i_1=2
        poly = independence_poly(2, [[1], [0]])
        self.assertEqual(poly, [1, 2])

    def test_path_3(self):
        # P_3: 0-1-2
        # Independent sets: {}, {0}, {1}, {2}, {0,2}
        # i_0=1, i_1=3, i_2=1
        poly = independence_poly(3, [[1], [0, 2], [1]])
        self.assertEqual(poly, [1, 3, 1])

    def test_path_4(self):
        # P_4: 0-1-2-3
        # Independent sets: {}, {0},{1},{2},{3}, {0,2},{0,3},{1,3}, {0,2} wait...
        # Careful: {0,2}, {0,3}, {1,3} — size 2 = 3
        # No size-3 independent set (0,2 blocks 1,3 — no, {0,2} is valid but
        # we can't add anything). Actually {1,3} is valid size 2.
        # i_0=1, i_1=4, i_2=3
        poly = independence_poly(4, [[1], [0, 2], [1, 3], [2]])
        self.assertEqual(poly, [1, 4, 3])

    def test_path_5(self):
        # P_5: 0-1-2-3-4
        # i_0=1, i_1=5, i_2=6, i_3=1
        poly = independence_poly(5, [[1], [0, 2], [1, 3], [2, 4], [3]])
        self.assertEqual(poly, [1, 5, 6, 1])

    def test_star_k13(self):
        # K_{1,3}: centre 0, leaves 1,2,3
        # Independent sets: {}, {0}, {1},{2},{3}, {1,2},{1,3},{2,3}, {1,2,3}
        # i_0=1, i_1=4, i_2=3, i_3=1
        adj = [[1, 2, 3], [0], [0], [0]]
        poly = independence_poly(4, adj)
        self.assertEqual(poly, [1, 4, 3, 1])

    def test_star_k14(self):
        # K_{1,4}: centre 0, leaves 1,2,3,4
        # i_0=1, i_1=5, i_2=6, i_3=4, i_4=1
        adj = [[1, 2, 3, 4], [0], [0], [0], [0]]
        poly = independence_poly(5, adj)
        self.assertEqual(poly, [1, 5, 6, 4, 1])

    def test_empty_graph(self):
        poly = independence_poly(0, [])
        self.assertEqual(poly, [1])


class TestUnimodality(unittest.TestCase):
    """Test unimodality checker."""

    def test_constant(self):
        self.assertTrue(is_unimodal([5, 5, 5]))

    def test_increasing(self):
        self.assertTrue(is_unimodal([1, 2, 3, 4]))

    def test_decreasing(self):
        self.assertTrue(is_unimodal([4, 3, 2, 1]))

    def test_peak(self):
        self.assertTrue(is_unimodal([1, 3, 5, 4, 2]))

    def test_plateau(self):
        self.assertTrue(is_unimodal([1, 3, 3, 3, 2]))

    def test_valley(self):
        # 1, 3, 2, 4 has a valley at index 2
        self.assertFalse(is_unimodal([1, 3, 2, 4]))

    def test_double_peak(self):
        # 1, 3, 2, 3, 1 — valley at index 2
        self.assertFalse(is_unimodal([1, 3, 2, 3, 1]))

    def test_short_sequences(self):
        self.assertTrue(is_unimodal([]))
        self.assertTrue(is_unimodal([1]))
        self.assertTrue(is_unimodal([1, 2]))
        self.assertTrue(is_unimodal([2, 1]))

    def test_plateau_then_rise(self):
        # 3, 3, 2, 3 — decreasing then increasing
        self.assertFalse(is_unimodal([3, 3, 2, 3]))


class TestLogConcavityAndNearMiss(unittest.TestCase):
    """Test log-concavity and near-miss utilities."""

    def test_log_concave_true(self):
        self.assertTrue(is_log_concave([1, 4, 6, 4, 1]))

    def test_log_concave_false(self):
        self.assertFalse(is_log_concave([1, 1, 2]))

    def test_log_concavity_ratio(self):
        ratio, pos = log_concavity_ratio([1, 1, 2])
        self.assertEqual(pos, 1)
        self.assertAlmostEqual(ratio, 2.0)

    def test_near_miss_monotone(self):
        self.assertEqual(near_miss_ratio([1, 2, 3]), (0.0, -1))

    def test_near_miss_unimodal(self):
        ratio, pos = near_miss_ratio([1, 3, 2, 1])
        self.assertEqual(pos, 2)
        self.assertAlmostEqual(ratio, 0.5)

    def test_near_miss_violation(self):
        ratio, pos = near_miss_ratio([1, 3, 2, 4])
        self.assertEqual(pos, 2)
        self.assertAlmostEqual(ratio, 2.0)


class TestTrees(unittest.TestCase):
    """Test tree enumeration."""

    def test_n1(self):
        result = list(trees(1, backend="networkx"))
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (1, [[]]))

    def test_n2(self):
        result = list(trees(2, backend="networkx"))
        self.assertEqual(len(result), 1)

    def test_n5(self):
        # 3 non-isomorphic trees on 5 vertices
        result = list(trees(5, backend="networkx"))
        self.assertEqual(len(result), 3)

    def test_n8(self):
        # 23 non-isomorphic trees on 8 vertices
        result = list(trees(8, backend="networkx"))
        self.assertEqual(len(result), 23)

    def test_n0(self):
        result = list(trees(0, backend="networkx"))
        self.assertEqual(len(result), 0)


class TestIntegration(unittest.TestCase):
    """Integration: all trees on n=8 should be unimodal."""

    def test_all_n8_unimodal(self):
        for n, adj in trees(8, backend="networkx"):
            poly = independence_poly(n, adj)
            self.assertTrue(
                is_unimodal(poly),
                f"Non-unimodal tree found on n={n}: {poly}",
            )

    def test_all_n10_unimodal(self):
        """All 106 trees on n=10 should be unimodal."""
        count = 0
        for n, adj in trees(10, backend="networkx"):
            poly = independence_poly(n, adj)
            self.assertTrue(
                is_unimodal(poly),
                f"Non-unimodal tree found on n={n}: {poly}",
            )
            count += 1
        self.assertEqual(count, 106)


class TestCrossBackend(unittest.TestCase):
    """Cross-check geng vs networkx on small n (if geng available)."""

    @unittest.skipUnless(shutil.which("geng"), "geng not available")
    def test_geng_vs_networkx_polys_n8_n10(self):
        for n in [8, 10]:
            polys_geng = sorted(
                tuple(independence_poly(n, adj))
                for n, adj in trees(n, backend="geng")
            )
            polys_nx = sorted(
                tuple(independence_poly(n, adj))
                for n, adj in trees(n, backend="networkx")
            )
            self.assertEqual(polys_geng, polys_nx)


class TestPolyConsistency(unittest.TestCase):
    """Cross-check indpoly vs erdos993_hunt polynomial implementations."""

    def test_indpoly_vs_hunt_n8(self):
        for n, adj in trees(8, backend="networkx"):
            edges = []
            for u in range(n):
                for v in adj[u]:
                    if u < v:
                        edges.append((u, v))
            poly_hunt = independence_poly_tree(n, edges, root=0)
            poly_ind = independence_poly(n, adj)
            self.assertEqual(poly_hunt, poly_ind)


class TestNmOptimizer(unittest.TestCase):
    """Focused checks for evolutionary nm search helpers."""

    def test_selection_score_rewards_observed_upside(self):
        self.assertAlmostEqual(_selection_score(0.82, 0.90, 0.5), 0.86)
        self.assertAlmostEqual(_selection_score(0.82, 0.80, 0.5), 0.82)

    def test_archive_keeps_best_distinct_tree(self):
        base = {
            "fingerprint": "a",
            "label": "seed_a",
            "nm": 0.80,
            "score": 0.81,
            "generation": 0,
        }
        improved_same = {
            "fingerprint": "a",
            "label": "seed_a_better",
            "nm": 0.83,
            "score": 0.85,
            "generation": 1,
        }
        distinct = {
            "fingerprint": "b",
            "label": "seed_b",
            "nm": 0.81,
            "score": 0.82,
            "generation": 1,
        }

        archive = _update_archive([base], [improved_same, distinct], archive_size=2)
        self.assertEqual(len(archive), 2)
        self.assertEqual(archive[0]["label"], "seed_a_better")
        self.assertEqual({item["fingerprint"] for item in archive}, {"a", "b"})

    def test_run_optimizer_records_score_and_archive(self):
        result = run_optimizer(
            n=12,
            pop_size=8,
            generations=2,
            elite_frac=0.25,
            mutations_per_ind=1,
            archive_size=4,
            prospect_bonus=0.5,
            seed=7,
            verbose=False,
        )
        self.assertFalse(result["counterexample_found"])
        self.assertEqual(result["archive_size"], 4)
        self.assertEqual(result["prospect_bonus"], 0.5)
        self.assertIn("best_score", result["history"][0])
        self.assertIn("archive_best_nm", result["history"][0])
        self.assertGreaterEqual(result["history"][0]["best_score"], result["history"][0]["best_nm"])


class TestLcBreakerOptimizer(unittest.TestCase):
    """Focused checks for the LC-breaker archive/prospect logic."""

    def test_selection_score_rewards_observed_lc_upside(self):
        score = lc_breaker_optimizer._selection_score(1.20, 1.32, 0.5)
        self.assertAlmostEqual(score, 1.26)
        self.assertAlmostEqual(
            lc_breaker_optimizer._selection_score(1.20, 1.10, 0.5),
            1.20,
        )

    def test_archive_keeps_best_distinct_lc_tree(self):
        base = {
            "fingerprint": "a",
            "origin": "seed_a",
            "lc_ratio": 1.4,
            "score": 1.45,
            "nm_ratio": 0.82,
            "generation": 0,
        }
        improved_same = {
            "fingerprint": "a",
            "origin": "seed_a_better",
            "lc_ratio": 1.5,
            "score": 1.55,
            "nm_ratio": 0.83,
            "generation": 1,
        }
        distinct = {
            "fingerprint": "b",
            "origin": "seed_b",
            "lc_ratio": 1.45,
            "score": 1.46,
            "nm_ratio": 0.80,
            "generation": 1,
        }

        archive = lc_breaker_optimizer._update_archive(
            [base],
            [improved_same, distinct],
            archive_size=2,
        )
        self.assertEqual(len(archive), 2)
        self.assertEqual(archive[0]["origin"], "seed_a_better")
        self.assertEqual({item["fingerprint"] for item in archive}, {"a", "b"})

    def test_run_records_score_and_archive(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            analysis_path = os.path.join(tmpdir, "analysis.json")
            out_path = os.path.join(tmpdir, "out.json")
            with open(analysis_path, "w", encoding="utf-8") as f:
                f.write('{"lc_failures": [{"graph6": "Cs"}]}')

            with contextlib.redirect_stdout(io.StringIO()):
                result = lc_breaker_optimizer.run(
                    n=4,
                    analysis_path=analysis_path,
                    out_path=out_path,
                    pop_size=8,
                    generations=2,
                    elite_frac=0.25,
                    archive_size=4,
                    prospect_bonus=0.5,
                    seed=7,
                    verbose_every=1,
                )

            self.assertEqual(result["archive_size"], 4)
            self.assertEqual(result["prospect_bonus"], 0.5)
            self.assertIn("best_score", result["history"][0])
            self.assertIn("archive_best_lc_ratio", result["history"][0])
            self.assertGreaterEqual(
                result["history"][0]["best_score"],
                result["history"][0]["best_lc_ratio"],
            )
            self.assertTrue(os.path.exists(out_path))


if __name__ == "__main__":
    unittest.main()
