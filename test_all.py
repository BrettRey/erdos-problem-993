"""Tests for the Erdős Problem #993 search code."""

import shutil
import unittest

from graph6 import parse_graph6
from indpoly import (
    independence_poly,
    is_log_concave,
    is_unimodal,
    log_concavity_ratio,
    near_miss_ratio,
)
from trees import trees


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


if __name__ == "__main__":
    unittest.main()
