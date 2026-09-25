"""Small, deterministic package regressions; not a new proof search."""

from __future__ import annotations

import unittest

import networkx as nx

from replay import analyze, replay
from scripts.audit_b1_zero_forbidden_core_20260905 import extendable_defects, polynomial
from scripts.audit_pendant_p2_boundary_20260904 import independent_masks
from indpoly import is_unimodal


def brute_counts(graph: nx.Graph) -> tuple[list[int], list[int]]:
    graph = nx.convert_node_labels_to_integers(graph)
    sets = independent_masks(graph)
    alpha = max(s.bit_count() for s in sets)
    maxima = [s for s in sets if s.bit_count() == alpha]
    poly = [sum(s.bit_count() == k for s in sets) for k in range(alpha + 1)]
    e = [sum(s.bit_count() == alpha - d and any(s & ~m == 0 for m in maxima)
             for s in sets) for d in range(5)]
    return poly, e


class PacketTests(unittest.TestCase):
    def test_counts_on_48_small_trees(self) -> None:
        graphs = [nx.empty_graph(1)]
        for n in range(2, 9):
            graphs.extend(nx.generators.nonisomorphic_trees(n))
        self.assertEqual(len(graphs), 48)
        for graph in graphs:
            with self.subTest(n=len(graph), edges=tuple(graph.edges())):
                poly, e = brute_counts(graph)
                self.assertEqual(polynomial(graph), poly)
                self.assertEqual(extendable_defects(graph), e)

    def test_forest_and_empty_cases(self) -> None:
        graphs = [nx.empty_graph(0), nx.empty_graph(5),
                  nx.disjoint_union(nx.path_graph(3), nx.star_graph(3))]
        for graph in graphs:
            with self.subTest(n=len(graph)):
                poly, e = brute_counts(graph)
                self.assertEqual(polynomial(graph), poly)
                self.assertEqual(extendable_defects(graph), e)

    def test_three_exact_witnesses(self) -> None:
        result = replay()
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["fixed_witnesses_recomputed"], 3)
        self.assertTrue(all(row["in_target_window"] for row in result["witnesses"]))
        self.assertTrue(all(row["full_margin"] > 0 for row in result["witnesses"]))

    def test_relabelling(self) -> None:
        graph = nx.path_graph(9)
        relabelled = nx.relabel_nodes(graph, {v: 20 - 2 * v for v in graph})
        left, right = analyze(graph), analyze(relabelled)
        for key in ("independence_polynomial", "e_0_to_4", "b_0_to_4", "full_margin"):
            self.assertEqual(left[key], right[key])

    def test_reject_non_tree(self) -> None:
        for graph in (nx.cycle_graph(8), nx.empty_graph(3), nx.empty_graph(0)):
            with self.assertRaises(ValueError):
                analyze(graph)

    def test_unimodality_plateaus(self) -> None:
        self.assertTrue(is_unimodal([1, 3, 3, 2, 2, 1]))
        self.assertFalse(is_unimodal([1, 3, 2, 2, 3]))
        self.assertFalse(is_unimodal([1, 4, 3, 5, 2]))


if __name__ == "__main__":
    unittest.main()
