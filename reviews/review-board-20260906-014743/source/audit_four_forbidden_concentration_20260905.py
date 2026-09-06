#!/usr/bin/env python3
"""Exact finite certificate for the four-forbidden b1-zero tree family.

The mathematical domination proof is in the companion note. This script
independently enumerates its 1,842 rooted-forest representatives, checks their
polynomials by graph DP and subset enumeration, and tests the reductions.
It is not a verification of the full depth-three window or Erdős #993.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import product
from math import comb
from pathlib import Path
from time import perf_counter

import networkx as nx

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from indpoly import _polymul_python as multiply
from scripts.audit_b1_zero_forbidden_core_20260905 import (
    audit_graph,
    coefficient,
    extendable_defects,
    polynomial,
)
from scripts.audit_pendant_p2_boundary_20260904 import independent_masks

RootCode = tuple  # Sorted tuple of the recursively encoded children.
BOUND = Fraction(52513, 217404)
RESERVE = Fraction(7, 27)
CORE_HYPEREDGES = {
    "shared_four": ((0, 1, 2, 3),),
    "shared_three": ((0, 1, 2), (2, 3)),
    "pair_path": ((0, 1), (1, 2), (2, 3)),
    "pair_star": ((0, 1), (0, 2), (0, 3)),
}


def canonical_root(graph: nx.Graph, vertex: int, parent: int | None = None) -> RootCode:
    return tuple(sorted(canonical_root(graph, child, vertex)
                        for child in graph[vertex] if child != parent))


@cache
def rooted_forests(order: int) -> tuple[RootCode, ...]:
    """Enumerate unordered rooted forests by multisets, independently of nx."""
    if order == 0:
        return ((),)
    types = [(size, code) for size in range(1, order + 1)
             for code in rooted_forests(size - 1)]
    result = []

    def extend(remaining: int, first: int, chosen: tuple) -> None:
        if remaining == 0:
            result.append(tuple(sorted(chosen)))
            return
        for index in range(first, len(types)):
            size, code = types[index]
            if size > remaining:
                break
            extend(remaining - size, index, chosen + (code,))

    extend(order, 0, ())
    assert len(set(result)) == len(result)
    return tuple(sorted(result))


def encode(code: RootCode) -> str:
    return "(" + "".join(encode(child) for child in code) + ")"


def forest_from_code(code: RootCode) -> tuple[nx.Graph, list[int]]:
    graph = nx.Graph()

    def add(branch: RootCode, parent: int | None = None) -> int:
        vertex = len(graph)
        graph.add_node(vertex)
        if parent is not None:
            graph.add_edge(vertex, parent)
        for child in branch:
            add(child, vertex)
        return vertex

    roots = [add(branch) for branch in code]
    return graph, roots


def corona(base: nx.Graph) -> nx.Graph:
    graph = base.copy()
    order = len(base)
    graph.add_edges_from((vertex, vertex + order) for vertex in range(order))
    return graph


def corona_subset_polynomials(base: nx.Graph, roots: list[int]) -> tuple[list[int], list[int]]:
    """P and root-deleted Q by all subsets of the ten base vertices."""
    order = len(base)
    root_mask = sum(1 << vertex for vertex in roots)
    p = [0] * (order + 1)
    q = [0] * (order + 1)
    for independent in independent_masks(base):
        size = independent.bit_count()
        for leaves in range(order - size + 1):
            count = comb(order - size, leaves)
            p[size + leaves] += count
            if independent & root_mask == 0:
                q[size + leaves] += count
    return p, q


def plus(left: list[int], right: list[int], sign: int = 1) -> list[int]:
    return [coefficient(left, j) + sign * coefficient(right, j)
            for j in range(max(len(left), len(right)))]


def concentrated_polynomials(p: list[int], q: list[int]) -> tuple[list[int], list[int]]:
    a = [1, 2, 1]  # (1+x)^2
    g = [1, 3, 1]  # a+x
    g3 = multiply(multiply(g, g), g)
    a3 = multiply(multiply(a, a), a)
    e = multiply([comb(9, j) for j in range(10)], p)
    b = plus(multiply(multiply(a, p), plus(g3, a3, -1)),
             [0] + multiply(q, g3))
    return e, b


def core_graph(kind: str) -> nx.Graph:
    graph = nx.empty_graph(4)
    for hyperedge in CORE_HYPEREDGES[kind]:
        forced = len(graph)
        graph.add_edges_from((vertex, forced) for vertex in hyperedge)
    for vertex in range(4):
        while graph.degree(vertex) < 3:
            graph.add_edge(vertex, len(graph))
    assert len(graph) == 13 and nx.is_tree(graph)
    return graph


def decorated_tree(base: nx.Graph, roots: list[int], kind: str = "shared_four",
                   owners: tuple[int, ...] | None = None,
                   leaf_roots: tuple[int, ...] | None = None) -> nx.Graph:
    owners = owners if owners is not None else (0,) * len(roots)
    leaf_roots = leaf_roots if leaf_roots is not None else (0,) * len(roots)
    graph = nx.disjoint_union(core_graph(kind), corona(base))
    for root, owner, leaf in zip(roots, owners, leaf_roots):
        graph.add_edge(owner, 13 + root + leaf * len(base))
    assert nx.is_tree(graph)
    return graph


def profiles(e: list[int], b: list[int], alpha: int) -> tuple[list[int], list[int]]:
    return ([coefficient(e, alpha - d) for d in range(5)],
            [coefficient(b, alpha - d) for d in range(5)])


def test_decoration(base: nx.Graph, roots: list[int], e: list[int], upper_b: list[int],
                    kind: str, owners: tuple, leaves: tuple) -> None:
    graph = decorated_tree(base, roots, kind, owners, leaves)
    full = polynomial(graph)
    b = plus(full, e, -1)
    alpha = 9 + len(base)
    assert len(full) - 1 == alpha
    assert all(0 <= value <= coefficient(upper_b, j) for j, value in enumerate(b))
    ed, bd = profiles(e, b, alpha)
    assert bd[0] == bd[1] == 0
    assert 2 * bd[3] >= (12 + len(base)) * bd[2]
    if len(base) == 10:
        assert Fraction(bd[4], ed[4]) <= BOUND
        assert (ed[3] + bd[3]) ** 2 > (ed[2] + bd[2]) * (ed[4] + bd[4])


def run(random_cases: int, small_order: int) -> dict:
    start = perf_counter()
    counters = Counter()
    recursive = rooted_forests(10)
    via_unrooted = {canonical_root(graph, vertex)
                    for graph in nx.nonisomorphic_trees(11) for vertex in graph}
    assert set(recursive) == via_unrooted and len(recursive) == 1842
    rows, data = [], []
    worst = None
    for code in recursive:
        base, roots = forest_from_code(code)
        assert len(base) == 10
        p, q = corona_subset_polynomials(base, roots)
        k = corona(base)
        assert p == polynomial(k)
        assert q == polynomial(k.subgraph(set(k) - set(roots)))
        e, b = concentrated_polynomials(p, q)
        graph = decorated_tree(base, roots)
        assert plus(e, b) == polynomial(graph)
        ed, bd = profiles(e, b, 19)
        assert ed == extendable_defects(graph)
        assert bd[0] == bd[1] == 0
        assert bd[3] >= 11 * bd[2]
        assert BOUND.denominator * bd[4] <= BOUND.numerator * ed[4]
        assert 7 * ed[4] - 27 * bd[4] > 0
        assert (ed[3] + bd[3]) ** 2 > (ed[2] + bd[2]) * (ed[4] + bd[4])
        load = (Fraction(bd[2], ed[2]) + Fraction(bd[4], ed[4])
                + Fraction(bd[2] * bd[4], ed[2] * ed[4]))
        if load >= RESERVE:
            counters["stronger_joint_endpoint_bound_failures"] += 1
        ratio = Fraction(bd[4], ed[4])
        if worst is None or ratio > worst[0]:
            worst = (ratio, graph, code, ed, bd)
        rows.append([encode(code), *ed[2:5], *bd[2:5],
                     7 * ed[4] - 27 * bd[4]])
        data.append((base, roots, e, b))
        counters["rooted_forest_representatives"] += 1
        counters["subset_and_graph_dp_crosschecks"] += 1
        counters["matching_bag_profile_crosschecks"] += 1
    assert worst is not None and worst[0] == BOUND
    assert counters["stronger_joint_endpoint_bound_failures"] == 17
    for order in range(small_order + 1):
        for code in rooted_forests(order):
            base, roots = forest_from_code(code)
            p, q = corona_subset_polynomials(base, roots)
            e, b = concentrated_polynomials(p, q)
            for kind in CORE_HYPEREDGES:
                for owners in product(range(4), repeat=len(roots)):
                    for leaves in product(range(2), repeat=len(roots)):
                        test_decoration(base, roots, e, b, kind, owners, leaves)
                        counters["small_decoration_domination_checks"] += 1
    rng = random.Random(20260905)
    for _ in range(random_cases):
        base, roots, e, b = rng.choice(data)
        owners = tuple(rng.randrange(4) for _ in roots)
        leaves = tuple(rng.randrange(2) for _ in roots)
        kind = rng.choice(tuple(CORE_HYPEREDGES))
        test_decoration(base, roots, e, b, kind, owners, leaves)
        counters["order33_decoration_domination_checks"] += 1
    worst_profile = audit_graph(worst[1])
    worst_profile["refined_depth3_case_closed"] = True
    assert len(worst_profile["forbidden_vertices"]) == 4
    assert len(worst_profile["forced_vertices"]) == 9
    slack = RESERVE - BOUND
    assert slack > 0
    return {
        "kind": "four_forbidden_concentration_certificate",
        "certificate_date": "2026-09-05",
        "claim_status": "written domination proof plus exhaustive exact finite certificate; not independently reviewed, Lean-verified, or a solution of #993",
        "scope": {"n": 33, "alpha": 19, "deficiency": 5,
                  "forbidden_count": 4, "b1": 0, "flexible_matching_pairs": 10,
                  "small_decoration_max_pairs": small_order,
                  "random_decoration_cases": random_cases, "seed": 20260905},
        "counters": dict(counters), "bound_b4_over_e4": str(BOUND),
        "reserve": str(RESERVE), "reserve_slack": str(slack),
        "bound_failures": 0, "domination_failures": 0,
        "worst_rooted_forest_code": encode(worst[2]),
        "worst_representative": worst_profile,
        "certificate_columns": ["rooted_forest_code", "e2", "e3", "e4",
                                "b2", "b3", "b4", "7e4_minus_27b4"],
        "certificate_rows": rows,
        "elapsed_seconds": round(perf_counter() - start, 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-cases", type=int, default=1000)
    parser.add_argument("--small-order", type=int, default=3)
    parser.add_argument("--output", type=Path,
                        default=REPO / "results/four_forbidden_concentration_20260905.json")
    args = parser.parse_args()
    result = run(args.random_cases, args.small_order)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                     ("counters", "bound_b4_over_e4", "reserve_slack",
                      "bound_failures", "domination_failures", "elapsed_seconds")}, indent=2))


if __name__ == "__main__":
    main()
