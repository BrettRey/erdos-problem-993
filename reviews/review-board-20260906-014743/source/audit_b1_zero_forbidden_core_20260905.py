#!/usr/bin/env python3
"""Exact audit of the b1=0 forbidden-core reduction; not a proof of #993."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb
from pathlib import Path
from time import perf_counter

import networkx as nx

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from indpoly import independence_poly
from scripts.audit_pendant_p2_boundary_20260904 import independent_masks
from scripts.probe_blocked_profile_depth3_20260828 import matching_bags


def polynomial(graph: nx.Graph) -> list[int]:
    graph = nx.convert_node_labels_to_integers(graph)
    return independence_poly(len(graph), [list(graph.neighbors(v)) for v in graph])


def coefficient(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def choose(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def density_bounds(alpha: int, delta: int, r: int) -> dict:
    """Exact cone-ray certificate for the b1=0 endpoint bounds (m >= 4)."""
    c = r + delta
    m = alpha - c
    assert 1 <= r <= delta - 1 and m >= 4
    beta = [Fraction(comb(m, j), 2 ** j) for j in range(5)]
    numerator = [r * choose(c - 3, 2) + choose(r, 2) * (c - 5) + choose(r, 3),
                 r * (c - 3) + choose(r, 2), r, 0, 0]
    denominator = [choose(c, 4 - j) for j in range(5)]
    ray_numerators = [sum(numerator[j] * beta[j] for j in range(k, 5))
                      for k in range(5)]
    ray_denominators = [sum(denominator[j] * beta[j] for j in range(k, 5))
                        for k in range(5)]
    ratios = [a / b for a, b in zip(ray_numerators, ray_denominators)]
    v = max(ratios)
    u = r / (beta[2] + c * beta[1] + choose(c, 2))
    load = u + v + u * v
    reserve = Fraction(5 * alpha + 17, 27 * (alpha - 3))
    slacks = [v * b - a for a, b in zip(ray_numerators, ray_denominators)]
    assert min(slacks) >= 0
    return {"n": 2 * alpha - delta, "alpha": alpha, "deficiency": delta,
            "forbidden_count": r, "forced_count": c, "flexible_pairs": m,
            "beta": list(map(str, beta)), "b4_coefficient_upper_bound": numerator,
            "e4_coefficients": denominator, "cone_ray_ratios": list(map(str, ratios)),
            "cone_ray_slacks": list(map(str, slacks)),
            "x2_upper": str(u), "x4_upper": str(v), "joint_load_upper": str(load),
            "reserve": str(reserve), "reserve_slack": str(reserve - load),
            "load_over_reserve": str(load / reserve), "closes_case": load < reserve}


def density_certificate() -> dict:
    rows = [density_bounds(alpha, delta, r)
            for alpha in (17, 18, 19) for delta in range(6)
            if 33 <= 2 * alpha - delta <= 38
            for r in range(1, delta)]
    closed = [row for row in rows if row["closes_case"]]
    residual = [row for row in rows if not row["closes_case"]]
    assert len(rows) == 13 and len(closed) == 12 and len(residual) == 1
    assert all(row["forbidden_count"] <= 3 for row in closed)
    assert (residual[0]["alpha"], residual[0]["deficiency"],
            residual[0]["forbidden_count"]) == (19, 5, 4)
    worst = max(closed, key=lambda row: Fraction(row["joint_load_upper"]))
    assert Fraction(worst["joint_load_upper"]) == Fraction(170838, 691067)
    assert Fraction(worst["joint_load_upper"]) < Fraction(7, 27)
    return {"parameter_cases": rows, "closed_cases": len(closed),
            "worst_closed_case": worst,
            "residual_relaxation_failures": residual,
            "warning": "A failed upper-bound relaxation is not a graph counterexample."}


def multiply(left: list[int], right: list[int], depth: int) -> list[int]:
    return [sum(left[j] * right[k - j] for j in range(k + 1))
            for k in range(depth + 1)]


def extendable_defects(graph: nx.Graph, depth: int = 4) -> list[int]:
    """Count extendable sets by deterministic feasible-domain DP on matching bags.

    A partial bag assignment has a unique set of possible root endpoints in
    its completions. Group by that set, not by the number of completions, so
    each partial independent set is counted once. Degrees count empty bags.
    """
    bags = matching_bags(graph)
    quotient = nx.Graph()
    quotient.add_nodes_from(range(len(bags)))
    location = {v: (i, j) for i, bag in enumerate(bags) for j, v in enumerate(bag)}
    ports = {}
    for u, v in graph.edges():
        i, p = location[u]
        j, q = location[v]
        if i != j:
            quotient.add_edge(i, j)
            ports[i, j] = (1 << p, 1 << q)
            ports[j, i] = (1 << q, 1 << p)
    assert len(bags) == len(polynomial(graph)) - 1
    total = [1] + [0] * depth
    for component in nx.connected_components(quotient):
        root = min(component)
        parent = {root: None}
        order = [root]
        for i in order:
            for j in quotient.neighbors(i):
                if j != parent[i]:
                    assert j not in parent
                    parent[j] = i
                    order.append(j)
        tables = {}
        for i in reversed(order):
            full = (1 << len(bags[i])) - 1
            states = {}
            for j in range(len(bags[i])):
                states[1 << j] = [1] + [0] * depth
            if depth:
                states.setdefault(full, [0] * (depth + 1))[1] += 1
            for child in quotient.neighbors(i):
                if parent.get(child) != i:
                    continue
                parent_port, child_port = ports[i, child]
                merged = {}
                for mask, values in states.items():
                    for child_mask, child_values in tables[child].items():
                        allowed = full ^ parent_port if child_mask == child_port else full
                        new_mask = mask & allowed
                        if new_mask:
                            contribution = multiply(values, child_values, depth)
                            target = merged.setdefault(new_mask, [0] * (depth + 1))
                            for d in range(depth + 1):
                                target[d] += contribution[d]
                states = merged
            tables[i] = states
        component_values = [sum(row[d] for row in tables[root].values())
                            for d in range(depth + 1)]
        total = multiply(total, component_values, depth)
    return total


def roles(graph: nx.Graph, alpha: int) -> tuple[set[int], set[int]]:
    allowed, forced = set(), set()
    for v in graph:
        closed = {v, *graph.neighbors(v)}
        if len(polynomial(graph.subgraph(set(graph) - closed))) == alpha:
            allowed.add(v)
        if len(polynomial(graph.subgraph(set(graph) - {v}))) - 1 < alpha:
            forced.add(v)
    return allowed, forced


def audit_graph(graph: nx.Graph, brute: bool = False) -> dict:
    poly = polynomial(graph)
    alpha = len(poly) - 1
    e = extendable_defects(graph)
    b = [coefficient(poly, alpha - d) - e[d] for d in range(5)]
    assert min(b) >= 0
    if alpha >= 4:
        assert 6 * b[3] >= (alpha - 4) * b[2]
        assert (alpha - 3) * e[3] >= 4 * e[4]
    if brute:
        sets = independent_masks(graph)
        maxima = [s for s in sets if s.bit_count() == alpha]
        brute_e = [sum(s.bit_count() == alpha - d and
                       any(s & ~m == 0 for m in maxima) for s in sets)
                   for d in range(5)]
        assert e == brute_e, (list(graph.edges()), e, brute_e)
        assert poly == [sum(s.bit_count() == k for s in sets) for k in range(alpha + 1)]
    delta = 2 * alpha - len(graph)
    row = {"graph6": nx.to_graph6_bytes(graph, header=False).decode().strip(),
           "n": len(graph), "alpha": alpha, "deficiency": delta,
           "e_0_to_4": e, "b_0_to_4": b}
    correction = b[3] ** 2 - b[2] * b[4] + 2 * e[3] * b[3] - e[2] * b[4] - e[4] * b[2]
    row["combined_correction"] = correction
    row["joint_margin_numerator"] = ((5 * alpha + 17) * e[2] * e[4]
        - 27 * (alpha - 3) * (e[2] * b[4] + e[4] * b[2] + b[2] * b[4]))
    if b[1]:
        return row
    allowed, forced = roles(graph, alpha)
    forbidden = sorted(set(graph) - allowed)
    core = graph.subgraph(allowed)
    core_poly = polynomial(core)
    assert e == [coefficient(core_poly, alpha - d) for d in range(5)]
    assert len(forced) == len(forbidden) + delta
    assert all(core.degree(v) == 0 for v in forced)
    assert all(len(set(graph.neighbors(v)) & forced) >= 3 for v in forbidden)
    assert not forbidden or len(forbidden) <= delta - 1
    for bag in matching_bags(graph):
        if len(bag) == 2 and set(bag) <= allowed:
            assert min(core.degree(v) for v in bag) == 1
    terms = []
    earlier = set()
    for v in forbidden:
        deleted = {v, *graph.neighbors(v)} | earlier
        residual = polynomial(graph.subgraph(set(graph) - deleted))
        terms.append([coefficient(residual, alpha - d - 1) for d in range(5)])
        earlier.add(v)
    assert b == [sum(term[d] for term in terms) for d in range(5)]
    row.update(forbidden_vertices=forbidden, forced_vertices=sorted(forced),
               ordered_blocked_terms=terms)
    if forbidden:
        r, c = len(forbidden), len(forced)
        flexible = graph.subgraph(allowed - forced)
        flexible_poly = polynomial(flexible)
        m = len(flexible_poly) - 1
        assert len(flexible) == 2 * m and m == alpha - c
        p = [coefficient(flexible_poly, m - j) for j in range(5)]
        upper_b2 = r * p[0]
        upper_b4 = (r * (p[2] + (c - 3) * p[1] + choose(c - 3, 2) * p[0])
                    + choose(r, 2) * (p[1] + (c - 5) * p[0])
                    + choose(r, 3) * p[0])
        assert b[2] <= upper_b2 and b[4] <= upper_b4
        for j in range(min(4, m)):
            assert (m - j) * p[j] <= 2 * (j + 1) * p[j + 1]
        row.update(b2_coefficient_bound=upper_b2, b4_coefficient_bound=upper_b4)
        if alpha in (17, 18, 19) and delta <= 5 and 33 <= len(graph) <= 38:
            bounds = density_bounds(alpha, delta, r)
            assert Fraction(b[2], e[2]) <= Fraction(bounds["x2_upper"])
            assert Fraction(b[4], e[4]) <= Fraction(bounds["x4_upper"])
            if bounds["closes_case"]:
                assert row["joint_margin_numerator"] > 0
            row["live_density_case_closed"] = bounds["closes_case"]
            if r == 4 and nx.is_connected(graph):
                rigid_core = graph.subgraph(set(forbidden) | forced)
                assert nx.is_tree(rigid_core) and len(rigid_core) == 13
                assert graph.subgraph(forbidden).number_of_edges() == 0
                assert all(len(set(graph.neighbors(v)) & forced) == 3 for v in forbidden)
                for component in nx.connected_components(flexible):
                    attachments = [(v, w) for v in component for w in graph.neighbors(v)
                                   if w not in component]
                    assert len(attachments) == 1 and attachments[0][1] in forbidden
    return row


def sharp_tree(r: int, alpha: int = 19) -> nx.Graph:
    """r forbidden degree-three vertices, joined by forced vertices, plus P2s."""
    graph = nx.empty_graph(r)
    for v in range(r - 1):
        connector = len(graph)
        graph.add_edges_from([(v, connector), (connector, v + 1)])
    for v in range(r):
        while graph.degree(v) < 3:
            graph.add_edge(v, len(graph))
    for j in range(alpha - (2 * r + 1)):
        inner = len(graph)
        graph.add_edges_from([(j % r, inner), (inner, inner + 1)])
    assert nx.is_tree(graph)
    return graph


def archived_rows(obj: object):
    if isinstance(obj, dict):
        if all(k in obj for k in ("graph6", "e_0_to_4", "b_0_to_4")):
            yield obj
        for value in obj.values():
            yield from archived_rows(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from archived_rows(value)


def run(max_n: int) -> dict:
    start = perf_counter()
    counters = Counter()
    by_delta = Counter()
    for n in range(2, max_n + 1):
        for graph in nx.nonisomorphic_trees(n):
            row = audit_graph(graph, brute=True)
            counters["small_trees_bruteforce_crosschecked"] += 1
            if row["alpha"] >= 4:
                counters["small_trees_general_shadow_bounds_checked"] += 1
            if row["b_0_to_4"][1] == 0:
                counters["small_b1_zero_trees"] += 1
                by_delta[row["deficiency"]] += 1
    components = [nx.empty_graph(1)] + [graph for n in range(2, 7)
                                       for graph in nx.nonisomorphic_trees(n)]
    for left, right in combinations_with_replacement(components, 2):
        row = audit_graph(nx.disjoint_union(left, right), brute=True)
        counters["disjoint_forests_bruteforce_crosschecked"] += 1
        if row["b_0_to_4"][1] == 0:
            counters["disjoint_b1_zero_forests"] += 1
    archived = {}
    sources = ["cross_reserve_witness_lifts_20260904.json",
               "cross_reserve_multicell_lifts_20260904.json",
               "negative_correction_seeds_n18_20260904.json"]
    for name in sources:
        for row in archived_rows(json.loads((REPO / "results" / name).read_text())):
            archived[row["graph6"]] = row
    rejected_split = None
    for old in archived.values():
        graph = nx.from_graph6_bytes(old["graph6"].encode())
        row = audit_graph(graph)
        assert row["e_0_to_4"] == old["e_0_to_4"]
        assert row["b_0_to_4"] == old["b_0_to_4"]
        counters["archived_profiles_independently_replayed"] += 1
        if row["b_0_to_4"][1] == 0:
            counters["archived_b1_zero_trees"] += 1
            if row["alpha"] >= 17 and row["joint_margin_numerator"] < 0:
                counters["b1_zero_joint_bound_failures"] += 1
        elif row["alpha"] >= 5 and row["combined_correction"] < 0:
            counters["b1_positive_nonadversity_counterexamples"] += 1
            if rejected_split is None:
                rejected_split = row
    sharp = [audit_graph(sharp_tree(r)) for r in range(1, 5)]
    for r, row in enumerate(sharp, 1):
        assert row["b_0_to_4"][1] == 0 and len(row["forbidden_vertices"]) == r
        assert row["deficiency"] == r + 1 and row["alpha"] == 19
    return {"kind": "b1_zero_forbidden_core_audit", "certificate_date": "2026-09-05",
            "claim_status": "exact audit accompanying a written structural proof; not Lean-verified or a solution of #993",
            "scope": {"small_tree_orders": [2, max_n], "archive_sources": sources,
                      "arithmetic": "exact integers; independent-set brute force and feasible-domain matching-bag DP"},
            "counters": dict(counters), "small_eligible_by_deficiency": dict(sorted(by_delta.items())),
            "structural_failures": 0, "sharp_forbidden_vertex_examples": sharp,
            "density_certificate": density_certificate(),
            "rejected_b1_positive_nonadversity_split": rejected_split,
            "elapsed_seconds": round(perf_counter() - start, 3)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=14)
    parser.add_argument("--output", type=Path,
                        default=REPO / "results/b1_zero_forbidden_core_20260905.json")
    args = parser.parse_args()
    result = run(args.max_n)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in ("counters", "structural_failures", "elapsed_seconds")}, indent=2))


if __name__ == "__main__":
    main()
