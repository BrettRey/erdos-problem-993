#!/usr/bin/env python3
"""Independently verify the order-60 counterexample to EADS.

EADS is the proposed auxiliary statement that every tree T has a vertex v
for which the modal intervals of

    I(T-v)  and  x I(T-N[v])

are at distance at most one.  The tree below has distance exactly two at
every vertex.  This does *not* disprove independence-polynomial unimodality:
the full polynomial of this tree is unimodal (indeed, log-concave).

Two exact computations are compared:

1. the directed-message tree DP used by the EADS search; and
2. direct subset summation over the smaller side of each bipartite induced
   forest.  The latter has at most 11 enumerated vertices here, so it checks
   all relevant polynomials without using the tree recurrence.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from math import comb
from pathlib import Path

from indpoly import (
    independence_poly,
    is_log_concave,
    is_unimodal,
    log_concavity_ratio,
    near_miss_ratio,
)
from nm_optimizer import _validate_tree
from scratch_eads_local_constraint_breaker_20260808 import signed_orientation
from scratch_existential_split_mode_20260807 import (
    interval_distance,
    mode_interval,
    tree_vertex_splits,
)


EDGES = (
    (0, 1),
    (0, 14), (0, 15), (0, 16), (0, 17), (0, 18),
    (1, 2), (1, 3), (1, 4), (1, 5), (1, 7), (1, 9), (1, 10),
    (2, 19), (2, 20), (2, 21), (2, 22), (2, 23),
    (3, 24), (3, 25), (3, 26), (3, 27),
    (4, 12), (4, 28), (4, 29), (4, 30),
    (5, 6), (5, 31), (5, 32), (5, 33),
    (6, 8), (6, 11), (6, 13),
    (7, 34), (7, 35), (7, 36), (7, 37), (7, 38),
    (8, 39), (8, 40), (8, 41), (8, 42),
    (9, 43), (9, 44), (9, 45), (9, 46),
    (10, 47), (10, 48), (10, 49), (10, 50),
    (11, 51), (11, 52), (11, 53), (11, 54),
    (13, 55), (13, 56), (13, 57), (13, 58), (13, 59),
)


def adjacency(n: int = 60) -> list[list[int]]:
    adj: list[list[int]] = [[] for _ in range(n)]
    for u, v in EDGES:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def induced_adjacency(
    adj: list[list[int]], removed: set[int],
) -> list[list[int]]:
    kept = [vertex for vertex in range(len(adj)) if vertex not in removed]
    relabel = {old: new for new, old in enumerate(kept)}
    return [
        [relabel[v] for v in adj[u] if v not in removed]
        for u in kept
    ]


def bipartite_independence_poly(adj: list[list[int]]) -> list[int]:
    """Enumerate subsets of one bipartition side, component by component."""
    n = len(adj)
    color = [-1] * n
    chosen_side: list[int] = []
    for start in range(n):
        if color[start] >= 0:
            continue
        sides = [[], []]
        color[start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            sides[color[u]].append(u)
            for v in adj[u]:
                if color[v] < 0:
                    color[v] = 1 - color[u]
                    queue.append(v)
                elif color[v] == color[u]:
                    raise AssertionError("certificate graph is not bipartite")
        chosen_side.extend(min(sides, key=len))

    other_side = [v for v in range(n) if v not in set(chosen_side)]
    other_index = {vertex: bit for bit, vertex in enumerate(other_side)}
    neighbor_masks = []
    for vertex in chosen_side:
        mask = 0
        for neighbor in adj[vertex]:
            mask |= 1 << other_index[neighbor]
        neighbor_masks.append(mask)

    poly = [0] * (n + 1)
    for subset in range(1 << len(chosen_side)):
        selected = subset.bit_count()
        forbidden = 0
        bits = subset
        while bits:
            low_bit = bits & -bits
            index = low_bit.bit_length() - 1
            forbidden |= neighbor_masks[index]
            bits -= low_bit
        available = len(other_side) - forbidden.bit_count()
        for extra in range(available + 1):
            poly[selected + extra] += comb(available, extra)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_counterexample_n60_20260808.json"),
    )
    args = parser.parse_args()
    adj = adjacency()
    if not _validate_tree(len(adj), adj):
        raise AssertionError("hard-coded certificate is not a tree")

    full_message, splits = tree_vertex_splits(adj)
    full_standard = tuple(independence_poly(len(adj), adj))
    full_bipartite = tuple(bipartite_independence_poly(adj))
    if not full_message == full_standard == full_bipartite:
        raise AssertionError("independent full-polynomial computations disagree")

    records: list[dict[str, object]] = []
    orientation_counts: Counter[str] = Counter()
    distance_counts: Counter[int] = Counter()
    for vertex, split in enumerate(splits):
        removed_a = {vertex}
        removed_b = {vertex, *adj[vertex]}
        a_check = tuple(bipartite_independence_poly(
            induced_adjacency(adj, removed_a),
        ))
        closed_check = bipartite_independence_poly(
            induced_adjacency(adj, removed_b),
        )
        b_check = (0, *closed_check)
        if a_check != split.a_poly or b_check != split.b_poly:
            raise AssertionError(f"independent split check failed at {vertex}")
        a_modes = mode_interval(a_check)
        b_modes = mode_interval(b_check)
        distance = interval_distance(a_modes, b_modes)
        orientation, signed_distance = signed_orientation(split)
        if distance != split.distance or signed_distance != distance:
            raise AssertionError(f"mode check failed at {vertex}")
        orientation_counts[orientation] += 1
        distance_counts[distance] += 1
        records.append({
            "vertex": vertex,
            "degree": len(adj[vertex]),
            "orientation": orientation,
            "distance": distance,
            "a_modes": list(a_modes),
            "b_modes": list(b_modes),
            "a_poly": list(a_check),
            "b_poly": list(b_check),
        })

    minimum_distance = min(record["distance"] for record in records)
    if minimum_distance != 2:
        raise AssertionError(f"expected minimum distance 2, got {minimum_distance}")
    unimodal = is_unimodal(list(full_message))
    log_concave = is_log_concave(list(full_message))
    if not unimodal:
        raise AssertionError("unexpectedly found a main-conjecture counterexample")

    report = {
        "claim_scope": "exact counterexample to EADS; not to Erdos Problem 993",
        "statement_refuted": (
            "Every tree has a vertex v such that the modal intervals of "
            "I(T-v) and x I(T-N[v]) are at distance at most one."
        ),
        "order": len(adj),
        "edges": [list(edge) for edge in EDGES],
        "degrees": [len(row) for row in adj],
        "tree_validated": True,
        "verification_methods": [
            "directed-message tree DP",
            "iterative rooted tree DP",
            "direct subset summation over the smaller bipartition side",
        ],
        "full_polynomial": list(full_message),
        "full_polynomial_unimodal": unimodal,
        "full_polynomial_log_concave": log_concave,
        "worst_log_concavity_ratio": list(
            log_concavity_ratio(list(full_message))
        ),
        "near_miss_ratio": list(near_miss_ratio(list(full_message))),
        "minimum_split_distance": minimum_distance,
        "distance_counts": {
            str(key): value for key, value in sorted(distance_counts.items())
        },
        "orientation_counts": dict(sorted(orientation_counts.items())),
        "splits": records,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "event": "verified",
        "order": len(adj),
        "minimum_split_distance": minimum_distance,
        "distance_counts": report["distance_counts"],
        "orientation_counts": report["orientation_counts"],
        "full_polynomial_unimodal": unimodal,
        "full_polynomial_log_concave": log_concave,
        "near_miss_ratio": report["near_miss_ratio"],
        "out": str(args.out),
    }), flush=True)


if __name__ == "__main__":
    main()
