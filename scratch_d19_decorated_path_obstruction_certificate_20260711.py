#!/usr/bin/env python3
"""Exact audit and compression obstruction for D19 decorated paths.

For an aggregated prefix rooted at its last block, let (A_j,B_j) be the
polynomials conditional on that boundary root being excluded/selected and
F_j=A_j+B_j.  Adding the next rooted block (E,S) gives

    F_(j+1) = A_j(E+S)+B_j E = F_j(E+S)-B_j S.

Thus every longer decorated path is again a single-edge product-minus-corner
bridge.  The several negative terms visible after full inclusion--exclusion
expansion are not independent perturbations.  Moreover deletion of a selected
boundary root injects B_{j,k} into A_{j,k-1}; applying this on both sides gives

    [x^k](B_j S) <= [x^(k-2)](A_j E).

This script exactly replays the recurrence against tree DP, checks the
compression identity in ZZ[x], and records hard spherical regression rows.
"""

from __future__ import annotations

import json
import random

from sympy.polys.domains import ZZ
from sympy.polys.rings import ring

from indpoly import independence_poly
from scratch_d17_checkpointed_search_20260711 import exact_valley
from scratch_d17_decorated_path_search_20260711 import (
    build_decorated_path,
    exact_path_polynomial,
    treehood_certificate,
)


def exact_states(polynomial_ring, x, branching: tuple[int, ...]):
    excluded = polynomial_ring.one
    selected = x
    for branch in branching:
        old_excluded = excluded
        excluded = (excluded + selected) ** branch
        selected = x * old_excluded**branch
    return excluded, selected


def post_descent_lc_bumps(coefficients: list[int]) -> list[int]:
    first = next(
        k for k in range(len(coefficients) - 1)
        if coefficients[k] > coefficients[k + 1]
    )
    return [
        k
        for k in range(max(1, first), len(coefficients) - 1)
        if coefficients[k - 1] * coefficients[k + 1] > coefficients[k] ** 2
    ]


def adjacency_from_blocks(blocks: tuple[tuple[int, ...], ...]) -> list[list[int]]:
    edges, order = build_decorated_path(blocks)
    adjacency = [[] for _ in range(order)]
    for u, v in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    return adjacency


def main() -> None:
    polynomial_ring, x = ring("x", ZZ)
    hard_rows = (
        (1, 4, 3),   # Galvin T_(3,4,1), 28 vertices.
        (1, 6, 6),   # Galvin prefix-GSB stress row, 79 vertices.
        (1, 8, 14),  # Galvin global-GSB stress row, 239 vertices.
        (1, 11, 21), # Galvin near-miss row, 484 vertices.
    )

    # Verify the exact prefix compression identity for heterogeneous blocks.
    compression_blocks = (
        (1, 4, 3),
        (2, 3, 2, 1),
        (1, 6, 6),
        (3, 2, 1, 2),
        (1, 8, 14),
        (2, 2, 3),
    )
    prefix_out = None
    prefix_on = None
    prefix_total = None
    compression_steps = 0
    for index, block in enumerate(compression_blocks):
        excluded, selected = exact_states(polynomial_ring, x, block)
        if index == 0:
            current_out = excluded
            current_on = selected
        else:
            direct_out = excluded * prefix_total
            direct_on = selected * prefix_out
            direct_total = direct_out + direct_on
            compressed_total = prefix_total * (excluded + selected) - prefix_on * selected
            assert direct_total == compressed_total
            current_out, current_on = direct_out, direct_on
            compression_steps += 1
        prefix_out = current_out
        prefix_on = current_on
        prefix_total = current_out + current_on
    compressed_coefficients = [int(value) for value in reversed(prefix_total.to_dense())]
    assert compressed_coefficients == exact_path_polynomial(compression_blocks)

    # Independent tree-DP replay on 32 deterministic m=3..6 paths.
    rng = random.Random(20_260_711)
    tree_dp_cases = 0
    for block_count in range(3, 7):
        for _ in range(8):
            blocks = tuple(
                tuple(rng.randint(1, 3) for _ in range(rng.randint(2, 4)))
                for _ in range(block_count)
            )
            adjacency = adjacency_from_blocks(blocks)
            exact = exact_path_polynomial(blocks)
            assert exact == independence_poly(len(adjacency), adjacency)
            certificate = treehood_certificate(blocks)
            assert certificate["edges"] == certificate["vertices"] - 1
            tree_dp_cases += 1

    # Exact hard-family regression: repeated Galvin blocks are smoothed, not
    # turned into valleys, by path joins through m=6.
    hard_regressions = []
    for row in hard_rows:
        for block_count in range(3, 7):
            blocks = (row,) * block_count
            coefficients = exact_path_polynomial(blocks)
            certificate = treehood_certificate(blocks)
            valley = exact_valley(coefficients)
            bumps = post_descent_lc_bumps(coefficients)
            assert valley is None
            assert not bumps
            hard_regressions.append(
                {
                    "branching_leaves_to_root": row,
                    "blocks": block_count,
                    "vertices": certificate["vertices"],
                    "degree": len(coefficients) - 1,
                    "post_descent_lc_bumps": len(bumps),
                    "valley": valley,
                }
            )

    result = {
        "compression_identity": "F_next=F_prefix*(E+S)-S_prefix*S",
        "coefficient_injection": "[x^k](S_prefix*S) <= [x^(k-2)](E_prefix*E)",
        "compression_steps_checked_in_ZZx": compression_steps,
        "tree_dp_replay_cases": tree_dp_cases,
        "hard_regressions": hard_regressions,
        "certificate": "passed",
    }
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
