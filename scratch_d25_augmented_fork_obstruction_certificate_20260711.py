#!/usr/bin/env python3
"""Exact all-lambda obstruction to the D25 GSB-augmented fork family.

Join a new center to five copies of Galvin's rooted tree ``G(60,18)``.  At
the last required prefix rank, the center inequality fails already at
``lambda=0`` and its gap decreases with ``lambda``.  Hence every fixed
``lambda in [0,1]`` fails on this one tree/center/rank instance.

The calculation uses exact ``ZZ[x]`` arithmetic.  It independently compares
the closed-form branch polynomials with the repository's rooted tree DP and
checks that the full tree still satisfies GSB and is unimodal.  Compact
fingerprints avoid printing integers with thousands of digits.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from decimal import Decimal, localcontext
from fractions import Fraction

from sympy.polys.domains import ZZ
from sympy.polys.rings import ring

from indpoly import is_unimodal
from scratch_d25_hard_branch_composition_search_20260711 import make_state
from scripts.analyze_prufer_corpus import make_galvin_tree


M = 60
T = 18
BRANCHES = 5
RANK = 3799


def compact_integer(value: int) -> dict[str, int | str]:
    text = str(abs(value))
    return {
        "sign": (value > 0) - (value < 0),
        "digits": len(text),
        "sha256_20": hashlib.sha256(str(value).encode()).hexdigest()[:20],
    }


def decimal_fraction(value: Fraction) -> str:
    with localcontext() as context:
        context.prec = 40
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def connected(adj: list[list[int]]) -> bool:
    seen = {0}
    queue = [0]
    for vertex in queue:
        for neighbor in adj[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    return len(seen) == len(adj)


def main() -> None:
    sys.set_int_max_str_digits(30_000)

    branch_adj = make_galvin_tree(M, T)
    branch_edges = sum(map(len, branch_adj)) // 2
    assert len(branch_adj) == 1 + M + 2 * M * T == 2221
    assert branch_edges == len(branch_adj) - 1
    assert connected(branch_adj) and len(branch_adj[0]) == M

    state = make_state("G(60,18)", branch_adj, 0)
    polynomial_ring, x = ring("x", ZZ)

    def from_coefficients(coefficients: tuple[int, ...]):
        return polynomial_ring.from_dict(
            {
                (degree,): coefficient
                for degree, coefficient in enumerate(coefficients)
                if coefficient
            }
        )

    reserve = (1 + 2 * x) ** (M * T)
    excluded = ((1 + 2 * x) ** T + x * (1 + x) ** T) ** M
    total = excluded + x * reserve
    assert total == from_coefficients(state.total)
    assert excluded == from_coefficients(state.excluded)
    assert reserve == from_coefficients(state.reserve)

    whole = total**BRANCHES + x * excluded**BRANCHES
    alpha = whole.degree()
    limit = math.ceil((2 * alpha - 1) / 3)
    assert (alpha, limit, RANK) == (5701, 3801, limit - 2)

    coefficient = lambda poly, rank: int(poly[(rank,)])
    n_r = coefficient(whole, RANK)
    a_c = coefficient(excluded**BRANCHES, RANK)
    a = coefficient(reserve * total ** (BRANCHES - 1), RANK)
    joint = coefficient(
        reserve**2 * total ** (BRANCHES - 2), RANK
    )

    # There are C(5,2)=10 unordered forks, hence the ordered covariance term.
    q_ordered = 20 * (n_r * joint - a**2)

    # Each branch root has full-tree degree M+1=61.  Clearing that denominator
    # from RHS-Q gives gap(lambda)=gap_zero+lambda*slope.
    cleared_gap_zero = (
        BRANCHES * a**2
        + BRANCHES * (M + 1) * a_c * a
        + BRANCHES * n_r * a
        - (M + 1) * q_ordered
    )
    cleared_slope = n_r * ((M + 1) * a_c - BRANCHES * a)
    assert cleared_gap_zero < 0 and cleared_slope < 0
    equality_lambda = Fraction(-cleared_gap_zero, cleared_slope)
    assert equality_lambda < 0

    i_r1 = coefficient(whole, RANK + 1)
    i_r2 = coefficient(whole, RANK + 2)
    gsb_gap = (
        (RANK + 1) * i_r1**2
        + n_r * i_r1
        - (RANK + 2) * n_r * i_r2
    )
    assert gsb_gap > 0
    coefficients = [int(value) for value in reversed(whole.to_dense())]
    assert len(coefficients) == alpha + 1 and is_unimodal(coefficients)

    full_order = 1 + BRANCHES * len(branch_adj)
    full_edges = BRANCHES * (branch_edges + 1)
    assert (full_order, full_edges) == (11106, 11105)

    report = {
        "kind": "d25_augmented_fork_all_lambda_obstruction",
        "tree": {
            "description": "new center joined to five copies of G(60,18)",
            "n": full_order,
            "edges": full_edges,
            "branch_root_degree": M + 1,
        },
        "prefix": {"alpha": alpha, "limit": limit, "rank": RANK},
        "candidate": {
            "lambda_domain": "[0,1]",
            "cleared_gap_at_lambda_zero": compact_integer(cleared_gap_zero),
            "cleared_slope": compact_integer(cleared_slope),
            "equality_lambda": decimal_fraction(equality_lambda),
            "conclusion": "gap is negative at zero and strictly decreases",
        },
        "global_checks": {
            "gsb_gap": compact_integer(gsb_gap),
            "independence_sequence_unimodal": True,
        },
        "exact_inputs": {
            "i_r": compact_integer(n_r),
            "a_center": compact_integer(a_c),
            "a_branch_root": compact_integer(a),
            "joint_two_branch_roots": compact_integer(joint),
            "q_ordered": compact_integer(q_ordered),
        },
        "certificate": "passed",
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
