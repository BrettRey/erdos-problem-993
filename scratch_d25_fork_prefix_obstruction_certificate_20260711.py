#!/usr/bin/env python3
"""Exact compact certificate for two prefix-local fork obstructions.

Let ``c`` be a new hub joined to the root ``p`` of ``G(40,20)`` and to
the roots ``u_1,...,u_48`` of 48 copies of ``TG(4,6)``.  At rank

    r = ceil((2 alpha(T)-1)/3) - 2 = 3279,

this script certifies failures of both local allocations considered in D25:

* the parent-oriented allocation; and
* the symmetric inverse-degree allocation.

All decisions are made with exact integers (or ``Fraction``).  The printed
record is kept compact by reporting decimal ratios, digit counts, and SHA-256
prefixes rather than thousands of decimal digits.  Re-running the file
reconstructs every exact integer from the two explicit rooted trees.
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

from scratch_d25_hard_branch_composition_search_20260711 import make_state
from scripts.analyze_prufer_corpus import (
    make_bautista_ramos_tree,
    make_galvin_tree,
)


sys.set_int_max_str_digits(20_000)


def attach(adj: list[list[int]], branch: list[list[int]]) -> None:
    """Attach a fresh copy of ``branch`` at its vertex 0 to hub 0."""
    offset = len(adj)
    adj.extend([] for _ in branch)
    for vertex, neighbors in enumerate(branch):
        for neighbor in neighbors:
            if vertex < neighbor:
                left, right = offset + vertex, offset + neighbor
                adj[left].append(right)
                adj[right].append(left)
    adj[0].append(offset)
    adj[offset].append(0)


def composed_tree(parent: list[list[int]], child: list[list[int]], m: int):
    adj: list[list[int]] = [[]]
    attach(adj, parent)
    for _ in range(m):
        attach(adj, child)
    return adj


def compact_integer(value: int) -> dict[str, int | str]:
    text = str(abs(value))
    return {
        "sign": (value > 0) - (value < 0),
        "digits": len(text),
        "sha256_20": hashlib.sha256(str(value).encode()).hexdigest()[:20],
    }


def decimal_ratio(numerator: int, denominator: int) -> str:
    with localcontext() as context:
        context.prec = 24
        return str(Decimal(numerator) / Decimal(denominator))


def main() -> None:
    multiplicity = 48
    parent_adj = make_galvin_tree(40, 20)
    child_adj = make_bautista_ramos_tree(4, 6)
    parent = make_state("G(40,20)", parent_adj, 0)
    child = make_state("TG(4,6)", child_adj, 0)

    # Audit that the advertised composition really is a tree of order 9418.
    adj = composed_tree(parent_adj, child_adj, multiplicity)
    edge_count = sum(map(len, adj)) // 2
    seen = {0}
    queue = [0]
    for vertex in queue:
        for neighbor in adj[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    assert len(adj) == 1 + len(parent_adj) + multiplicity * len(child_adj) == 9418
    assert len(seen) == len(adj) and edge_count == len(adj) - 1

    polynomial_ring, x = ring("x", ZZ)

    def as_poly(coefficients: tuple[int, ...]):
        return polynomial_ring.from_dict(
            {(degree,): coefficient for degree, coefficient in enumerate(coefficients)}
        )

    p0, e0, r0 = map(as_poly, (parent.total, parent.excluded, parent.reserve))
    p1, e1, r1 = map(as_poly, (child.total, child.excluded, child.reserve))
    assert p0 == e0 + x * r0
    assert p1 == e1 + x * r1

    # Reuse the adjacent powers needed by every marked coefficient.
    p46 = p1 ** (multiplicity - 2)
    p47 = p46 * p1
    p48 = p47 * p1
    e48 = e1**multiplicity
    whole = p0 * p48 + x * e0 * e48
    alpha = whole.degree()
    limit = math.ceil((2 * alpha - 1) / 3)
    rank = limit - 2
    assert (alpha, limit, rank) == (4921, 3281, 3279)

    coefficient = lambda poly: int(poly[(rank,)])
    n_r = coefficient(p0 * p48 + x * e0 * e48)
    a_c = coefficient(e0 * e48)
    a_p = coefficient(r0 * p48)
    a_u = coefficient(r1 * p0 * p47)
    j_pc = coefficient(r0 * r1 * p47)
    j_cc = coefficient(r1 * r1 * p0 * p46)
    j_ordered = (
        2 * multiplicity * j_pc
        + multiplicity * (multiplicity - 1) * j_cc
    )

    # Parent-oriented center allocation:
    # N J_ord <= C(C+2(a_p+a_c)), C=sum_i a_{u_i}=m a_u.
    child_sum = multiplicity * a_u
    oriented_lhs = n_r * j_ordered
    oriented_rhs = child_sum * (child_sum + 2 * (a_p + a_c))
    oriented_gap = oriented_rhs - oriented_lhs
    assert oriented_gap < 0

    # Symmetric inverse-degree center allocation:
    # q_ord <= sum_{v in N(c)} a_v^2/deg(v) + a_c sum_{v in N(c)}a_v.
    # Here deg(p)=41 and deg(u_i)=6 in the composed tree.
    product_ordered = (
        (a_p + child_sum) ** 2 - a_p**2 - multiplicity * a_u**2
    )
    q_ordered = n_r * j_ordered - product_ordered
    inverse_budget = (
        Fraction(a_p**2, 41)
        + Fraction(multiplicity * a_u**2, 6)
        + a_c * (a_p + child_sum)
    )
    inverse_gap = inverse_budget - q_ordered
    assert q_ordered > 0 and inverse_gap < 0

    # The actual ordered log-concavity target for the whole tree still has
    # positive slack.  This confirms that the two failures are failures of
    # the proposed local allocations, not counterexamples to the theorem.
    i_r1 = int(whole[(rank + 1,)])
    i_r2 = int(whole[(rank + 2,)])
    ordered_lc_lhs = (rank + 2) * n_r * i_r2
    ordered_lc_rhs = (rank + 1) * i_r1**2
    ordered_lc_gap = ordered_lc_rhs - ordered_lc_lhs
    assert ordered_lc_gap > 0
    gsb_rhs = ordered_lc_rhs + n_r * i_r1
    gsb_gap = gsb_rhs - ordered_lc_lhs
    assert gsb_gap > ordered_lc_gap > 0

    # Clear the rational denominator before compactly fingerprinting the gap.
    inverse_gap_numerator = inverse_gap.numerator
    report = {
        "tree": {
            "description": "hub + G(40,20) + 48 copies of TG(4,6)",
            "n": len(adj),
            "edges": edge_count,
            "parent_root_degree_in_T": len(parent_adj[0]) + 1,
            "child_root_degree_in_T": len(child_adj[0]) + 1,
        },
        "prefix": {"alpha": alpha, "limit": limit, "rank": rank},
        "oriented": {
            "lhs_over_rhs": decimal_ratio(oriented_lhs, oriented_rhs),
            "rhs_minus_lhs": compact_integer(oriented_gap),
        },
        "inverse_degree": {
            "q_ordered_over_budget": decimal_ratio(
                q_ordered * inverse_budget.denominator,
                inverse_budget.numerator,
            ),
            "budget_minus_q_denominator": inverse_gap.denominator,
            "budget_minus_q_numerator": compact_integer(inverse_gap_numerator),
        },
        "global_ordered_log_concavity": {
            "lhs_over_rhs": decimal_ratio(ordered_lc_lhs, ordered_lc_rhs),
            "rhs_minus_lhs": compact_integer(ordered_lc_gap),
        },
        "global_gsb": {
            "lhs_over_rhs": decimal_ratio(ordered_lc_lhs, gsb_rhs),
            "rhs_minus_lhs": compact_integer(gsb_gap),
        },
        "exact_inputs": {
            "i_r": compact_integer(n_r),
            "a_c": compact_integer(a_c),
            "a_p": compact_integer(a_p),
            "a_u": compact_integer(a_u),
            "j_pc": compact_integer(j_pc),
            "j_cc": compact_integer(j_cc),
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
