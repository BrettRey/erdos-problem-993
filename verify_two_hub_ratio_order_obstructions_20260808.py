#!/usr/bin/env python3
"""Exact obstructions to two naive ratio-order proofs for two-hub trees.

The witnesses refute sufficient-condition routes only.  Every polynomial
certified here is log-concave and unimodal, so nothing here refutes Erdős
Problem #993 or the conjecture that all two-branch-vertex trees are
log-concave.
"""

from __future__ import annotations

import json
from math import comb
from pathlib import Path

from indpoly import independence_poly, is_log_concave, is_unimodal


def poly_add(left: list[int], right: list[int]) -> list[int]:
    """Add coefficient lists and remove trailing zero coefficients."""
    out = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def shift(poly: list[int]) -> list[int]:
    return [0] + poly


def trim(poly: list[int]) -> list[int]:
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def partial_ratio_failure(
    left: list[int], right: list[int],
) -> dict[str, int | str] | None:
    """Return an exact witness that ``left`` is not below ``right``."""
    left = trim(left)
    right = trim(right)
    degree = max(len(left), len(right)) - 1
    a = left + [0] * (degree + 1 - len(left))
    b = right + [0] * (degree + 1 - len(right))
    delta_left = next(k for k, value in enumerate(a) if value)
    degree_right = max(k for k, value in enumerate(b) if value)
    if delta_left > degree_right + 1:
        return {
            "type": "support",
            "delta_left": delta_left,
            "degree_right": degree_right,
        }
    for k in range(1, degree + 1):
        lhs = a[k] * b[k - 1]
        rhs = b[k] * a[k - 1]
        if lhs > rhs:
            return {"type": "ratio", "k": k, "lhs": lhs, "rhs": rhs}
    return None


def synchronization_failures(
    left: list[int], right: list[int],
) -> list[dict[str, object]]:
    failures = []
    first = partial_ratio_failure(left, shift(right))
    if first is not None:
        failures.append({"relation": "left <= x*right", "witness": first})
    second = partial_ratio_failure(right, shift(left))
    if second is not None:
        failures.append({"relation": "right <= x*left", "witness": second})
    return failures


def partial_synchronization_failure(
    left: list[int], right: list[int],
) -> dict[str, int] | None:
    degree = max(len(left), len(right)) - 1
    a = left + [0] * (degree + 3 - len(left))
    b = right + [0] * (degree + 3 - len(right))

    def coefficient(poly: list[int], index: int) -> int:
        return 0 if index < 0 else poly[index]

    for lower in range(degree + 2):
        for upper in range(lower, degree + 2):
            lhs = (
                coefficient(a, upper + 1) * coefficient(b, lower - 1)
                + coefficient(a, lower - 1) * coefficient(b, upper + 1)
            )
            rhs = (
                coefficient(a, upper) * coefficient(b, lower)
                + coefficient(a, lower) * coefficient(b, upper)
            )
            if lhs > rhs:
                return {
                    "lower": lower,
                    "upper": upper,
                    "lhs": lhs,
                    "rhs": rhs,
                }
    return None


def lc_gaps(poly: list[int]) -> list[int]:
    return [
        poly[k] * poly[k] - poly[k - 1] * poly[k + 1]
        for k in range(1, len(poly) - 1)
    ]


def double_star(a: int, b: int) -> list[list[int]]:
    adjacency = [[] for _ in range(a + b + 2)]

    def add_edge(u: int, v: int) -> None:
        adjacency[u].append(v)
        adjacency[v].append(u)

    add_edge(0, 1)
    for leaf in range(2, 2 + a):
        add_edge(0, leaf)
    for leaf in range(2 + a, 2 + a + b):
        add_edge(1, leaf)
    return adjacency


def contract_hub_edge(adjacency: list[list[int]]) -> list[list[int]]:
    """Contract hubs 0 and 1 in the explicit double-star labelling."""
    leaves = len(adjacency) - 2
    return [list(range(1, leaves + 1))] + [[0] for _ in range(leaves)]


def subdivide_hub_edge(adjacency: list[list[int]]) -> list[list[int]]:
    out = [list(neighbors) for neighbors in adjacency] + [[]]
    new_vertex = len(adjacency)
    out[0].remove(1)
    out[1].remove(0)
    out[0].append(new_vertex)
    out[1].append(new_vertex)
    out[new_vertex] = [0, 1]
    return out


def main() -> None:
    # First obstruction: the natural two-term decomposition of DS(3,3).
    a_term = [comb(6, k) for k in range(7)]
    d_term = [0] + [2 * comb(3, k) for k in range(4)]
    full_33 = poly_add(a_term, d_term)
    sync_33 = synchronization_failures(a_term, d_term)
    partial_sync_33 = partial_synchronization_failure(a_term, d_term)
    assert sync_33
    assert partial_sync_33 is not None
    assert is_log_concave(a_term) and is_log_concave(d_term)
    assert is_log_concave(full_33) and is_unimodal(full_33)

    # Second obstruction: subdivision-contraction summands for DS(2,2).
    tree_22 = double_star(2, 2)
    contraction = contract_hub_edge(tree_22)
    subdivision = subdivide_hub_edge(tree_22)
    poly_tree = independence_poly(len(tree_22), tree_22)
    poly_contraction = independence_poly(len(contraction), contraction)
    poly_subdivision = independence_poly(len(subdivision), subdivision)
    shifted_contraction = shift(poly_contraction)
    expected_subdivision = poly_add(poly_tree, shifted_contraction)
    subdivision_sync = synchronization_failures(poly_tree, shifted_contraction)
    assert poly_subdivision == expected_subdivision
    assert subdivision_sync
    assert all(
        is_log_concave(poly) and is_unimodal(poly)
        for poly in (poly_tree, poly_contraction, poly_subdivision)
    )

    report = {
        "claim_scope": (
            "Exact counterexamples to two sufficient-condition proof routes "
            "for two-hub trees; not counterexamples to log-concavity, "
            "unimodality, or Erdos Problem 993."
        ),
        "adjacent_double_star_decomposition": {
            "tree": "double star with three leaves at each hub",
            "order": 8,
            "A=(1+x)^6": a_term,
            "D=2x(1+x)^3": d_term,
            "I=A+D": full_33,
            "A_D_synchronization_failures": sync_33,
            "A_D_partial_synchronization_failure": partial_sync_33,
            "I_log_concavity_gaps": lc_gaps(full_33),
            "I_log_concave": True,
            "I_unimodal": True,
        },
        "subdivision_contraction_decomposition": {
            "tree": "double star with two leaves at each hub",
            "order": 6,
            "subdivided_edge": [0, 1],
            "I(T)": poly_tree,
            "I(T/e)": poly_contraction,
            "xI(T/e)": shifted_contraction,
            "I(T_e)": poly_subdivision,
            "identity_verified": poly_subdivision == expected_subdivision,
            "summand_synchronization_failures": subdivision_sync,
            "I(T_e)_log_concavity_gaps": lc_gaps(poly_subdivision),
            "all_three_polynomials_log_concave_and_unimodal": True,
        },
    }
    out = Path("results/two_hub_ratio_order_obstructions_20260808.json")
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "event": "verified",
        "out": str(out),
        "double_star_33_partial_sync_witness": partial_sync_33,
        "double_star_22_subdivision_identity": True,
    }))


if __name__ == "__main__":
    main()
