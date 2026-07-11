#!/usr/bin/env python3
"""Exact D20 state-grammar audit and leaf-smoothing obstruction.

The audit freezes three facts.

1. At their natural roots, Galvin states have no post-descent LC bump in E or
   S; their bump is created only in P=E+S.  Bautista--Ramos TG_(m,t) puts its
   m bumps in E (and P), while S has none.  A parent-selected state can thus
   inherit TG bumps but not Galvin bumps.

2. Direct hub leaves never occur in the selected state and impose the exact
   mass penalty 2^{-d}.  They can flatten the excluded branch only by making
   the correcting selected branch exponentially irrelevant.

3. Nested stars preserve selected mass, but force

       S_half=x (1+x)^q R(x),   q=b*s,

   for a hard-core polynomial R.  Exact regressions show the TG bumps being
   successively erased by this binomial convolution.  For any bounded finite
   hard library the fixed-core near-miss expansion applies uniformly, giving
   a positive eventual 1/q reserve.  Co-scaling the hard TG parameter is the
   sole pivot tested by the companion driver.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb

from indpoly import independence_poly
from scratch_d20_coscaling_search_20260711 import (
    CoCandidate,
    CoSpec,
    build_tree as build_coscale_tree,
    exact_half_states as exact_coscale_half_states,
)
from scratch_d20_state_grammar_search_20260711 import (
    Candidate,
    Feature,
    HalfSpec,
    build_tree,
    exact_half_states,
    poly_add,
    poly_mul,
    rooted_states,
    species_libraries,
)
from scripts.analyze_prufer_corpus import (
    make_bautista_ramos_tree,
    make_galvin_tree,
)
from sympy.polys.domains import ZZ
from sympy.polys.rings import ring


def post_descent_bumps(coefficients: list[int]) -> list[int]:
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    first = next(
        k for k in range(len(coefficients) - 1)
        if coefficients[k] > coefficients[k + 1]
    )
    return [
        k
        for k in range(max(1, first), len(coefficients) - 1)
        if coefficients[k - 1] * coefficients[k + 1] > coefficients[k] ** 2
    ]


def adjacency_from_edges(edges: list[tuple[int, int]], order: int) -> list[list[int]]:
    adjacency = [[] for _ in range(order)]
    for u, v in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    assert len(edges) == order - 1
    return adjacency


def main() -> None:
    tg_parameters = ((2, 5), (4, 6), (6, 7), (8, 8))
    galvin_parameters = ((3, 4), (6, 6), (14, 8), (21, 11))
    state_profiles = []

    for m, t in galvin_parameters:
        excluded, selected = rooted_states(make_galvin_tree(m, t), 0)
        total = poly_add(excluded, selected)
        profile = {
            "family": "Galvin",
            "parameters": [m, t],
            "excluded_bumps": len(post_descent_bumps(excluded.copy())),
            "selected_bumps": len(post_descent_bumps(selected.copy())),
            "total_bumps": len(post_descent_bumps(total.copy())),
        }
        assert profile == {
            **profile,
            "excluded_bumps": 0,
            "selected_bumps": 0,
            "total_bumps": 1,
        }
        state_profiles.append(profile)

    tg_excluded: dict[tuple[int, int], list[int]] = {}
    for m, t in tg_parameters:
        excluded, selected = rooted_states(make_bautista_ramos_tree(m, t), 0)
        total = poly_add(excluded, selected)
        tg_excluded[(m, t)] = excluded
        profile = {
            "family": "Bautista-Ramos",
            "parameters": [m, t],
            "excluded_bumps": len(post_descent_bumps(excluded.copy())),
            "selected_bumps": len(post_descent_bumps(selected.copy())),
            "total_bumps": len(post_descent_bumps(total.copy())),
        }
        assert profile["excluded_bumps"] == m
        assert profile["selected_bumps"] == 0
        assert profile["total_bumps"] == m
        state_profiles.append(profile)

    # Exact selected-mass penalty for one mixed direct/nested half.
    tg_library, galvin_library = species_libraries()
    spec = HalfSpec(
        tg_index=0,
        tg_count=2,
        galvin_index=0,
        galvin_count=1,
        nested_star_leaves=7,
        nested_star_count=2,
        direct_leaves=9,
    )
    polynomial_ring, x = ring("x", ZZ)
    half_e, half_s = exact_half_states(
        spec, tg_library, galvin_library, polynomial_ring, x
    )
    actual_mass_ratio = Fraction(int(half_s(1)), int(half_e(1)))
    tg_q = Fraction(
        sum(tg_library[0].excluded_exact),
        sum(tg_library[0].excluded_exact) + sum(tg_library[0].selected_exact),
    )
    galvin_q = Fraction(
        sum(galvin_library[0].excluded_exact),
        sum(galvin_library[0].excluded_exact)
        + sum(galvin_library[0].selected_exact),
    )
    star_q = Fraction(2**spec.nested_star_leaves, 2**spec.nested_star_leaves + 1)
    predicted_mass_ratio = (
        Fraction(1, 2**spec.direct_leaves)
        * tg_q**spec.tg_count
        * galvin_q**spec.galvin_count
        * star_q**spec.nested_star_count
    )
    assert actual_mass_ratio == predicted_mass_ratio
    assert actual_mass_ratio <= Fraction(1, 2**spec.direct_leaves)

    # Exact formula versus independent tree DP for the mixed grammar.
    dummy_feature = Feature(None, None, None, 0.0, 0.0, None)
    other = HalfSpec(1, 1, 1, 0, 5, 1, 2)
    candidate = Candidate(
        0.0,
        spec,
        other,
        0,
        0,
        0.0,
        0,
        0.0,
        dummy_feature,
        dummy_feature,
        0.0,
        0.0,
    )
    left_e, left_s = exact_half_states(
        spec, tg_library, galvin_library, polynomial_ring, x
    )
    right_e, right_s = exact_half_states(
        other, tg_library, galvin_library, polynomial_ring, x
    )
    formula = left_e * right_e + left_s * right_e + left_e * right_s
    formula_coefficients = [int(value) for value in reversed(formula.to_dense())]
    edges, order = build_tree(candidate, tg_library, galvin_library)
    assert formula_coefficients == independence_poly(
        order, adjacency_from_edges(edges, order)
    )

    # Exact bump-erasure regressions for S=x E_TG (1+x)^q.
    smoothing_cases = ((2, 5, 256), (4, 6, 1024), (6, 7, 4096))
    smoothing = []
    for m, t, q in smoothing_cases:
        before = len(post_descent_bumps(tg_excluded[(m, t)].copy()))
        smoothed = poly_mul(
            tg_excluded[(m, t)], [comb(q, k) for k in range(q + 1)]
        )
        after = len(post_descent_bumps(smoothed))
        assert before == m
        assert after == 0
        smoothing.append(
            {"parameters": [m, t], "binomial_exponent": q, "before": before, "after": after}
        )

    # Co-scaling formula versus independent tree DP on a small heterogeneous pair.
    co_left = CoSpec(2, 3, 1, 4, 1)
    co_right = CoSpec(3, 4, 1, 5, 1)
    co_candidate = CoCandidate(
        0.0,
        co_left,
        co_right,
        0,
        0,
        0.0,
        0,
        0.0,
        dummy_feature,
        dummy_feature,
        0.0,
        0.0,
    )
    co_left_e, co_left_s = exact_coscale_half_states(
        co_left, polynomial_ring, x
    )
    co_right_e, co_right_s = exact_coscale_half_states(
        co_right, polynomial_ring, x
    )
    co_formula = (
        co_left_e * co_right_e
        + co_left_s * co_right_e
        + co_left_e * co_right_s
    )
    co_coefficients = [int(value) for value in reversed(co_formula.to_dense())]
    co_edges, co_order = build_coscale_tree(co_candidate)
    assert co_coefficients == independence_poly(
        co_order, adjacency_from_edges(co_edges, co_order)
    )

    result = {
        "state_profiles": state_profiles,
        "direct_leaf_mass_identity": "S(1)/E(1)=2^(-d) q_TG^a q_G^g (2^s/(2^s+1))^b",
        "direct_leaf_mass_ratio_exact": [
            actual_mass_ratio.numerator,
            actual_mass_ratio.denominator,
        ],
        "nested_selected_state": "S=x E_TG^a E_G^g (1+x)^(bs)",
        "smoothing_regressions": smoothing,
        "mixed_formula_tree_dp_order": order,
        "coscaling_formula_tree_dp_order": co_order,
        "certificate": "passed",
    }
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
