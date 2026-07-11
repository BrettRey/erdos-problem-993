#!/usr/bin/env python3
"""Exact certificate for the D23 recursive phase-grammar obstruction.

For a layer with b recursive children and rooted decoration D,

    E' = P^b P_D,               S' = x E^b E_D,

where E_D counts sets excluding the decoration root.  At x=1 put
rho=S/E, p=rho/(1+rho), delta=mu(S)-mu(E), and V=Var(P).  Direct logarithmic
differentiation gives the exact recurrences

    rho'   = w_D (1+rho)^(-b),             w_D=E_D(1)/P_D(1),
    delta' = c_D - b p delta,              c_D=1+mu(E_D)-mu(P_D),
    V'     >= b V / 2.

The last inequality holds because rho'<=1, so the excluded component has
mixture weight at least 1/2 and variance bV+Var(P_D).

Consequently, for periodic layer types with b>=3, any stable phase cycle has
bounded delta (at exact marginality, at most linear growth), whereas
sqrt(V) grows at least (b_min/2)^(h/2).  A periodic orbit with response
multiplier greater than one is repelling.  Thus the dichotomy is exact:
stable/marginal root-occupancy phase coexistence cannot produce cardinality
phases separated on the standard-deviation scale.

This script audits the identities in QQ against expanded ZZ[x] polynomials
and independent tree DP, and gives an exact rational trapping/contraction
certificate for the genuine near-critical 5-ary two-cycle.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from indpoly import independence_poly
from scratch_d23_phase_grammar_search_20260711 import (
    Moments,
    mixture,
    moments,
    product,
)
from scratch_d23_two_type_phase_search_20260711 import (
    GRAMMARS,
    LayerSpec,
    build_half,
    decoration_data,
    exact_states,
    order,
    specs,
)
from scratch_d17_checkpointed_search_20260711 import append_record, utc_now


def state_moments(layers, data):
    excluded = Moments(1, Fraction(0), Fraction(0))
    selected = Moments(1, Fraction(1), Fraction(0))
    audit = []
    for layer in layers:
        old_e = excluded
        old_s = selected
        total = mixture(old_e, old_s)
        decor = data[layer.decor]
        total_power = product(total, layer.branch)
        excluded_power = product(old_e, layer.branch)
        excluded = Moments(
            total_power.mass * decor.total_moments.mass,
            total_power.mean + decor.total_moments.mean,
            total_power.variance + decor.total_moments.variance,
        )
        selected = Moments(
            excluded_power.mass * decor.excluded_moments.mass,
            Fraction(1) + excluded_power.mean + decor.excluded_moments.mean,
            excluded_power.variance + decor.excluded_moments.variance,
        )
        rho = Fraction(selected.mass, excluded.mass)
        old_rho = Fraction(old_s.mass, old_e.mass)
        old_p = Fraction(old_rho, 1 + old_rho)
        w = Fraction(decor.excluded_moments.mass, decor.total_moments.mass)
        assert rho == w / (1 + old_rho) ** layer.branch
        delta = selected.mean - excluded.mean
        old_delta = old_s.mean - old_e.mean
        c = Fraction(1) + decor.excluded_moments.mean - decor.total_moments.mean
        assert delta == c - layer.branch * old_p * old_delta
        new_total = mixture(excluded, selected)
        assert new_total.variance >= Fraction(layer.branch, 2) * total.variance
        audit.append(
            {
                "branch": layer.branch,
                "decor": layer.decor,
                "selected_probability": float(Fraction(selected.mass, new_total.mass)),
                "mean_delta": float(delta),
                "variance": float(new_total.variance),
                "standardized_delta": float(delta) / math.sqrt(float(new_total.variance)),
            }
        )
    return excluded, selected, mixture(excluded, selected), audit


def adjacency_from_edges(edges, n):
    adjacency = [[] for _ in range(n)]
    for u, v in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    return adjacency


def critical_cycle_certificate(data):
    critical = data["critical_280_367"]
    w = Fraction(critical.excluded_moments.mass, critical.total_moments.mass)
    threshold = Fraction(5**5, 4**6)
    assert w == Fraction(280, 367)
    assert w - threshold == Fraction(5, 367 * 4096)

    def phase_map(r):
        return w / (1 + r) ** 5

    # Disjoint rational intervals around the two phase values.  They map into
    # one another exactly, and f^2 is a contraction on the first interval.
    first = (
        Fraction(248834950528794879, 10**18),
        Fraction(248834970528794879, 10**18),
    )
    second = (
        Fraction(251169389085038573, 10**18),
        Fraction(251169409197338573, 10**18),
    )
    assert first[1] < second[0]
    assert phase_map(first[1]) >= second[0]
    assert phase_map(first[0]) <= second[1]
    assert phase_map(second[1]) >= first[0]
    assert phase_map(second[0]) <= first[1]
    contraction_upper = (
        25
        * first[1]
        / (1 + first[0])
        * second[1]
        / (1 + second[1])
    )
    assert contraction_upper < 1

    # Interval evaluation also bounds both selected-state masses away from 0.
    selected_probability_interval = (
        first[0] / (1 + first[0]),
        second[1] / (1 + second[1]),
    )
    assert selected_probability_interval[0] > Fraction(19, 100)
    assert selected_probability_interval[1] < Fraction(21, 100)
    return {
        "activity": [w.numerator, w.denominator],
        "critical_activity": [threshold.numerator, threshold.denominator],
        "activity_excess": [5, 367 * 4096],
        "phase_interval_1": [str(first[0]), str(first[1])],
        "phase_interval_2": [str(second[0]), str(second[1])],
        "f2_contraction_upper": float(contraction_upper),
        "selected_probability_bounds": [
            float(selected_probability_interval[0]),
            float(selected_probability_interval[1]),
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/d23_phase_obstruction_certificate_20260711.jsonl"),
    )
    args = parser.parse_args()
    data = decoration_data()
    profiles = specs(7, data, 500000)

    # Exact recurrence and variance-bound audit on every frozen profile.
    maximum_standardized_delta = 0.0
    variance_steps = 0
    for _, layers in profiles:
        _, _, _, audit = state_moments(layers, data)
        variance_steps += len(audit)
        maximum_standardized_delta = max(
            maximum_standardized_delta,
            *(abs(row["standardized_delta"]) for row in audit),
        )

    # Expanded polynomial moment replay and independent tree-DP replay at
    # depths that remain deliberately small enough for a standalone audit.
    expanded_replays = 0
    tree_dp_replays = 0
    for grammar in GRAMMARS:
        for depth in (1, 2, 3):
            layers = tuple(
                grammar.period[k % len(grammar.period)] for k in range(depth)
            )
            if order(layers, data) > 800:
                continue
            e_poly, s_poly = exact_states(layers, data)
            e_m, s_m, _, _ = state_moments(layers, data)
            assert moments(e_poly) == e_m
            assert moments(s_poly) == s_m
            expanded_replays += 1

            edges, n, _ = build_half(layers, data)
            assert len(edges) == n - 1
            adjacency = adjacency_from_edges(edges, n)
            total_poly = [0] * max(len(e_poly), len(s_poly))
            for k, value in enumerate(e_poly):
                total_poly[k] += value
            for k, value in enumerate(s_poly):
                total_poly[k] += value
            assert independence_poly(n, adjacency) == total_poly
            tree_dp_replays += 1

    cycle = critical_cycle_certificate(data)
    theorem = {
        "mass_recurrence": "rho_h=w_j/(1+rho_(h-1))^b_j",
        "mean_recurrence": "delta_h=c_j-b_j p_(h-1) delta_(h-1)",
        "variance_bound": "V_h >= (b_j/2) V_(h-1)",
        "periodic_response": "delta_(h+L)=C+(-1)^L M delta_h, M=product_j b_j p_(j-1)",
        "stable_case": "M<1 implies delta=O(1)",
        "marginal_case": "M=1 implies delta=O(h)",
        "repelling_case": "M>1 is a repelling phase orbit",
        "standardized_conclusion": "for b_min>2, |delta_h|/sqrt(V_h)->0 on every stable or marginal periodic phase orbit",
    }
    result = {
        "at": utc_now(),
        "kind": "d23_phase_obstruction_certificate",
        "certificate": "passed",
        "profiles": len(profiles),
        "variance_steps_checked_exactly": variance_steps,
        "expanded_ZZx_replays": expanded_replays,
        "independent_tree_dp_replays": tree_dp_replays,
        "maximum_finite_standardized_delta": maximum_standardized_delta,
        "critical_phase_cycle": cycle,
        "theorem": theorem,
    }
    append_record(args.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
