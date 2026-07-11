#!/usr/bin/env python3
"""Tilted scan of GSB-augmented fork lambda constraints on hard hubs."""

from __future__ import annotations

import math

from scratch_d17_spherical_search_20260711 import add, multiply, power
from scratch_d25_hard_branch_composition_search_20260711 import make_state
from scratch_d25_two_species_asymptotic_20260711 import (
    centering_log_z,
    shift,
    tilted,
    weighted,
)
from scripts.analyze_prufer_corpus import (
    make_bautista_ramos_tree,
    make_galvin_tree,
)


def constraint(parent, child, multiplicity, rank, log_z, parent_degree, child_degree):
    p0, e0, r0 = (tilted(poly, log_z) for poly in
                  (parent.total, parent.excluded, parent.reserve))
    p1, e1, r1 = (tilted(poly, log_z) for poly in
                  (child.total, child.excluded, child.reserve))
    p_m2 = power(p1, multiplicity - 2)
    p_m1 = multiply(p_m2, p1)
    p_m = multiply(p_m1, p1)
    e_m = power(e1, multiplicity)
    whole = add(multiply(p0, p_m), shift(multiply(e0, e_m), log_z))
    a_c = multiply(e0, e_m)
    a_p = multiply(r0, p_m)
    a_u = multiply(multiply(r1, p0), p_m1)
    j_pc = multiply(multiply(r0, r1), p_m1)
    j_cc = multiply(multiply(multiply(r1, r1), p0), p_m2)
    joint = add(
        weighted(j_pc, 2 * multiplicity),
        weighted(j_cc, multiplicity * (multiplicity - 1)),
    )
    items = (whole, joint, a_c, a_p, a_u)
    if any(rank >= len(item.values) or item.values[rank] <= 0 for item in items):
        return None
    log_n, log_j, log_ac, log_ap, log_au = [
        math.log(item.values[rank]) + item.log_scale for item in items
    ]
    top = max(
        log_n + log_j,
        2 * log_ac,
        2 * log_ap,
        2 * log_au,
        log_n + log_ac,
        log_n + log_ap,
        log_n + log_au,
    )
    n_joint = math.exp(log_n + log_j - top)
    ac = math.exp(log_ac - top / 2)
    ap = math.exp(log_ap - top / 2)
    au = math.exp(log_au - top / 2)
    product_ordered = (
        (ap + multiplicity * au) ** 2
        - ap**2
        - multiplicity * au**2
    )
    q_ordered = n_joint - product_ordered
    base = (
        ap**2 / parent_degree
        + multiplicity * au**2 / child_degree
        + ac * (ap + multiplicity * au)
    )
    neighbor_allocation = (
        math.exp(log_n + log_ap - top) / parent_degree
        + multiplicity * math.exp(log_n + log_au - top) / child_degree
    )
    center_allocation = math.exp(log_n + log_ac - top)
    residual_at_zero = q_ordered - base - neighbor_allocation
    slope = center_allocation - neighbor_allocation
    if slope > 0:
        kind, bound = "lower", residual_at_zero / slope
    elif slope < 0:
        kind, bound = "upper", residual_at_zero / slope
    elif residual_at_zero <= 0:
        kind, bound = "none", 0.0
    else:
        kind, bound = "impossible", math.inf
    half_gap = base + (neighbor_allocation + center_allocation) / 2 - q_ordered
    scale = max(abs(q_ordered), base + (neighbor_allocation + center_allocation) / 2, 1e-300)
    return {
        "kind": kind,
        "bound": bound,
        "half_relative_gap": half_gap / scale,
        "inverse_ratio": q_ordered / base if base else 0.0,
    }


def main() -> None:
    parent_parameters = (
        (28, 14), (32, 16), (36, 18), (40, 18), (40, 20),
        (44, 20), (44, 22), (48, 22), (48, 24), (52, 24), (52, 26),
    )
    child_parameters = tuple(
        (m, t) for m in (3, 4, 5, 6) for t in (5, 6, 7, 8)
    )
    multiplicities = (24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160)
    parents = [
        (m, t, make_state(f"G_{m}_{t}", make_galvin_tree(m, t), 0))
        for m, t in parent_parameters
    ]
    children = [
        (m, t, make_state(f"TG_{m}_{t}", make_bautista_ramos_tree(m, t), 0))
        for m, t in child_parameters
    ]
    lower, upper = 0.0, 1.0
    lower_witness = upper_witness = None
    min_half = None
    tested = 0
    for parent_m, parent_t, parent in parents:
        for child_m, child_t, child in children:
            for multiplicity in multiplicities:
                alpha = max(
                    len(parent.total) - 1
                    + multiplicity * (len(child.total) - 1),
                    1 + len(parent.excluded) - 1
                    + multiplicity * (len(child.excluded) - 1),
                )
                limit = math.ceil((2 * alpha - 1) / 3)
                rank = limit - 2
                log_z = centering_log_z(parent, child, multiplicity, rank)
                row = constraint(
                    parent, child, multiplicity, rank, log_z,
                    parent_m + 1, child_m + 2,
                )
                if row is None:
                    continue
                tested += 1
                payload = {
                    "parent": (parent_m, parent_t),
                    "child": (child_m, child_t),
                    "multiplicity": multiplicity,
                    "alpha": alpha,
                    "limit": limit,
                    "rank": rank,
                    **row,
                }
                if row["kind"] == "lower" and row["bound"] > lower:
                    lower, lower_witness = row["bound"], payload
                elif row["kind"] == "upper" and row["bound"] < upper:
                    upper, upper_witness = row["bound"], payload
                if min_half is None or row["half_relative_gap"] < min_half[0]:
                    min_half = (row["half_relative_gap"], payload)
                if lower > upper or row["half_relative_gap"] < -1e-7:
                    print({
                        "failure": payload,
                        "interval": (lower, upper),
                        "lower_witness": lower_witness,
                        "upper_witness": upper_witness,
                    }, flush=True)
                    raise SystemExit(1)
    print({
        "tested": tested,
        "lambda_interval": (lower, upper),
        "lower_witness": lower_witness,
        "upper_witness": upper_witness,
        "min_half": min_half,
    }, flush=True)


if __name__ == "__main__":
    main()
