#!/usr/bin/env python3
"""Tilted boundary scan around the sharp G(40,20)/TG(4,6) fork family."""

from __future__ import annotations

import heapq
import math

import numpy as np

from scratch_d25_hard_branch_composition_search_20260711 import make_state
from scratch_d25_two_species_asymptotic_20260711 import (
    centering_log_z,
    ratio_rows,
    shift,
    tilted,
    weighted,
)
from scratch_d17_spherical_search_20260711 import add, multiply, power
from scripts.analyze_prufer_corpus import make_bautista_ramos_tree, make_galvin_tree


def inverse_ratio(parent, child, multiplicity, log_z, rank, parent_degree, child_degree):
    p0, e0, r0 = (tilted(poly, log_z) for poly in
                  (parent.total, parent.excluded, parent.reserve))
    p1, e1, r1 = (tilted(poly, log_z) for poly in
                  (child.total, child.excluded, child.reserve))
    p1_m = power(p1, multiplicity)
    e1_m = power(e1, multiplicity)
    whole = add(multiply(p0, p1_m), shift(multiply(e0, e1_m), log_z))
    ac = multiply(e0, e1_m)
    ap = multiply(r0, p1_m)
    p1_m1 = power(p1, multiplicity - 1)
    au = multiply(multiply(r1, p0), p1_m1)
    jpc = multiply(multiply(r0, r1), p1_m1)
    jcc = multiply(
        multiply(multiply(r1, r1), p0), power(p1, multiplicity - 2)
    )
    joint = add(
        weighted(jpc, 2 * multiplicity),
        weighted(jcc, multiplicity * (multiplicity - 1)),
    )
    items = (whole, joint, ac, ap, au)
    if any(rank >= len(item.values) or item.values[rank] <= 0 for item in items):
        return None
    logs = [math.log(item.values[rank]) + item.log_scale for item in items]
    log_n, log_j, log_ac, log_ap, log_au = logs
    top = max(
        log_n + log_j,
        2 * log_ap,
        2 * log_au + 2 * math.log(multiplicity),
        log_ac + max(log_ap, log_au + math.log(multiplicity)),
    )
    n_j = math.exp(log_n + log_j - top)
    ap_value = math.exp(log_ap - top / 2)
    au_value = math.exp(log_au - top / 2)
    ac_value = math.exp(log_ac - top / 2)
    sum_a = ap_value + multiplicity * au_value
    product_ordered = (
        sum_a * sum_a - ap_value * ap_value
        - multiplicity * au_value * au_value
    )
    q_ordered = n_j - product_ordered
    budget = (
        ap_value * ap_value / parent_degree
        + multiplicity * au_value * au_value / child_degree
        + ac_value * sum_a
    )
    return q_ordered / budget


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
        make_state(f"G_{m}_{t}@0", make_galvin_tree(m, t), 0)
        for m, t in parent_parameters
    ]
    children = [
        make_state(f"TG_{m}_{t}@0", make_bautista_ramos_tree(m, t), 0)
        for m, t in child_parameters
    ]
    top: list[tuple[float, tuple]] = []
    crossings = []
    inverse_crossings = []
    inverse_top: list[tuple[float, tuple]] = []
    tested = 0
    for parent in parents:
        for child in children:
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
                rows = dict(
                    (r, ratio)
                    for ratio, r in ratio_rows(
                        parent, child, multiplicity, log_z, rank
                    )
                )
                ratio = rows.get(rank)
                if ratio is None:
                    continue
                tested += 1
                payload = (
                    parent.label, child.label, multiplicity,
                    alpha, limit, rank, log_z,
                )
                if len(top) < 20:
                    heapq.heappush(top, (ratio, payload))
                elif ratio > top[0][0]:
                    heapq.heapreplace(top, (ratio, payload))
                if ratio > 1.0 + 1e-6:
                    crossings.append((ratio, payload))
                parent_m = int(parent.label.split("_")[1])
                child_m = int(child.label.split("_")[1])
                inverse = inverse_ratio(
                    parent,
                    child,
                    multiplicity,
                    log_z,
                    rank,
                    parent_m + 1,
                    child_m + 2,
                )
                if inverse is not None:
                    if len(inverse_top) < 20:
                        heapq.heappush(inverse_top, (inverse, payload))
                    elif inverse > inverse_top[0][0]:
                        heapq.heapreplace(inverse_top, (inverse, payload))
                    if inverse > 1.0 + 1e-6:
                        inverse_crossings.append((inverse, payload))
    print(
        {
            "tested": tested,
            "crossings": len(crossings),
            "strongest_crossings": sorted(crossings, reverse=True)[:20],
            "inverse_crossings": len(inverse_crossings),
            "strongest_inverse_crossings": sorted(inverse_crossings, reverse=True)[:20],
            "inverse_top": sorted(inverse_top, reverse=True),
            "top": sorted(top, reverse=True),
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
