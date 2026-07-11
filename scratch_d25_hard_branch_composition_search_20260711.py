#!/usr/bin/env python3
"""Adversarial marked-ULC/fork scan on repeated hard rooted gadgets.

This targets the exact oriented local inequality at a hub c with one parent
branch and m child branches.  The children are repeated copies of a hard
rooted state; the parent may be a different hard state.  Float convolution
ranks candidates only.  Any apparent crossing is replayed with exact Python
integers before being reported.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

import numpy as np

from graph6 import parse_graph6
from indpoly import _polyadd, _polymul
from scratch_d17_spherical_search_20260711 import Scaled, add, multiply, normalize, power
from scratch_d20_state_grammar_search_20260711 import rooted_states
from scripts.analyze_prufer_corpus import make_bautista_ramos_tree, make_galvin_tree


@dataclass(frozen=True)
class HardState:
    label: str
    total: tuple[int, ...]
    excluded: tuple[int, ...]
    reserve: tuple[int, ...]


def make_state(label: str, adj: list[list[int]], root: int) -> HardState:
    excluded, selected = rooted_states(adj, root)
    assert selected[0] == 0
    reserve = selected[1:]
    total = _polyadd(excluded, selected)
    def trim(poly: list[int]) -> list[int]:
        while len(poly) > 1 and poly[-1] == 0:
            poly.pop()
        return poly
    total, excluded, reserve = map(trim, (total, excluded, reserve))
    return HardState(label, tuple(total), tuple(excluded), tuple(reserve))


def scaled(poly: tuple[int, ...] | list[int]) -> Scaled:
    total = sum(poly)
    values = np.array([value / total for value in poly], dtype=float)
    return normalize(Scaled(values, math.log(total)))


def shift(item: Scaled) -> Scaled:
    return Scaled(np.pad(item.values, (1, 0)), item.log_scale)


def weighted(item: Scaled, weight: int) -> Scaled:
    return Scaled(item.values, item.log_scale + math.log(weight))


def best_pointwise_ratio(
    n_poly: Scaled,
    joint: Scaled,
    child_sum: Scaled,
    augmented_sum: Scaled,
    max_rank: int,
) -> tuple[float, int]:
    limit = min(
        len(n_poly.values),
        len(joint.values),
        len(child_sum.values),
        len(augmented_sum.values),
    )
    best = 0.0
    best_rank = -1
    log_scale_ratio = (
        n_poly.log_scale
        + joint.log_scale
        - child_sum.log_scale
        - augmented_sum.log_scale
    )
    floors = tuple(
        np.max(item.values) * 1e-13
        for item in (n_poly, joint, child_sum, augmented_sum)
    )
    for rank in range(min(limit, max_rank + 1)):
        values = (
            n_poly.values[rank],
            joint.values[rank],
            child_sum.values[rank],
            augmented_sum.values[rank],
        )
        if any(value <= floor for value, floor in zip(values, floors)):
            continue
        log_ratio = log_scale_ratio + sum(map(math.log, values[:2])) - sum(
            map(math.log, values[2:])
        )
        ratio = math.exp(log_ratio)
        if ratio > best:
            best, best_rank = float(ratio), rank
    return best, best_rank


def float_candidate(parent: HardState, child: HardState, multiplicity: int):
    p0, e0, r0 = map(scaled, (parent.total, parent.excluded, parent.reserve))
    p1, e1, r1 = map(scaled, (child.total, child.excluded, child.reserve))
    p1_m = power(p1, multiplicity)
    e1_m = power(e1, multiplicity)
    whole = add(multiply(p0, p1_m), shift(multiply(e0, e1_m)))
    a_c = multiply(e0, e1_m)
    a_p = multiply(r0, p1_m)
    a_child = multiply(multiply(r1, p0), power(p1, multiplicity - 1))
    j_pc = multiply(multiply(r0, r1), power(p1, multiplicity - 1))
    if multiplicity >= 2:
        j_cc = multiply(
            multiply(multiply(r1, r1), p0),
            power(p1, multiplicity - 2),
        )
        joint_ordered = add(
            weighted(j_pc, 2 * multiplicity),
            weighted(j_cc, multiplicity * (multiplicity - 1)),
        )
    else:
        joint_ordered = weighted(j_pc, 2)
    child_sum = weighted(a_child, multiplicity)
    boundary = add(a_p, a_c)
    augmented_sum = add(child_sum, weighted(boundary, 2))
    alpha = max(
        len(parent.total) - 1 + multiplicity * (len(child.total) - 1),
        1
        + len(parent.excluded)
        - 1
        + multiplicity * (len(child.excluded) - 1),
    )
    prefix_limit = math.ceil((2 * alpha - 1) / 3)
    ratio, rank = best_pointwise_ratio(
        whole,
        joint_ordered,
        child_sum,
        augmented_sum,
        max_rank=prefix_limit - 2,
    )
    return ratio, rank, alpha, prefix_limit


def exact_mul(left: list[int], right: list[int], cap: int) -> list[int]:
    out = [0] * min(cap + 1, len(left) + len(right) - 1)
    for i, a in enumerate(left[: cap + 1]):
        if not a:
            continue
        for j, b in enumerate(right[: cap + 1 - i]):
            if b:
                out[i + j] += a * b
    return out


def exact_power(poly: list[int], exponent: int, cap: int) -> list[int]:
    out = [1]
    while exponent:
        if exponent & 1:
            out = exact_mul(out, poly, cap)
        exponent //= 2
        if exponent:
            poly = exact_mul(poly, poly, cap)
    return out


def exact_replay(parent: HardState, child: HardState, multiplicity: int, rank: int):
    p0, e0, r0 = map(list, (parent.total, parent.excluded, parent.reserve))
    p1, e1, r1 = map(list, (child.total, child.excluded, child.reserve))
    p1_m = exact_power(p1, multiplicity, rank)
    e1_m = exact_power(e1, multiplicity, rank)
    excluded_product = exact_mul(e0, e1_m, rank)
    whole = _polyadd(exact_mul(p0, p1_m, rank), [0] + excluded_product[:rank])
    ac = excluded_product[rank] if rank < len(excluded_product) else 0
    ap_poly = exact_mul(r0, p1_m, rank)
    au_poly = exact_mul(
        exact_mul(r1, p0, rank),
        exact_power(p1, multiplicity - 1, rank),
        rank,
    )
    jpc_poly = exact_mul(
        exact_mul(r0, r1, rank),
        exact_power(p1, multiplicity - 1, rank),
        rank,
    )
    j_ordered = 2 * multiplicity * (jpc_poly[rank] if rank < len(jpc_poly) else 0)
    if multiplicity >= 2:
        jcc_poly = exact_mul(
            exact_mul(exact_mul(r1, r1, rank), p0, rank),
            exact_power(p1, multiplicity - 2, rank),
            rank,
        )
        j_ordered += multiplicity * (multiplicity - 1) * (
            jcc_poly[rank] if rank < len(jcc_poly) else 0
        )
    n_r = whole[rank] if rank < len(whole) else 0
    ap = ap_poly[rank] if rank < len(ap_poly) else 0
    au = au_poly[rank] if rank < len(au_poly) else 0
    child_sum = multiplicity * au
    lhs = n_r * j_ordered
    rhs = child_sum * (child_sum + 2 * (ap + ac))
    return {"rank": rank, "lhs": lhs, "rhs": rhs, "gap": rhs - lhs}


def library() -> list[HardState]:
    out: list[HardState] = []
    data = json.load(open("results/analysis_n26.json"))
    for index, row in enumerate(data["lc_failures"]):
        _, adj = parse_graph6(row["graph6"].encode())
        for root in range(len(adj)):
            out.append(make_state(f"LC26_{index}@{root}", adj, root))
    for m, t in ((2, 5), (4, 6), (6, 7), (8, 8), (12, 10)):
        adj = make_bautista_ramos_tree(m, t)
        out.append(make_state(f"TG_{m}_{t}@0", adj, 0))
    for m, t in ((3, 4), (6, 6), (14, 8), (21, 11), (40, 20)):
        adj = make_galvin_tree(m, t)
        out.append(make_state(f"G_{m}_{t}@0", adj, 0))
    return out


def main() -> None:
    states = library()
    best = (0.0, None)
    tested = 0
    for parent in states:
        for child in states:
            for multiplicity in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32):
                ratio, rank, alpha, prefix_limit = float_candidate(
                    parent, child, multiplicity
                )
                tested += 1
                payload = (
                    parent.label,
                    child.label,
                    multiplicity,
                    rank,
                    alpha,
                    prefix_limit,
                )
                if ratio > best[0]:
                    best = (ratio, payload)
                if ratio > 1.0 + 1e-7:
                    print({"prefix_float_crossing": payload, "ratio": ratio}, flush=True)
                    exact = exact_replay(parent, child, multiplicity, rank)
                    if exact["gap"] < 0:
                        print({"failure": payload, "float_ratio": ratio, "exact": exact})
                        raise SystemExit(1)
    parent_label, child_label, multiplicity, rank, _, _ = best[1]
    parent = next(state for state in states if state.label == parent_label)
    child = next(state for state in states if state.label == child_label)
    exact = (
        exact_replay(parent, child, multiplicity, rank)
        if best[1][4] <= 1200
        else {"skipped": "float ranking only above alpha 1200"}
    )
    print(
        {
            "states": len(states),
            "candidates": tested,
            "best_float_ratio": best[0],
            "best": best[1],
            "best_exact": exact,
            "failure": None,
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
