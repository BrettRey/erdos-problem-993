#!/usr/bin/env python3
"""Search mixed two-level hard-hub compositions at the last prefix rank.

The outer hub has one G(40,20) branch, ``k`` copies of the inner hard hub
H_m(TG(4,6)), and ``b-k`` ordinary TG(4,6) branches.  This interpolates
between the original D25 obstruction and the regular two-level composition.
Large coefficient products are evaluated by saddle-point tilting and FFT.
"""

from __future__ import annotations

import argparse
import json
import math

from scratch_d25_hard_branch_composition_search_20260711 import make_state
from scratch_d25_nested_hub_scalar_search_20260711 import (
    DegreeState,
    MomentState,
    ScaledState,
    add,
    base_degree_state,
    base_moment_state,
    base_scaled_state,
    bottom_inverse_ratio,
    compose_degree,
    compose_moment,
    compose_scaled,
    moment_add,
    moment_power,
    moment_product,
    moment_shift,
    multiply,
    power,
    shift,
)
from scripts.analyze_prufer_corpus import (
    make_bautista_ramos_tree,
    make_galvin_tree,
)


def multi_moment(
    parent: MomentState,
    children: list[tuple[MomentState, int]],
    log_z: float,
) -> MomentState:
    excluded = parent.total
    reserve = parent.excluded
    for child, count in children:
        excluded = moment_product(excluded, moment_power(child.total, count))
        reserve = moment_product(reserve, moment_power(child.excluded, count))
    return MomentState(
        moment_add(excluded, moment_shift(reserve, log_z)),
        excluded,
        reserve,
    )


def multi_scaled(
    parent: ScaledState,
    children: list[tuple[ScaledState, int]],
    log_z: float,
) -> ScaledState:
    excluded = parent.total
    reserve = parent.excluded
    for child, count in children:
        excluded = multiply(excluded, power(child.total, count))
        reserve = multiply(reserve, power(child.excluded, count))
    return ScaledState(add(excluded, shift(reserve, log_z)), excluded, reserve)


def multi_degree(
    parent: DegreeState,
    children: list[tuple[DegreeState, int]],
) -> DegreeState:
    excluded = parent.total + sum(count * child.total for child, count in children)
    reserve = parent.excluded + sum(
        count * child.excluded for child, count in children
    )
    return DegreeState(
        max(excluded, 1 + reserve),
        excluded,
        reserve,
        1 + parent.order + sum(count * child.order for child, count in children),
        1 + sum(count * child.hubs for child, count in children),
    )


def evaluate(inner_m: int, hard_copies: int, outer_branches: int, parent, base):
    if not 1 <= hard_copies <= outer_branches:
        raise ValueError((hard_copies, outer_branches))
    ordinary_copies = outer_branches - hard_copies
    parent_degree = base_degree_state(parent, 1641)
    base_degree = base_degree_state(base, 162)
    hard_degree = compose_degree(parent_degree, base_degree, inner_m)
    degree = multi_degree(
        parent_degree,
        [(hard_degree, hard_copies), (base_degree, ordinary_copies)],
    )
    alpha = degree.total
    rank = math.ceil((2 * alpha - 1) / 3) - 2

    def state_at(log_z: float, scaled: bool = False):
        if scaled:
            parent_s = base_scaled_state(parent, log_z)
            base_s = base_scaled_state(base, log_z)
            hard_s = compose_scaled(parent_s, base_s, inner_m, log_z)
            whole_s = multi_scaled(
                parent_s,
                [(hard_s, hard_copies), (base_s, ordinary_copies)],
                log_z,
            )
            return parent_s, base_s, hard_s, whole_s
        parent_s = base_moment_state(parent, log_z)
        base_s = base_moment_state(base, log_z)
        hard_s = compose_moment(parent_s, base_s, inner_m, log_z)
        return multi_moment(
            parent_s,
            [(hard_s, hard_copies), (base_s, ordinary_copies)],
            log_z,
        )

    lo, hi = -8.0, 8.0
    for _ in range(70):
        mid = (lo + hi) / 2
        if state_at(mid).total.mean < rank + 1:
            lo = mid
        else:
            hi = mid
    log_z = (lo + hi) / 2
    z = math.exp(log_z)
    parent_s, base_s, hard_s, whole_s = state_at(log_z, scaled=True)
    w0, w1, w2 = map(float, whole_s.total.values[rank : rank + 3])
    ordered = (rank + 2) * w0 * w2 / ((rank + 1) * w1 * w1)
    gsb = (rank + 2) * w0 * w2 / (
        (rank + 1) * w1 * w1 + z * w0 * w1
    )

    # Audit a bottom hard hub.  Removing its edge to the outer hub leaves an
    # exterior rooted at the outer hub with one fewer hard branch.
    exterior = multi_scaled(
        parent_s,
        [(hard_s, hard_copies - 1), (base_s, ordinary_copies)],
        log_z,
    )
    local_inverse = bottom_inverse_ratio(
        [
            (parent_s, 1, 41),
            (base_s, inner_m, 6),
            (exterior, 1, outer_branches + 1),
        ],
        rank,
        log_z,
    )
    return {
        "inner_m": inner_m,
        "hard_copies": hard_copies,
        "outer_branches": outer_branches,
        "n": degree.order,
        "alpha": alpha,
        "hubs": degree.hubs,
        "rank": rank,
        "z": z,
        "ordered_lc_ratio": ordered,
        "gsb_ratio": gsb,
        "bottom_inverse_ratio": local_inverse,
        "bottom_inverse_fails": local_inverse > 1.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-branches", type=int, default=48)
    parser.add_argument("--inner", default="32,48,80,112,160")
    parser.add_argument("--hard", default="1,2,4,8,16,24,32,40,48")
    args = parser.parse_args()
    parent = make_state("G", make_galvin_tree(40, 20), 0)
    base = make_state("TG", make_bautista_ramos_tree(4, 6), 0)
    rows = []
    for inner_m in map(int, args.inner.split(",")):
        for hard_copies in map(int, args.hard.split(",")):
            if hard_copies > args.outer_branches:
                continue
            row = evaluate(inner_m, hard_copies, args.outer_branches, parent, base)
            rows.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)
    print(json.dumps({
        "tested": len(rows),
        "max_ordered_lc": max(rows, key=lambda row: row["ordered_lc_ratio"]),
        "max_gsb": max(rows, key=lambda row: row["gsb_ratio"]),
        "local_inverse_failures": sum(row["bottom_inverse_fails"] for row in rows),
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
