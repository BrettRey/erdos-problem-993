#!/usr/bin/env python3
"""Focused two-type recursive phase-grammar search for D23.

For a rooted decoration D and a previous rooted state (E,S), a layer with b
recursive children has

    E'=(E+S)^b P_D,       S'=x E^b E_D.

The frozen library contains exact rooted decorations selected by mass/mean/
variance diagnostics: one nearly critical 5-ary type, two alternating pairs
with a long near-neutral transient, and one finite-transient extremizer.  Thus
the search tests phase susceptibility, rather than another terminal leaf
plateau.

Floating point polynomials rank candidates.  A reported crossing is replayed
with integer coefficients and accompanied by a materialized tree certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from scratch_d17_checkpointed_search_20260711 import (
    append_record,
    strict_post_descent_pressure,
    utc_now,
)
from scratch_d17_spherical_search_20260711 import Scaled, add, multiply, normalize, power
from scratch_d20_state_grammar_search_20260711 import feature, rooted_states
from scratch_d23_phase_grammar_search_20260711 import (
    Moments,
    mixture,
    moments,
    poly_add,
    poly_mul,
    poly_pow,
    valleys,
)


Edge = tuple[int, int]


@dataclass(frozen=True)
class DecorSpec:
    label: str
    order: int
    root: int
    edges: tuple[Edge, ...]


@dataclass(frozen=True)
class LayerSpec:
    branch: int
    decor: str


@dataclass(frozen=True)
class GrammarSpec:
    label: str
    period: tuple[LayerSpec, ...]


DECORATIONS = (
    DecorSpec("empty", 0, 0, ()),
    # q=280/367, only 3.326e-6 above the regular 5-ary bifurcation.
    DecorSpec(
        "critical_280_367",
        13,
        11,
        ((0, 9), (1, 0), (1, 2), (2, 3), (2, 6), (2, 8),
         (3, 4), (3, 5), (6, 7), (9, 10), (10, 11), (11, 12)),
    ),
    # First alternating pair: a long transient has response product near one;
    # the eventual stable two-step orbit has multiplier about 0.959965.
    DecorSpec(
        "marginal_a",
        12,
        11,
        ((0, 3), (0, 5), (0, 7), (0, 9), (0, 11), (1, 0),
         (1, 2), (3, 4), (5, 6), (7, 8), (9, 10)),
    ),
    DecorSpec(
        "marginal_b",
        12,
        5,
        ((0, 5), (1, 0), (1, 2), (1, 4), (2, 3), (5, 6),
         (5, 7), (5, 8), (5, 9), (5, 10), (5, 11)),
    ),
    # Second alternating pair: another long near-neutral transient; the
    # eventual stable two-step orbit has multiplier about 0.964581.
    DecorSpec(
        "near_a",
        12,
        1,
        ((0, 5), (0, 9), (0, 11), (1, 0), (1, 2), (2, 3),
         (2, 4), (5, 6), (5, 7), (5, 8), (9, 10)),
    ),
    DecorSpec(
        "near_b",
        9,
        0,
        ((0, 4), (0, 7), (1, 0), (1, 2), (2, 3), (4, 5),
         (4, 6), (7, 8)),
    ),
    # Best finite-transient decoration in the scalar diagnostic round.
    DecorSpec(
        "transient_b",
        12,
        11,
        ((0, 6), (1, 0), (1, 2), (1, 5), (2, 3), (3, 4),
         (6, 7), (6, 9), (6, 11), (7, 8), (9, 10)),
    ),
)


GRAMMARS = (
    GrammarSpec("critical5", (LayerSpec(5, "critical_280_367"),)),
    GrammarSpec("marginal47", (LayerSpec(4, "marginal_a"), LayerSpec(7, "marginal_b"))),
    GrammarSpec("near56", (LayerSpec(5, "near_a"), LayerSpec(6, "near_b"))),
    GrammarSpec("transient27", (LayerSpec(2, "empty"), LayerSpec(7, "transient_b"))),
)


@dataclass(frozen=True)
class DecorData:
    spec: DecorSpec
    excluded_exact: list[int]
    total_exact: list[int]
    excluded_float: Scaled
    total_float: Scaled
    excluded_moments: Moments
    total_moments: Moments


def exact_to_scaled(coefficients: list[int]) -> Scaled:
    logs = np.array([math.log(value) if value else -math.inf for value in coefficients])
    scale = float(np.max(logs))
    return normalize(Scaled(np.exp(logs - scale), scale))


def decoration_data() -> dict[str, DecorData]:
    result: dict[str, DecorData] = {}
    for spec in DECORATIONS:
        if spec.order == 0:
            excluded, selected = [1], [0]
        else:
            adjacency = [[] for _ in range(spec.order)]
            for u, v in spec.edges:
                adjacency[u].append(v)
                adjacency[v].append(u)
            excluded, selected = rooted_states(adjacency, spec.root)
        total = poly_add(excluded, selected)
        result[spec.label] = DecorData(
            spec,
            excluded,
            total,
            exact_to_scaled(excluded),
            exact_to_scaled(total),
            moments(excluded),
            moments(total),
        )
    return result


def shift_x(polynomial: Scaled) -> Scaled:
    return Scaled(np.pad(polynomial.values, (1, 0)), polynomial.log_scale)


def float_states(
    layers: tuple[LayerSpec, ...], data: dict[str, DecorData]
) -> tuple[Scaled, Scaled]:
    excluded = Scaled(np.array([1.0]), 0.0)
    selected = Scaled(np.array([0.0, 1.0]), 0.0)
    for layer in layers:
        old_excluded = excluded
        total = add(excluded, selected)
        decor = data[layer.decor]
        excluded = multiply(power(total, layer.branch), decor.total_float)
        selected = shift_x(multiply(power(old_excluded, layer.branch), decor.excluded_float))
    return excluded, selected


def exact_states(
    layers: tuple[LayerSpec, ...], data: dict[str, DecorData]
) -> tuple[list[int], list[int]]:
    excluded, selected = [1], [0, 1]
    for layer in layers:
        old_excluded = excluded
        total = poly_add(excluded, selected)
        decor = data[layer.decor]
        excluded = poly_mul(poly_pow(total, layer.branch), decor.total_exact)
        selected = poly_mul([0, 1], poly_mul(poly_pow(old_excluded, layer.branch), decor.excluded_exact))
    return excluded, selected


def order(layers: tuple[LayerSpec, ...], data: dict[str, DecorData]) -> int:
    value = 1
    for layer in layers:
        value = layer.branch * value + 1 + data[layer.decor].spec.order
    return value


def specs(max_depth: int, data: dict[str, DecorData], max_order: int):
    result = []
    for grammar in GRAMMARS:
        for offset in range(len(grammar.period)):
            rotated = grammar.period[offset:] + grammar.period[:offset]
            for depth in range(1, max_depth + 1):
                layers = tuple(rotated[k % len(rotated)] for k in range(depth))
                if order(layers, data) <= max_order:
                    result.append((f"{grammar.label}_o{offset}_d{depth}", layers))
    unique = {layers: label for label, layers in result}
    return sorted(((label, layers) for layers, label in unique.items()), key=lambda z: order(z[1], data))


def build_half(
    layers: tuple[LayerSpec, ...], data: dict[str, DecorData]
) -> tuple[list[Edge], int, int]:
    edges: list[Edge] = []
    next_vertex = [0]

    def copy_decoration(spec: DecorSpec) -> int:
        start = next_vertex[0]
        next_vertex[0] += spec.order
        edges.extend((start + u, start + v) for u, v in spec.edges)
        return start + spec.root

    def recurse(depth: int) -> int:
        root = next_vertex[0]
        next_vertex[0] += 1
        if depth < 0:
            return root
        layer = layers[depth]
        for _ in range(layer.branch):
            child = recurse(depth - 1)
            edges.append((root, child))
        decor = data[layer.decor].spec
        if decor.order:
            decor_root = copy_decoration(decor)
            edges.append((root, decor_root))
        return root

    root = recurse(len(layers) - 1)
    return edges, next_vertex[0], root


def joined_tree_certificate(
    left: tuple[LayerSpec, ...],
    right: tuple[LayerSpec, ...],
    data: dict[str, DecorData],
) -> dict[str, int | bool]:
    le, ln, lr = build_half(left, data)
    re, rn, rr = build_half(right, data)
    edges = le + [(ln + u, ln + v) for u, v in re] + [(lr, ln + rr)]
    n = ln + rn
    adjacency = [[] for _ in range(n)]
    for u, v in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    seen = {0}
    queue = [0]
    for vertex in queue:
        for neighbor in adjacency[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    return {"order": n, "edges": len(edges), "connected": len(seen) == n, "acyclic": len(edges) == n - 1}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-depth", type=int, default=6)
    parser.add_argument("--max-order", type=int, default=120000)
    parser.add_argument("--output", type=Path, default=Path("results/d23_two_type_phase_search_20260711.jsonl"))
    parser.add_argument("--skip-exact", action="store_true")
    parser.add_argument("--stop-first-crossing", action="store_true")
    parser.add_argument("--relative-floor", type=float, default=1e-10)
    args = parser.parse_args()

    data = decoration_data()
    profiles = specs(args.max_depth, data, args.max_order)
    states = {label: float_states(layers, data) for label, layers in profiles}
    rows = []
    crossings = []
    for left_label, left_layers in profiles:
        for right_label, right_layers in profiles:
            n = order(left_layers, data) + order(right_layers, data)
            if n > args.max_order:
                continue
            le, ls = states[left_label]
            re, rs = states[right_label]
            full = add(multiply(le, re), add(multiply(ls, re), multiply(le, rs)))
            values = np.maximum(full.values, 0.0)
            values /= values.sum()
            pressure = strict_post_descent_pressure(values, relative_floor=args.relative_floor, descent_tolerance=1e-10, min_separation=1)
            if pressure is None:
                continue
            lf = feature(ls, relative_floor=args.relative_floor, descent_tolerance=1e-10, min_separation=1)
            rf = feature(rs, relative_floor=args.relative_floor, descent_tolerance=1e-10, min_separation=1)
            row = {"left": left_label, "right": right_label, "order": n, "left_selected": lf.__dict__, "right_selected": rf.__dict__, **pressure}
            rows.append(row)
            if float(pressure["best_later_ratio"]) > 1.0 + 1e-10:
                append_record(args.output, {"at": utc_now(), "kind": "d23_two_type_float_crossing", **row})
                print(json.dumps({"float_crossing": row}, sort_keys=True), flush=True)
                if args.stop_first_crossing:
                    return
                if args.skip_exact:
                    crossings.append(row)
                    continue
                lee, les = exact_states(left_layers, data)
                ree, res = exact_states(right_layers, data)
                coefficients = poly_add(poly_mul(lee, ree), poly_add(poly_mul(les, ree), poly_mul(lee, res)))
                row["exact_valleys"] = valleys(coefficients)
                row["tree"] = joined_tree_certificate(left_layers, right_layers, data)
                crossings.append(row)
                append_record(args.output, {"at": utc_now(), "kind": "d23_two_type_exact_crossing", **row})

    best = max(rows, key=lambda row: float(row["best_later_ratio"]), default=None)
    append_record(args.output, {"at": utc_now(), "kind": "d23_two_type_complete", "profiles": len(profiles), "pairs_with_descent": len(rows), "best": best, "exact_crossings": crossings})
    print(json.dumps({"profiles": len(profiles), "pairs_with_descent": len(rows), "best": best, "crossings": crossings}, sort_keys=True))


if __name__ == "__main__":
    main()
