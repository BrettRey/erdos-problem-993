#!/usr/bin/env python3
"""Exact/FFT audit of recursively decorated hard-core phase trees.

A layer ``(b,s)`` consists of a new root, ``b`` copies of the previous
rooted state, and one additional child which is the centre of an ``s``-leaf
star.  If ``(E,P)`` are the root-excluded and total states below, then

    E_new = P^b ((1+x)^s+x),
    S_new = x E^b (1+x)^s,
    P_new = E_new+S_new.

At x=1 the effective activity is 2^s/(2^s+1).  In particular (b,s)=(5,2)
lies just above the regular 5-ary hard-core bifurcation and has a genuine
stable two-cycle, while neither root state has exponentially small mass.

The script first computes exact mass/mean/variance diagnostics using integer
arithmetic.  Polynomial floats rank joined pairs only.  Any valley is replayed
in ZZ[x], checked against an independently built tree, and emitted with a
treehood certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import numpy as np

from scratch_d17_checkpointed_search_20260711 import (
    append_record,
    strict_post_descent_pressure,
    utc_now,
)
from scratch_d17_spherical_search_20260711 import Scaled, add, multiply, normalize, power


def poly_add(a: list[int], b: list[int]) -> list[int]:
    out = [0] * max(len(a), len(b))
    for k, value in enumerate(a):
        out[k] += value
    for k, value in enumerate(b):
        out[k] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        if left:
            for j, right in enumerate(b):
                if right:
                    out[i + j] += left * right
    return out


def poly_pow(base: list[int], exponent: int) -> list[int]:
    result = [1]
    while exponent:
        if exponent & 1:
            result = poly_mul(result, base)
        exponent //= 2
        if exponent:
            base = poly_mul(base, base)
    return result


def binomial(exponent: int) -> list[int]:
    return [math.comb(exponent, k) for k in range(exponent + 1)]


def exact_states(layers: tuple[tuple[int, int], ...]) -> tuple[list[int], list[int]]:
    """Return (E,S) for layers ordered from leaves toward the root."""
    excluded = [1]
    selected = [0, 1]
    for branch, leaves in layers:
        old_excluded = excluded
        total = poly_add(excluded, selected)
        decoration_e = binomial(leaves)
        decoration_p = poly_add(decoration_e, [0, 1])
        excluded = poly_mul(poly_pow(total, branch), decoration_p)
        selected = poly_mul([0, 1], poly_mul(poly_pow(old_excluded, branch), decoration_e))
    return excluded, selected


def exact_to_scaled(coefficients: list[int]) -> Scaled:
    logs = np.array([math.log(value) if value else -math.inf for value in coefficients])
    scale = float(np.max(logs))
    values = np.exp(logs - scale)
    return normalize(Scaled(values, scale))


def float_states(layers: tuple[tuple[int, int], ...]) -> tuple[Scaled, Scaled]:
    excluded = Scaled(np.array([1.0]), 0.0)
    selected = Scaled(np.array([0.0, 1.0]), 0.0)
    binomial_cache: dict[int, Scaled] = {}
    for branch, leaves in layers:
        old_excluded = excluded
        total = add(excluded, selected)
        if leaves not in binomial_cache:
            binomial_cache[leaves] = exact_to_scaled(binomial(leaves))
        decoration_e = binomial_cache[leaves]
        decoration_p = add(decoration_e, Scaled(np.array([0.0, 1.0]), 0.0))
        excluded = multiply(power(total, branch), decoration_p)
        selected = multiply(
            Scaled(np.array([0.0, 1.0]), 0.0),
            multiply(power(old_excluded, branch), decoration_e),
        )
    return excluded, selected


@dataclass(frozen=True)
class Moments:
    mass: int
    mean: Fraction
    variance: Fraction


def moments(coefficients: list[int]) -> Moments:
    mass = sum(coefficients)
    first = sum(k * value for k, value in enumerate(coefficients))
    second = sum(k * k * value for k, value in enumerate(coefficients))
    mean = Fraction(first, mass)
    return Moments(mass, mean, Fraction(second, mass) - mean * mean)


def mixture(left: Moments, right: Moments) -> Moments:
    mass = left.mass + right.mass
    mean = Fraction(left.mass * left.mean + right.mass * right.mean, mass)
    second = Fraction(
        left.mass * (left.variance + left.mean * left.mean)
        + right.mass * (right.variance + right.mean * right.mean),
        mass,
    )
    return Moments(mass, mean, second - mean * mean)


def product(moment: Moments, exponent: int) -> Moments:
    return Moments(
        moment.mass**exponent,
        exponent * moment.mean,
        exponent * moment.variance,
    )


def moment_states(layers: tuple[tuple[int, int], ...]) -> tuple[Moments, Moments]:
    """Exact x=1 mass, mean and variance without expanding coefficients."""
    excluded = Moments(1, Fraction(0), Fraction(0))
    selected = Moments(1, Fraction(1), Fraction(0))
    for branch, leaves in layers:
        old_excluded = excluded
        total = mixture(excluded, selected)
        decoration_e = Moments(2**leaves, Fraction(leaves, 2), Fraction(leaves, 4))
        decoration_s = Moments(1, Fraction(1), Fraction(0))
        decoration_p = mixture(decoration_e, decoration_s)
        total_power = product(total, branch)
        excluded = Moments(
            total_power.mass * decoration_p.mass,
            total_power.mean + decoration_p.mean,
            total_power.variance + decoration_p.variance,
        )
        excluded_power = product(old_excluded, branch)
        selected = Moments(
            excluded_power.mass * decoration_e.mass,
            Fraction(1) + excluded_power.mean + decoration_e.mean,
            excluded_power.variance + decoration_e.variance,
        )
    return excluded, selected


def layer_order(layers: tuple[tuple[int, int], ...]) -> int:
    order = 1
    for branch, leaves in layers:
        order = branch * order + leaves + 2
    return order


def joined_float(
    left: tuple[tuple[int, int], ...],
    right: tuple[tuple[int, int], ...],
) -> np.ndarray:
    left_e, left_s = float_states(left)
    right_e, right_s = float_states(right)
    full = add(
        multiply(left_e, right_e),
        add(multiply(left_s, right_e), multiply(left_e, right_s)),
    )
    values = np.maximum(full.values, 0.0)
    return values / values.sum()


def joined_exact(
    left: tuple[tuple[int, int], ...],
    right: tuple[tuple[int, int], ...],
) -> list[int]:
    left_e, left_s = exact_states(left)
    right_e, right_s = exact_states(right)
    return poly_add(
        poly_mul(left_e, right_e),
        poly_add(poly_mul(left_s, right_e), poly_mul(left_e, right_s)),
    )


def valleys(coefficients: list[int]) -> list[int]:
    return [
        k
        for k in range(1, len(coefficients) - 1)
        if coefficients[k - 1] > coefficients[k] < coefficients[k + 1]
    ]


def build_tree(layers: tuple[tuple[int, int], ...]) -> tuple[list[tuple[int, int]], int, int]:
    """Materialize one rooted half; returns edges, order, root."""
    edges: list[tuple[int, int]] = []

    def copy_layer(depth: int) -> int:
        root = next_vertex[0]
        next_vertex[0] += 1
        if depth < 0:
            return root
        branch, leaves = layers[depth]
        for _ in range(branch):
            child = copy_layer(depth - 1)
            edges.append((root, child))
        centre = next_vertex[0]
        next_vertex[0] += 1
        edges.append((root, centre))
        for _ in range(leaves):
            leaf = next_vertex[0]
            next_vertex[0] += 1
            edges.append((centre, leaf))
        return root

    next_vertex = [0]
    root = copy_layer(len(layers) - 1)
    return edges, next_vertex[0], root


def tree_certificate(
    left: tuple[tuple[int, int], ...],
    right: tuple[tuple[int, int], ...],
) -> dict[str, int | bool]:
    left_edges, left_order, left_root = build_tree(left)
    right_edges, right_order, right_root = build_tree(right)
    shifted = [(u + left_order, v + left_order) for u, v in right_edges]
    edges = left_edges + shifted + [(left_root, right_root + left_order)]
    order = left_order + right_order
    seen = {0}
    frontier = [0]
    adjacency = [[] for _ in range(order)]
    for u, v in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    for vertex in frontier:
        for neighbor in adjacency[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                frontier.append(neighbor)
    return {
        "order": order,
        "edges": len(edges),
        "connected": len(seen) == order,
        "acyclic": len(edges) == order - 1,
    }


def layer_sequences(max_depth: int) -> Iterable[tuple[tuple[int, int], ...]]:
    phase = ((5, 2),)
    alternates = (
        ((5, 2), (5, 2)),
        ((4, 0), (6, 1)),
        ((6, 1), (4, 0)),
        ((5, 1), (5, 3)),
        ((5, 3), (5, 1)),
    )
    for depth in range(1, max_depth + 1):
        yield phase * depth
        for period in alternates:
            yield tuple(period[k % len(period)] for k in range(depth))


def diagnostic(layers: tuple[tuple[int, int], ...]) -> dict[str, object]:
    e, s = moment_states(layers)
    total = mixture(e, s)
    delta = s.mean - e.mean
    return {
        "layers": [list(item) for item in layers],
        "order": layer_order(layers),
        "selected_probability": float(Fraction(s.mass, e.mass + s.mass)),
        "mean_excluded": float(e.mean),
        "mean_selected": float(s.mean),
        "mean_delta": float(delta),
        "variance_total": float(total.variance),
        "standardized_delta": float(abs(delta) / math.sqrt(float(total.variance))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-depth", type=int, default=5)
    parser.add_argument("--max-order", type=int, default=40000)
    parser.add_argument("--output", type=Path, default=Path("results/d23_phase_grammar_search_20260711.jsonl"))
    args = parser.parse_args()

    sequences = sorted({
        layers
        for layers in layer_sequences(args.max_depth)
        if layer_order(layers) <= args.max_order
    }, key=lambda layers: (layer_order(layers), layers))
    diagnostics = [diagnostic(layers) for layers in sequences]
    append_record(
        args.output,
        {
            "at": utc_now(),
            "kind": "d23_phase_diagnostics",
            "diagnostics": diagnostics,
        },
    )

    best: dict[str, object] | None = None
    crossings: list[dict[str, object]] = []
    state_cache = {layers: float_states(layers) for layers in sequences}
    for left in sequences:
        for right in sequences:
            order = layer_order(left) + layer_order(right)
            if order > args.max_order:
                continue
            left_e, left_s = state_cache[left]
            right_e, right_s = state_cache[right]
            full = add(
                multiply(left_e, right_e),
                add(multiply(left_s, right_e), multiply(left_e, right_s)),
            )
            values = np.maximum(full.values, 0.0)
            values /= values.sum()
            pressure = strict_post_descent_pressure(
                values,
                relative_floor=1e-14,
                descent_tolerance=1e-10,
                min_separation=1,
            )
            if pressure is None:
                continue
            row = {
                "left": [list(item) for item in left],
                "right": [list(item) for item in right],
                "order": order,
                **pressure,
            }
            if best is None or float(row["best_later_ratio"]) > float(best["best_later_ratio"]):
                best = row
            if float(row["best_later_ratio"]) > 1.0 + 1e-10:
                coefficients = joined_exact(left, right)
                exact_valleys = valleys(coefficients)
                row["exact_valleys"] = exact_valleys
                row["tree"] = tree_certificate(left, right)
                crossings.append(row)
                append_record(args.output, {"at": utc_now(), "kind": "d23_exact_crossing", **row})

    append_record(
        args.output,
        {
            "at": utc_now(),
            "kind": "d23_phase_search_complete",
            "profiles": len(sequences),
            "pairs": len(sequences) ** 2,
            "best": best,
            "exact_crossings": crossings,
        },
    )
    print(json.dumps({"profiles": len(sequences), "best": best, "crossings": crossings}, sort_keys=True))


if __name__ == "__main__":
    main()
