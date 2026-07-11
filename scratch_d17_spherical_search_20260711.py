#!/usr/bin/env python3
"""Search balanced double spherical trees for a genuine coefficient valley.

A branching sequence is read from the leaves toward the root.  Starting with
the rooted leaf states E=1 and S=x, one parent layer of branching b gives

    E <- (E+S)^b,       S <- x E_old^b.

Two copies are joined at their roots, so I=E^2+2ES.  Floating arithmetic is
used only to rank candidates; any crossing must be replayed with exact integer
polynomials before it is a certificate.
"""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass

import numpy as np
from scipy.signal import fftconvolve


@dataclass(frozen=True)
class Scaled:
    values: np.ndarray
    log_scale: float


def normalize(item: Scaled) -> Scaled:
    maximum = float(np.max(item.values))
    if maximum <= 0:
        return item
    return Scaled(item.values / maximum, item.log_scale + math.log(maximum))


def add(a: Scaled, b: Scaled, a_weight: float = 1.0, b_weight: float = 1.0) -> Scaled:
    a_log = a.log_scale + math.log(a_weight)
    b_log = b.log_scale + math.log(b_weight)
    common = max(a_log, b_log)
    out = np.zeros(max(len(a.values), len(b.values)))
    out[: len(a.values)] += a.values * math.exp(a_log - common)
    out[: len(b.values)] += b.values * math.exp(b_log - common)
    return normalize(Scaled(out, common))


def multiply(a: Scaled, b: Scaled) -> Scaled:
    values = fftconvolve(a.values, b.values)
    values = np.maximum(values, 0.0)
    return normalize(Scaled(values, a.log_scale + b.log_scale))


def power(a: Scaled, exponent: int) -> Scaled:
    result = Scaled(np.array([1.0]), 0.0)
    while exponent:
        if exponent & 1:
            result = multiply(result, a)
        exponent //= 2
        if exponent:
            a = multiply(a, a)
    return result


def spherical_states(branching: tuple[int, ...]) -> tuple[Scaled, Scaled, int]:
    excluded = Scaled(np.array([1.0]), 0.0)
    selected = Scaled(np.array([0.0, 1.0]), 0.0)
    rooted_order = 1
    for branch in branching:
        arbitrary = add(excluded, selected)
        old_excluded = excluded
        excluded = power(arbitrary, branch)
        selected_base = power(old_excluded, branch)
        selected = Scaled(np.concatenate((np.array([0.0]), selected_base.values)), selected_base.log_scale)
        rooted_order = 1 + branch * rooted_order
    return excluded, selected, rooted_order


def joined_polynomial(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[np.ndarray, int]:
    left_excluded, left_selected, left_order = spherical_states(left)
    right_excluded, right_selected, right_order = spherical_states(right)
    full = add(
        multiply(left_excluded, right_excluded),
        add(
            multiply(left_selected, right_excluded),
            multiply(left_excluded, right_selected),
        ),
    )
    probability = np.maximum(full.values, 0.0)
    probability /= probability.sum()
    return probability, left_order + right_order


def spherical_polynomial(branching: tuple[int, ...]) -> tuple[np.ndarray, int]:
    return joined_polynomial(branching, branching)


def spherical_order(branching: tuple[int, ...]) -> int:
    rooted_order = 1
    for branch in branching:
        rooted_order = 1 + branch * rooted_order
    return 2 * rooted_order


def pressure(
    poly: np.ndarray,
    relative_floor: float = 1e-12,
    min_separation: int = 1,
) -> dict[str, float | int | None]:
    support = np.flatnonzero(poly > float(np.max(poly)) * relative_floor)
    if len(support) < 3:
        return {"first": None, "best_at": None, "best": 0.0, "first_ratio": None}
    first = None
    first_ratio = None
    best = 0.0
    best_at = None
    for k in range(int(support[0]), int(support[-1])):
        ratio = float(poly[k + 1] / poly[k])
        if first is None and ratio < 1.0 - 1e-9:
            first = k
            first_ratio = ratio
        elif first is not None and k >= first + min_separation and ratio > best:
            best = ratio
            best_at = k
    return {"first": first, "best_at": best_at, "best": best, "first_ratio": first_ratio}


def candidates(max_depth: int, max_branch: int) -> set[tuple[int, ...]]:
    out: set[tuple[int, ...]] = set()
    for depth in range(4, max_depth + 1):
        for split in range(1, depth):
            for left in range(1, max_branch + 1):
                for right in range(1, max_branch + 1):
                    out.add((left,) * split + (right,) * (depth - split))
                    out.add(tuple(left if j % 2 == 0 else right for j in range(depth)))
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-depth", type=int, default=12)
    parser.add_argument("--max-branch", type=int, default=8)
    parser.add_argument("--max-order", type=int, default=100_000)
    parser.add_argument("--random", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--min-separation", type=int, default=1)
    parser.add_argument("--asymmetric-pairs", type=int, default=0)
    args = parser.parse_args()

    rows = candidates(args.max_depth, args.max_branch)
    rng = random.Random(args.seed)
    for _ in range(args.random):
        depth = rng.randint(4, args.max_depth)
        rows.add(tuple(rng.randint(1, args.max_branch) for _ in range(depth)))

    if args.asymmetric_pairs:
        eligible = [row for row in rows if spherical_order(row) // 2 <= args.max_order]
        best_pairs = []
        for _ in range(args.asymmetric_pairs):
            left = rng.choice(eligible)
            right = rng.choice(eligible)
            if spherical_order(left) // 2 + spherical_order(right) // 2 > args.max_order:
                continue
            poly, order = joined_polynomial(left, right)
            result = pressure(poly, min_separation=args.min_separation)
            score = float(result["best"])
            best_pairs.append((score, left, right, order, result))
            if score > 1.0 + 1e-7:
                print(
                    {
                        "candidate_crossing": [left, right],
                        "order": order,
                        **result,
                    },
                    flush=True,
                )
                return
        best_pairs.sort(reverse=True)
        for row in best_pairs[:30]:
            print(row)
        return

    best: list[tuple[float, tuple[int, ...], int, dict[str, float | int | None]]] = []
    for branching in rows:
        if spherical_order(branching) > args.max_order:
            continue
        poly, order = spherical_polynomial(branching)
        result = pressure(poly, min_separation=args.min_separation)
        score = float(result["best"])
        best.append((score, branching, order, result))
        if score > 1.0 + 1e-7:
            print({"candidate_crossing": branching, "order": order, **result}, flush=True)
            return
    best.sort(reverse=True)
    for row in best[:30]:
        print(row)


if __name__ == "__main__":
    main()
