#!/usr/bin/env python3
"""Tilted-convolution frontier for the sharp D25 two-species family.

The hub has one parent branch G(40,20) and m identical child branches
TG(4,6).  For each m this evaluates the oriented fork ratio at every prefix
rank visible under a grid of coefficient-centering exponential tilts.  The
tilt multiplies every rank-r coefficient by the same z^r and therefore leaves
the target ratio unchanged while preventing floating underflow.
"""

from __future__ import annotations

import math

import numpy as np

from scratch_d17_spherical_search_20260711 import Scaled, add, multiply, normalize, power
from scratch_d25_hard_branch_composition_search_20260711 import library


def tilted(poly: tuple[int, ...], log_z: float) -> Scaled:
    logs = np.array(
        [math.log(value) + rank * log_z if value else -math.inf
         for rank, value in enumerate(poly)],
        dtype=float,
    )
    maximum = float(np.max(logs))
    values = np.exp(logs - maximum)
    values[~np.isfinite(logs)] = 0.0
    return normalize(Scaled(values, maximum))


def shift(item: Scaled, log_z: float) -> Scaled:
    # Multiplication by x after the z^rank tilt contributes one extra z.
    return Scaled(np.pad(item.values, (1, 0)), item.log_scale + log_z)


def weighted(item: Scaled, weight: int) -> Scaled:
    return Scaled(item.values, item.log_scale + math.log(weight))


def mean(poly: tuple[int, ...], log_z: float) -> float:
    logs = np.array(
        [math.log(value) + rank * log_z if value else -math.inf
         for rank, value in enumerate(poly)],
        dtype=float,
    )
    maximum = float(np.max(logs))
    weights = np.exp(logs - maximum)
    ranks = np.arange(len(poly), dtype=float)
    return float(np.dot(ranks, weights) / np.sum(weights))


def centering_log_z(parent, child, multiplicity: int, target: float) -> float:
    lo, hi = -12.0, 12.0
    for _ in range(70):
        mid = (lo + hi) / 2
        value = mean(parent.total, mid) + multiplicity * mean(child.total, mid)
        if value < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def ratio_rows(parent, child, multiplicity: int, log_z: float, max_rank: int):
    p0, e0, r0 = (tilted(poly, log_z) for poly in
                  (parent.total, parent.excluded, parent.reserve))
    p1, e1, r1 = (tilted(poly, log_z) for poly in
                  (child.total, child.excluded, child.reserve))
    p1_m = power(p1, multiplicity)
    e1_m = power(e1, multiplicity)
    whole = add(multiply(p0, p1_m), shift(multiply(e0, e1_m), log_z))
    a_c = multiply(e0, e1_m)
    a_p = multiply(r0, p1_m)
    p1_m1 = power(p1, multiplicity - 1)
    a_u = multiply(multiply(r1, p0), p1_m1)
    j_pc = multiply(multiply(r0, r1), p1_m1)
    j_cc = multiply(
        multiply(multiply(r1, r1), p0),
        power(p1, multiplicity - 2),
    )
    joint = add(
        weighted(j_pc, 2 * multiplicity),
        weighted(j_cc, multiplicity * (multiplicity - 1)),
    )
    child_sum = weighted(a_u, multiplicity)
    augmented = add(child_sum, weighted(add(a_p, a_c), 2))
    arrays = (whole, joint, child_sum, augmented)
    limit = min(max_rank + 1, *(len(item.values) for item in arrays))
    floors = [float(np.max(item.values)) * 1e-11 for item in arrays]
    log_scale = (
        whole.log_scale + joint.log_scale
        - child_sum.log_scale - augmented.log_scale
    )
    out = []
    for rank in range(limit):
        values = [item.values[rank] for item in arrays]
        if any(value <= floor for value, floor in zip(values, floors)):
            continue
        log_ratio = log_scale + math.log(values[0]) + math.log(values[1]) \
            - math.log(values[2]) - math.log(values[3])
        out.append((math.exp(log_ratio), rank))
    return out


def derivative_stats(poly: tuple[int, ...], reserve: tuple[int, ...], log_z: float):
    z = math.exp(log_z)
    ranks_p = np.arange(len(poly), dtype=float)
    coeff_p = np.array([float(v / sum(poly)) for v in poly])
    weights = coeff_p * np.exp(ranks_p * log_z - np.max(ranks_p * log_z))
    weights /= weights.sum()
    mu = float(np.dot(ranks_p, weights))
    variance = float(np.dot((ranks_p - mu) ** 2, weights))

    def log_eval(seq: tuple[int, ...], t: float) -> float:
        ranks = np.arange(len(seq), dtype=float)
        logs = np.array([
            math.log(v) + rank * t if v else -math.inf
            for rank, v in enumerate(seq)
        ])
        top = float(np.max(logs))
        return top + math.log(float(np.sum(np.exp(logs - top))))

    eps = 1e-5
    s_plus = log_eval(reserve, log_z + eps) - log_eval(poly, log_z + eps)
    s_minus = log_eval(reserve, log_z - eps) - log_eval(poly, log_z - eps)
    slope = (s_plus - s_minus) / (2 * eps)
    leading = -1.0 - slope * slope / variance
    return {"z": z, "mu_child": mu, "variance_child": variance,
            "log_q_slope": slope, "predicted_m_times_ratio_minus_1": leading}


def main() -> None:
    states = library()
    parent = next(item for item in states if item.label == "G_40_20@0")
    child = next(item for item in states if item.label == "TG_4_6@0")
    multiplicities = (16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024)
    for multiplicity in multiplicities:
        alpha = max(
            len(parent.total) - 1 + multiplicity * (len(child.total) - 1),
            1 + len(parent.excluded) - 1
            + multiplicity * (len(child.excluded) - 1),
        )
        limit = math.ceil((2 * alpha - 1) / 3)
        targets = np.linspace(max(4, 0.42 * alpha), limit - 2, 28)
        log_tilts = {
            round(centering_log_z(parent, child, multiplicity, float(target)), 10)
            for target in targets
        }
        best = (0.0, -1, 0.0)
        seen: set[int] = set()
        for log_z in sorted(log_tilts):
            for ratio, rank in ratio_rows(
                parent, child, multiplicity, log_z, limit - 2
            ):
                if rank in seen and ratio <= best[0]:
                    continue
                seen.add(rank)
                if ratio > best[0]:
                    best = (ratio, rank, log_z)
        stats = derivative_stats(child.total, child.reserve, best[2])
        print(
            {
                "m": multiplicity,
                "alpha": alpha,
                "limit": limit,
                "best_ratio": best[0],
                "best_rank": best[1],
                "rank_over_alpha": best[1] / alpha,
                "m_times_ratio_minus_1": multiplicity * (best[0] - 1.0),
                "saddle": stats,
            },
            flush=True,
        )


if __name__ == "__main__":
    main()
