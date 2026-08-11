#!/usr/bin/env python3
"""Exact effective threshold for fixed-arm multi-arm stars.

For an arm vector ``a = (a_1, ..., a_m)``, the multi-arm star with ``s``
additional hub leaves has independence polynomial

    F_s(x) = (1+x)^s Q(x) + x R(x),
    Q(x) = product I(P_{a_j}; x),
    R(x) = product I(P_{a_j-1}; x).

This module computes an explicit threshold above which ``F_s`` is proved
unimodal.  The proof uses only integer inequalities; the finite-range scan is
reported separately as computational evidence.
"""

from __future__ import annotations

import argparse
import json
import sys
from math import comb
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from indpoly import independence_poly, is_unimodal


def poly_add(a: list[int], b: list[int]) -> list[int]:
    out = [0] * max(len(a), len(b))
    for index, value in enumerate(a):
        out[index] += value
    for index, value in enumerate(b):
        out[index] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def path_independence_poly(order: int) -> list[int]:
    """Return ``I(P_order; x)``; the order-zero path has polynomial one."""
    if order < 0:
        raise ValueError("path order must be nonnegative")
    if order == 0:
        return [1]
    if order == 1:
        return [1, 1]
    previous_two = [1]
    previous_one = [1, 1]
    for _ in range(2, order + 1):
        current = poly_add(previous_one, [0] + previous_two)
        previous_two, previous_one = previous_one, current
    return previous_one


def fixed_arm_terms(arms: tuple[int, ...]) -> tuple[list[int], list[int]]:
    """Return ``Q`` and ``e=xR`` for a positive arm vector."""
    if not arms or any(arm < 1 for arm in arms):
        raise ValueError("arms must be a nonempty tuple of positive integers")
    q = [1]
    r_poly = [1]
    for arm in arms:
        q = poly_mul(q, path_independence_poly(arm))
        r_poly = poly_mul(r_poly, path_independence_poly(arm - 1))
    return q, [0] + r_poly


def fixed_arm_polynomial(s: int, arms: tuple[int, ...]) -> list[int]:
    """Compute ``F_s`` exactly from its include/exclude decomposition."""
    if s < 0:
        raise ValueError("s must be nonnegative")
    q, e = fixed_arm_terms(arms)
    binomial = [comb(s, k) for k in range(s + 1)]
    return poly_add(poly_mul(binomial, q), e)


def multi_arm_adjacency(
    s: int, arms: tuple[int, ...],
) -> list[list[int]]:
    """Construct ``M(s; arms)`` for independent tree-DP replay."""
    n = 1 + s + sum(arms)
    adjacency: list[list[int]] = [[] for _ in range(n)]

    def add_edge(u: int, v: int) -> None:
        adjacency[u].append(v)
        adjacency[v].append(u)

    next_vertex = 1
    for _ in range(s):
        add_edge(0, next_vertex)
        next_vertex += 1
    for arm in arms:
        previous = 0
        for _ in range(arm):
            add_edge(previous, next_vertex)
            previous = next_vertex
            next_vertex += 1
    return adjacency


def d_coefficient(s: int, k: int, q: list[int]) -> int:
    """Coefficient ``[x^k](1+x)^s Q`` without building the full row."""
    total = 0
    for j, qj in enumerate(q):
        binomial_index = k - j
        if 0 <= binomial_index <= s:
            total += qj * comb(s, binomial_index)
    return total


def early_difference(
    s: int, k: int, q: list[int], e: list[int],
) -> int:
    """Return ``[x^(k+1)]F_s - [x^k]F_s`` exactly."""
    e_k = e[k] if k < len(e) else 0
    e_next = e[k + 1] if k + 1 < len(e) else 0
    return d_coefficient(s, k + 1, q) - d_coefficient(s, k, q) + e_next - e_k


def binomial_difference_lower_bound(s: int, k: int) -> int:
    """The ``q_0=1`` contribution to the kth difference of ``(1+x)^sQ``."""
    return comb(s, k + 1) - comb(s, k)


def certified_threshold(e: list[int]) -> tuple[int, int, list[int]]:
    """Return a sufficient threshold, a crude cap, and correction drops.

    Let ``r=deg(e)``.  For ``s >= 2r+1``, every contribution of a
    nonconstant coefficient of ``Q`` to each difference through index ``r``
    is nonnegative.  It therefore suffices that the binomial contribution
    dominate ``max(0, e_k-e_{k+1})`` at every such index.
    """
    r = len(e) - 1
    drops = [
        max(0, e[k] - (e[k + 1] if k + 1 < len(e) else 0))
        for k in range(r + 1)
    ]
    max_drop = max(drops)
    lower = 2 * r + 1
    crude = lower + (r + 1) * max_drop

    def sufficient(s: int) -> bool:
        return all(
            binomial_difference_lower_bound(s, k) >= drops[k]
            for k in range(r + 1)
        )

    lo, hi = lower, crude
    while lo < hi:
        middle = (lo + hi) // 2
        if sufficient(middle):
            hi = middle
        else:
            lo = middle + 1
    if not sufficient(lo):
        raise AssertionError("the proved crude threshold did not suffice")
    return lo, crude, drops


def analyze_arms(
    arms: tuple[int, ...], finite_scan_cap: int,
) -> dict[str, object]:
    q, e = fixed_arm_terms(arms)
    r = len(e) - 1
    threshold, crude, drops = certified_threshold(e)

    replay_s = [threshold, threshold + 1, threshold + 17]
    replay = []
    for s in replay_s:
        differences = [early_difference(s, k, q, e) for k in range(r + 1)]
        lower_margins = [
            binomial_difference_lower_bound(s, k) - drops[k]
            for k in range(r + 1)
        ]
        actual_minus_lower = [
            differences[k] - lower_margins[k]
            for k in range(r + 1)
        ]
        if (
            min(differences) < 0
            or min(lower_margins) < 0
            or min(actual_minus_lower) < 0
        ):
            raise AssertionError("effective-threshold replay failed")
        replay.append({
            "s": s,
            "minimum_actual_early_difference": min(differences),
            "minimum_proof_lower_margin": min(lower_margins),
            "minimum_actual_minus_proof_lower_bound": min(actual_minus_lower),
        })

    dp_crosschecks = []
    for s in sorted({0, 1, 2, 7}):
        formula = fixed_arm_polynomial(s, arms)
        adjacency = multi_arm_adjacency(s, arms)
        direct = independence_poly(len(adjacency), adjacency)
        if formula != direct:
            raise AssertionError(f"tree-DP mismatch for arms={arms}, s={s}")
        dp_crosschecks.append({"s": s, "coefficient_count": len(formula)})

    finite_limit = min(threshold - 1, finite_scan_cap)
    nonunimodal = []
    early_failure_count = 0
    last_early_failure = -1
    for s in range(finite_limit + 1):
        coefficients = fixed_arm_polynomial(s, arms)
        if not is_unimodal(coefficients):
            nonunimodal.append({"s": s, "coefficients": coefficients})
        if any(
            early_difference(s, k, q, e) < 0
            for k in range(r + 1)
        ):
            early_failure_count += 1
            last_early_failure = s

    finite_complete = finite_limit == threshold - 1
    return {
        "arms": list(arms),
        "q": q,
        "x_r": e,
        "correction_degree": r,
        "correction_drops": drops,
        "certified_threshold": threshold,
        "crude_threshold": crude,
        "proof_replay": replay,
        "tree_dp_crosschecks": dp_crosschecks,
        "finite_scan": {
            "range": [0, finite_limit],
            "complete_below_threshold": finite_complete,
            "nonunimodal_count": len(nonunimodal),
            "nonunimodal_witnesses": nonunimodal[:5],
            "early_difference_failure_count": early_failure_count,
            "last_early_difference_failure": last_early_failure,
        },
        "all_s_unimodal_certified": finite_complete and not nonunimodal,
    }


def parse_arms(raw: str) -> tuple[int, ...]:
    try:
        arms = tuple(int(part) for part in raw.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error
    if not arms or any(arm < 1 for arm in arms):
        raise argparse.ArgumentTypeError("use comma-separated positive arm lengths")
    return arms


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--arms",
        action="append",
        type=parse_arms,
        help="comma-separated arm vector; repeat for multiple vectors",
    )
    parser.add_argument("--finite-scan-cap", type=int, default=100_000)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/fixed_arm_unimodality_threshold_20260808.json"),
    )
    args = parser.parse_args()
    if args.finite_scan_cap < 0:
        parser.error("--finite-scan-cap must be nonnegative")

    arm_vectors = args.arms or [(5, 5, 4, 2), (2, 3, 6), (6, 6, 1), (8,)]
    records = [analyze_arms(arms, args.finite_scan_cap) for arms in arm_vectors]
    report = {
        "claim_scope": (
            "Theorem: each fixed arm vector is unimodal for all s at or above "
            "its certified threshold. Finite scans below the threshold are "
            "computational certificates, reported separately. This does not "
            "resolve Erdos Problem 993."
        ),
        "method": "exact integer arithmetic; independent tree-DP spot checks",
        "records": records,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "event": "complete",
        "out": str(args.out),
        "records": [
            {
                "arms": record["arms"],
                "certified_threshold": record["certified_threshold"],
                "all_s_unimodal_certified": record["all_s_unimodal_certified"],
            }
            for record in records
        ],
    }))


if __name__ == "__main__":
    main()
