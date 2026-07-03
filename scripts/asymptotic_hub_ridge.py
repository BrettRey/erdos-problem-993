#!/usr/bin/env python3
"""Closed-form scan of hub-bouquet crossing pressure.

For a hub with ``s`` pendant leaves and path arms of lengths
``a_1, ..., a_m`` beyond the hub, the independence polynomial is

    (1+x)^s prod_j I(P_{a_j}) + x prod_j I(P_{a_j-1}).

This lets us scan the broom / finite-arm hub ridge at larger star sizes
without running tree DP on every vertex.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from decimal import Decimal, getcontext
from functools import lru_cache
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import independence_poly  # noqa: E402
from scripts.scan_crossing_families import make_multiarm_star, parse_int_values  # noqa: E402
from targeted import make_broom  # noqa: E402


getcontext().prec = 60


def poly_add(a: list[int], b: list[int]) -> list[int]:
    size = max(len(a), len(b))
    out = [0] * size
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return trim(out)


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    if len(a) < len(b):
        short, long = a, b
    else:
        short, long = b, a
    for i, c in enumerate(short):
        if c == 0:
            continue
        for j, d in enumerate(long):
            out[i + j] += c * d
    return trim(out)


def trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


@lru_cache(maxsize=None)
def path_poly(length: int) -> tuple[int, ...]:
    if length < 0:
        return (1,)
    if length == 0:
        return (1,)
    if length == 1:
        return (1, 1)
    prev2 = [1]
    prev1 = [1, 1]
    for _ in range(2, length + 1):
        current = poly_add(prev1, [0, *prev2])
        prev2, prev1 = prev1, current
    return tuple(prev1)


@lru_cache(maxsize=None)
def arm_products(arms: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    exclude_hub = [1]
    include_hub = [1]
    for arm in arms:
        exclude_hub = poly_mul(exclude_hub, list(path_poly(arm)))
        include_hub = poly_mul(include_hub, list(path_poly(arm - 1)))
    return tuple(exclude_hub), tuple(include_hub)


@lru_cache(maxsize=None)
def binomial_coeffs(s: int) -> tuple[int, ...]:
    if s < 0:
        raise ValueError("s must be nonnegative")
    coeffs = [1] * (s + 1)
    value = 1
    for k in range(1, s + 1):
        value = value * (s - k + 1) // k
        coeffs[k] = value
    return tuple(coeffs)


def hub_bouquet_poly(s: int, arms: tuple[int, ...]) -> list[int]:
    exclude_hub, include_hub = arm_products(arms)
    leaf_term = poly_mul(list(binomial_coeffs(s)), list(exclude_hub))
    hub_term = [0, *include_hub]
    return poly_add(leaf_term, hub_term)


def first_descent(poly: list[int]) -> int | None:
    for k in range(1, len(poly)):
        if poly[k] < poly[k - 1]:
            return k
    return None


def max_post_descent_ratio(poly: list[int]) -> tuple[int, int, int | None, int, int | None]:
    descent = first_descent(poly)
    if descent is None:
        return 0, 1, None, 0, None
    best_num = 0
    best_den = 1
    best_pos: int | None = None
    upward_count = 0
    first_upward: int | None = None
    for j in range(descent, len(poly) - 1):
        den = poly[j]
        num = poly[j + 1]
        if den == 0:
            continue
        if num > den:
            upward_count += 1
            if first_upward is None:
                first_upward = j
        if num * best_den > best_num * den:
            best_num = num
            best_den = den
            best_pos = j
    return best_num, best_den, best_pos, upward_count, first_upward


def decimal_ratio(num: int, den: int) -> Decimal:
    if den == 0:
        return Decimal("Infinity")
    return Decimal(num) / Decimal(den)


def evaluate(family: str, s: int, arms: tuple[int, ...]) -> dict[str, Any]:
    poly = hub_bouquet_poly(s, arms)
    num, den, pos, upward_count, first_upward = max_post_descent_ratio(poly)
    ratio = decimal_ratio(num, den)
    reserve = Decimal(1) - ratio
    n = 1 + s + sum(arms)
    return {
        "family": family,
        "label": label_for(family, s, arms),
        "params": {"s": s, "arms": list(arms)},
        "n": n,
        "alpha": len(poly) - 1,
        "first_descent": first_descent(poly),
        "crossing_pos": pos,
        "crossing_pressure": float(ratio),
        "crossing_pressure_decimal": format(ratio, "f"),
        "crossing_reserve": float(reserve),
        "crossing_reserve_decimal": format(reserve, "f"),
        "s_times_reserve": float(Decimal(s) * reserve),
        "s_times_reserve_decimal": format(Decimal(s) * reserve, "f"),
        "n_times_reserve": float(Decimal(n) * reserve),
        "n_times_reserve_decimal": format(Decimal(n) * reserve, "f"),
        "post_descent_upward_count": upward_count,
        "first_upward_pos": first_upward,
    }


def label_for(family: str, s: int, arms: tuple[int, ...]) -> str:
    if family == "broom":
        return f"broom({arms[0] + 1},{s})"
    return f"M({s};{','.join(map(str, arms))})"


def arm_tuples(values: list[int], max_arms: int) -> list[tuple[int, ...]]:
    tuples: list[tuple[int, ...]] = []
    for size in range(2, max_arms + 1):
        tuples.extend(itertools.combinations_with_replacement(values, size))
    return tuples


def parse_s_values(spec: str) -> list[int]:
    values = parse_int_values(spec)
    if any(value < 1 for value in values):
        raise ValueError("all s values must be positive")
    return values


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "family",
        "label",
        "params",
        "n",
        "alpha",
        "first_descent",
        "crossing_pos",
        "crossing_pressure",
        "crossing_reserve",
        "s_times_reserve",
        "n_times_reserve",
        "post_descent_upward_count",
        "first_upward_pos",
    ]
    return {key: row[key] for key in keys}


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    best = max(rows, key=lambda row: row["crossing_pressure"])
    return {
        "processed": len(rows),
        "counts": {
            "post_descent_upward_rows": sum(
                1 for row in rows if row["post_descent_upward_count"] > 0
            ),
            "broom_rows": sum(1 for row in rows if row["family"] == "broom"),
            "multiarm_star_rows": sum(1 for row in rows if row["family"] == "multiarm_star"),
        },
        "numeric": {
            "max_crossing_pressure": best["crossing_pressure"],
            "min_crossing_reserve": best["crossing_reserve"],
            "min_s_times_reserve": min(row["s_times_reserve"] for row in rows),
            "max_s_times_reserve": max(row["s_times_reserve"] for row in rows),
        },
    }


def build_rows(
    *,
    s_values: list[int],
    broom_arms: list[int],
    arm_values: list[int],
    max_arms: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    multiarm_grid = arm_tuples(arm_values, max_arms)
    for s in s_values:
        for arm in broom_arms:
            rows.append(evaluate("broom", s, (arm,)))
        for arms in multiarm_grid:
            rows.append(evaluate("multiarm_star", s, arms))
    return rows


def self_check() -> None:
    cases = [
        (1, (1,)),
        (5, (4,)),
        (7, (2, 3, 6)),
        (12, (2, 4, 5)),
    ]
    for s, arms in cases:
        formula = hub_bouquet_poly(s, arms)
        if len(arms) == 1:
            n, adj = make_broom(arms[0] + 1, s)
        else:
            n, adj = make_multiarm_star(s, arms)
        dp_poly = independence_poly(n, adj)
        if formula != dp_poly:
            raise AssertionError(f"closed form mismatch for s={s}, arms={arms}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s-values", default="100,200,500,1000,2000,5000")
    parser.add_argument("--broom-max-arm", type=int, default=120)
    parser.add_argument(
        "--broom-arm-values",
        default=None,
        help="optional comma/range spec overriding --broom-max-arm",
    )
    parser.add_argument("--arm-values", default="2-10")
    parser.add_argument("--max-arms", type=int, default=3)
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/asymptotic_hub_ridge_2026-07-03.json"),
    )
    args = parser.parse_args()
    if args.broom_max_arm < 1:
        parser.error("--broom-max-arm must be positive")
    if args.max_arms < 2:
        parser.error("--max-arms must be at least 2")

    if args.self_check:
        self_check()

    s_values = parse_s_values(args.s_values)
    broom_arms = (
        parse_int_values(args.broom_arm_values)
        if args.broom_arm_values is not None
        else list(range(1, args.broom_max_arm + 1))
    )
    if any(value < 1 for value in broom_arms):
        parser.error("all broom arm values must be positive")
    arm_values = parse_int_values(args.arm_values)
    rows = build_rows(
        s_values=s_values,
        broom_arms=broom_arms,
        arm_values=arm_values,
        max_arms=args.max_arms,
    )
    top_rows = sorted(rows, key=lambda row: row["crossing_pressure"], reverse=True)[: args.top]
    best_by_s = [
        compact(max((row for row in rows if row["params"]["s"] == s), key=lambda row: row["crossing_pressure"]))
        for s in s_values
    ]
    tracked_labels = {"M(488;2,3,6)", "M(5000;2,3,6)"}
    tracked_rows = [
        compact(row)
        for row in rows
        if tuple(row["params"]["arms"]) == (2, 3, 6) or row["label"] in tracked_labels
    ]
    summary = {
        "source": {
            "kind": "asymptotic_hub_ridge",
            "s_values": s_values,
            "broom_max_arm": args.broom_max_arm,
            "broom_arm_values": broom_arms,
            "arm_values": arm_values,
            "max_arms": args.max_arms,
            "formula": "(1+x)^s prod I(P_a) + x prod I(P_{a-1})",
        },
        **summarize(rows),
        "best_by_s": best_by_s,
        "top_by_crossing_pressure": [compact(row) for row in top_rows],
        "tracked_M_s_2_3_6": tracked_rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Scanned {summary['processed']} hub-bouquet rows; "
        f"post-descent-upward={summary['counts']['post_descent_upward_rows']}, "
        f"best-pressure={summary['numeric']['max_crossing_pressure']}, "
        f"wrote {args.out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
