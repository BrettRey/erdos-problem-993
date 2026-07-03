#!/usr/bin/env python3
"""Exact rational checks for the binomial variance reserve lemma."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


def first_descent(n: int, p: Fraction) -> int | None:
    """First k with b_k < b_{k-1} for Bin(n,p), or None if none."""
    a = (n + 1) * p
    descent = a.numerator // a.denominator + 1
    if descent > n:
        return None
    return descent


def evaluate(n: int, p: Fraction) -> dict[str, Any] | None:
    q = 1 - p
    if not (0 < p < 1):
        return None
    descent = first_descent(n, p)
    if descent is None or descent >= n:
        return None

    ratio = Fraction(n - descent, descent + 1) * p / q
    reserve = 1 - ratio
    variance = n * p * q
    a = (n + 1) * p
    theta = descent - a
    exact_reserve = (theta + 1) / ((descent + 1) * q)
    exact_v_reserve = n * p * (theta + 1) / (descent + 1)
    lower_v_over_v_plus_3 = variance / (variance + 3)
    return {
        "n": n,
        "p": p,
        "variance": variance,
        "first_descent": descent,
        "ratio": ratio,
        "reserve": reserve,
        "theta": theta,
        "v_reserve": variance * reserve,
        "exact_reserve_ok": reserve == exact_reserve,
        "exact_v_reserve_ok": variance * reserve == exact_v_reserve,
        "v_over_v_plus_3_ok": variance * reserve >= lower_v_over_v_plus_3,
        "quarter_ok_when_v_ge_1": variance < 1 or variance * reserve >= Fraction(1, 4),
    }


def frac_to_json(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def compact(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "n": row["n"],
        "p": frac_to_json(row["p"]),
        "variance": frac_to_json(row["variance"]),
        "first_descent": row["first_descent"],
        "ratio": frac_to_json(row["ratio"]),
        "reserve": frac_to_json(row["reserve"]),
        "v_reserve": frac_to_json(row["v_reserve"]),
        "v_reserve_float": float(row["v_reserve"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=200)
    parser.add_argument("--max-den", type=int, default=100)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/binomial_variance_reserve_check_2026-07-03.json"),
    )
    args = parser.parse_args()
    if args.max_n < 2:
        parser.error("--max-n must be at least 2")
    if args.max_den < 2:
        parser.error("--max-den must be at least 2")

    rows = []
    failures: list[dict[str, Any]] = []
    for n in range(2, args.max_n + 1):
        for den in range(2, args.max_den + 1):
            for num in range(1, den):
                row = evaluate(n, Fraction(num, den))
                if row is None:
                    continue
                rows.append(row)
                if not (
                    row["exact_reserve_ok"]
                    and row["exact_v_reserve_ok"]
                    and row["v_over_v_plus_3_ok"]
                    and row["quarter_ok_when_v_ge_1"]
                ):
                    failures.append(row)

    v_ge_1 = [row for row in rows if row["variance"] >= 1]
    top_all = sorted(rows, key=lambda row: row["v_reserve"])[: args.top]
    top_v_ge_1 = sorted(v_ge_1, key=lambda row: row["v_reserve"])[: args.top]
    summary: dict[str, Any] = {
        "source": {
            "kind": "binomial_variance_reserve_exact_rational_check",
            "max_n": args.max_n,
            "max_den": args.max_den,
        },
        "processed": len(rows),
        "counts": {
            "variance_ge_1": len(v_ge_1),
            "failures": len(failures),
        },
        "checks": {
            "reserve_equals_theta_formula": all(row["exact_reserve_ok"] for row in rows),
            "v_reserve_equals_theta_formula": all(row["exact_v_reserve_ok"] for row in rows),
            "v_reserve_ge_v_over_v_plus_3": all(row["v_over_v_plus_3_ok"] for row in rows),
            "v_reserve_ge_one_quarter_when_v_ge_1": all(
                row["quarter_ok_when_v_ge_1"] for row in rows
            ),
        },
        "top_smallest_v_reserve": [compact(row) for row in top_all],
        "top_smallest_v_reserve_with_v_ge_1": [compact(row) for row in top_v_ge_1],
        "failures": [compact(row) for row in failures[: args.top]],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(rows)}, failures={len(failures)}",
        flush=True,
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
