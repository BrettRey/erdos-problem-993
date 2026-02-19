#!/usr/bin/env python3
"""Exact rational check for the five Sub-claim A helper exceptions.

Context:
  The surrogate helper inequality
    lb_u2(k,j) := mu_{k,j+2}(u2_{k,j}) - mu_{k,j}(lambda_{k,j}) >= 1
  fails on exactly five small pairs in scans:
    (6,0), (6,2), (7,1), (8,0), (8,2).

Sub-claim A target needs instead:
  Delta_total(k,j) :=
    [margin(k,j+2) + (mode(k,j+2)-1)] - [margin(k,j) + (mode(k,j)-1)] >= 1.

This script verifies Delta_total > 1 exactly (Fraction arithmetic)
on those five pairs.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from typing import Any

from conjecture_a_mixed_spider_exact_margin import compute_margin

EXCEPTIONS: list[tuple[int, int]] = [(6, 0), (6, 2), (7, 1), (8, 0), (8, 2)]


def main() -> None:
    ap = argparse.ArgumentParser(description="Exact check on five Sub-claim A exception pairs.")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    records: list[dict[str, Any]] = []
    bad_mode_step: list[dict[str, Any]] = []
    bad_delta: list[dict[str, Any]] = []

    min_delta_minus_one: Fraction | None = None
    min_witness: dict[str, int] | None = None

    for k, j in EXCEPTIONS:
        d0 = compute_margin(k, j)
        d2 = compute_margin(k, j + 2)

        mode_step = d2.mode - d0.mode
        delta_total = (d2.margin + Fraction(d2.mode - 1, 1)) - (
            d0.margin + Fraction(d0.mode - 1, 1)
        )
        delta_minus_one = delta_total - Fraction(1, 1)

        rec = {
            "k": k,
            "j": j,
            "mode_j": d0.mode,
            "mode_j2": d2.mode,
            "mode_step": mode_step,
            "lambda_j": str(d0.lam),
            "lambda_j2": str(d2.lam),
            "delta_total": str(delta_total),
            "delta_total_float": float(delta_total),
            "delta_minus_one": str(delta_minus_one),
            "delta_minus_one_float": float(delta_minus_one),
            "delta_gt_one": bool(delta_total > 1),
        }
        records.append(rec)

        if mode_step != 1:
            bad_mode_step.append(rec)
        if not (delta_total > 1):
            bad_delta.append(rec)
        if min_delta_minus_one is None or delta_minus_one < min_delta_minus_one:
            min_delta_minus_one = delta_minus_one
            min_witness = {"k": k, "j": j}

    payload = {
        "exceptions": EXCEPTIONS,
        "records": records,
        "summary": {
            "count": len(EXCEPTIONS),
            "mode_step_fail_count": len(bad_mode_step),
            "delta_gt_one_fail_count": len(bad_delta),
            "min_delta_minus_one": str(min_delta_minus_one),
            "min_delta_minus_one_float": float(min_delta_minus_one)
            if min_delta_minus_one is not None
            else None,
            "min_witness": min_witness,
        },
    }

    print("Sub-claim A: exact check on five helper-exception pairs")
    print(f"pairs={len(EXCEPTIONS)}")
    print(f"mode_step_fail_count={len(bad_mode_step)}")
    print(f"delta_gt_one_fail_count={len(bad_delta)}")
    if min_delta_minus_one is not None and min_witness is not None:
        print(
            "min(delta_total-1)="
            f"{float(min_delta_minus_one):.15f} at (k,j)=({min_witness['k']},{min_witness['j']})"
        )
    for rec in records:
        print(
            f"(k,j)=({rec['k']},{rec['j']}): "
            f"delta_total={rec['delta_total_float']:.15f}, "
            f"delta_minus_one={rec['delta_minus_one_float']:.15f}"
        )

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f_out:
            json.dump(payload, f_out, indent=2)
        print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
