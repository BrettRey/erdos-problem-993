#!/usr/bin/env python3
"""Diagnostics for Sub-claim A on mixed spiders.

Sub-claim A target:
  margin(k,j+2) >= margin(k,j)   for all k>=6, j>=0,
where margin(k,j) is the tie-fugacity margin on S(2^k,1^j).

This script scans (k,j) ranges and reports:
  1) direct Sub-claim A failures,
  2) mode-step behavior m(k,j+2)-m(k,j),
  3) tie-fugacity monotonicity lambda_{j+2} >= lambda_j,
  4) a decomposition:
       full_delta = margin_{j+2} - margin_j
                  = core_delta + lift_delta,
     where
       core_delta = mu(k,j+2; lambda_j) - mu(k,j; lambda_j) - (m_{j+2}-m_j),
       lift_delta = [mu(k,j+2; lambda_{j+2}) - mu(k,j+2; lambda_j)].

The decomposition isolates the compensation mechanism: the fixed-lambda core term
can be negative, while lambda-lift is positive and appears to dominate.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from typing import Any

from conjecture_a_mixed_spider_combined_scan import (
    build_g_coeffs_shift_binom,
    build_h_coeffs_2k,
    family_mu,
    mode_lambda_from_fg,
    next_f_times_1px,
)
from conjecture_a_mixed_spider_exact_margin import compute_margin


def maybe_record(items: list[dict[str, Any]], item: dict[str, Any], cap: int) -> None:
    if len(items) < cap:
        items.append(item)


def to_float(x: Any) -> float:
    return float(x)


def main() -> None:
    ap = argparse.ArgumentParser(description="Sub-claim A diagnostics on mixed spiders.")
    ap.add_argument("--k-min", type=int, default=6)
    ap.add_argument("--k-max", type=int, default=3000)
    ap.add_argument("--j-min", type=int, default=0)
    ap.add_argument("--j-max", type=int, default=120)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--store-cap", type=int, default=200)
    ap.add_argument(
        "--exact-top",
        type=int,
        default=8,
        help="Exact-rational check on the tightest top-N full deltas (0 disables).",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.k_min < 6:
        raise ValueError("k-min must be >= 6 for Sub-claim A scope")
    if args.k_max < args.k_min:
        raise ValueError("k-max must be >= k-min")
    if args.j_min < 0:
        raise ValueError("j-min must be >= 0")
    if args.j_max < args.j_min + 2:
        raise ValueError("j-max must be >= j-min+2")
    if args.exact_top < 0:
        raise ValueError("exact-top must be >= 0")

    t0 = time.time()
    print(
        f"Sub-claim A diagnostics: k={args.k_min}..{args.k_max}, j={args.j_min}..{args.j_max}",
        flush=True,
    )

    bad_subclaim_a: list[dict[str, Any]] = []
    bad_mode_step: list[dict[str, Any]] = []
    bad_lambda_mono: list[dict[str, Any]] = []

    # Tight witnesses (tracked as minima).
    min_full: dict[str, Any] | None = None
    min_core: dict[str, Any] | None = None
    min_lift: dict[str, Any] | None = None

    # Keep a small top list of tight full deltas for optional exact checks.
    tight_full: list[dict[str, Any]] = []

    total_pairs = 0

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)
        f = h[:]

        modes: list[int] = []
        lambdas: list[float] = []
        margins: list[float] = []

        for j in range(0, args.j_max + 1):
            m, lam, _, _ = mode_lambda_from_fg(f, g)
            mu = family_mu(k, j, lam)
            margin = mu - (m - 1)

            modes.append(m)
            lambdas.append(lam)
            margins.append(margin)

            if j < args.j_max:
                f = next_f_times_1px(f)

        j_lo = args.j_min
        j_hi = args.j_max - 2

        for j in range(j_lo, j_hi + 1):
            total_pairs += 1

            m0 = modes[j]
            m2 = modes[j + 2]
            lam0 = lambdas[j]
            lam2 = lambdas[j + 2]
            dm = m2 - m0

            full_delta = margins[j + 2] - margins[j]
            core_delta = family_mu(k, j + 2, lam0) - family_mu(k, j, lam0) - dm
            lift_delta = full_delta - core_delta

            rec = {
                "k": k,
                "j": j,
                "mode_j": m0,
                "mode_j2": m2,
                "mode_step": dm,
                "lambda_j": lam0,
                "lambda_j2": lam2,
                "margin_j": margins[j],
                "margin_j2": margins[j + 2],
                "full_delta": full_delta,
                "core_delta": core_delta,
                "lift_delta": lift_delta,
            }

            # Sub-claim A failure.
            if full_delta < -args.tol:
                maybe_record(bad_subclaim_a, rec, args.store_cap)

            # Structural diagnostics.
            if dm != 1:
                maybe_record(bad_mode_step, rec, args.store_cap)
            if lam2 + args.tol < lam0:
                maybe_record(bad_lambda_mono, rec, args.store_cap)

            # Min trackers.
            if min_full is None or full_delta < min_full["full_delta"]:
                min_full = rec
            if min_core is None or core_delta < min_core["core_delta"]:
                min_core = rec
            if min_lift is None or lift_delta < min_lift["lift_delta"]:
                min_lift = rec

            # Maintain a tiny list of tight full deltas.
            if len(tight_full) < max(1, args.exact_top):
                tight_full.append(rec)
                tight_full.sort(key=lambda x: x["full_delta"])
            else:
                if full_delta < tight_full[-1]["full_delta"]:
                    tight_full[-1] = rec
                    tight_full.sort(key=lambda x: x["full_delta"])

        if (k - args.k_min) % 500 == 0 or k == args.k_max:
            cur = min_full["full_delta"] if min_full is not None else float("nan")
            print(
                f"k={k:5d}: badA={len(bad_subclaim_a)} badMode={len(bad_mode_step)} "
                f"badLam={len(bad_lambda_mono)} min_full={cur:.12g}",
                flush=True,
            )

    # Optional exact checks on the tightest floating witnesses.
    exact_checks: list[dict[str, Any]] = []
    if args.exact_top > 0 and tight_full:
        seen: set[tuple[int, int]] = set()
        for rec in tight_full:
            k = int(rec["k"])
            j = int(rec["j"])
            if (k, j) in seen:
                continue
            seen.add((k, j))
            d0 = compute_margin(k, j)
            d2 = compute_margin(k, j + 2)
            exact_delta = d2.margin - d0.margin
            exact_checks.append(
                {
                    "k": k,
                    "j": j,
                    "float_full_delta": rec["full_delta"],
                    "exact_full_delta": str(exact_delta),
                    "exact_full_delta_float": to_float(exact_delta),
                    "mode_j": d0.mode,
                    "mode_j2": d2.mode,
                    "lambda_j": str(d0.lam),
                    "lambda_j2": str(d2.lam),
                }
            )

    summary: dict[str, Any] = {
        "params": vars(args),
        "overall": {
            "total_pairs": total_pairs,
            "wall_s": time.time() - t0,
        },
        "fails": {
            "subclaim_a": {"count": len(bad_subclaim_a), "examples": bad_subclaim_a},
            "mode_step_not_1": {"count": len(bad_mode_step), "examples": bad_mode_step},
            "lambda_j2_lt_lambda_j": {"count": len(bad_lambda_mono), "examples": bad_lambda_mono},
        },
        "extrema": {
            "min_full_delta": min_full,
            "min_core_delta": min_core,
            "min_lift_delta": min_lift,
        },
        "tight_float_witnesses": tight_full,
        "exact_checks_top": exact_checks,
    }

    print("-" * 96, flush=True)
    print(f"total_pairs={total_pairs}", flush=True)
    print(f"Sub-claim A failures={len(bad_subclaim_a)}", flush=True)
    print(f"Mode-step failures={len(bad_mode_step)}", flush=True)
    print(f"Lambda-monotonic failures={len(bad_lambda_mono)}", flush=True)
    if min_full is not None:
        print(f"min full delta={min_full['full_delta']:.15f} at (k,j)=({min_full['k']},{min_full['j']})", flush=True)
    if min_core is not None:
        print(f"min core delta={min_core['core_delta']:.15f} at (k,j)=({min_core['k']},{min_core['j']})", flush=True)
    if min_lift is not None:
        print(f"min lift delta={min_lift['lift_delta']:.15f} at (k,j)=({min_lift['k']},{min_lift['j']})", flush=True)
    if exact_checks:
        print("Exact checks on tight witnesses:", flush=True)
        for ex in exact_checks:
            print(
                f"  (k,j)=({ex['k']},{ex['j']}): exact full delta={ex['exact_full_delta_float']:.15f}",
                flush=True,
            )
    print(f"wall={summary['overall']['wall_s']:.2f}s", flush=True)

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f_out:
            json.dump(summary, f_out, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
