#!/usr/bin/env python3
"""Fixed-lambda + lane-2 proof harness for Sub-claim A.

Sub-claim A (mixed spiders):
  margin(k,j+2) >= margin(k,j),  for k>=6, j>=0.

Equivalent (when mode step is +1):
  mu_{k,j+2}(lambda_{j+2}) - mu_{k,j}(lambda_j) >= 1.

This harness centers the requested lane:
  total_step = fixed_step + lambda_gain
  fixed_step = mu_{k,j+2}(lambda_j) - mu_{k,j}(lambda_j)
  lambda_gain = mu_{k,j+2}(lambda_{j+2}) - mu_{k,j+2}(lambda_j)

and uses lane #2 data:
  lambda_{j+2} >= lambda_j (and stronger proxy u2 >= lambda_j).

It reports:
  1) total-step failures (<1),
  2) lane-2 failures,
  3) fixed-lambda identity residuals,
  4) the lower bound
       lb_u2 := mu_{k,j+2}(u2) - mu_{k,j}(lambda_j),
     where u2 = a'_{m_j}/a'_{m_j+1} from (1+2x)^k(1+x)^{j+2}.
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


def maybe_store(dst: list[dict[str, Any]], item: dict[str, Any], cap: int) -> None:
    if len(dst) < cap:
        dst.append(item)


def main() -> None:
    ap = argparse.ArgumentParser(description="Fixed-lambda + lane-2 Sub-claim A harness.")
    ap.add_argument("--k-min", type=int, default=6)
    ap.add_argument("--k-max", type=int, default=1200)
    ap.add_argument("--j-max", type=int, default=120)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--store-cap", type=int, default=200)
    ap.add_argument("--exact-top", type=int, default=8)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.k_min < 6:
        raise ValueError("k-min must be >= 6")
    if args.k_max < args.k_min:
        raise ValueError("k-max must be >= k-min")
    if args.j_max < 2:
        raise ValueError("j-max must be >= 2")

    t0 = time.time()
    print(
        f"Sub-claim A mean-step lane: k={args.k_min}..{args.k_max}, j=0..{args.j_max}",
        flush=True,
    )

    bad_mode_step: list[dict[str, Any]] = []
    bad_lambda_step: list[dict[str, Any]] = []
    bad_total_step: list[dict[str, Any]] = []
    bad_lb_u2: list[dict[str, Any]] = []

    min_total: dict[str, Any] | None = None
    min_fixed: dict[str, Any] | None = None
    min_gain: dict[str, Any] | None = None
    min_lb_u2: dict[str, Any] | None = None
    max_id_err: dict[str, Any] | None = None

    # Track tight total-step witnesses for exact rational confirmation.
    tight_total: list[dict[str, Any]] = []

    total_pairs = 0
    checked_pairs = 0

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)

        # Precompute rows for j=0..j_max+2
        rows: list[dict[str, Any]] = []
        f = h[:]
        for j in range(args.j_max + 3):
            m, lam, _, _ = mode_lambda_from_fg(f, g)
            mu = family_mu(k, j, lam)
            rows.append({"j": j, "m": m, "lam": lam, "mu": mu, "f": f[:]})
            if j < args.j_max + 2:
                f = next_f_times_1px(f)

        for j in range(args.j_max + 1):
            total_pairs += 1
            r0 = rows[j]
            r2 = rows[j + 2]

            m0 = int(r0["m"])
            m2 = int(r2["m"])
            lam0 = float(r0["lam"])
            lam2 = float(r2["lam"])
            mu0 = float(r0["mu"])
            mu2 = float(r2["mu"])
            mode_step = m2 - m0

            if mode_step != 1:
                maybe_store(
                    bad_mode_step,
                    {"k": k, "j": j, "m_j": m0, "m_j2": m2, "mode_step": mode_step},
                    args.store_cap,
                )
                continue

            checked_pairs += 1

            lam_step = lam2 - lam0
            if lam_step < -args.tol:
                maybe_store(
                    bad_lambda_step,
                    {"k": k, "j": j, "lambda_j": lam0, "lambda_j2": lam2, "lambda_step": lam_step},
                    args.store_cap,
                )

            mu2_at_lam0 = family_mu(k, j + 2, lam0)
            fixed_step = mu2_at_lam0 - mu0
            gain = mu2 - mu2_at_lam0
            total_step = mu2 - mu0
            id_err = (fixed_step + gain) - total_step

            # lane-2 helper ratio u2 from A_{j+2} coefficients
            f2 = r2["f"]
            u2 = None
            lb_u2 = None
            if m0 >= 1 and (m0 + 1) < len(f2):
                # Use integer division directly to avoid float(int) overflow on huge coefficients.
                u2 = f2[m0] / f2[m0 + 1]
                lb_u2 = family_mu(k, j + 2, u2) - mu0

            rec = {
                "k": k,
                "j": j,
                "mode_j": m0,
                "mode_j2": m2,
                "lambda_j": lam0,
                "lambda_j2": lam2,
                "lambda_step": lam_step,
                "fixed_step": fixed_step,
                "lambda_gain": gain,
                "total_step": total_step,
                "identity_err": id_err,
                "u2": u2,
                "lb_u2": lb_u2,
            }

            if min_total is None or total_step < min_total["total_step"]:
                min_total = rec
            if min_fixed is None or fixed_step < min_fixed["fixed_step"]:
                min_fixed = rec
            if min_gain is None or gain < min_gain["lambda_gain"]:
                min_gain = rec
            if lb_u2 is not None and (min_lb_u2 is None or lb_u2 < min_lb_u2["lb_u2"]):
                min_lb_u2 = rec
            if max_id_err is None or abs(id_err) > abs(max_id_err["identity_err"]):
                max_id_err = rec

            if total_step < 1.0 - args.tol:
                maybe_store(bad_total_step, rec, args.store_cap)

            if lb_u2 is not None and lb_u2 < 1.0 - args.tol:
                maybe_store(bad_lb_u2, rec, args.store_cap)

            # tight total-step list for exact checks
            if len(tight_total) < max(1, args.exact_top):
                tight_total.append(rec)
                tight_total.sort(key=lambda x: x["total_step"])
            elif total_step < tight_total[-1]["total_step"]:
                tight_total[-1] = rec
                tight_total.sort(key=lambda x: x["total_step"])

        if (k - args.k_min) % 300 == 0 or k == args.k_max:
            cur = float("nan") if min_total is None else min_total["total_step"]
            print(
                f"k={k:5d}: mode_bad={len(bad_mode_step)} lam_bad={len(bad_lambda_step)} "
                f"step_bad={len(bad_total_step)} min_total={cur:.12g}",
                flush=True,
            )

    exact_checks: list[dict[str, Any]] = []
    if args.exact_top > 0:
        seen: set[tuple[int, int]] = set()
        for rec in tight_total:
            k = int(rec["k"])
            j = int(rec["j"])
            if (k, j) in seen:
                continue
            seen.add((k, j))
            d0 = compute_margin(k, j)
            d2 = compute_margin(k, j + 2)
            # total mean-step exact
            exact_total = (d2.margin + (d2.mode - 1)) - (d0.margin + (d0.mode - 1))
            exact_checks.append(
                {
                    "k": k,
                    "j": j,
                    "float_total_step": rec["total_step"],
                    "exact_total_step": str(exact_total),
                    "exact_total_step_float": float(exact_total),
                    "mode_j": d0.mode,
                    "mode_j2": d2.mode,
                    "lambda_j": str(d0.lam),
                    "lambda_j2": str(d2.lam),
                }
            )

    # Exact checks on lb_u2<1 witnesses (if any).
    exact_lb_exception_checks: list[dict[str, Any]] = []
    if bad_lb_u2:
        seen_exc: set[tuple[int, int]] = set()
        for rec in bad_lb_u2:
            k = int(rec["k"])
            j = int(rec["j"])
            if (k, j) in seen_exc:
                continue
            seen_exc.add((k, j))
            d0 = compute_margin(k, j)
            d2 = compute_margin(k, j + 2)
            exact_total = (d2.margin + (d2.mode - 1)) - (d0.margin + (d0.mode - 1))
            exact_lb_exception_checks.append(
                {
                    "k": k,
                    "j": j,
                    "lb_u2_float": rec["lb_u2"],
                    "exact_total_step": str(exact_total),
                    "exact_total_step_float": float(exact_total),
                }
            )

    payload = {
        "params": vars(args),
        "overall": {
            "pairs_total": total_pairs,
            "pairs_mode_step_plus1": checked_pairs,
            "wall_s": time.time() - t0,
        },
        "failures": {
            "mode_step_plus1": {"count": len(bad_mode_step), "examples": bad_mode_step},
            "lane2_lambda_nondecreasing": {"count": len(bad_lambda_step), "examples": bad_lambda_step},
            "total_step_ge_1": {"count": len(bad_total_step), "examples": bad_total_step},
            "lb_u2_ge_1": {"count": len(bad_lb_u2), "examples": bad_lb_u2},
        },
        "extrema": {
            "min_total_step": min_total,
            "min_fixed_step": min_fixed,
            "min_lambda_gain": min_gain,
            "min_lb_u2": min_lb_u2,
            "max_identity_error_abs": max_id_err,
        },
        "tight_total_exact_checks": exact_checks,
        "lb_u2_exception_exact_checks": exact_lb_exception_checks,
    }

    print("-" * 96, flush=True)
    print(f"pairs_total={total_pairs}", flush=True)
    print(f"pairs_mode_step_plus1={checked_pairs}", flush=True)
    print(f"mode_step_plus1 fails={len(bad_mode_step)}", flush=True)
    print(f"lane2 lambda-step fails={len(bad_lambda_step)}", flush=True)
    print(f"total_step>=1 fails={len(bad_total_step)}", flush=True)
    print(f"lb_u2>=1 fails={len(bad_lb_u2)}", flush=True)
    if min_total is not None:
        print(f"min total_step={min_total['total_step']:.15f} at (k,j)=({min_total['k']},{min_total['j']})", flush=True)
    if min_fixed is not None:
        print(f"min fixed_step={min_fixed['fixed_step']:.15f} at (k,j)=({min_fixed['k']},{min_fixed['j']})", flush=True)
    if min_gain is not None:
        print(f"min lambda_gain={min_gain['lambda_gain']:.15f} at (k,j)=({min_gain['k']},{min_gain['j']})", flush=True)
    if min_lb_u2 is not None:
        print(f"min lb_u2={min_lb_u2['lb_u2']:.15f} at (k,j)=({min_lb_u2['k']},{min_lb_u2['j']})", flush=True)
    if max_id_err is not None:
        print(f"max |identity_err|={abs(max_id_err['identity_err']):.3e}", flush=True)
    if exact_checks:
        print("Exact checks on tight total-step witnesses:", flush=True)
        for ex in exact_checks:
            print(
                f"  (k,j)=({ex['k']},{ex['j']}): exact total_step={ex['exact_total_step_float']:.15f}",
                flush=True,
            )
    if exact_lb_exception_checks:
        print("Exact checks on lb_u2<1 exceptions:", flush=True)
        for ex in exact_lb_exception_checks:
            print(
                f"  (k,j)=({ex['k']},{ex['j']}): lb_u2={ex['lb_u2_float']:.12f}, "
                f"exact total_step={ex['exact_total_step_float']:.15f}",
                flush=True,
            )
    print(f"wall={payload['overall']['wall_s']:.2f}s", flush=True)

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f_out:
            json.dump(payload, f_out, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
