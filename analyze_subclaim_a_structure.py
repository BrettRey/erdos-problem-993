#!/usr/bin/env python3
"""Structural diagnostics for Sub-claim A on mixed spiders S(2^k,1^j).

Sub-claim A target:
  margin(k,j+2) >= margin(k,j) for k>=6, j>=0.

This script scans (k,j) and reports support for candidate intermediate lemmas:
  1) mode(k,j+2) - mode(k,j) == 1,
  2) lambda(k,j+2) - lambda(k,j) >= 0,
  3) decomposition of mean increment:
       mu_{j+2}(lam_{j+2}) - mu_j(lam_j)
       = [mu_{j+2}(lam_j)-mu_j(lam_j)] + [mu_{j+2}(lam_{j+2})-mu_{j+2}(lam_j)].
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


def main() -> None:
    ap = argparse.ArgumentParser(description="Scan structural lemmas for Sub-claim A.")
    ap.add_argument("--k-min", type=int, default=6)
    ap.add_argument("--k-max", type=int, default=2000)
    ap.add_argument("--j-max", type=int, default=200)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.k_min < 1:
        raise ValueError("k-min must be >= 1")
    if args.k_max < args.k_min:
        raise ValueError("k-max must be >= k-min")
    if args.j_max < 2:
        raise ValueError("j-max must be >= 2")

    t0 = time.time()

    bad_mode_step: list[dict[str, Any]] = []
    bad_lambda_step: list[dict[str, Any]] = []
    bad_margin_step: list[dict[str, Any]] = []

    min_mode_step = {"value": float("inf"), "k": -1, "j": -1}
    min_lambda_step = {"value": float("inf"), "k": -1, "j": -1}
    min_margin_step = {"value": float("inf"), "k": -1, "j": -1}
    min_mu_fixed_step = {"value": float("inf"), "k": -1, "j": -1}
    min_mu_lambda_gain = {"value": float("inf"), "k": -1, "j": -1}
    min_mu_total_step = {"value": float("inf"), "k": -1, "j": -1}

    total_pairs = 0

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)

        # Precompute mode/lambda/margin for j=0..j_max+2.
        f = h[:]
        arr: list[dict[str, float | int]] = []
        for j in range(args.j_max + 3):
            m, lam, _, _ = mode_lambda_from_fg(f, g)
            if m == 0:
                row = {"j": j, "mode": 0, "lam": 0.0, "margin": 0.0, "mu": 0.0}
            else:
                mu = family_mu(k, j, lam)
                margin = mu - (m - 1)
                row = {"j": j, "mode": m, "lam": lam, "margin": margin, "mu": mu}
            arr.append(row)
            f = next_f_times_1px(f)

        for j in range(args.j_max + 1):
            a = arr[j]
            b = arr[j + 2]

            mode_step = int(b["mode"]) - int(a["mode"])
            lam_step = float(b["lam"]) - float(a["lam"])
            margin_step = float(b["margin"]) - float(a["margin"])

            if mode_step < min_mode_step["value"]:
                min_mode_step = {"value": mode_step, "k": k, "j": j}
            if lam_step < min_lambda_step["value"]:
                min_lambda_step = {"value": lam_step, "k": k, "j": j}
            if margin_step < min_margin_step["value"]:
                min_margin_step = {"value": margin_step, "k": k, "j": j}

            if mode_step != 1 and len(bad_mode_step) < 100:
                bad_mode_step.append(
                    {
                        "k": k,
                        "j": j,
                        "mode_j": int(a["mode"]),
                        "mode_j2": int(b["mode"]),
                        "mode_step": mode_step,
                    }
                )

            if lam_step < -args.tol and len(bad_lambda_step) < 100:
                bad_lambda_step.append(
                    {
                        "k": k,
                        "j": j,
                        "lambda_j": float(a["lam"]),
                        "lambda_j2": float(b["lam"]),
                        "lambda_step": lam_step,
                    }
                )

            if margin_step < -args.tol and len(bad_margin_step) < 100:
                bad_margin_step.append(
                    {
                        "k": k,
                        "j": j,
                        "margin_j": float(a["margin"]),
                        "margin_j2": float(b["margin"]),
                        "margin_step": margin_step,
                    }
                )

            # Mean increment decomposition:
            # total = mu_{j+2}(lam_{j+2}) - mu_j(lam_j)
            # fixed = mu_{j+2}(lam_j) - mu_j(lam_j)
            # gain  = mu_{j+2}(lam_{j+2}) - mu_{j+2}(lam_j)
            lam_j = float(a["lam"])
            lam_j2 = float(b["lam"])
            mu_j_lj = float(a["mu"])
            mu_j2_lj = family_mu(k, j + 2, lam_j)
            mu_j2_lj2 = float(b["mu"])
            fixed_step = mu_j2_lj - mu_j_lj
            lambda_gain = mu_j2_lj2 - mu_j2_lj
            total_step = mu_j2_lj2 - mu_j_lj

            if fixed_step < min_mu_fixed_step["value"]:
                min_mu_fixed_step = {"value": fixed_step, "k": k, "j": j}
            if lambda_gain < min_mu_lambda_gain["value"]:
                min_mu_lambda_gain = {"value": lambda_gain, "k": k, "j": j}
            if total_step < min_mu_total_step["value"]:
                min_mu_total_step = {"value": total_step, "k": k, "j": j}

            total_pairs += 1

        if (k - args.k_min) % 100 == 0 or k == args.k_max:
            print(
                f"k={k:5d} mode_bad={len(bad_mode_step)} lam_bad={len(bad_lambda_step)} "
                f"margin_bad={len(bad_margin_step)}",
                flush=True,
            )

    payload = {
        "params": vars(args),
        "overall": {
            "pairs": total_pairs,
            "wall_s": time.time() - t0,
        },
        "failures": {
            "mode_step_eq_1": bad_mode_step,
            "lambda_step_nonneg": bad_lambda_step,
            "margin_step_nonneg": bad_margin_step,
        },
        "mins": {
            "mode_step": min_mode_step,
            "lambda_step": min_lambda_step,
            "margin_step": min_margin_step,
            "mu_fixed_step": min_mu_fixed_step,
            "mu_lambda_gain": min_mu_lambda_gain,
            "mu_total_step": min_mu_total_step,
        },
    }

    print("-" * 96, flush=True)
    print(f"pairs={total_pairs}", flush=True)
    print(f"mode_step==1 failures={len(bad_mode_step)}", flush=True)
    print(f"lambda_step>=0 failures={len(bad_lambda_step)}", flush=True)
    print(f"margin_step>=0 failures={len(bad_margin_step)}", flush=True)
    print(f"min mode_step={min_mode_step}", flush=True)
    print(f"min lambda_step={min_lambda_step}", flush=True)
    print(f"min margin_step={min_margin_step}", flush=True)
    print(f"min mu_fixed_step={min_mu_fixed_step}", flush=True)
    print(f"min mu_lambda_gain={min_mu_lambda_gain}", flush=True)
    print(f"min mu_total_step={min_mu_total_step}", flush=True)
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

