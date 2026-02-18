#!/usr/bin/env python3
"""Sub-claim A diagnostics for mixed spiders.

Sub-claim A target:
  margin(k,j+2) >= margin(k,j)
for mixed spiders T = S(2^k,1^j), with k >= 6 and j >= 0.

At each (k,j), margin is the tie-fugacity margin
  margin = mu(T, lambda_m) - (m-1),
where m is the leftmost mode at lambda=1 and lambda_m = i_{m-1}/i_m.

This script computes:
  - direct delta: margin(k,j+2)-margin(k,j),
  - decomposition deltaA + deltaE = delta_margin, where
      A = k*2λ/(1+2λ) + j*λ/(1+λ) - (m-1),
      E = r(B-A)/(1+r), r = λ(1+λ)^(k-j)/(1+2λ)^k,
      B = 1 + k*λ/(1+λ) - (m-1),
  - mode shift diagnostics m(k,j+2)-m(k,j).

Use this as a shared harness for analytic lanes from multiple agents.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from typing import Any

from conjecture_a_mixed_spider_combined_scan import (
    build_g_coeffs_shift_binom,
    build_h_coeffs_2k,
    family_log_r,
    family_mu,
    mode_lambda_from_fg,
    next_f_times_1px,
)


def exp_clamped(x: float) -> float:
    if x < -745.0:
        return 0.0
    if x > 709.0:
        return float("inf")
    return math.exp(x)


def margin_data(k: int, j: int, f: list[int], g: list[int]) -> dict[str, float | int]:
    m, lam, _, _ = mode_lambda_from_fg(f, g)
    if m == 0:
        return {"mode": 0, "lambda_mode": 0.0, "margin": 0.0, "A": 0.0, "E": 0.0, "r": 0.0}

    log_r = family_log_r(k, j, lam)
    r = exp_clamped(log_r)

    A = k * (2.0 * lam / (1.0 + 2.0 * lam)) + j * (lam / (1.0 + lam)) - (m - 1)
    B = 1.0 + k * (lam / (1.0 + lam)) - (m - 1)
    E = r * (B - A) / (1.0 + r)
    margin = family_mu(k, j, lam) - (m - 1)
    return {
        "mode": m,
        "lambda_mode": lam,
        "margin": margin,
        "A": A,
        "E": E,
        "r": r,
    }


def maybe_record(out: list[dict[str, Any]], item: dict[str, Any], cap: int) -> None:
    if len(out) < cap:
        out.append(item)


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify Sub-claim A on mixed spiders.")
    ap.add_argument("--k-min", type=int, default=6)
    ap.add_argument("--k-max", type=int, default=8000)
    ap.add_argument("--j-max", type=int, default=80)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--store-cap", type=int, default=300)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.k_min < 1:
        raise ValueError("k-min must be >= 1")
    if args.k_max < args.k_min:
        raise ValueError("k-max must be >= k-min")
    if args.j_max < 2:
        raise ValueError("j-max must be >= 2")

    t0 = time.time()
    print(
        f"Sub-claim A scan: k={args.k_min}..{args.k_max}, j=0..{args.j_max}",
        flush=True,
    )

    full_fail: list[dict[str, Any]] = []
    even_tail_fail: list[dict[str, Any]] = []
    odd_tail_fail: list[dict[str, Any]] = []

    mode_shift_hist: dict[int, int] = {}
    mode_shift_fail: list[dict[str, Any]] = []

    min_delta: dict[str, Any] | None = None
    min_even_tail_delta: dict[str, Any] | None = None
    min_odd_tail_delta: dict[str, Any] | None = None
    min_delta_A: dict[str, Any] | None = None
    min_delta_E: dict[str, Any] | None = None

    total_pairs = 0

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)
        f = h[:]

        # Need data at j and j+2 up to j_max, so compute 0..j_max+2.
        data: list[dict[str, float | int]] = []
        for j in range(0, args.j_max + 3):
            d = margin_data(k, j, f, g)
            data.append(d)
            if j < args.j_max + 2:
                f = next_f_times_1px(f)

        for j in range(0, args.j_max + 1):
            d0 = data[j]
            d2 = data[j + 2]

            delta = float(d2["margin"]) - float(d0["margin"])
            delta_A = float(d2["A"]) - float(d0["A"])
            delta_E = float(d2["E"]) - float(d0["E"])
            mode_shift = int(d2["mode"]) - int(d0["mode"])
            total_pairs += 1

            mode_shift_hist[mode_shift] = mode_shift_hist.get(mode_shift, 0) + 1
            if mode_shift not in (0, 1, 2):
                maybe_record(
                    mode_shift_fail,
                    {"k": k, "j": j, "mode_j": int(d0["mode"]), "mode_j2": int(d2["mode"])},
                    args.store_cap,
                )

            wit = {
                "k": k,
                "j": j,
                "delta_margin": delta,
                "delta_A": delta_A,
                "delta_E": delta_E,
                "mode_j": int(d0["mode"]),
                "mode_j2": int(d2["mode"]),
                "mode_shift": mode_shift,
                "lambda_j": float(d0["lambda_mode"]),
                "lambda_j2": float(d2["lambda_mode"]),
                "margin_j": float(d0["margin"]),
                "margin_j2": float(d2["margin"]),
            }

            if min_delta is None or delta < min_delta["delta_margin"]:
                min_delta = wit
            if min_delta_A is None or delta_A < min_delta_A["delta_A"]:
                min_delta_A = wit
            if min_delta_E is None or delta_E < min_delta_E["delta_E"]:
                min_delta_E = wit

            if delta < -args.tol:
                maybe_record(full_fail, wit, args.store_cap)

            if j >= 4 and j % 2 == 0:
                if min_even_tail_delta is None or delta < min_even_tail_delta["delta_margin"]:
                    min_even_tail_delta = wit
                if delta < -args.tol:
                    maybe_record(even_tail_fail, wit, args.store_cap)
            if j >= 5 and j % 2 == 1:
                if min_odd_tail_delta is None or delta < min_odd_tail_delta["delta_margin"]:
                    min_odd_tail_delta = wit
                if delta < -args.tol:
                    maybe_record(odd_tail_fail, wit, args.store_cap)

        if (k - args.k_min) % 100 == 0 or k == args.k_max:
            print(
                f"k={k:5d}: full_fail={len(full_fail)} even_tail_fail={len(even_tail_fail)} "
                f"odd_tail_fail={len(odd_tail_fail)}",
                flush=True,
            )

    summary = {
        "params": vars(args),
        "overall": {
            "total_pairs": total_pairs,
            "wall_s": time.time() - t0,
            "min_delta": min_delta,
            "min_even_tail_delta": min_even_tail_delta,
            "min_odd_tail_delta": min_odd_tail_delta,
            "min_delta_A": min_delta_A,
            "min_delta_E": min_delta_E,
            "mode_shift_hist": mode_shift_hist,
        },
        "claims": {
            "full_subclaim_A_j_ge_0": {"fail_count": len(full_fail), "fails": full_fail},
            "even_tail_j_ge_4": {"fail_count": len(even_tail_fail), "fails": even_tail_fail},
            "odd_tail_j_ge_5": {"fail_count": len(odd_tail_fail), "fails": odd_tail_fail},
            "mode_shift_outside_0_1_2": {
                "fail_count": len(mode_shift_fail),
                "fails": mode_shift_fail,
            },
        },
    }

    print("-" * 96, flush=True)
    print(f"total_pairs={total_pairs}", flush=True)
    print(f"min_delta={min_delta}", flush=True)
    print(f"min_even_tail_delta={min_even_tail_delta}", flush=True)
    print(f"min_odd_tail_delta={min_odd_tail_delta}", flush=True)
    print(f"min_delta_A={min_delta_A}", flush=True)
    print(f"min_delta_E={min_delta_E}", flush=True)
    print(f"mode_shift_hist={mode_shift_hist}", flush=True)
    print(f"full_fail={len(full_fail)}", flush=True)
    print(f"even_tail_fail={len(even_tail_fail)}", flush=True)
    print(f"odd_tail_fail={len(odd_tail_fail)}", flush=True)
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
