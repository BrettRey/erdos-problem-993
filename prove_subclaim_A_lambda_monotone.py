#!/usr/bin/env python3
"""Lane #2 diagnostics: lambda(k,j+2) >= lambda(k,j) for mixed spiders.

Family: T_{k,j} = S(2^k,1^j), with
  I_{k,j}(x) = (1+2x)^k (1+x)^j + x(1+x)^k.

For each (k,j), let:
  m_j   = leftmost mode index of I_{k,j},
  lam_j = i_{m_j-1}/i_{m_j}.

Target (Sub-claim A, lane #2):
  lam_{j+2} >= lam_j.

This script scans and reports:
  1) direct lambda monotonicity failures,
  2) mode-step failures (m_{j+2} != m_j+1),
  3) a sufficient condition check:
       u2 >= lam_j,
     where
       u2 = a'_{m_j} / a'_{m_j+1},
       a'_t = [x^t](1+2x)^k(1+x)^(j+2).

If mode-step is +1, then lam_{j+2} is a weighted average of u2 and g-ratio y>0,
so lam_{j+2} >= u2. Thus u2 >= lam_j implies lam_{j+2} >= lam_j.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from fractions import Fraction
from typing import Any

from conjecture_a_mixed_spider_combined_scan import (
    build_g_coeffs_shift_binom,
    build_h_coeffs_2k,
    mode_lambda_from_fg,
    next_f_times_1px,
)


def maybe_store(dst: list[dict[str, Any]], item: dict[str, Any], cap: int) -> None:
    if len(dst) < cap:
        dst.append(item)


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify lane #2 (lambda parity-step monotonicity).")
    ap.add_argument("--k-min", type=int, default=6)
    ap.add_argument("--k-max", type=int, default=4000)
    ap.add_argument("--j-max", type=int, default=120)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--store-cap", type=int, default=100)
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
        f"Lane #2 scan: k={args.k_min}..{args.k_max}, j=0..{args.j_max}",
        flush=True,
    )

    bad_mode_step: list[dict[str, Any]] = []
    bad_lambda_step: list[dict[str, Any]] = []
    bad_u2_ge_lamj: list[dict[str, Any]] = []

    # Decomposition pieces for the u2 >= lam_j inequality:
    # F = (A+2B+C)(A+G) - (D+2A+B)(B+H).
    # Track T1 + T2 = F with:
    #   T1 = (A^2 - B D) + (A C - B^2)
    #   T2 = G(A+2B+C) - H(2A+B+D)
    min_f = {"value": float("inf"), "k": -1, "j": -1}
    min_t1 = {"value": float("inf"), "k": -1, "j": -1}
    min_lambda_step: dict[str, Any] | None = None
    min_u2_minus_lamj: dict[str, Any] | None = None

    total_pairs = 0
    checked_pairs = 0

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)

        # Precompute f_j, m_j, lam_j for j=0..j_max+2.
        f = h[:]
        rows: list[dict[str, Any]] = []
        for j in range(args.j_max + 3):
            m, lam, im1, im = mode_lambda_from_fg(f, g)
            rows.append({"j": j, "m": m, "lam": lam, "im1": im1, "im": im, "f": f[:]})
            f = next_f_times_1px(f)

        for j in range(args.j_max + 1):
            total_pairs += 1
            a = rows[j]
            b = rows[j + 2]

            m = int(a["m"])
            m2 = int(b["m"])
            im1_j = int(a["im1"])
            im_j = int(a["im"])
            im1_j2 = int(b["im1"])
            im_j2 = int(b["im"])

            # lam_step >= 0  <=> im1_j2/im_j2 >= im1_j/im_j
            step_num = im1_j2 * im_j - im1_j * im_j2
            if min_lambda_step is None or step_num * min_lambda_step["den"] < min_lambda_step["num"] * (
                im_j * im_j2
            ):
                min_lambda_step = {
                    "num": step_num,
                    "den": im_j * im_j2,
                    "k": k,
                    "j": j,
                }

            if m2 != m + 1:
                maybe_store(
                    bad_mode_step,
                    {"k": k, "j": j, "m_j": m, "m_j2": m2, "mode_step": m2 - m},
                    args.store_cap,
                )
                continue

            checked_pairs += 1

            if step_num < 0:
                maybe_store(
                    bad_lambda_step,
                    {
                        "k": k,
                        "j": j,
                        "im1_j": im1_j,
                        "im_j": im_j,
                        "im1_j2": im1_j2,
                        "im_j2": im_j2,
                    },
                    args.store_cap,
                )

            f_j = a["f"]
            f_j2 = b["f"]

            # Need neighbors around m in f_j and m+1 in f_j2.
            if not (m >= 2 and (m + 1) < len(f_j) and (m + 1) < len(f_j2)):
                continue

            A = int(f_j[m])
            B = int(f_j[m - 1])
            C = int(f_j[m - 2])
            D = int(f_j[m + 1])
            G = int(g[m]) if m < len(g) else 0
            H = int(g[m - 1]) if (m - 1) < len(g) else 0

            # u2 = a'_{m}/a'_{m+1} from j+2 a-part.
            u2_num = int(f_j2[m])
            u2_den = int(f_j2[m + 1])
            # u2 - lam_j >= 0  <=> u2_num/u2_den >= im1_j/im_j
            u2_gap_num = u2_num * im_j - u2_den * im1_j
            u2_gap_den = u2_den * im_j

            if (
                min_u2_minus_lamj is None
                or u2_gap_num * min_u2_minus_lamj["den"] < min_u2_minus_lamj["num"] * u2_gap_den
            ):
                min_u2_minus_lamj = {"num": u2_gap_num, "den": u2_gap_den, "k": k, "j": j}

            if u2_gap_num < 0:
                maybe_store(
                    bad_u2_ge_lamj,
                    {
                        "k": k,
                        "j": j,
                        "u2_num": u2_num,
                        "u2_den": u2_den,
                        "im1_j": im1_j,
                        "im_j": im_j,
                    },
                    args.store_cap,
                )

            f_val = (A + 2 * B + C) * (A + G) - (D + 2 * A + B) * (B + H)
            t1 = (A * A - B * D) + (A * C - B * B)
            t2 = G * (A + 2 * B + C) - H * (2 * A + B + D)

            if f_val < min_f["value"]:
                min_f = {"value": f_val, "k": k, "j": j}
            if t1 < min_t1["value"]:
                min_t1 = {"value": t1, "k": k, "j": j}
        if (k - args.k_min) % 200 == 0 or k == args.k_max:
            print(
                f"k={k:5d} mode_bad={len(bad_mode_step)} "
                f"lam_bad={len(bad_lambda_step)} u2_bad={len(bad_u2_ge_lamj)}",
                flush=True,
            )

    payload = {
        "params": vars(args),
        "overall": {
            "pairs_total": total_pairs,
            "pairs_with_mode_step_plus1": checked_pairs,
            "wall_s": time.time() - t0,
        },
        "failures": {
            "mode_step_plus1": bad_mode_step,
            "lambda_step_nonnegative": bad_lambda_step,
            "u2_ge_lambda_j": bad_u2_ge_lamj,
        },
        "mins": {
            "lambda_step": min_lambda_step,
            "u2_minus_lambda_j": min_u2_minus_lamj,
            "F_value": min_f,
            "T1_value": min_t1,
        },
    }

    print("-" * 96, flush=True)
    print(f"pairs_total={total_pairs}", flush=True)
    print(f"pairs_mode_step_plus1={checked_pairs}", flush=True)
    print(f"mode_step_plus1 fails={len(bad_mode_step)}", flush=True)
    print(f"lambda_step_nonnegative fails={len(bad_lambda_step)}", flush=True)
    print(f"u2_ge_lambda_j fails={len(bad_u2_ge_lamj)}", flush=True)
    if min_lambda_step is not None:
        print(
            f"min lambda_step={Fraction(min_lambda_step['num'], min_lambda_step['den'])} "
            f"at (k={min_lambda_step['k']}, j={min_lambda_step['j']})",
            flush=True,
        )
    else:
        print("min lambda_step=N/A", flush=True)
    if min_u2_minus_lamj is not None:
        print(
            f"min u2-lambda_j={Fraction(min_u2_minus_lamj['num'], min_u2_minus_lamj['den'])} "
            f"at (k={min_u2_minus_lamj['k']}, j={min_u2_minus_lamj['j']})",
            flush=True,
        )
    else:
        print("min u2-lambda_j=N/A", flush=True)
    print(f"min F={min_f}", flush=True)
    print(f"min T1={min_t1}", flush=True)
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
