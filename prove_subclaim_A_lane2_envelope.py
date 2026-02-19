#!/usr/bin/env python3
"""Lane #2 envelope check for Sub-claim A (mixed spiders).

Target lane:
  lambda(k,j+2) >= lambda(k,j)  for k>=6, j>=0
on T_{k,j} = S(2^k,1^j), with
  I_{k,j}(x) = (1+2x)^k(1+x)^j + x(1+x)^k.

At m = mode(k,j), define
  A=f_m, B=f_{m-1}, C=f_{m-2}, D=f_{m+1},   f_t=[x^t](1+2x)^k(1+x)^j
  G=g_m, H=g_{m-1},                           g_t=[x^t]x(1+x)^k.

Lane-2 sufficient inequality (equivalent in this setup):
  F := (A+2B+C)(A+G) - (D+2A+B)(B+H) >= 0.

Write F = T1 + T2 with
  T1 = (A^2 - BD) + (AC - B^2),
  T2 = G(A+2B+C) - H(2A+B+D).

When T2 < 0, use s=G/A and
  F/A^2 = T1/A^2 + s * (T2/(AG)).
If s <= s_hat, then
  F/A^2 >= T1/A^2 + s_hat * (T2/(AG)).

This script uses the explicit coefficient lower bound on A:
  - if m <= k:
      A >= 2^m*C(k,m) + j*2^(m-1)*C(k,m-1),
      so s <= m / (2^m*(k-m+1) + j*2^(m-1)*m).
  - if m = k+1:
      A >= 2^k*j + 2^(k-1)*k*C(j,2),
      so s <= 1 / (2^k*j + 2^(k-1)*k*C(j,2)).

The script reports direct F-failures and envelope-failures.
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
    mode_lambda_from_fg,
    next_f_times_1px,
)


def maybe_store(dst: list[dict[str, Any]], item: dict[str, Any], cap: int) -> None:
    if len(dst) < cap:
        dst.append(item)


def leftmost_mode(coeffs: list[int]) -> int:
    m = 0
    best = coeffs[0]
    for i in range(1, len(coeffs)):
        v = coeffs[i]
        if v > best:
            best = v
            m = i
    return m


def main() -> None:
    ap = argparse.ArgumentParser(description="Lane-2 envelope verifier for Sub-claim A.")
    ap.add_argument("--k-min", type=int, default=6)
    ap.add_argument("--k-max", type=int, default=3000)
    ap.add_argument("--j-max", type=int, default=120)
    ap.add_argument("--store-cap", type=int, default=100)
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
        f"Lane-2 envelope scan: k={args.k_min}..{args.k_max}, j=0..{args.j_max}",
        flush=True,
    )

    bad_mode_step: list[dict[str, Any]] = []
    bad_mode_match_f: list[dict[str, Any]] = []
    bad_direct_f: list[dict[str, Any]] = []
    bad_s_upper: list[dict[str, Any]] = []
    bad_envelope: list[dict[str, Any]] = []
    bad_t1_when_t2_nonneg: list[dict[str, Any]] = []

    total_pairs = 0
    checked_pairs = 0
    t2_nonneg_count = 0
    t2_neg_count = 0
    g_zero_count = 0

    min_f_ratio_num: int | None = None
    min_f_ratio_den: int | None = None
    min_f_wit: dict[str, int] | None = None

    min_lb_ratio_num: int | None = None
    min_lb_ratio_den: int | None = None
    min_lb_wit: dict[str, int] | None = None

    min_env_ratio_num: int | None = None
    min_env_ratio_den: int | None = None
    min_env_ratio_wit: dict[str, int] | None = None

    max_s_tight_num: int | None = None
    max_s_tight_den: int | None = None
    max_s_tight_wit: dict[str, int] | None = None

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)

        rows: list[dict[str, Any]] = []
        f = h[:]
        for j in range(args.j_max + 3):
            m_i, _, _, _ = mode_lambda_from_fg(f, g)
            m_f = leftmost_mode(f)
            rows.append({"j": j, "m_i": m_i, "m_f": m_f, "f": f[:]})
            f = next_f_times_1px(f)

        for j in range(args.j_max + 1):
            total_pairs += 1
            r0 = rows[j]
            r2 = rows[j + 2]

            m = int(r0["m_i"])
            m2 = int(r2["m_i"])
            mf = int(r0["m_f"])

            if m != mf:
                maybe_store(
                    bad_mode_match_f,
                    {"k": k, "j": j, "mode_i": m, "mode_f": mf},
                    args.store_cap,
                )

            if m2 != m + 1:
                maybe_store(
                    bad_mode_step,
                    {"k": k, "j": j, "mode_j": m, "mode_j2": m2, "mode_step": m2 - m},
                    args.store_cap,
                )
                continue

            f0 = r0["f"]
            if not (m >= 2 and (m + 1) < len(f0)):
                continue

            checked_pairs += 1

            A = int(f0[m])
            B = int(f0[m - 1])
            C = int(f0[m - 2])
            D = int(f0[m + 1])
            G = int(g[m]) if m < len(g) else 0
            H = int(g[m - 1]) if (m - 1) < len(g) else 0

            F = (A + 2 * B + C) * (A + G) - (D + 2 * A + B) * (B + H)
            T1 = (A * A - B * D) + (A * C - B * B)
            T2 = G * (A + 2 * B + C) - H * (2 * A + B + D)

            if F < 0:
                maybe_store(
                    bad_direct_f,
                    {"k": k, "j": j, "mode": m, "F": F},
                    args.store_cap,
                )

            f_den = A * A
            if min_f_ratio_num is None or F * min_f_ratio_den < min_f_ratio_num * f_den:
                min_f_ratio_num = F
                min_f_ratio_den = f_den
                min_f_wit = {"k": k, "j": j, "mode": m}

            if G == 0:
                g_zero_count += 1
                continue

            # Explicit upper bound s=G/A <= num_s/den_s from two-term lower bound on A.
            if m <= k:
                num_s = m
                den_s = (1 << m) * (k - m + 1) + j * (1 << (m - 1)) * m
            elif m == k + 1:
                # Here g_{k+1}=1 and j>=1 on every observed pair.
                if j <= 0:
                    maybe_store(
                        bad_s_upper,
                        {"k": k, "j": j, "mode": m, "reason": "m=k+1 with j<=0"},
                        args.store_cap,
                    )
                    continue
                num_s = 1
                den_s = (1 << k) * j + (1 << (k - 1)) * k * (j * (j - 1) // 2)
            else:
                # Should not happen when G>0.
                maybe_store(
                    bad_s_upper,
                    {"k": k, "j": j, "mode": m, "reason": "G>0 but m>k+1"},
                    args.store_cap,
                )
                continue

            # Verify the s-upper bound itself: G/A <= num_s/den_s.
            # Cross-multiplied: G*den_s <= A*num_s.
            if G * den_s > A * num_s:
                maybe_store(
                    bad_s_upper,
                    {
                        "k": k,
                        "j": j,
                        "mode": m,
                        "lhs": G * den_s,
                        "rhs": A * num_s,
                    },
                    args.store_cap,
                )

            # Track tightness max of (G/A)/(num_s/den_s) = G*den_s/(A*num_s).
            tight_num = G * den_s
            tight_den = A * num_s
            if (
                max_s_tight_num is None
                or tight_num * max_s_tight_den > max_s_tight_num * tight_den
            ):
                max_s_tight_num = tight_num
                max_s_tight_den = tight_den
                max_s_tight_wit = {"k": k, "j": j, "mode": m}

            if T2 >= 0:
                t2_nonneg_count += 1
                if T1 < 0:
                    maybe_store(
                        bad_t1_when_t2_nonneg,
                        {"k": k, "j": j, "mode": m, "T1": T1, "T2": T2},
                        args.store_cap,
                    )
                continue

            t2_neg_count += 1

            # Ratio form of envelope condition:
            #   T1*den_s*G >= |T2|*num_s*A
            # i.e.
            #   env_ratio := (T1*den_s*G) / (|T2|*num_s*A) >= 1.
            env_ratio_num = T1 * den_s * G
            env_ratio_den = (-T2) * num_s * A
            if (
                min_env_ratio_num is None
                or env_ratio_num * min_env_ratio_den < min_env_ratio_num * env_ratio_den
            ):
                min_env_ratio_num = env_ratio_num
                min_env_ratio_den = env_ratio_den
                min_env_ratio_wit = {"k": k, "j": j, "mode": m}

            # Envelope lower bound for F/A^2 in the T2<0 regime:
            #   LB = T1/A^2 + (num_s/den_s) * T2/(A*G).
            # Equivalent integer numerator:
            #   lb_num = T1*den_s*G + T2*num_s*A.
            lb_num = T1 * den_s * G + T2 * num_s * A
            lb_den = A * A * den_s * G

            if min_lb_ratio_num is None or lb_num * min_lb_ratio_den < min_lb_ratio_num * lb_den:
                min_lb_ratio_num = lb_num
                min_lb_ratio_den = lb_den
                min_lb_wit = {"k": k, "j": j, "mode": m}

            if lb_num < 0:
                maybe_store(
                    bad_envelope,
                    {
                        "k": k,
                        "j": j,
                        "mode": m,
                        "lb_num": lb_num,
                    },
                    args.store_cap,
                )

        if (k - args.k_min) % 500 == 0 or k == args.k_max:
            print(
                f"k={k:5d} direct_bad={len(bad_direct_f)} env_bad={len(bad_envelope)} "
                f"mode_bad={len(bad_mode_step)}",
                flush=True,
            )

    def frac_float(num: int | None, den: int | None) -> float | None:
        if num is None or den is None:
            return None
        return num / den

    payload: dict[str, Any] = {
        "params": vars(args),
        "overall": {
            "pairs_total": total_pairs,
            "pairs_mode_step_plus1": checked_pairs,
            "g_zero_pairs": g_zero_count,
            "t2_negative_pairs": t2_neg_count,
            "t2_nonnegative_pairs": t2_nonneg_count,
            "wall_s": time.time() - t0,
        },
        "failures": {
            "mode_step_plus1": bad_mode_step,
            "mode_equals_product_mode": bad_mode_match_f,
            "direct_F_nonnegative": bad_direct_f,
            "s_upper_bound": bad_s_upper,
            "envelope_nonnegative": bad_envelope,
            "t1_nonnegative_when_t2_nonnegative": bad_t1_when_t2_nonneg,
        },
        "extrema": {
            "min_F_over_A2": {
                "value": frac_float(min_f_ratio_num, min_f_ratio_den),
                "witness": min_f_wit,
            },
            "min_envelope_lower_bound": {
                "value": frac_float(min_lb_ratio_num, min_lb_ratio_den),
                "witness": min_lb_wit,
            },
            "min_envelope_ratio": {
                "value": frac_float(min_env_ratio_num, min_env_ratio_den),
                "witness": min_env_ratio_wit,
            },
            "max_s_over_s_hat": {
                "value": frac_float(max_s_tight_num, max_s_tight_den),
                "witness": max_s_tight_wit,
            },
        },
    }

    print("-" * 96, flush=True)
    print(f"pairs_total={total_pairs}", flush=True)
    print(f"pairs_mode_step_plus1={checked_pairs}", flush=True)
    print(f"direct_F_nonnegative fails={len(bad_direct_f)}", flush=True)
    print(f"envelope_nonnegative fails={len(bad_envelope)}", flush=True)
    print(f"s_upper_bound fails={len(bad_s_upper)}", flush=True)
    print(f"mode_step_plus1 fails={len(bad_mode_step)}", flush=True)
    print(f"mode_equals_product_mode fails={len(bad_mode_match_f)}", flush=True)
    print(f"T2<0 pairs={t2_neg_count}, T2>=0 pairs={t2_nonneg_count}", flush=True)
    if min_f_ratio_num is not None and min_f_ratio_den is not None and min_f_wit is not None:
        print(
            f"min F/A^2={min_f_ratio_num/min_f_ratio_den:.18g} at "
            f"(k,j,m)=({min_f_wit['k']},{min_f_wit['j']},{min_f_wit['mode']})",
            flush=True,
        )
    if min_lb_ratio_num is not None and min_lb_ratio_den is not None and min_lb_wit is not None:
        print(
            f"min envelope LB={min_lb_ratio_num/min_lb_ratio_den:.18g} at "
            f"(k,j,m)=({min_lb_wit['k']},{min_lb_wit['j']},{min_lb_wit['mode']})",
            flush=True,
        )
    if (
        min_env_ratio_num is not None
        and min_env_ratio_den is not None
        and min_env_ratio_wit is not None
    ):
        print(
            f"min envelope ratio={min_env_ratio_num/min_env_ratio_den:.18g} at "
            f"(k,j,m)=({min_env_ratio_wit['k']},{min_env_ratio_wit['j']},{min_env_ratio_wit['mode']})",
            flush=True,
        )
    if max_s_tight_num is not None and max_s_tight_den is not None and max_s_tight_wit is not None:
        print(
            f"max (s/s_hat)={max_s_tight_num/max_s_tight_den:.18g} at "
            f"(k,j,m)=({max_s_tight_wit['k']},{max_s_tight_wit['j']},{max_s_tight_wit['mode']})",
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
