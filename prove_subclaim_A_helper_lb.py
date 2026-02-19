#!/usr/bin/env python3
"""Check/attempt closure of the helper lower bound in Sub-claim A.

Helper bound (lane-2 surrogate):
  lb_u2(k,j) := mu_{k,j+2}(u2_{k,j}) - mu_{k,j}(lambda_{k,j}) >= 1,
where
  lambda_{k,j} is tie fugacity at mode m_{k,j},
  u2_{k,j} = a'_{m_{k,j}} / a'_{m_{k,j}+1},
  a'_t = [x^t](1+2x)^k(1+x)^(j+2).

This script:
  1) scans for failures of lb_u2 >= 1,
  2) checks whether failures are confined to tiny finite bases,
  3) computes exact-rational lb_u2 on the failure set.
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
    family_mu,
    mode_lambda_from_fg,
    next_f_times_1px,
)


def maybe_store(dst: list[dict[str, Any]], item: dict[str, Any], cap: int) -> None:
    if len(dst) < cap:
        dst.append(item)


def exact_mu(k: int, j: int, lam: Fraction) -> Fraction:
    """Exact mu(S(2^k,1^j); lam) as a rational."""
    one = Fraction(1, 1)
    u = lam / (one + lam)
    q = (2 * lam) / (one + 2 * lam)
    a = Fraction(k) * q
    b = one + Fraction(k) * u

    if k >= j:
        r = lam * (one + lam) ** (k - j) / (one + 2 * lam) ** k
    else:
        r = lam / ((one + 2 * lam) ** k * (one + lam) ** (j - k))

    return (a + Fraction(j) * u + r * b) / (one + r)


def main() -> None:
    ap = argparse.ArgumentParser(description="Helper lower bound checker for Sub-claim A.")
    ap.add_argument("--k-min", type=int, default=6)
    ap.add_argument("--k-max", type=int, default=4000)
    ap.add_argument("--j-max", type=int, default=80)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--store-cap", type=int, default=200)
    ap.add_argument("--exact-cap", type=int, default=20)
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
        f"Helper-lb scan: k={args.k_min}..{args.k_max}, j=0..{args.j_max}",
        flush=True,
    )

    bad_mode_step: list[dict[str, Any]] = []
    bad_lb: list[dict[str, Any]] = []

    min_lb: dict[str, Any] | None = None
    min_lb_k_ge_9: dict[str, Any] | None = None

    total_pairs = 0
    checked_pairs = 0

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)
        f = h[:]

        rows: list[dict[str, Any]] = []
        for j in range(args.j_max + 3):
            m, lam, im1, im = mode_lambda_from_fg(f, g)
            rows.append({"j": j, "m": m, "lam": lam, "im1": im1, "im": im, "f": f[:]})
            if j < args.j_max + 2:
                f = next_f_times_1px(f)

        for j in range(args.j_max + 1):
            total_pairs += 1
            r0 = rows[j]
            r2 = rows[j + 2]

            m0 = int(r0["m"])
            m2 = int(r2["m"])
            if m2 != m0 + 1:
                maybe_store(
                    bad_mode_step,
                    {"k": k, "j": j, "m_j": m0, "m_j2": m2, "mode_step": m2 - m0},
                    args.store_cap,
                )
                continue

            checked_pairs += 1

            f2 = r2["f"]
            if not (m0 >= 1 and (m0 + 1) < len(f2)):
                continue

            lam_j = float(r0["lam"])
            u2 = f2[m0] / f2[m0 + 1]
            lb = family_mu(k, j + 2, u2) - family_mu(k, j, lam_j)

            rec = {
                "k": k,
                "j": j,
                "mode_j": m0,
                "mode_j2": m2,
                "lambda_j": lam_j,
                "u2": u2,
                "lb_u2": lb,
            }

            if min_lb is None or lb < min_lb["lb_u2"]:
                min_lb = rec
            if k >= 9 and (min_lb_k_ge_9 is None or lb < min_lb_k_ge_9["lb_u2"]):
                min_lb_k_ge_9 = rec

            if lb < 1.0 - args.tol:
                maybe_store(bad_lb, rec, args.store_cap)

        if (k - args.k_min) % 300 == 0 or k == args.k_max:
            cur = float("nan") if min_lb is None else min_lb["lb_u2"]
            print(
                f"k={k:5d}: mode_bad={len(bad_mode_step)} lb_bad={len(bad_lb)} min_lb={cur:.12g}",
                flush=True,
            )

    # Exact-rational check of lb_u2 on found failures.
    exact_bad: list[dict[str, Any]] = []
    for rec in bad_lb[: args.exact_cap]:
        k = int(rec["k"])
        j = int(rec["j"])

        # Rebuild exact data for this (k,j).
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)
        f = h[:]
        rows: list[dict[str, Any]] = []
        for jj in range(j + 3):
            m, _lam, im1, im = mode_lambda_from_fg(f, g)
            rows.append({"m": m, "im1": im1, "im": im, "f": f[:]})
            if jj < j + 2:
                f = next_f_times_1px(f)

        m0 = int(rows[j]["m"])
        lam_j = Fraction(int(rows[j]["im1"]), int(rows[j]["im"]))
        f2 = rows[j + 2]["f"]
        u2 = Fraction(int(f2[m0]), int(f2[m0 + 1]))

        lb_exact = exact_mu(k, j + 2, u2) - exact_mu(k, j, lam_j)
        exact_bad.append(
            {
                "k": k,
                "j": j,
                "lb_float": rec["lb_u2"],
                "lb_exact": str(lb_exact),
                "lb_exact_float": float(lb_exact),
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
            "lb_u2_ge_1": {"count": len(bad_lb), "examples": bad_lb},
        },
        "extrema": {
            "min_lb_all": min_lb,
            "min_lb_k_ge_9": min_lb_k_ge_9,
        },
        "exact_bad_lb": exact_bad,
    }

    print("-" * 96, flush=True)
    print(f"pairs_total={total_pairs}", flush=True)
    print(f"pairs_mode_step_plus1={checked_pairs}", flush=True)
    print(f"mode_step_plus1 fails={len(bad_mode_step)}", flush=True)
    print(f"lb_u2>=1 fails={len(bad_lb)}", flush=True)
    if min_lb is not None:
        print(f"min lb_u2 (all)={min_lb['lb_u2']:.15f} at (k,j)=({min_lb['k']},{min_lb['j']})", flush=True)
    if min_lb_k_ge_9 is not None:
        print(
            f"min lb_u2 (k>=9)={min_lb_k_ge_9['lb_u2']:.15f} at (k,j)=({min_lb_k_ge_9['k']},{min_lb_k_ge_9['j']})",
            flush=True,
        )
    if exact_bad:
        print("Exact checks on lb_u2<1 witnesses:", flush=True)
        for ex in exact_bad:
            print(
                f"  (k,j)=({ex['k']},{ex['j']}): lb_exact={ex['lb_exact_float']:.15f}",
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
