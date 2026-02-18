#!/usr/bin/env python3
"""Check minimizer-reduction structure for mixed spiders S(2^k,1^j).

We analyze the combined bound (equal to tie margin):
  margin(k,j) = mu(T_{k,j}, lambda_m) - (m-1),
where T_{k,j} = S(2^k,1^j), m is leftmost mode at lambda=1, and
lambda_m = i_{m-1} / i_m.

For each k we compute:
  - argmin over j in [0, j_max],
  - argmin over j in [2, j_max],
  - parity-tail monotonicity for j>=4:
      margin(k, j+2) >= margin(k, j),
  - residue comparisons:
      k mod 3 = 1  -> expect j=0 minimal,
      k mod 3 = 0,2 -> expect j=1 minimal (with small-k exceptions).
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
    family_log_i,
    family_mu,
    mode_lambda_from_fg,
    next_f_times_1px,
    product_log_i,
    product_mu,
    weights_from_log_i,
)


def best_margin_for_j(k: int, j: int, f: list[int], g: list[int]) -> tuple[float, str, int, float]:
    """Return best combined margin over tip/unit leaf choices at fixed (k,j)."""
    m, lam, _, _ = mode_lambda_from_fg(f, g)
    if m == 0:
        return 0.0, "none", m, lam

    mu_a_tip = family_mu(k - 1, j + 1, lam)
    mu_b_tip = family_mu(k - 1, j, lam)
    c1_tip = mu_a_tip - (m - 1)
    c2_tip = mu_b_tip - (m - 2)
    w1_tip, w2_tip = weights_from_log_i(
        family_log_i(k - 1, j + 1, lam),
        family_log_i(k - 1, j, lam),
        lam,
    )
    tip = w1_tip * c1_tip + w2_tip * c2_tip

    if j == 0:
        return tip, "tip", m, lam

    mu_a_unit = family_mu(k, j - 1, lam)
    mu_b_unit = product_mu(k, j - 1, lam)
    c1_unit = mu_a_unit - (m - 1)
    c2_unit = mu_b_unit - (m - 2)
    w1_unit, w2_unit = weights_from_log_i(
        family_log_i(k, j - 1, lam),
        product_log_i(k, j - 1, lam),
        lam,
    )
    unit = w1_unit * c1_unit + w2_unit * c2_unit

    if tip <= unit:
        return tip, "tip", m, lam
    return unit, "unit", m, lam


def maybe_record(out: list[dict[str, Any]], item: dict[str, Any], cap: int) -> None:
    if len(out) < cap:
        out.append(item)


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify mixed-spider minimizer-reduction structure.")
    ap.add_argument("--k-min", type=int, default=1)
    ap.add_argument("--k-max", type=int, default=800)
    ap.add_argument("--j-max", type=int, default=80)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--store-cap", type=int, default=200)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.k_min < 1:
        raise ValueError("k-min must be >= 1")
    if args.k_max < args.k_min:
        raise ValueError("k-max must be >= k-min")
    if args.j_max < 5:
        raise ValueError("j-max must be >= 5")

    t0 = time.time()
    print(
        f"Minimizer-reduction scan: k={args.k_min}..{args.k_max}, j=0..{args.j_max}",
        flush=True,
    )

    global_min: dict[str, Any] | None = None
    total_pairs = 0

    # Sub-claim A diagnostics.
    bad_even_tail: list[dict[str, Any]] = []
    bad_odd_tail: list[dict[str, Any]] = []
    bad_min_ge2_not_23: list[dict[str, Any]] = []

    # Sub-claim B/C diagnostics.
    bad_mod1_j0: list[dict[str, Any]] = []
    bad_mod02_j1: list[dict[str, Any]] = []
    bad_global_j_gt1: list[dict[str, Any]] = []

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)
        g = build_g_coeffs_shift_binom(k)
        f = h[:]

        vals: list[float] = []
        meta: list[dict[str, Any]] = []

        for j in range(0, args.j_max + 1):
            val, leaf_type, m, lam = best_margin_for_j(k, j, f, g)
            vals.append(val)
            meta.append(
                {
                    "k": k,
                    "j": j,
                    "leaf_type": leaf_type,
                    "mode": m,
                    "lambda_mode": lam,
                    "margin": val,
                }
            )
            total_pairs += 1

            if global_min is None or val < global_min["margin"]:
                global_min = dict(meta[-1])

            if j < args.j_max:
                f = next_f_times_1px(f)

        # A: parity-tail monotonicity (j>=4 by parity).
        for j in range(4, args.j_max - 1, 2):
            if vals[j + 2] < vals[j] - args.tol:
                maybe_record(
                    bad_even_tail,
                    {
                        "k": k,
                        "j": j,
                        "margin_j": vals[j],
                        "margin_j_plus_2": vals[j + 2],
                    },
                    args.store_cap,
                )
                break
        for j in range(5, args.j_max - 1, 2):
            if vals[j + 2] < vals[j] - args.tol:
                maybe_record(
                    bad_odd_tail,
                    {
                        "k": k,
                        "j": j,
                        "margin_j": vals[j],
                        "margin_j_plus_2": vals[j + 2],
                    },
                    args.store_cap,
                )
                break

        # A': min over j>=2 is at j=2 or j=3.
        j_star_ge2 = min(range(2, args.j_max + 1), key=lambda jj: vals[jj])
        if j_star_ge2 not in (2, 3):
            maybe_record(
                bad_min_ge2_not_23,
                {
                    "k": k,
                    "j_star_ge2": j_star_ge2,
                    "margin_star": vals[j_star_ge2],
                    "margin_j2": vals[2],
                    "margin_j3": vals[3],
                },
                args.store_cap,
            )

        # B/C: residue-wise j=0/j=1 comparison.
        if k % 3 == 1:
            if vals[0] > vals[1] + args.tol:
                maybe_record(
                    bad_mod1_j0,
                    {
                        "k": k,
                        "margin_j0": vals[0],
                        "margin_j1": vals[1],
                    },
                    args.store_cap,
                )
        else:
            if vals[1] > vals[0] + args.tol:
                maybe_record(
                    bad_mod02_j1,
                    {
                        "k": k,
                        "margin_j0": vals[0],
                        "margin_j1": vals[1],
                    },
                    args.store_cap,
                )

        # Global j-min over scanned range.
        j_star = min(range(0, args.j_max + 1), key=lambda jj: vals[jj])
        if j_star > 1:
            maybe_record(
                bad_global_j_gt1,
                {
                    "k": k,
                    "j_star": j_star,
                    "margin_star": vals[j_star],
                    "margin_j0": vals[0],
                    "margin_j1": vals[1],
                },
                args.store_cap,
            )

        if (k - args.k_min) % 50 == 0 or k == args.k_max:
            print(
                f"k={k:5d}: bad_even={len(bad_even_tail)} bad_odd={len(bad_odd_tail)} "
                f"bad_ge2={len(bad_min_ge2_not_23)} bad_global={len(bad_global_j_gt1)}",
                flush=True,
            )

    summary = {
        "params": vars(args),
        "overall": {
            "total_pairs": total_pairs,
            "wall_s": time.time() - t0,
            "global_min": global_min,
        },
        "claims": {
            "A_even_tail_monotone_from_j4": {
                "fail_count": len(bad_even_tail),
                "fails": bad_even_tail,
            },
            "A_odd_tail_monotone_from_j5": {
                "fail_count": len(bad_odd_tail),
                "fails": bad_odd_tail,
            },
            "A_prime_min_ge2_in_{2,3}": {
                "fail_count": len(bad_min_ge2_not_23),
                "fails": bad_min_ge2_not_23,
            },
            "B_mod1_prefers_j0": {
                "fail_count": len(bad_mod1_j0),
                "fails": bad_mod1_j0,
            },
            "C_mod0or2_prefers_j1": {
                "fail_count": len(bad_mod02_j1),
                "fails": bad_mod02_j1,
            },
            "global_min_over_j0_jmax_in_{0,1}": {
                "fail_count": len(bad_global_j_gt1),
                "fails": bad_global_j_gt1,
            },
        },
    }

    print("-" * 96, flush=True)
    print(f"total_pairs={total_pairs}", flush=True)
    print(f"global_min={global_min}", flush=True)
    print(f"A_even_tail fails={len(bad_even_tail)}", flush=True)
    print(f"A_odd_tail fails={len(bad_odd_tail)}", flush=True)
    print(f"A' min_ge2_in{{2,3}} fails={len(bad_min_ge2_not_23)}", flush=True)
    print(f"B mod1->j0 fails={len(bad_mod1_j0)}", flush=True)
    print(f"C mod0/2->j1 fails={len(bad_mod02_j1)}", flush=True)
    print(f"global j-min in {{0,1}} fails={len(bad_global_j_gt1)}", flush=True)
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
