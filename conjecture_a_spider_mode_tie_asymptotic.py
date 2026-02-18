#!/usr/bin/env python3
"""Stable large-k scan for focused mode-tie quantities on S(2^k).

For T = S(2^k), choose a tip leaf l and its degree-2 support s, and define:
  A = T - l = S(2^(k-1), 1)
  B = T - {l,s} = S(2^(k-1))

At lambda_m(T) = i_{m-1}(T) / i_m(T), with m = mode(T), compute:
  margin = mu(T, lambda_m) - (m - 1)
  c1     = mu(A, lambda_m) - (m - 1)
  c2     = mu(B, lambda_m) - (m - 2)
  dmu    = mu(A, lambda_m) - mu(B, lambda_m)

This script uses closed forms and log-scaled ratios to avoid overflow.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from typing import Any


def coeff_s2k(k: int, j: int) -> int:
    """Coefficient i_j of I(S(2^k); x) = (1+2x)^k + x(1+x)^k."""
    t1 = math.comb(k, j) * (1 << j) if 0 <= j <= k else 0
    t2 = math.comb(k, j - 1) if 1 <= j <= k + 1 else 0
    return t1 + t2


def exact_mode_s2k(k: int) -> int:
    """Leftmost mode by exact coefficient scan (O(k))."""
    best_j = 0
    best_v = coeff_s2k(k, 0)
    for j in range(1, k + 2):
        v = coeff_s2k(k, j)
        if v > best_v:
            best_v = v
            best_j = j
    return best_j


def mode_formula_s2k(k: int) -> int:
    """Empirical exact formula for leftmost mode (verified in scans)."""
    if k <= 2:
        return exact_mode_s2k(k)
    return (2 * k + 1) // 3


def exp_clamped(x: float) -> float:
    """exp(x) with hard clamps for float-safe tails."""
    if x < -745.0:
        return 0.0
    if x > 709.0:
        return float("inf")
    return math.exp(x)


def stable_quantities_for_k(k: int) -> dict[str, float | int]:
    """Compute mode-tie quantities for a single k using stable formulas."""
    m = mode_formula_s2k(k)
    i_prev = coeff_s2k(k, m - 1)
    i_mode = coeff_s2k(k, m)
    lam = exp_clamped(math.log(i_prev) - math.log(i_mode))

    # r_k = lambda * ((1+lambda)/(1+2lambda))^k
    log_q = math.log1p(lam) - math.log1p(2.0 * lam)
    r_k = lam * exp_clamped(k * log_q)
    r_km1 = lam * exp_clamped((k - 1) * log_q)

    # Mean for T = S(2^k).
    mu_t = (
        2.0 * k * lam / (1.0 + 2.0 * lam)
        + r_k
        + r_k * k * lam / (1.0 + lam)
    ) / (1.0 + r_k)

    # Shared term for A/B with k-1.
    b0 = 1.0 + lam * (k - 1) / (1.0 + lam)

    # Mean for B = S(2^(k-1)).
    mu_b = (2.0 * (k - 1) * lam / (1.0 + 2.0 * lam) + r_km1 * b0) / (1.0 + r_km1)

    # Mean for A = S(2^(k-1), 1).
    a0 = lam + 2.0 * (k - 1) * lam * (1.0 + lam) / (1.0 + 2.0 * lam)
    mu_a = (a0 + r_km1 * b0) / ((1.0 + lam) + r_km1)

    c1 = mu_a - (m - 1)
    c2 = mu_b - (m - 2)
    margin = mu_t - (m - 1)
    dmu = mu_a - mu_b

    # Partition-function weights in mu(T)-decomposition:
    # I_T = I_A + lambda I_B, so alpha + beta = 1.
    den = (1.0 + 2.0 * lam) + (1.0 + lam) * r_km1
    alpha = ((1.0 + lam) + r_km1) / den
    beta = lam * (1.0 + r_km1) / den
    combined = alpha * c1 + beta * c2

    return {
        "k": k,
        "n": 2 * k + 1,
        "mode": m,
        "lambda_mode": lam,
        "mu_t": mu_t,
        "mu_a": mu_a,
        "mu_b": mu_b,
        "margin": margin,
        "c1": c1,
        "c2": c2,
        "dmu": dmu,
        "alpha": alpha,
        "beta": beta,
        "alpha_plus_beta_minus_1": alpha + beta - 1.0,
        "combined_minus_margin": combined - margin,
        "r_k": r_k,
        "r_km1": r_km1,
    }


def bump_min(store: dict[str, Any], key: str, value: float, k: int) -> None:
    cur = store.get(key)
    if cur is None or value < cur["value"]:
        store[key] = {"value": value, "k": k}


def main() -> None:
    ap = argparse.ArgumentParser(description="Stable S(2^k) mode-tie asymptotic scan.")
    ap.add_argument("--min-k", type=int, default=2)
    ap.add_argument("--max-k", type=int, default=1000)
    ap.add_argument(
        "--verify-mode-max",
        type=int,
        default=1000,
        help="Exact-check mode formula for k<=this bound (0 disables).",
    )
    ap.add_argument("--sample-every", type=int, default=100)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.min_k < 2:
        raise ValueError("min-k must be >= 2")
    if args.max_k < args.min_k:
        raise ValueError("max-k must be >= min-k")

    t0 = time.time()
    print(
        f"S(2^k) mode-tie asymptotic scan: k={args.min_k}..{args.max_k}",
        flush=True,
    )

    mode_formula_fail = 0
    mode_formula_first_fail: dict[str, int] | None = None

    mins: dict[str, Any] = {}
    by_mod3: dict[str, dict[str, Any]] = {
        "0": {},
        "1": {},
        "2": {},
    }

    max_abs_alpha_sum_err = 0.0
    max_abs_combined_err = 0.0

    samples: list[dict[str, Any]] = []

    for k in range(args.min_k, args.max_k + 1):
        if args.verify_mode_max > 0 and k <= args.verify_mode_max:
            m_exact = exact_mode_s2k(k)
            m_formula = mode_formula_s2k(k)
            if m_exact != m_formula:
                mode_formula_fail += 1
                if mode_formula_first_fail is None:
                    mode_formula_first_fail = {
                        "k": k,
                        "m_exact": m_exact,
                        "m_formula": m_formula,
                    }

        row = stable_quantities_for_k(k)
        mod = str(k % 3)

        bump_min(mins, "margin", row["margin"], k)  # type: ignore[arg-type]
        bump_min(mins, "c1", row["c1"], k)  # type: ignore[arg-type]
        bump_min(mins, "c2", row["c2"], k)  # type: ignore[arg-type]
        bump_min(mins, "dmu", row["dmu"], k)  # type: ignore[arg-type]
        bump_min(by_mod3[mod], "margin", row["margin"], k)  # type: ignore[arg-type]
        bump_min(by_mod3[mod], "c1", row["c1"], k)  # type: ignore[arg-type]
        bump_min(by_mod3[mod], "c2", row["c2"], k)  # type: ignore[arg-type]
        bump_min(by_mod3[mod], "dmu", row["dmu"], k)  # type: ignore[arg-type]

        max_abs_alpha_sum_err = max(
            max_abs_alpha_sum_err,
            abs(row["alpha_plus_beta_minus_1"]),  # type: ignore[arg-type]
        )
        max_abs_combined_err = max(
            max_abs_combined_err,
            abs(row["combined_minus_margin"]),  # type: ignore[arg-type]
        )

        if (
            k == args.min_k
            or k == args.max_k
            or (args.sample_every > 0 and (k - args.min_k) % args.sample_every == 0)
        ):
            samples.append(
                {
                    "k": row["k"],
                    "n": row["n"],
                    "mode": row["mode"],
                    "lambda_mode": row["lambda_mode"],
                    "margin": row["margin"],
                    "c1": row["c1"],
                    "c2": row["c2"],
                    "dmu": row["dmu"],
                }
            )

        if k == args.min_k or k == args.max_k or k % 500 == 0:
            print(
                f"k={k:6d}: margin={row['margin']:.12f} "
                f"c1={row['c1']:.12f} c2={row['c2']:.12f} dmu={row['dmu']:.12f}",
                flush=True,
            )

    payload = {
        "params": {
            "min_k": args.min_k,
            "max_k": args.max_k,
            "verify_mode_max": args.verify_mode_max,
            "sample_every": args.sample_every,
        },
        "summary": {
            "checked_k": args.max_k - args.min_k + 1,
            "mode_formula_fail": mode_formula_fail,
            "mode_formula_first_fail": mode_formula_first_fail,
            "mins": mins,
            "mins_by_k_mod_3": by_mod3,
            "max_abs_alpha_plus_beta_minus_1": max_abs_alpha_sum_err,
            "max_abs_combined_minus_margin": max_abs_combined_err,
            "wall_s": time.time() - t0,
        },
        "samples": samples,
    }

    print("-" * 88, flush=True)
    print(f"mode formula mismatches: {mode_formula_fail}", flush=True)
    print(
        f"global mins: margin={mins['margin']} c1={mins['c1']} "
        f"c2={mins['c2']} dmu={mins['dmu']}",
        flush=True,
    )
    print(
        f"max errors: |alpha+beta-1|={max_abs_alpha_sum_err:.3e}, "
        f"|combined-margin|={max_abs_combined_err:.3e}",
        flush=True,
    )
    print(f"wall: {time.time()-t0:.2f}s", flush=True)

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
