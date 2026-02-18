#!/usr/bin/env python3
"""Analyze combined-bound structure on mixed spiders S(2^k,1^j).

Tree family:
  T = S(2^k,1^j) with independence polynomial
    I_T(x) = (1+2x)^k (1+x)^j + x(1+x)^k.

At m = mode(I_T) (leftmost) and lambda_m = i_{m-1}/i_m, define:
  margin = mu(T,lambda_m) - (m-1).

For a chosen leaf l with support s:
  A = T-l, B = T-{l,s},
  c1 = mu(A,lambda_m) - (m-1),
  c2 = mu(B,lambda_m) - (m-2),
  w1 = I_A(lambda_m)/I_T(lambda_m),
  w2 = lambda_m I_B(lambda_m)/I_T(lambda_m),
  combined = w1*c1 + w2*c2.

This script computes these quantities analytically for both leaf types:
  - tip leaf from a length-2 arm (requires k>=1),
  - unit leaf directly attached to hub (requires j>=1),
and scans ranges of (k,j) for positivity and anti-correlation patterns.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from typing import Any


def exp_clamped(x: float) -> float:
    if x < -745.0:
        return 0.0
    if x > 709.0:
        return float("inf")
    return math.exp(x)


def log1p_exp(x: float) -> float:
    """Stable log(1+exp(x))."""
    if x > 0:
        return x + math.log1p(math.exp(-x))
    return math.log1p(math.exp(x))


def mode_lambda_from_fg(f: list[int], g: list[int]) -> tuple[int, float, int, int]:
    """Leftmost mode and tie ratio from P = F + G."""
    deg = max(len(f), len(g)) - 1

    def coeff(i: int) -> int:
        if i < 0:
            return 0
        return (f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)

    m = 0
    best = coeff(0)
    for i in range(1, deg + 1):
        v = coeff(i)
        if v > best:
            best = v
            m = i

    if m == 0:
        return m, 0.0, 0, coeff(0)

    im1 = coeff(m - 1)
    im = coeff(m)
    lam = im1 / im
    return m, lam, im1, im


def next_f_times_1px(f: list[int]) -> list[int]:
    out = [0] * (len(f) + 1)
    for i in range(len(out)):
        a = f[i] if i < len(f) else 0
        b = f[i - 1] if i - 1 >= 0 else 0
        out[i] = a + b
    return out


def build_h_coeffs_2k(k: int) -> list[int]:
    """Coeffs of (1+2x)^k."""
    out = [0] * (k + 1)
    c = 1
    p2 = 1
    out[0] = 1
    for a in range(1, k + 1):
        c = c * (k - a + 1) // a
        p2 <<= 1
        out[a] = c * p2
    return out


def build_g_coeffs_shift_binom(k: int) -> list[int]:
    """Coeffs of x(1+x)^k."""
    out = [0] * (k + 2)
    c = 1
    out[1] = 1
    for t in range(2, k + 2):
        r = t - 1
        c = c * (k - (r - 1)) // r
        out[t] = c
    return out


def family_log_r(k: int, j: int, lam: float) -> float:
    # r = lam * (1+lam)^(k-j) / (1+2lam)^k
    return (
        math.log(lam)
        + (k - j) * math.log1p(lam)
        - k * math.log1p(2.0 * lam)
    )


def family_log_i(k: int, j: int, lam: float) -> float:
    log_a = k * math.log1p(2.0 * lam) + j * math.log1p(lam)
    log_r = family_log_r(k, j, lam)
    return log_a + log1p_exp(log_r)


def family_mu(k: int, j: int, lam: float) -> float:
    # I = A + B, with B/A = r.
    log_r = family_log_r(k, j, lam)
    r = exp_clamped(log_r)
    mu_a = k * (2.0 * lam / (1.0 + 2.0 * lam)) + j * (lam / (1.0 + lam))
    mu_b = 1.0 + k * (lam / (1.0 + lam))
    return (mu_a + r * mu_b) / (1.0 + r)


def product_mu(k: int, j: int, lam: float) -> float:
    # (1+2x)^k (1+x)^j
    return k * (2.0 * lam / (1.0 + 2.0 * lam)) + j * (lam / (1.0 + lam))


def product_log_i(k: int, j: int, lam: float) -> float:
    return k * math.log1p(2.0 * lam) + j * math.log1p(lam)


def weights_from_log_i(log_i_a: float, log_i_b: float, lam: float) -> tuple[float, float]:
    # w2/w1 = lam * I_B / I_A.
    log_q = math.log(lam) + log_i_b - log_i_a
    if log_q >= 50.0:
        w2 = 1.0
    elif log_q <= -50.0:
        w2 = 0.0
    else:
        q = math.exp(log_q)
        w2 = q / (1.0 + q)
    w1 = 1.0 - w2
    return w1, w2


def init_case_stats() -> dict[str, Any]:
    return {
        "count": 0,
        "c1_neg_count": 0,
        "min_combined": None,
        "min_c1": None,
        "min_c2": None,
        "min_margin": None,
        "min_combined_when_c1_neg": None,
        "max_w1_when_c1_neg": None,  # weight on c1 term
        "min_beta_guard_when_c1_neg": None,  # w2 - needed_w2
        "max_abs_identity_err": 0.0,  # |combined-margin|
    }


def update_min(
    stats: dict[str, Any],
    key: str,
    value: float,
    witness: dict[str, Any],
) -> None:
    cur = stats[key]
    if cur is None or value < cur["value"]:
        stats[key] = {"value": value, "witness": witness}


def main() -> None:
    ap = argparse.ArgumentParser(description="Mixed-spider combined-bound scan.")
    ap.add_argument("--k-min", type=int, default=1)
    ap.add_argument("--k-max", type=int, default=300)
    ap.add_argument("--j-min", type=int, default=0)
    ap.add_argument("--j-max", type=int, default=300)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    if args.k_min < 1:
        raise ValueError("k-min must be >= 1")
    if args.k_max < args.k_min:
        raise ValueError("k-max must be >= k-min")
    if args.j_min < 0:
        raise ValueError("j-min must be >= 0")
    if args.j_max < args.j_min:
        raise ValueError("j-max must be >= j-min")

    t0 = time.time()
    print(
        f"Mixed spider scan S(2^k,1^j): k={args.k_min}..{args.k_max}, "
        f"j={args.j_min}..{args.j_max}",
        flush=True,
    )

    tip = init_case_stats()
    unit = init_case_stats()
    global_min_combined: dict[str, Any] | None = None
    combined_fail = 0

    for k in range(args.k_min, args.k_max + 1):
        h = build_h_coeffs_2k(k)  # (1+2x)^k
        g = build_g_coeffs_shift_binom(k)  # x(1+x)^k
        f = h[:]  # j=0: (1+2x)^k

        for j in range(0, args.j_max + 1):
            if j >= args.j_min:
                m, lam, _, _ = mode_lambda_from_fg(f, g)
                if m == 0:
                    continue

                mu_t = family_mu(k, j, lam)
                margin = mu_t - (m - 1)

                # Tip leaf case: A=S(2^(k-1),1^(j+1)), B=S(2^(k-1),1^j)
                mu_a_tip = family_mu(k - 1, j + 1, lam)
                mu_b_tip = family_mu(k - 1, j, lam)
                c1_tip = mu_a_tip - (m - 1)
                c2_tip = mu_b_tip - (m - 2)
                log_i_a_tip = family_log_i(k - 1, j + 1, lam)
                log_i_b_tip = family_log_i(k - 1, j, lam)
                w1_tip, w2_tip = weights_from_log_i(log_i_a_tip, log_i_b_tip, lam)
                combined_tip = w1_tip * c1_tip + w2_tip * c2_tip
                err_tip = abs(combined_tip - margin)

                wit_tip = {
                    "k": k,
                    "j": j,
                    "n": 1 + 2 * k + j,
                    "mode": m,
                    "lambda_mode": lam,
                    "margin": margin,
                    "c1": c1_tip,
                    "c2": c2_tip,
                    "w1": w1_tip,  # weight on c1
                    "w2": w2_tip,  # weight on c2
                    "combined": combined_tip,
                    "identity_err": combined_tip - margin,
                    "leaf_type": "tip",
                }

                tip["count"] += 1
                tip["max_abs_identity_err"] = max(tip["max_abs_identity_err"], err_tip)
                update_min(tip, "min_combined", combined_tip, wit_tip)
                update_min(tip, "min_c1", c1_tip, wit_tip)
                update_min(tip, "min_c2", c2_tip, wit_tip)
                update_min(tip, "min_margin", margin, wit_tip)

                if c1_tip < -args.tol:
                    tip["c1_neg_count"] += 1
                    update_min(tip, "min_combined_when_c1_neg", combined_tip, wit_tip)
                    cur = tip["max_w1_when_c1_neg"]
                    if cur is None or w1_tip > cur["value"]:
                        tip["max_w1_when_c1_neg"] = {"value": w1_tip, "witness": wit_tip}
                    need_w2 = (-c1_tip) / (c2_tip - c1_tip) if c2_tip > c1_tip else float("inf")
                    beta_guard = w2_tip - need_w2
                    update_min(tip, "min_beta_guard_when_c1_neg", beta_guard, wit_tip)

                if combined_tip < -args.tol:
                    combined_fail += 1

                if global_min_combined is None or combined_tip < global_min_combined["combined"]:
                    global_min_combined = wit_tip

                # Unit leaf case: A=S(2^k,1^(j-1)), B=(1+2x)^k(1+x)^(j-1)
                if j >= 1:
                    mu_a_unit = family_mu(k, j - 1, lam)
                    mu_b_unit = product_mu(k, j - 1, lam)
                    c1_unit = mu_a_unit - (m - 1)
                    c2_unit = mu_b_unit - (m - 2)
                    log_i_a_unit = family_log_i(k, j - 1, lam)
                    log_i_b_unit = product_log_i(k, j - 1, lam)
                    w1_unit, w2_unit = weights_from_log_i(log_i_a_unit, log_i_b_unit, lam)
                    combined_unit = w1_unit * c1_unit + w2_unit * c2_unit
                    err_unit = abs(combined_unit - margin)

                    wit_unit = {
                        "k": k,
                        "j": j,
                        "n": 1 + 2 * k + j,
                        "mode": m,
                        "lambda_mode": lam,
                        "margin": margin,
                        "c1": c1_unit,
                        "c2": c2_unit,
                        "w1": w1_unit,  # weight on c1
                        "w2": w2_unit,  # weight on c2
                        "combined": combined_unit,
                        "identity_err": combined_unit - margin,
                        "leaf_type": "unit",
                    }

                    unit["count"] += 1
                    unit["max_abs_identity_err"] = max(unit["max_abs_identity_err"], err_unit)
                    update_min(unit, "min_combined", combined_unit, wit_unit)
                    update_min(unit, "min_c1", c1_unit, wit_unit)
                    update_min(unit, "min_c2", c2_unit, wit_unit)
                    update_min(unit, "min_margin", margin, wit_unit)

                    if c1_unit < -args.tol:
                        unit["c1_neg_count"] += 1
                        update_min(unit, "min_combined_when_c1_neg", combined_unit, wit_unit)
                        cur = unit["max_w1_when_c1_neg"]
                        if cur is None or w1_unit > cur["value"]:
                            unit["max_w1_when_c1_neg"] = {"value": w1_unit, "witness": wit_unit}
                        need_w2 = (
                            (-c1_unit) / (c2_unit - c1_unit)
                            if c2_unit > c1_unit
                            else float("inf")
                        )
                        beta_guard = w2_unit - need_w2
                        update_min(unit, "min_beta_guard_when_c1_neg", beta_guard, wit_unit)

                    if combined_unit < -args.tol:
                        combined_fail += 1

                    if (
                        global_min_combined is None
                        or combined_unit < global_min_combined["combined"]
                    ):
                        global_min_combined = wit_unit

            # Advance j -> j+1 for F_j.
            if j < args.j_max:
                f = next_f_times_1px(f)

        if (k - args.k_min) % 25 == 0 or k == args.k_max:
            tip_min = None if tip["min_combined"] is None else tip["min_combined"]["value"]
            unit_min = None if unit["min_combined"] is None else unit["min_combined"]["value"]
            print(
                f"k={k:4d}: tip_min_combined={tip_min} unit_min_combined={unit_min} "
                f"fails={combined_fail}",
                flush=True,
            )

    summary = {
        "params": vars(args),
        "overall": {
            "combined_fail": combined_fail,
            "global_min_combined": global_min_combined,
            "wall_s": time.time() - t0,
        },
        "tip_leaf": tip,
        "unit_leaf": unit,
    }

    print("-" * 96, flush=True)
    print(f"combined_fail={combined_fail}", flush=True)
    print(f"global_min_combined={global_min_combined}", flush=True)
    print(
        f"tip: min_combined={tip['min_combined']} c1_neg={tip['c1_neg_count']} "
        f"max|id_err|={tip['max_abs_identity_err']}",
        flush=True,
    )
    print(
        f"unit: min_combined={unit['min_combined']} c1_neg={unit['c1_neg_count']} "
        f"max|id_err|={unit['max_abs_identity_err']}",
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
