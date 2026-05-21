"""Closed-form Route-2 diagnostics for the pure spider lane S(2^a).

This isolates the large-a pattern behind the spider-lane scan.  For the spider
with ``a`` arms of length 2,

    I_a(x) = (1 + 2x)^a + x(1 + x)^a.

For Route-2, remove one length-2 arm, so B = S(2^(a-1)).  The script checks
the exact slack

    Delta_a = mu_B(lambda_a) - (m_a - 3/2),

where ``m_a`` is the leftmost mode of I_a and
``lambda_a = [x^(m_a-1)] I_a / [x^m_a] I_a``.

The binomial-spine approximation drops the exponentially small hub-on term.
It predicts

    Delta_a^0 = 3/2 - 2 m_a / (a + 1),

which tends to 1/6 from above in every residue class.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def coeff(a: int, k: int) -> int:
    """Coefficient of x^k in I(S(2^a))."""
    off = comb(a, k) * (2**k) if 0 <= k <= a else 0
    on = comb(a, k - 1) if 0 <= k - 1 <= a else 0
    return off + on


def mode_formula(a: int) -> int:
    """Leftmost mode formula proved in the accompanying note."""
    if a == 2:
        return 2
    return (2 * a + 1) // 3


def mode_exact(a: int) -> int:
    vals = [coeff(a, k) for k in range(a + 2)]
    mx = max(vals)
    return next(k for k, v in enumerate(vals) if v == mx)


def mean_spider_fraction(a: int, lam: Fraction) -> Fraction:
    """Exact mean independent-set size for S(2^a) at fugacity lam."""
    off_weight = (1 + 2 * lam) ** a
    on_weight = lam * (1 + lam) ** a
    off_mean = Fraction(2 * a) * lam / (1 + 2 * lam)
    on_mean = 1 + Fraction(a) * lam / (1 + lam)
    return (off_weight * off_mean + on_weight * on_mean) / (off_weight + on_weight)


def exact_record(a: int) -> dict:
    m = mode_exact(a)
    lam = Fraction(coeff(a, m - 1), coeff(a, m))
    slack = mean_spider_fraction(a - 1, lam) - Fraction(2 * m - 3, 2)
    spine = Fraction(3, 2) - Fraction(2 * m, a + 1)
    return {
        "a": a,
        "m": m,
        "lambda": lam,
        "slack": slack,
        "spine_slack": spine,
        "slack_minus_one_sixth": slack - Fraction(1, 6),
        "error_vs_spine": slack - spine,
    }


def coeff_ratio_epsilon(a: int, k: int) -> float:
    """Hub-on/off coefficient ratio C(a,k-1)/(2^k C(a,k))."""
    if not (1 <= k <= a):
        return 0.0
    log_eps = math.log(k) - math.log(a - k + 1) - k * math.log(2.0)
    if log_eps < -745.0:
        return 0.0
    return math.exp(log_eps)


def lambda_float(a: int, m: int) -> float:
    """Stable float lambda from coefficient ratios without large integers."""
    lam0 = m / (2.0 * (a - m + 1))
    return lam0 * (1.0 + coeff_ratio_epsilon(a, m - 1)) / (
        1.0 + coeff_ratio_epsilon(a, m)
    )


def mean_spider_float(a: int, lam: float) -> float:
    """Stable float mean for S(2^a) at fugacity lam."""
    off_mean = 2.0 * a * lam / (1.0 + 2.0 * lam)
    on_mean = 1.0 + a * lam / (1.0 + lam)
    log_ratio = math.log(lam) + a * (
        math.log1p(lam) - math.log1p(2.0 * lam)
    )
    ratio = 0.0 if log_ratio < -745.0 else math.exp(log_ratio)
    return (off_mean + ratio * on_mean) / (1.0 + ratio)


def float_record(a: int) -> dict:
    m = mode_formula(a)
    lam = lambda_float(a, m)
    slack = mean_spider_float(a - 1, lam) - (m - 1.5)
    spine = 1.5 - 2.0 * m / (a + 1)
    return {
        "a": a,
        "m": m,
        "lambda": lam,
        "slack": slack,
        "spine_slack": spine,
        "slack_minus_one_sixth": slack - 1.0 / 6.0,
        "error_vs_spine": slack - spine,
    }


def as_jsonable(rec: dict) -> dict:
    out = {}
    for key, value in rec.items():
        if isinstance(value, Fraction):
            out[key] = {
                "num": value.numerator,
                "den": value.denominator,
                "float": float(value),
            }
        else:
            out[key] = value
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Pure S(2^a) Route-2 asymptotic check.")
    ap.add_argument("--exact-max", type=int, default=250)
    ap.add_argument("--float-max", type=int, default=5000)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    mode_failures = [
        (a, mode_formula(a), mode_exact(a))
        for a in range(2, args.exact_max + 1)
        if mode_formula(a) != mode_exact(a)
    ]

    exact_rows = [exact_record(a) for a in range(2, args.exact_max + 1)]
    min_exact = min(exact_rows, key=lambda r: r["slack"])
    min_margin = min(exact_rows, key=lambda r: r["slack_minus_one_sixth"])
    min_error = min(exact_rows, key=lambda r: r["error_vs_spine"])

    sample_as = sorted(
        {
            2,
            3,
            5,
            8,
            12,
            20,
            50,
            100,
            199,
            500,
            1000,
            args.float_max,
        }
    )
    samples = [float_record(a) for a in sample_as if a <= args.float_max]

    print(f"checked exact a=2..{args.exact_max}")
    print(f"mode formula failures: {len(mode_failures)}")
    if mode_failures:
        print(mode_failures[:10])
    print(
        "min exact slack: "
        f"a={min_exact['a']} m={min_exact['m']} "
        f"slack={float(min_exact['slack']):.12g}"
    )
    print(
        "min exact slack-1/6: "
        f"a={min_margin['a']} margin={float(min_margin['slack_minus_one_sixth']):.12g}"
    )
    print(
        "min exact error vs spine: "
        f"a={min_error['a']} error={float(min_error['error_vs_spine']):.12g}"
    )
    print("\nsamples")
    for rec in samples:
        print(
            f"a={rec['a']:5d} m={rec['m']:5d} "
            f"lambda={rec['lambda']:.12g} "
            f"slack={rec['slack']:.12g} "
            f"spine={rec['spine_slack']:.12g} "
            f"slack-1/6={rec['slack_minus_one_sixth']:.12g}"
        )

    if args.out:
        payload = {
            "params": {"exact_max": args.exact_max, "float_max": args.float_max},
            "mode_failures": mode_failures,
            "min_exact_slack": as_jsonable(min_exact),
            "min_exact_margin_over_one_sixth": as_jsonable(min_margin),
            "min_exact_error_vs_spine": as_jsonable(min_error),
            "samples": [as_jsonable(r) for r in samples],
        }
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
