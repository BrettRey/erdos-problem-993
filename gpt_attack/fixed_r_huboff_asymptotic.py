"""Asymptotic hub-off reserve expansion for fixed-r spider lanes.

For F_a=P_r(x)(1+2x)^a and m=2a/3+delta, this computes the
first-order reserve constant

    C_{r,q} = lim_{a -> infinity, a=q mod 3}
              a * (mu_{P_r(1+2x)^(a-1)}(lambda_0) - m + 4/3),

where lambda_0=F_{m-1}/F_m and delta=D_{r,q}/3.

The computation uses truncated series in h=1/a and never forms the dense
rational functions in the residue parameter t.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from fixed_r_huboff_certificate import stabilized_shifts
from route2_spider_lane_scan import path_polys


Series = list[Fraction]


def trim_order(poly: Series, order: int) -> Series:
    out = poly[: order + 1]
    out.extend([Fraction(0)] * (order + 1 - len(out)))
    return out


def add(lhs: Series, rhs: Series, order: int) -> Series:
    return [lhs[i] + rhs[i] for i in range(order + 1)]


def sub(lhs: Series, rhs: Series, order: int) -> Series:
    return [lhs[i] - rhs[i] for i in range(order + 1)]


def scale(poly: Series, scalar: Fraction, order: int) -> Series:
    return [scalar * poly[i] for i in range(order + 1)]


def mul(lhs: Series, rhs: Series, order: int) -> Series:
    out = [Fraction(0)] * (order + 1)
    for i in range(order + 1):
        if lhs[i] == 0:
            continue
        for j in range(order + 1 - i):
            out[i + j] += lhs[i] * rhs[j]
    return out


def inv(poly: Series, order: int) -> Series:
    if poly[0] == 0:
        raise ZeroDivisionError("series has zero constant term")
    out = [Fraction(0)] * (order + 1)
    out[0] = 1 / poly[0]
    for n in range(1, order + 1):
        total = Fraction(0)
        for i in range(1, n + 1):
            total += poly[i] * out[n - i]
        out[n] = -total / poly[0]
    return out


def div(lhs: Series, rhs: Series, order: int) -> Series:
    return mul(lhs, inv(rhs, order), order)


def ratio_series(u: int, delta: Fraction, order: int) -> Series:
    """Return 2^-u binom(a,m-u)/binom(a,m) as a series in h=1/a."""
    out = [Fraction(1)] + [Fraction(0)] * order
    for i in range(u):
        # Normalized factor:
        # ((2/3)+h(delta-i)) / (2*((1/3)+h(-delta+i+1)))
        # = (1 + (3/2)(delta-i)h) / (1 + 3(-delta+i+1)h).
        num = [Fraction(1), Fraction(3, 2) * (delta - i)] + [Fraction(0)] * (order - 1)
        den = [Fraction(1), 3 * (-delta + i + 1)] + [Fraction(0)] * (order - 1)
        out = mul(out, div(num, den, order), order)
    return out


def coeff_ratio_sums(path_poly: list[int], delta: Fraction, order: int) -> tuple[Series, Series]:
    s_zero = [Fraction(0)] * (order + 1)
    s_minus = [Fraction(0)] * (order + 1)
    for s, coeff in enumerate(path_poly):
        if coeff == 0:
            continue
        s_zero = add(s_zero, scale(ratio_series(s, delta, order), Fraction(coeff), order), order)
        s_minus = add(
            s_minus,
            scale(ratio_series(s + 1, delta, order), Fraction(coeff), order),
            order,
        )
    return s_zero, s_minus


def eval_poly_at_lambda_series(path_poly: list[int], lam: Series, order: int) -> Series:
    out = [Fraction(0)] * (order + 1)
    power = [Fraction(1)] + [Fraction(0)] * order
    for coeff in path_poly:
        if coeff:
            out = add(out, scale(power, Fraction(coeff), order), order)
        power = mul(power, lam, order)
    return out


def path_mean_series(path_poly: list[int], lam: Series, order: int) -> Series:
    deriv = [s * coeff for s, coeff in enumerate(path_poly)][1:]
    p_val = eval_poly_at_lambda_series(path_poly, lam, order)
    dp_val = eval_poly_at_lambda_series(deriv, lam, order)
    return div(mul(lam, dp_val, order), p_val, order)


def raw_moments_at_one(path_poly: list[int]) -> tuple[Fraction, Fraction, Fraction]:
    total = sum(path_poly)
    m1 = sum(s * coeff for s, coeff in enumerate(path_poly))
    m2 = sum(s * s * coeff for s, coeff in enumerate(path_poly))
    m3 = sum(s * s * s * coeff for s, coeff in enumerate(path_poly))
    return Fraction(m1, total), Fraction(m2, total), Fraction(m3, total)


def central_moments_at_one(path_poly: list[int]) -> tuple[Fraction, Fraction, Fraction]:
    m1, m2, m3 = raw_moments_at_one(path_poly)
    variance = m2 - m1**2
    third_central = m3 - 3 * m1 * m2 + 2 * m1**3
    return m1, variance, third_central


def first_order_shift_candidates(path_poly: list[int]) -> dict[int, list[int]]:
    """Return shifts allowed by the first-order F-mode inequalities.

    The asymptotic adjacent checks for F=P_r(1+2x)^a give

        M_1 - 1/3 <= delta <= M_1 + 2/3,

    where delta=D/3.  For residue q, D must be congruent to q modulo 3.
    Boundary cases may leave two candidates at first order.
    """
    m1, _, _ = raw_moments_at_one(path_poly)
    lo = 3 * m1 - 1
    hi = 3 * m1 + 2
    candidates: dict[int, list[int]] = {}
    start = lo.numerator // lo.denominator - 6
    stop = hi.numerator // hi.denominator + 7
    for residue in [0, 1, 2]:
        candidates[residue] = [
            shift
            for shift in range(start, stop + 1)
            if shift % 3 == residue and lo <= shift <= hi
        ]
    return candidates


def reserve_constant_moment_formula(path_poly: list[int], shift: int) -> Fraction:
    """Return C_{r,q}=lim a*reserve using path moments at x=1.

    If delta=shift/3 and M_j=E[S^j] for S distributed by coefficients of P_r
    at x=1, then

        C = (-24 delta + 54 M_1^3 - 9 M_1^2 - 81 M_1 M_2
             + 24 M_1 + 9 M_2 + 27 M_3 + 16) / 12.

    Equivalently, with variance V and third central moment K_3,

        C = (9V + 27K_3 + 24(M_1-delta) + 16) / 12.
    """
    delta = Fraction(shift, 3)
    m1, m2, m3 = raw_moments_at_one(path_poly)
    return (
        -24 * delta
        + 54 * m1**3
        - 9 * m1**2
        - 81 * m1 * m2
        + 24 * m1
        + 9 * m2
        + 27 * m3
        + 16
    ) / 12


def huboff_asymptotic_for_shift(path_poly: list[int], shift: int) -> dict:
    order = 2
    delta = Fraction(shift, 3)
    s_zero, s_minus = coeff_ratio_sums(path_poly, delta, order)
    lam = div(s_minus, s_zero, order)
    two_lam = scale(lam, Fraction(2), order)
    denom = add([Fraction(1)] + [Fraction(0)] * order, two_lam, order)
    arm_rate = div(two_lam, denom, order)  # 2 lambda / (1+2 lambda)
    p_mean = path_mean_series(path_poly, lam, order)

    # (a-1) arm_rate = h^-1 arm_rate - arm_rate.
    # reserve = p_mean + (a-1)arm_rate - (2/(3h)+delta) + 4/3.
    h_minus_one_coeff = arm_rate[0] - Fraction(2, 3)
    constant = p_mean[0] + arm_rate[1] - arm_rate[0] - delta + Fraction(4, 3)
    h_coeff = p_mean[1] + arm_rate[2] - arm_rate[1]
    moment_formula = reserve_constant_moment_formula(path_poly, shift)
    if h_coeff != moment_formula:
        raise AssertionError((h_coeff, moment_formula))
    m1, m2, m3 = raw_moments_at_one(path_poly)
    _, variance, third_central = central_moments_at_one(path_poly)
    return {
        "shift": shift,
        "delta": str(delta),
        "path_moment_1": str(m1),
        "path_moment_2": str(m2),
        "path_moment_3": str(m3),
        "path_variance": str(variance),
        "path_third_central_moment": str(third_central),
        "m1_minus_delta": str(m1 - delta),
        "lambda_minus_1_times_a": str(lam[1]),
        "lambda_h2_coeff": str(lam[2]),
        "h_minus_one_coeff": str(h_minus_one_coeff),
        "constant_coeff": str(constant),
        "reserve_times_a_limit": str(h_coeff),
        "reserve_times_a_limit_moment_formula": str(moment_formula),
        "reserve_times_a_limit_float": float(h_coeff),
    }


def parse_r_values(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description="Fixed-r hub-off reserve asymptotics.")
    ap.add_argument("--r-values", default="8,20,80,120,160,200,240,280,320")
    ap.add_argument("--threshold", type=int, default=500)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    r_values = parse_r_values(args.r_values)
    paths = path_polys(max(r_values))
    records = []
    for r in r_values:
        shifts = stabilized_shifts(r, args.threshold, paths)
        shift_candidates = first_order_shift_candidates(paths[r])
        shift_matches = {
            residue: shifts[residue] in shift_candidates[residue] for residue in [0, 1, 2]
        }
        residues = []
        for residue in [0, 1, 2]:
            row = huboff_asymptotic_for_shift(paths[r], shifts[residue])
            row["residue"] = residue
            residues.append(row)
        record = {
            "r": r,
            "threshold_for_shift": args.threshold,
            "shifts": shifts,
            "first_order_shift_candidates": shift_candidates,
            "shift_matches_first_order": shift_matches,
            "residues": residues,
        }
        records.append(record)
        print(f"r={r}, shifts={shifts}")
        print(f"  first-order candidates={shift_candidates}, matches={shift_matches}")
        for row in residues:
            print(
                "  q={residue}: a(lambda-1)->{lam}, a*reserve->{res:.12g}".format(
                    residue=row["residue"],
                    lam=row["lambda_minus_1_times_a"],
                    res=row["reserve_times_a_limit_float"],
                )
            )

    result = {"params": {"r_values": r_values, "threshold": args.threshold}, "records": records}
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
