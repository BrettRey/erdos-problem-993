"""Effective-threshold probe for the fixed-r eventual Route-2 schema.

This script tests the theorem-schema route:

* choose the first-order stabilized shift D for each residue;
* build exact rational functions for the F adjacent margins, lambda_0 interval,
  and hub-off reserve;
* choose half of the positive leading constants as the margin/reserve targets;
* use exact real-root isolation, shifted-coefficient positivity, or Cauchy
  root bounds to find a threshold after which each rational function is
  positive.

It is a proof-engineering probe, not a replacement for the existing large-lane
coefficient certificates.  The direct rational functions are useful for small
and moderate r; for large r the recurrence/GMP path remains more practical.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import floor
from pathlib import Path

import sympy as sp

from fixed_r_huboff_asymptotic import (
    first_order_shift_candidates,
    reserve_constant_moment_formula,
)
from fixed_r_huboff_certificate import (
    add_int_poly,
    mul_int_poly,
    normalized_coeff_ratio,
    poly_pair_common_integer_coeffs,
    scale_int_poly,
    sub_int_poly,
    t,
)
from route2_spider_lane_scan import path_polys


def rational(frac: Fraction) -> sp.Rational:
    return sp.Rational(frac.numerator, frac.denominator)


def parse_r_values(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part.strip()]


def canonical_shift(path_poly: list[int], r: int, residue: int) -> int:
    """Return the shift predicted by the first-order rule.

    For r>=4 the candidate is unique.  The first-order boundary cases r=2,3
    are handled by the exact lower-candidate choices already identified.
    """
    if r == 2:
        return {0: 3, 1: 1, 2: 2}[residue]
    if r == 3:
        return {0: 3, 1: 4, 2: 2}[residue]
    candidates = first_order_shift_candidates(path_poly)[residue]
    if len(candidates) != 1:
        raise ValueError((r, residue, candidates))
    return candidates[0]


def cauchy_threshold(poly: sp.Poly, min_t: int) -> int:
    """Cheap Cauchy upper bound for positive roots."""
    coeffs = poly.all_coeffs()
    leading = coeffs[0]
    if leading <= 0:
        raise ValueError(f"nonpositive leading coefficient: {leading}")
    if len(coeffs) == 1:
        return min_t
    bound = sp.Rational(1) + max(sp.Rational(abs(coeff), leading) for coeff in coeffs[1:])
    return max(min_t, bound.p // bound.q + 1)


def shifted_coefficients_positive(poly: sp.Poly, shift: int) -> bool:
    """Return whether poly(u+shift) has strictly positive coefficients."""
    shift_int = sp.Integer(shift)
    # Horner-Taylor shift.  If result is the low-to-high coefficient list for
    # the processed high-degree suffix, then multiplying by (u+shift) sends
    # result[j] to shift*result[j] + result[j-1].
    shifted = [poly.LC()]
    for coeff in poly.all_coeffs()[1:]:
        nxt = [sp.Integer(0)] * (len(shifted) + 1)
        for i, value in enumerate(shifted):
            nxt[i] += shift_int * value
            nxt[i + 1] += value
        nxt[0] += coeff
        shifted = nxt
    return all(coeff > 0 for coeff in shifted)


def shifted_coefficient_threshold(poly: sp.Poly, min_t: int) -> int:
    """Find T such that poly(u+T) has strictly positive coefficients."""
    if shifted_coefficients_positive(poly, min_t):
        return min_t

    hi = max(1, min_t + 1)
    while not shifted_coefficients_positive(poly, hi):
        hi *= 2

    lo = min_t
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if shifted_coefficients_positive(poly, mid):
            hi = mid
        else:
            lo = mid
    return hi


def positive_root_threshold(poly: sp.Poly, min_t: int, method: str) -> int:
    """Return T such that poly(n)>0 for every integer n>=T.

    The polynomial must have positive leading coefficient.  `isolate` uses
    exact real-root isolation and verifies by evaluation.  `shift` searches
    for shifted-coefficient positivity.  `cauchy` uses the cheaper Cauchy
    upper bound for all positive roots.
    """
    if poly.is_zero:
        raise ValueError("zero polynomial")
    if poly.LC() <= 0:
        raise ValueError(f"nonpositive leading coefficient: {poly.LC()}")
    if method == "shift":
        return shifted_coefficient_threshold(poly, min_t)
    if method == "cauchy":
        return cauchy_threshold(poly, min_t)
    if method != "isolate":
        raise ValueError(method)

    threshold = min_t
    for interval, _mult in poly.intervals():
        lo, hi = interval
        if hi < min_t:
            continue
        if hi >= 0:
            hi_rat = sp.Rational(hi)
            threshold = max(threshold, hi_rat.p // hi_rat.q + 1)

    while poly.eval(threshold) <= 0:
        threshold += 1
    return threshold


def rational_positive_threshold_from_polys(
    num_poly: sp.Poly,
    den_poly: sp.Poly,
    min_t: int,
    method: str,
) -> tuple[int, dict]:
    """Return an eventual positivity threshold for numerator/denominator polys."""
    if den_poly.LC() < 0:
        num_poly = -num_poly
        den_poly = -den_poly

    num_t = positive_root_threshold(num_poly, min_t, method)
    den_t = positive_root_threshold(den_poly, min_t, method)
    threshold = max(num_t, den_t)
    if method == "isolate":
        while num_poly.eval(threshold) * den_poly.eval(threshold) <= 0:
            threshold += 1
    return threshold, {
        "method": method,
        "num_degree": num_poly.degree(),
        "den_degree": den_poly.degree(),
        "num_lc": str(num_poly.LC()),
        "den_lc": str(den_poly.LC()),
        "num_threshold_t": num_t,
        "den_threshold_t": den_t,
    }


def rational_positive_threshold(expr, min_t: int, method: str) -> tuple[int, dict]:
    """Return an eventual positivity threshold for a rational expression."""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    return rational_positive_threshold_from_polys(
        sp.Poly(num, t),
        sp.Poly(den, t),
        min_t,
        method,
    )


def leading_constant_for_a_times(expr, a) -> sp.Rational:
    return sp.cancel(sp.limit(a * expr, t, sp.oo))


def poly_from_low_coeffs(coeffs: list[int]) -> sp.Poly:
    return sp.Poly.from_list(list(reversed(coeffs)), gens=t)


def linear_int_poly(expr) -> list[int]:
    poly = sp.Poly(expr, t)
    return [int(poly.nth(i)) for i in range(poly.degree() + 1)]


def reserve_polys_from_recurrence(
    r: int,
    lam0,
    a,
    m,
    eta: sp.Rational,
) -> tuple[sp.Poly, sp.Poly]:
    """Return numerator/denominator for reserve - eta/a.

    This mirrors the cleared recurrence used by the hub-off certificate, but
    allows an arbitrary rational eta instead of only 1/(R a).
    """
    lam_num, lam_den = sp.fraction(sp.cancel(lam0))
    l_poly = sp.Poly(lam_num, t)
    z_poly = sp.Poly(lam_den, t)
    l_coeffs, z_coeffs = poly_pair_common_integer_coeffs(l_poly, z_poly)

    eval_prev2 = [1]  # E_0
    deriv_prev2 = [0]  # N_0
    if r == 0:
        path_eval = eval_prev2
        log_deriv_num = deriv_prev2
    else:
        eval_prev1 = add_int_poly(z_coeffs, l_coeffs)  # E_1
        deriv_prev1 = l_coeffs[:]  # N_1
        if r == 1:
            path_eval = eval_prev1
            log_deriv_num = deriv_prev1
        else:
            for n in range(2, r + 1):
                if n % 2 == 0:
                    path_eval = add_int_poly(
                        eval_prev1,
                        mul_int_poly(l_coeffs, eval_prev2),
                    )
                    log_deriv_num = add_int_poly(
                        deriv_prev1,
                        mul_int_poly(l_coeffs, add_int_poly(eval_prev2, deriv_prev2)),
                    )
                else:
                    path_eval = add_int_poly(
                        mul_int_poly(z_coeffs, eval_prev1),
                        mul_int_poly(l_coeffs, eval_prev2),
                    )
                    log_deriv_num = add_int_poly(
                        mul_int_poly(z_coeffs, deriv_prev1),
                        mul_int_poly(l_coeffs, add_int_poly(eval_prev2, deriv_prev2)),
                    )
                eval_prev2, eval_prev1 = eval_prev1, path_eval
                deriv_prev2, deriv_prev1 = deriv_prev1, log_deriv_num

    eta_num = int(eta.p)
    eta_den = int(eta.q)
    a_poly = linear_int_poly(a)
    m_poly = linear_int_poly(m)
    arm_mean_num = scale_int_poly(
        mul_int_poly(sub_int_poly(a_poly, [1]), l_coeffs),
        2,
    )
    arm_mean_den = add_int_poly(z_coeffs, scale_int_poly(l_coeffs, 2))
    target_den = scale_int_poly(a_poly, 3 * eta_den)
    target_num = sub_int_poly(
        scale_int_poly(
            mul_int_poly(a_poly, sub_int_poly([4], scale_int_poly(m_poly, 3))),
            eta_den,
        ),
        [3 * eta_num],
    )

    numerator = add_int_poly(
        add_int_poly(
            mul_int_poly(mul_int_poly(arm_mean_num, path_eval), target_den),
            mul_int_poly(mul_int_poly(log_deriv_num, arm_mean_den), target_den),
        ),
        mul_int_poly(mul_int_poly(target_num, arm_mean_den), path_eval),
    )
    denominator = mul_int_poly(mul_int_poly(arm_mean_den, path_eval), target_den)
    return poly_from_low_coeffs(numerator), poly_from_low_coeffs(denominator)


def residue_record(
    r: int,
    residue: int,
    path_poly: list[int],
    min_t: int,
    method: str,
) -> dict:
    shift = canonical_shift(path_poly, r, residue)
    a = 3 * t + residue
    m = sp.simplify((2 * a + shift) / 3)
    if not m.is_integer:
        raise ValueError((r, residue, shift, m))

    f_minus = normalized_coeff_ratio(path_poly, a, m, -1)
    f_zero = normalized_coeff_ratio(path_poly, a, m, 0)
    f_plus = normalized_coeff_ratio(path_poly, a, m, 1)
    left_margin = sp.cancel((f_zero - f_minus) / f_zero)
    right_margin = sp.cancel((f_zero - f_plus) / f_zero)
    lam0 = sp.cancel(f_minus / f_zero)

    left_const = leading_constant_for_a_times(left_margin, a)
    right_const = leading_constant_for_a_times(right_margin, a)
    if left_const <= 0 or right_const <= 0:
        raise ValueError((r, residue, shift, left_const, right_const))
    left_eta = left_const / 2
    right_eta = right_const / 2

    reserve_const = rational(reserve_constant_moment_formula(path_poly, shift))
    if reserve_const <= 0:
        raise ValueError((r, residue, shift, reserve_const))
    reserve_eta = reserve_const / 2
    checks = {
        "left_margin_ge_eta_over_a": left_margin - left_eta / a,
        "right_margin_ge_eta_over_a": right_margin - right_eta / a,
        "lambda0_ge_3_over_4": lam0 - sp.Rational(3, 4),
        "lambda0_le_2": sp.Rational(2) - lam0,
    }

    check_records = {}
    threshold_t = min_t
    for name, expr in checks.items():
        check_t, meta = rational_positive_threshold(expr, min_t, method)
        check_records[name] = {"threshold_t": check_t, **meta}
        threshold_t = max(threshold_t, check_t)

    reserve_num, reserve_den = reserve_polys_from_recurrence(r, lam0, a, m, reserve_eta)
    reserve_t, reserve_meta = rational_positive_threshold_from_polys(
        reserve_num,
        reserve_den,
        min_t,
        method,
    )
    check_records["reserve_ge_eta_over_a"] = {"threshold_t": reserve_t, **reserve_meta}
    threshold_t = max(threshold_t, reserve_t)

    return {
        "residue": residue,
        "shift": shift,
        "mode": str(m),
        "left_margin_constant": str(left_const),
        "right_margin_constant": str(right_const),
        "reserve_constant": str(reserve_const),
        "left_eta": str(left_eta),
        "right_eta": str(right_eta),
        "reserve_eta": str(reserve_eta),
        "threshold_t": threshold_t,
        "threshold_a": 3 * threshold_t + residue,
        "checks": check_records,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Fixed-r exact root-threshold probe.")
    ap.add_argument("--r-values", default="4,8,20")
    ap.add_argument("--min-t", type=int, default=0)
    ap.add_argument("--method", choices=["isolate", "shift", "cauchy"], default="isolate")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    r_values = parse_r_values(args.r_values)
    paths = path_polys(max(r_values))
    records = []
    for r in r_values:
        print(f"r={r}")
        residues = []
        for residue in [0, 1, 2]:
            row = residue_record(r, residue, paths[r], args.min_t, args.method)
            residues.append(row)
            print(
                "  q={q}: D={D}, A={A}, C={C}".format(
                    q=residue,
                    D=row["shift"],
                    A=row["threshold_a"],
                    C=row["reserve_constant"],
                )
            )
        record = {
            "r": r,
            "min_t": args.min_t,
            "threshold_a": max(row["threshold_a"] for row in residues),
            "residues": residues,
        }
        records.append(record)
        print(f"  combined threshold A={record['threshold_a']}")

    result = {
        "params": {"r_values": r_values, "min_t": args.min_t, "method": args.method},
        "records": records,
    }
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
