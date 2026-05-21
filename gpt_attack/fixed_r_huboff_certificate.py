"""General fixed-r hub-off certificate for lanes S(2^a,r).

For fixed r, write

    I_{a,r}(x) = P_r(x)(1+2x)^a + x P_{r-1}(x)(1+x)^a.

This script extracts the stabilized full-polynomial mode shifts and then
checks the symbolic hub-off inequalities needed by the fixed-r proof template:

* left/right margins for F_a=P_r(x)(1+2x)^a;
* hub-off Route-2 reserve for B's hub-off approximation;
* a compact fugacity interval for lambda_0.

It does not prove the hub-on perturbation bound.  Use it to identify which
fixed r lanes have the same proof shape as r=8.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import ceil
from math import gcd
from math import lcm

import sympy as sp

from route2_spider_lane_scan import (
    mode_left,
    path_polys,
    route2_float_for_removed_arm,
    spider_poly_from_counts,
)


t = sp.symbols("t", integer=True, positive=True)


def stabilized_shifts(r: int, threshold: int, paths: list[list[int]]) -> dict[int, int]:
    shifts: dict[int, int] = {}
    for residue in [0, 1, 2]:
        a = threshold
        while a % 3 != residue:
            a += 1
        counts = {2: a}
        if r == 2:
            counts[2] += 1
        else:
            counts[r] = 1
        m = mode_left(spider_poly_from_counts(counts, paths))
        shifts[residue] = 3 * m - 2 * a
    return shifts


def binom_offset_ratio(a, m, offset: int):
    """Return binom(a, m+offset) / binom(a, m) as a rational function."""
    if offset == 0:
        return sp.Integer(1)
    out = sp.Integer(1)
    if offset > 0:
        for i in range(1, offset + 1):
            out *= (a - m - i + 1) / (m + i)
    else:
        for i in range(1, -offset + 1):
            out *= (m - i + 1) / (a - m + i)
    return sp.cancel(out)


def normalized_coeff_ratio(path_poly: list[int], a, m, shift: int):
    """Return F_{m+shift}/(2^m binom(a,m)) for F=P_r(1+2x)^a."""
    ratio = binom_offset_ratio(a, m, shift)
    total = sp.Integer(0)
    for j, coeff in enumerate(path_poly):
        if coeff:
            total += sp.Rational(2) ** (shift - j) * coeff * ratio
        if j + 1 < len(path_poly):
            offset = shift - j
            ratio *= (m + offset) / (a - m - offset + 1)
    return sp.cancel(sp.together(total))


def shifted_poly(poly: sp.Poly, threshold_t: int) -> sp.Poly:
    """Return poly(u+threshold_t), represented using the global symbol t."""
    degree = poly.degree()
    coeffs = [poly.nth(i) for i in range(degree + 1)]
    threshold = sp.Integer(threshold_t)
    shifted_coeffs = []
    for power in range(degree + 1):
        shifted_coeff = sp.Integer(0)
        for i in range(power, degree + 1):
            shifted_coeff += coeffs[i] * sp.binomial(i, power) * threshold ** (i - power)
        shifted_coeffs.append(shifted_coeff)
    return sp.Poly.from_list(list(reversed(shifted_coeffs)), gens=t)


def poly_has_positive_coefficients(poly: sp.Poly) -> bool:
    return all(poly.nth(i) > 0 for i in range(poly.degree() + 1))


def shifted_poly_positive(poly: sp.Poly, threshold_t: int) -> bool:
    """Check that poly(u+threshold_t) has strictly positive coefficients."""
    return poly_has_positive_coefficients(shifted_poly(poly, threshold_t))


def shifted_positive(expr, threshold_t: int) -> bool:
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    return shifted_poly_positive(sp.Poly(num, t), threshold_t) and shifted_poly_positive(
        sp.Poly(den, t), threshold_t
    )


def trim_int_poly(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def add_int_poly(lhs: list[int], rhs: list[int]) -> list[int]:
    degree = max(len(lhs), len(rhs))
    out = [0] * degree
    for i, coeff in enumerate(lhs):
        out[i] += coeff
    for i, coeff in enumerate(rhs):
        out[i] += coeff
    return trim_int_poly(out)


def sub_int_poly(lhs: list[int], rhs: list[int]) -> list[int]:
    degree = max(len(lhs), len(rhs))
    out = [0] * degree
    for i, coeff in enumerate(lhs):
        out[i] += coeff
    for i, coeff in enumerate(rhs):
        out[i] -= coeff
    return trim_int_poly(out)


def scale_int_poly(poly: list[int], scalar: int) -> list[int]:
    if scalar == 0:
        return [0]
    return trim_int_poly([scalar * coeff for coeff in poly])


def mul_int_poly(lhs: list[int], rhs: list[int]) -> list[int]:
    if lhs == [0] or rhs == [0]:
        return [0]
    if len(lhs) > len(rhs):
        lhs, rhs = rhs, lhs
    out = [0] * (len(lhs) + len(rhs) - 1)
    for i, left_coeff in enumerate(lhs):
        if left_coeff == 0:
            continue
        for j, right_coeff in enumerate(rhs):
            if right_coeff:
                out[i + j] += left_coeff * right_coeff
    return trim_int_poly(out)


def poly_pair_common_integer_coeffs(lhs: sp.Poly, rhs: sp.Poly) -> tuple[list[int], list[int]]:
    """Return low-to-high integer coefficient lists after one common scaling."""
    common_den = 1
    for poly in [lhs, rhs]:
        for i in range(poly.degree() + 1):
            common_den = lcm(common_den, int(sp.denom(poly.nth(i))))

    def convert(poly: sp.Poly) -> list[int]:
        return trim_int_poly(
            [int(poly.nth(i) * common_den) for i in range(poly.degree() + 1)]
        )

    lhs_coeffs = convert(lhs)
    rhs_coeffs = convert(rhs)
    common_gcd = 0
    for coeff in lhs_coeffs + rhs_coeffs:
        common_gcd = gcd(common_gcd, abs(coeff))
    if common_gcd > 1:
        lhs_coeffs = [coeff // common_gcd for coeff in lhs_coeffs]
        rhs_coeffs = [coeff // common_gcd for coeff in rhs_coeffs]
    return lhs_coeffs, rhs_coeffs


def shifted_linear_int_poly(expr) -> list[int]:
    poly = shifted_poly(sp.Poly(expr, t), 0)
    return [int(poly.nth(i)) for i in range(poly.degree() + 1)]


def reserve_margin_positive(
    path_poly: list[int],
    a,
    m,
    lam0,
    reserve_denom: int,
    threshold_t: int,
) -> bool:
    """Certify huboff_mean - m + 4/3 >= 1/(reserve_denom*a).

    Direct substitution of lambda_0 into P'/P creates large nested rational
    expressions for larger r.  This keeps the expression as polynomial
    arithmetic.  If lambda=L/Z and E_n=Z^deg(P_n) P_n(L/Z), then the path
    recurrence P_n=P_{n-1}+xP_{n-2} gives E_n without expanding all powers of
    L and Z.  A companion recurrence gives N_n=Z^deg(P_n) lambda P_n'(lambda).
    """
    lam_num, lam_den = sp.fraction(sp.cancel(lam0))
    l_poly = shifted_poly(sp.Poly(lam_num, t), threshold_t)
    z_poly = shifted_poly(sp.Poly(lam_den, t), threshold_t)
    l_coeffs, z_coeffs = poly_pair_common_integer_coeffs(l_poly, z_poly)
    r = path_poly[1] if len(path_poly) > 1 else 0

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

    a_poly = shifted_linear_int_poly(a.subs(t, t + threshold_t))
    m_poly = shifted_linear_int_poly(m.subs(t, t + threshold_t))
    arm_mean_num = scale_int_poly(
        mul_int_poly(sub_int_poly(a_poly, [1]), l_coeffs),
        2,
    )
    arm_mean_den = add_int_poly(z_coeffs, scale_int_poly(l_coeffs, 2))
    target_den3 = scale_int_poly(a_poly, 3 * reserve_denom)
    target_num3 = sub_int_poly(
        scale_int_poly(
            mul_int_poly(a_poly, sub_int_poly([4], scale_int_poly(m_poly, 3))),
            reserve_denom,
        ),
        [3],
    )

    numerator = (
        add_int_poly(
            add_int_poly(
                mul_int_poly(mul_int_poly(arm_mean_num, path_eval), target_den3),
                mul_int_poly(mul_int_poly(log_deriv_num, arm_mean_den), target_den3),
            ),
            mul_int_poly(mul_int_poly(target_num3, arm_mean_den), path_eval),
        )
    )
    denominator = mul_int_poly(mul_int_poly(arm_mean_den, path_eval), target_den3)
    return all(coeff > 0 for coeff in numerator) and all(coeff > 0 for coeff in denominator)


def check_r(
    r: int,
    threshold: int,
    margin_denom: int,
    reserve_denom: int,
    exact_max: int,
) -> bool:
    paths = path_polys(max(2, r))
    pr = paths[r]
    shifts = stabilized_shifts(r, threshold, paths)
    print(f"r={r}, threshold a>={threshold}")
    print(f"stabilized shifts D_q=3m-2a: {shifts}")

    if exact_max:
        min_slack = None
        failures = []
        for a0 in range(1, exact_max + 1):
            counts = {2: a0}
            if r == 2:
                counts[2] += 1
            else:
                counts[r] = 1
            rec = route2_float_for_removed_arm(counts, 2, paths)
            if rec["route2_slack"] <= 0:
                failures.append(a0)
            if min_slack is None or rec["route2_slack"] < min_slack[1]:
                min_slack = (a0, rec["route2_slack"], rec["m"])
        print(
            f"float finite check a=1..{exact_max}: failures={len(failures)}, "
            f"min a={min_slack[0]} slack={min_slack[1]:.12g} m={min_slack[2]}"
        )

    all_ok = True
    for residue, shift in shifts.items():
        threshold_t = ceil((threshold - residue) / 3)
        a = 3 * t + residue
        m = sp.simplify((2 * a + shift) / 3)
        if not m.is_integer:
            print(f"  q={residue}: nonintegral mode expression {m}")
            all_ok = False
            continue
        f_minus = normalized_coeff_ratio(pr, a, m, -1)
        f_zero = normalized_coeff_ratio(pr, a, m, 0)
        f_plus = normalized_coeff_ratio(pr, a, m, 1)
        left_margin = sp.cancel((f_zero - f_minus) / f_zero)
        right_margin = sp.cancel((f_zero - f_plus) / f_zero)
        lam0 = sp.cancel(f_minus / f_zero)
        checks = {
            f"left_margin_ge_1_over_{margin_denom}a": left_margin
            - Fraction(1, margin_denom) / a,
            f"right_margin_ge_1_over_{margin_denom}a": right_margin
            - Fraction(1, margin_denom) / a,
            "lambda0_ge_3_over_4": lam0 - sp.Rational(3, 4),
            "lambda0_ge_1_over_2": lam0 - sp.Rational(1, 2),
            "lambda0_le_2": sp.Rational(2) - lam0,
        }
        print(f"  residue q={residue}, t>={threshold_t}, m={m}")
        for name, expr in checks.items():
            ok = shifted_positive(expr, threshold_t)
            all_ok = all_ok and ok
            print(f"    {name}: {ok}")
        reserve_ok = reserve_margin_positive(pr, a, m, lam0, reserve_denom, threshold_t)
        all_ok = all_ok and reserve_ok
        print(f"    reserve_ge_1_over_{reserve_denom}a: {reserve_ok}")
    return all_ok


def main() -> None:
    ap = argparse.ArgumentParser(description="Fixed-r hub-off symbolic certificate.")
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--threshold", type=int, default=200)
    ap.add_argument("--margin-denom", type=int, default=1000)
    ap.add_argument("--reserve-denom", type=int, default=1000)
    ap.add_argument("--exact-max", type=int, default=0)
    args = ap.parse_args()

    ok = check_r(args.r, args.threshold, args.margin_denom, args.reserve_denom, args.exact_max)
    assert ok


if __name__ == "__main__":
    main()
