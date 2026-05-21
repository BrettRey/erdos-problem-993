"""Exact finite Route-2 checks for fixed-r spider lanes S(2^a,r).

For each lane, this verifies the Route-2 inequality after removing a length-2
arm:

    mu_B(lambda) >= m - 3/2,

where ``m`` is the leftmost mode of ``S(2^a,r)`` and
``lambda=i_{m-1}(T)/i_m(T)``.

The sign check is exact integer arithmetic.  If ``lambda=L/M`` and
``B(x)=sum b_k x^k``, then the inequality is equivalent to

    sum_k (2k-2m+3) b_k L^k M^(d-k) >= 0,

where ``d=deg B``.
"""

from __future__ import annotations

import argparse
import json
import math
from math import comb
from pathlib import Path

from route2_spider_lane_scan import mode_left, path_polys, spider_poly_from_counts


DEFAULT_R_VALUES = [4, 8, 12, 16, 20, 24, 32, 40, 60, 80]


def b_counts_after_removing_length_2(counts: dict[int, int]) -> dict[int, int]:
    out = dict(counts)
    out[2] -= 1
    if out[2] == 0:
        del out[2]
    return out


def fixed_lane_poly(r: int, a: int, paths: list[list[int]]) -> list[int]:
    """Return I(S(2^a,r)) by the direct spider formula.

    This avoids generic polynomial exponentiation:

        P_r(x)(1+2x)^a + x P_{r-1}(x)(1+x)^a.
    """
    if r == 2:
        # S(2^a,2) is the pure S(2^(a+1)) lane.
        return fixed_lane_poly_no_special(2, a + 1, paths)
    return fixed_lane_poly_no_special(r, a, paths)


def fixed_lane_poly_no_special(r: int, a: int, paths: list[list[int]]) -> list[int]:
    pr = paths[r]
    prm1 = paths[r - 1]
    degree = max(a + len(pr) - 1, 1 + a + len(prm1) - 1)
    out = [0] * (degree + 1)

    for j, coeff in enumerate(pr):
        if coeff == 0:
            continue
        for kk in range(a + 1):
            out[j + kk] += coeff * (2**kk) * comb(a, kk)

    for j, coeff in enumerate(prm1):
        if coeff == 0:
            continue
        for kk in range(a + 1):
            out[1 + j + kk] += coeff * comb(a, kk)

    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def scaled_sums(poly: list[int], numerator: int, denominator: int) -> tuple[int, int]:
    """Return Z and first moment numerator after clearing denominators."""
    degree = len(poly) - 1
    powers_num = [1]
    powers_den = [1]
    for _ in range(degree):
        powers_num.append(powers_num[-1] * numerator)
        powers_den.append(powers_den[-1] * denominator)

    z = 0
    first = 0
    for k, coeff in enumerate(poly):
        if coeff == 0:
            continue
        weight = coeff * powers_num[k] * powers_den[degree - k]
        z += weight
        first += k * weight
    return z, first


def path_eval_and_deriv_scaled(r: int, numerator: int, denominator: int) -> tuple[int, int]:
    """Return M^d P_r(L/M) and M^d (L/M)P'_r(L/M)."""
    eval_prev2 = 1
    deriv_prev2 = 0
    if r == 0:
        return eval_prev2, deriv_prev2

    eval_prev1 = denominator + numerator
    deriv_prev1 = numerator
    if r == 1:
        return eval_prev1, deriv_prev1

    for n in range(2, r + 1):
        if n % 2 == 0:
            path_eval = eval_prev1 + numerator * eval_prev2
            log_deriv_num = deriv_prev1 + numerator * (eval_prev2 + deriv_prev2)
        else:
            path_eval = denominator * eval_prev1 + numerator * eval_prev2
            log_deriv_num = (
                denominator * deriv_prev1
                + numerator * (eval_prev2 + deriv_prev2)
            )
        eval_prev2, eval_prev1 = eval_prev1, path_eval
        deriv_prev2, deriv_prev1 = deriv_prev1, log_deriv_num
    return path_eval, log_deriv_num


def fixed_lane_b_scaled_sums(r: int, a: int, numerator: int, denominator: int) -> tuple[int, int]:
    """Return scaled ``B(lambda)`` and ``lambda B'(lambda)`` for B=S(2^(a-1),r).

    For r != 2, after removing one length-2 arm from S(2^a,r),

        B(x) = P_r(x)(1+2x)^(a-1) + xP_{r-1}(x)(1+x)^(a-1).

    This evaluates that expression directly at lambda=L/M and clears a single
    common denominator, avoiding a degree-by-degree coefficient sum.
    """
    l_val = numerator
    m_val = denominator
    a_minus = a - 1
    p_eval, p_deriv = path_eval_and_deriv_scaled(r, l_val, m_val)
    q_eval, q_deriv = path_eval_and_deriv_scaled(r - 1, l_val, m_val)
    p_degree = (r + 1) // 2
    q_degree = r // 2
    common_degree = max(p_degree + a_minus, q_degree + a)
    arm_base = m_val + 2 * l_val
    tail_base = m_val + l_val
    arm_pow = arm_base**a_minus
    tail_pow = tail_base**a_minus
    arm_scale = m_val ** (common_degree - (p_degree + a_minus))
    tail_scale = m_val ** (common_degree - (q_degree + a))

    z_arm = p_eval * arm_pow * arm_scale
    z_tail = l_val * q_eval * tail_pow * tail_scale

    first_arm_inner = p_deriv * arm_pow
    first_tail_inner = l_val * (q_eval + q_deriv) * tail_pow
    if a_minus:
        first_arm_inner += p_eval * 2 * a_minus * l_val * arm_base ** (a_minus - 1)
        first_tail_inner += (
            a_minus * l_val * l_val * q_eval * tail_base ** (a_minus - 1)
        )

    first_arm = first_arm_inner * arm_scale
    first_tail = first_tail_inner * tail_scale
    return z_arm + z_tail, first_arm + first_tail


def ratio_float(numerator: int, denominator: int) -> float:
    if numerator == 0:
        return 0.0
    sign = -1.0 if numerator < 0 else 1.0
    return sign * math.exp(math.log(abs(numerator)) - math.log(denominator))


def fraction_less(num_a: int, den_a: int, num_b: int, den_b: int) -> bool:
    return num_a * den_b < num_b * den_a


def route2_record(r: int, a: int, paths: list[list[int]]) -> dict | None:
    counts = {2: a}
    if r == 2:
        counts[2] += 1
    else:
        counts[r] = 1

    poly_t = fixed_lane_poly(r, a, paths)
    m = mode_left(poly_t)
    if m < 2 or poly_t[m - 1] <= 0 or poly_t[m] <= 0:
        return None

    lam_num = poly_t[m - 1]
    lam_den = poly_t[m]
    if r == 2:
        b_counts = b_counts_after_removing_length_2(counts)
        poly_b = spider_poly_from_counts(b_counts, paths)
        if m - 1 >= len(poly_b):
            return None
        z, first = scaled_sums(poly_b, lam_num, lam_den)
    else:
        z, first = fixed_lane_b_scaled_sums(r, a, lam_num, lam_den)
    route2_num = 2 * first - (2 * m - 3) * z
    route2_den = 2 * z
    stronger_num = (first - (m - 1) * z) * (lam_num + lam_den) + lam_num * z
    stronger_den = (lam_num + lam_den) * z
    slack_float = ratio_float(route2_num, route2_den)
    stronger_float = ratio_float(stronger_num, stronger_den)

    return {
        "r": r,
        "a": a,
        "n": 1 + 2 * a + (0 if r == 2 else r),
        "m": m,
        "lambda_num": lam_num,
        "lambda_den": lam_den,
        "route2_num": route2_num,
        "route2_den": route2_den,
        "stronger_num": stronger_num,
        "stronger_den": stronger_den,
        "route2_slack_float": slack_float,
        "stronger_slack_float": stronger_float,
    }


def compact_record(record: dict) -> dict:
    return {
        "r": record["r"],
        "a": record["a"],
        "n": record["n"],
        "m": record["m"],
        "route2_num_sign": (record["route2_num"] > 0) - (record["route2_num"] < 0),
        "stronger_num_sign": (record["stronger_num"] > 0) - (record["stronger_num"] < 0),
        "route2_slack_float": record["route2_slack_float"],
        "stronger_slack_float": record["stronger_slack_float"],
        "lambda_float": record["lambda_num"] / record["lambda_den"],
    }


def check_lanes(r_values: list[int], a_max: int) -> dict:
    paths = path_polys(max(2, max(r_values)))
    lane_summaries = []
    global_min_route2 = None
    global_min_stronger = None
    total_checked = 0
    skipped = 0

    for r in r_values:
        lane_min_route2 = None
        lane_min_stronger = None
        failures = []
        stronger_failures = []
        lane_checked = 0
        for a in range(1, a_max + 1):
            record = route2_record(r, a, paths)
            if record is None:
                skipped += 1
                continue
            total_checked += 1
            lane_checked += 1
            if record["route2_num"] < 0:
                failures.append(compact_record(record))
            if record["stronger_num"] < 0:
                stronger_failures.append(compact_record(record))
            if lane_min_route2 is None or fraction_less(
                record["route2_num"],
                record["route2_den"],
                lane_min_route2["_route2_num"],
                lane_min_route2["_route2_den"],
            ):
                lane_min_route2 = compact_record(record) | {
                    "_route2_num": record["route2_num"],
                    "_route2_den": record["route2_den"],
                }
            if lane_min_stronger is None or fraction_less(
                record["stronger_num"],
                record["stronger_den"],
                lane_min_stronger["_stronger_num"],
                lane_min_stronger["_stronger_den"],
            ):
                lane_min_stronger = compact_record(record) | {
                    "_stronger_num": record["stronger_num"],
                    "_stronger_den": record["stronger_den"],
                }

        if lane_min_route2 is not None and (
            global_min_route2 is None
            or fraction_less(
                lane_min_route2["_route2_num"],
                lane_min_route2["_route2_den"],
                global_min_route2["_route2_num"],
                global_min_route2["_route2_den"],
            )
        ):
            global_min_route2 = lane_min_route2
        if lane_min_stronger is not None and (
            global_min_stronger is None
            or fraction_less(
                lane_min_stronger["_stronger_num"],
                lane_min_stronger["_stronger_den"],
                global_min_stronger["_stronger_num"],
                global_min_stronger["_stronger_den"],
            )
        ):
            global_min_stronger = lane_min_stronger

        lane_summaries.append(
            {
                "r": r,
                "checked": lane_checked,
                "failures": failures,
                "stronger_failures": stronger_failures,
                "min_route2": strip_internal_min(lane_min_route2),
                "min_stronger": strip_internal_min(lane_min_stronger),
            }
        )

    return {
        "params": {"r_values": r_values, "a_max": a_max},
        "total_checked": total_checked,
        "skipped": skipped,
        "total_failures": sum(len(row["failures"]) for row in lane_summaries),
        "total_stronger_failures": sum(len(row["stronger_failures"]) for row in lane_summaries),
        "global_min_route2": strip_internal_min(global_min_route2),
        "global_min_stronger": strip_internal_min(global_min_stronger),
        "lanes": lane_summaries,
    }


def strip_internal_min(record: dict | None) -> dict | None:
    if record is None:
        return None
    return {k: v for k, v in record.items() if not k.startswith("_")}


def parse_r_values(raw: str) -> list[int]:
    if raw.strip().lower() == "selected":
        return DEFAULT_R_VALUES[:]
    return [int(part) for part in raw.split(",") if part.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description="Exact finite Route-2 fixed-r lane checks.")
    ap.add_argument("--r-values", default="selected")
    ap.add_argument("--a-max", type=int, default=199)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    r_values = parse_r_values(args.r_values)
    result = check_lanes(r_values, args.a_max)

    print(f"checked records: {result['total_checked']}")
    print(f"skipped records: {result['skipped']}")
    print(f"Route-2 failures: {result['total_failures']}")
    print(f"stronger-threshold failures: {result['total_stronger_failures']}")
    print(f"global min Route-2: {result['global_min_route2']}")
    print(f"global min stronger: {result['global_min_stronger']}")
    for lane in result["lanes"]:
        print(
            f"r={lane['r']:>3} "
            f"min_route2={lane['min_route2']['route2_slack_float']:.12g} "
            f"at a={lane['min_route2']['a']} "
            f"min_stronger={lane['min_stronger']['stronger_slack_float']:.12g} "
            f"at a={lane['min_stronger']['a']} "
            f"failures={len(lane['failures'])} "
            f"stronger_failures={len(lane['stronger_failures'])}"
        )

    if args.out:
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2))
        print(f"wrote {path}")

    assert result["total_failures"] == 0
    assert result["total_stronger_failures"] == 0


if __name__ == "__main__":
    main()
