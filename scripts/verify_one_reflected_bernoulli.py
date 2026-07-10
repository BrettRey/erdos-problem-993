#!/usr/bin/env python3
"""Exact verification harness for the one-reflected-Bernoulli lemma.

Let X be a Poisson-binomial sum whose parameters lie in (0, 1/2], let
Y ~ Bernoulli(q), q in [0, 1/2], and put Z = X - Y.  The proof shows that,
when Var(X) >= 1, the first strict descent D of Z is one of S-1 and S,
where S is the first strict descent of X, and

    (Var(X) + q(1-q)) * Delta_D(Z) >= 1/4.

This script is a deterministic exact-arithmetic verification harness, not the
proof.  It checks the descent-selection threshold, the curvature identity and
its three Newton lower bounds, the final proof lower bound, and the quarter
inequality on exhaustive finite rational grids.
"""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence


CORE_PROBABILITIES = (
    Fraction(1, 10),
    Fraction(1, 6),
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(2, 5),
    Fraction(1, 2),
)
HEAVY_PROBABILITIES = (
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(2, 5),
    Fraction(1, 2),
)
DUST_PROBABILITIES = (
    Fraction(1, 100),
    Fraction(1, 50),
    Fraction(1, 20),
    Fraction(1, 10),
)
Q_VALUES = (
    Fraction(0),
    Fraction(1, 1000),
    Fraction(1, 100),
    Fraction(1, 20),
    Fraction(1, 10),
    Fraction(20, 119),
    Fraction(1, 4),
    Fraction(1, 3),
    Fraction(2, 5),
    Fraction(1, 2),
)


def convolve_bernoulli(pmf: Sequence[Fraction], p: Fraction) -> list[Fraction]:
    out = [Fraction(0) for _ in range(len(pmf) + 1)]
    for k, value in enumerate(pmf):
        out[k] += value * (1 - p)
        out[k + 1] += value * p
    return out


def poisson_binomial_pmf(parameters: Sequence[Fraction]) -> list[Fraction]:
    pmf = [Fraction(1)]
    for p in parameters:
        if not 0 < p <= Fraction(1, 2):
            raise AssertionError(f"parameter outside (0, 1/2]: {p}")
        pmf = convolve_bernoulli(pmf, p)
    return pmf


def first_strict_descent(pmf: Sequence[Fraction], support_start: int = 0) -> int:
    for index in range(1, len(pmf)):
        if pmf[index] < pmf[index - 1]:
            return support_start + index
    raise AssertionError("finite positive pmf has no strict descent")


def signed_pmf(a: Sequence[Fraction], q: Fraction) -> list[Fraction]:
    """Return the pmf of X-Y on support -1, 0, ..., len(a)-1."""

    if not 0 <= q <= Fraction(1, 2):
        raise AssertionError(f"q outside [0, 1/2]: {q}")
    out = [Fraction(0) for _ in range(len(a) + 1)]
    for k, value in enumerate(a):
        out[k + 1] += (1 - q) * value
        out[k] += q * value
    return out


def at(seq: Sequence[Fraction], index: int) -> Fraction:
    return seq[index] if 0 <= index < len(seq) else Fraction(0)


def a_at(a: Sequence[Fraction], k: int) -> Fraction:
    return at(a, k)


def c_at(c: Sequence[Fraction], k: int) -> Fraction:
    return at(c, k + 1)


def effective_drop(c: Sequence[Fraction], descent: int) -> Fraction:
    center = c_at(c, descent)
    if center <= 0:
        raise AssertionError("effective drop has zero center coefficient")
    return 1 - c_at(c, descent - 1) * c_at(c, descent + 1) / center**2


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "exact": fraction_text(value),
        "decimal": f"{float(value):.17g}",
    }


def parameter_blocks(parameters: Sequence[Fraction]) -> list[dict[str, Any]]:
    counts = Counter(parameters)
    return [
        {"count": counts[p], "p": fraction_text(p)}
        for p in sorted(counts)
    ]


def selection_threshold(a: Sequence[Fraction], s: int) -> Fraction:
    if s < 2:
        raise AssertionError(f"selection threshold requires S >= 2, got {s}")
    u = a[s - 1] / a[s - 2]
    v = a[s] / a[s - 1]
    if u < 1 or v >= 1:
        raise AssertionError(f"invalid threshold ratios u={u}, v={v}, S={s}")
    return (u - 1) / (u * (1 - v))


def verify_curvature(
    a: Sequence[Fraction],
    c: Sequence[Fraction],
    q: Fraction,
    k: int,
) -> tuple[Fraction, Fraction]:
    """Check the exact decomposition and lower bounds at index k.

    Returns (actual Delta_k(c), the combined explicit curvature lower bound).
    """

    akm1 = a_at(a, k - 1)
    ak = a_at(a, k)
    akp1 = a_at(a, k + 1)
    akp2 = a_at(a, k + 2)
    first = ak**2 - akm1 * akp1
    second = akp1**2 - ak * akp2
    cross = ak * akp1 - akm1 * akp2

    lhs = c_at(c, k) ** 2 - c_at(c, k - 1) * c_at(c, k + 1)
    rhs = (1 - q) ** 2 * first + q**2 * second + q * (1 - q) * cross
    if lhs != rhs:
        raise AssertionError(
            f"curvature decomposition failed at k={k}, q={q}: {lhs} != {rhs}"
        )

    if first < ak**2 / (k + 1):
        raise AssertionError(f"first Newton bound failed at k={k}")
    if second < akp1**2 / (k + 2):
        raise AssertionError(f"second Newton bound failed at k={k}")
    if cross < 2 * ak * akp1 / (k + 2):
        raise AssertionError(f"cross Newton bound failed at k={k}")

    center = c_at(c, k)
    if center <= 0 or ak <= 0:
        raise AssertionError(f"nonpositive center at curvature index k={k}")
    delta = lhs / center**2
    t = q / (1 - q)
    r = akp1 / ak
    lower = Fraction(1, k + 2) + Fraction(1, (k + 1) * (k + 2)) / (
        1 + t * r
    ) ** 2
    if delta < lower:
        raise AssertionError(
            f"combined curvature lower bound failed at k={k}, q={q}: "
            f"Delta={delta}, lower={lower}"
        )
    return delta, lower


def verify_row(
    parameters: Sequence[Fraction],
    a: Sequence[Fraction],
    variance_x: Fraction,
    s: int,
    q: Fraction,
) -> dict[str, Any]:
    c = signed_pmf(a, q)
    d = first_strict_descent(c, support_start=-1)
    threshold = selection_threshold(a, s)
    t = q / (1 - q)
    expected = s - 1 if t > threshold else s
    if d != expected or d not in (s - 1, s):
        raise AssertionError(
            f"selection failed: S={s}, D={d}, expected={expected}, "
            f"q={q}, threshold={threshold}, params={parameters}"
        )

    delta, curvature_lower = verify_curvature(a, c, q, d)
    direct_delta = effective_drop(c, d)
    if delta != direct_delta:
        raise AssertionError(f"curvature/direct Delta mismatch: {delta} != {direct_delta}")

    if d == s - 1:
        proof_lower = Fraction(1, s + 1)
        if delta < proof_lower:
            raise AssertionError(
                f"shifted-PB Newton bound failed: Delta={delta}, lower={proof_lower}"
            )
    else:
        proof_lower = curvature_lower

    variance_total = variance_x + q * (1 - q)
    proof_scaled = variance_total * proof_lower
    scaled = variance_total * delta
    if proof_scaled < Fraction(1, 4):
        raise AssertionError(
            f"proof-chain quarter bound failed: {proof_scaled}, S={s}, D={d}, q={q}"
        )
    if scaled < Fraction(1, 4):
        raise AssertionError(
            f"quarter bound failed: {scaled}, S={s}, D={d}, q={q}"
        )

    return {
        "descent": d,
        "selection": "S-1" if d == s - 1 else "S",
        "threshold_equal": t == threshold,
        "delta": delta,
        "curvature_lower": curvature_lower,
        "proof_lower": proof_lower,
        "variance_total": variance_total,
        "proof_scaled": proof_scaled,
        "scaled": scaled,
    }


def core_profiles(min_m: int, max_m: int) -> Iterable[tuple[Fraction, ...]]:
    for m in range(min_m, max_m + 1):
        yield from itertools.combinations_with_replacement(CORE_PROBABILITIES, m)


def grouped_profiles() -> Iterable[tuple[Fraction, ...]]:
    for heavy_p, dust_p, heavy_count, dust_count in itertools.product(
        HEAVY_PROBABILITIES,
        DUST_PROBABILITIES,
        range(2, 9),
        (1, 2, 4, 8, 16),
    ):
        yield tuple(sorted((heavy_p,) * heavy_count + (dust_p,) * dust_count))


def anchored_checks() -> dict[str, Any]:
    p = Fraction(11, 31)
    parameters = (p,) * 5
    a = poisson_binomial_pmf(parameters)
    s = first_strict_descent(a)
    if s != 3:
        raise AssertionError(f"anchor S mismatch: {s}")
    u = a[s - 1] / a[s - 2]
    v = a[s] / a[s - 1]
    threshold = selection_threshold(a, s)
    if (u, v, threshold) != (
        Fraction(11, 10),
        Fraction(11, 20),
        Fraction(20, 99),
    ):
        raise AssertionError(f"anchor threshold mismatch: {(u, v, threshold)}")

    rows: dict[str, Any] = {}
    variance_x = sum((x * (1 - x) for x in parameters), Fraction(0))
    for label, q, expected_d in (
        ("below", Fraction(1, 10), 3),
        ("equal", Fraction(20, 119), 3),
        ("above", Fraction(1, 2), 2),
    ):
        row = verify_row(parameters, a, variance_x, s, q)
        if row["descent"] != expected_d:
            raise AssertionError(
                f"anchor {label} descent mismatch: {row['descent']} != {expected_d}"
            )
        c = signed_pmf(a, q)
        if label == "equal" and c_at(c, s - 1) != c_at(c, s - 2):
            raise AssertionError("threshold equality does not create the expected plateau")
        rows[label] = {
            "q": fraction_text(q),
            "t": fraction_text(q / (1 - q)),
            "descent": row["descent"],
            "scaled": fraction_record(row["scaled"]),
        }

    fair_parameters = (Fraction(1, 2),) * 5
    fair_a = poisson_binomial_pmf(fair_parameters)
    fair_s = first_strict_descent(fair_a)
    if fair_s != 4:
        raise AssertionError(f"fair-binomial S mismatch: {fair_s}")
    fair_variance = Fraction(5, 4)
    fair_zero = verify_row(
        fair_parameters, fair_a, fair_variance, fair_s, Fraction(0)
    )
    fair_perturbed = verify_row(
        fair_parameters, fair_a, fair_variance, fair_s, Fraction(1, 1000)
    )
    if fair_zero["descent"] != 4 or fair_perturbed["descent"] != 3:
        raise AssertionError("fair-binomial strict-descent discontinuity mismatch")

    return {
        "binomial_5_p_11_over_31": {
            "S": s,
            "u": fraction_text(u),
            "v": fraction_text(v),
            "t_threshold": fraction_text(threshold),
            "q_threshold": fraction_text(Fraction(20, 119)),
            "rows": rows,
        },
        "fair_binomial_plateau": {
            "S": fair_s,
            "q_0_descent": fair_zero["descent"],
            "q_1_over_1000_descent": fair_perturbed["descent"],
            "q_0_delta": fraction_text(fair_zero["delta"]),
            "q_1_over_1000_delta": fraction_text(fair_perturbed["delta"]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-m", type=int, default=4)
    parser.add_argument("--max-m", type=int, default=9)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/one_reflected_bernoulli_verify_2026-07-10.json"
        ),
    )
    args = parser.parse_args()
    if args.min_m < 1 or args.max_m < args.min_m:
        parser.error("require 1 <= min-m <= max-m")

    anchors = anchored_checks()
    seen: set[tuple[Fraction, ...]] = set()
    generated = Counter()
    eligible = Counter()
    selections = Counter()
    rows = 0
    threshold_equal_rows = 0
    minimum: tuple[Fraction, dict[str, Any]] | None = None

    for family, profiles in (
        ("core_multisets", core_profiles(args.min_m, args.max_m)),
        ("heavy_dust_blocks", grouped_profiles()),
    ):
        for parameters in profiles:
            generated[family] += 1
            if parameters in seen:
                continue
            seen.add(parameters)
            a = poisson_binomial_pmf(parameters)
            variance_x = sum(
                (p * (1 - p) for p in parameters),
                Fraction(0),
            )
            if variance_x < 1:
                continue
            eligible[family] += 1
            s = first_strict_descent(a)
            if s + 1 > 4 * variance_x:
                raise AssertionError(
                    f"one-sided localization failed: S={s}, V_X={variance_x}"
                )
            for q in Q_VALUES:
                row = verify_row(parameters, a, variance_x, s, q)
                rows += 1
                selections[row["selection"]] += 1
                threshold_equal_rows += int(row["threshold_equal"])
                if minimum is None or row["scaled"] < minimum[0]:
                    minimum = (
                        row["scaled"],
                        {
                            "family_first_seen": family,
                            "parameters": parameter_blocks(parameters),
                            "variance_x": fraction_record(variance_x),
                            "q": fraction_text(q),
                            "S": s,
                            "D": row["descent"],
                            "selection": row["selection"],
                            "delta": fraction_record(row["delta"]),
                            "proof_lower": fraction_record(row["proof_lower"]),
                            "variance_total": fraction_record(row["variance_total"]),
                            "scaled": fraction_record(row["scaled"]),
                            "proof_scaled": fraction_record(row["proof_scaled"]),
                        },
                    )

    if minimum is None:
        raise AssertionError("grid produced no eligible profiles")

    output = {
        "status": "passed",
        "arithmetic": "fractions.Fraction (exact rational)",
        "conditions": {
            "x_parameters": "0 < p_i <= 1/2",
            "variance_x": ">= 1",
            "q": "0 <= q <= 1/2",
            "q_values": [fraction_text(q) for q in Q_VALUES],
        },
        "grid": {
            "core_probability_values": [
                fraction_text(p) for p in CORE_PROBABILITIES
            ],
            "core_m_range": [args.min_m, args.max_m],
            "heavy_probability_values": [
                fraction_text(p) for p in HEAVY_PROBABILITIES
            ],
            "dust_probability_values": [
                fraction_text(p) for p in DUST_PROBABILITIES
            ],
            "heavy_count_range": [2, 8],
            "dust_counts": [1, 2, 4, 8, 16],
        },
        "counts": {
            "profiles_generated": dict(generated),
            "profiles_eligible_and_unique": dict(eligible),
            "unique_profiles_seen": len(seen),
            "rows": rows,
            "selection_cases": dict(selections),
            "threshold_equal_rows": threshold_equal_rows,
            "curvature_decomposition_checks": rows,
            "individual_newton_bound_checks": 3 * rows,
            "proof_chain_quarter_checks": rows,
            "actual_quarter_checks": rows,
            "failures": 0,
        },
        "minimum_scaled_effective_drop": minimum[1],
        "anchored_checks": anchors,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

    print(f"passed {rows:,} exact rows across {sum(eligible.values()):,} profiles")
    print(f"selection cases: {dict(selections)}")
    print(
        "minimum (Var(X)+q(1-q))*Delta = "
        f"{fraction_text(minimum[0])} = {float(minimum[0]):.17g}"
    )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
