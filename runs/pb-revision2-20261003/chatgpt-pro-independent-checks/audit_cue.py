#!/usr/bin/env python3
"""Exact numerical checks for the semicircle-count example.

Uses only the Python standard library.  Run with ``python audit_cue.py``.
Writes audit_cue.json and audit_cue.log beside this script.

Scope: the CUE variance formula and its strict monotonicity were checked
algebraically separately.  They are inputs to this numerical verification,
not conclusions of the program.  The program uses rational alternating-series
bounds and Machin's identity to enclose pi, then certifies V_1996 < 1 < V_1998
and the displayed numbers for N = 10000.  Together with the separately proved
monotonicity, the first pair establishes the claimed first even N = 1998.
No floating-point arithmetic or third-party packages are used.
"""

from fractions import Fraction
import json
from pathlib import Path


def arctan_interval(q, terms):
    """Enclose atan(1/q) by two consecutive alternating-series partial sums."""
    assert q > 1 and terms > 0
    partial = sum(
        (Fraction((-1) ** k, (2 * k + 1) * q ** (2 * k + 1))
         for k in range(terms)),
        Fraction(0),
    )
    next_partial = partial + Fraction(
        (-1) ** terms, (2 * terms + 1) * q ** (2 * terms + 1)
    )
    return min(partial, next_partial), max(partial, next_partial)


def fixed_decimal(integer, digits):
    """Render integer / 10**digits without rounding."""
    sign = "-" if integer < 0 else ""
    integer = abs(integer)
    scale = 10 ** digits
    whole, fractional = divmod(integer, scale)
    return f"{sign}{whole}.{fractional:0{digits}d}"


def decimal_enclosure(lower, upper, digits=60):
    """Return outward-rounded decimal endpoints, checked by exact comparison."""
    assert lower <= upper
    scale = 10 ** digits
    floor_lower = (lower.numerator * scale) // lower.denominator
    ceil_upper = -((-upper.numerator * scale) // upper.denominator)
    assert Fraction(floor_lower, scale) <= lower
    assert Fraction(ceil_upper, scale) >= upper
    return {
        "lower": fixed_decimal(floor_lower, digits),
        "upper": fixed_decimal(ceil_upper, digits),
        "decimal_places": digits,
    }


def variance_interval(n, pi_squared_lower, pi_squared_upper):
    """Enclose N/4 - 2/pi**2 * sum_{odd d < N} (N-d)/d**2."""
    assert n > 0 and n % 2 == 0
    total = sum(
        (Fraction(n - d, d * d) for d in range(1, n, 2)),
        Fraction(0),
    )
    lower = Fraction(n, 4) - 2 * total / pi_squared_lower
    upper = Fraction(n, 4) - 2 * total / pi_squared_upper
    assert lower < upper
    return lower, upper


def assert_strict_decimal_enclosure(interval, displayed_lower, displayed_upper):
    """Certify that the exact interval lies inside the given decimal interval."""
    lower, upper = interval
    assert Fraction(displayed_lower) < lower
    assert upper < Fraction(displayed_upper)


def main():
    # Machin's identity: pi = 16 atan(1/5) - 4 atan(1/239).
    # The alternating series terms strictly decrease, so consecutive partial
    # sums strictly enclose each arctangent.  Interval subtraction reverses
    # the second interval's endpoints.
    a_lower, a_upper = arctan_interval(5, 50)
    b_lower, b_upper = arctan_interval(239, 15)
    pi_lower = 16 * a_lower - 4 * b_upper
    pi_upper = 16 * a_upper - 4 * b_lower
    assert Fraction(3) < pi_lower < pi_upper < Fraction(22, 7)
    assert pi_upper - pi_lower < Fraction("4.02e-72")
    assert_strict_decimal_enclosure(
        (pi_lower, pi_upper),
        "3.14159265358979323846264338327950288419716939937510",
        "3.14159265358979323846264338327950288419716939937511",
    )
    pi_squared_lower = pi_lower * pi_lower
    pi_squared_upper = pi_upper * pi_upper

    intervals = {
        n: variance_interval(n, pi_squared_lower, pi_squared_upper)
        for n in (1996, 1998, 10000)
    }

    # These comparisons use the exact rational endpoints, not rounded output.
    assert intervals[1996][1] < 1
    assert intervals[1998][0] > 1
    displayed_variance_intervals = {
        1996: (
            "0.99996543523160400122908765",
            "0.99996543523160400122908766",
        ),
        1998: (
            "1.00006690864230001837548359",
            "1.00006690864230001837548360",
        ),
        10000: (
            "1.16323843886831071519829808",
            "1.16323843886831071519829809",
        ),
    }
    for n, endpoints in displayed_variance_intervals.items():
        assert_strict_decimal_enclosure(intervals[n], *endpoints)

    # Reciprocal is decreasing on the positive interval.
    variance_lower, variance_upper = intervals[10000]
    curvature_interval = (
        Fraction(1, 4) / variance_upper,
        Fraction(1, 4) / variance_lower,
    )
    assert_strict_decimal_enclosure(
        curvature_interval,
        "0.21491724451886196499928543",
        "0.21491724451886196499928544",
    )
    n = 10000
    newton = Fraction(n + 1, (n // 2 + 2) * (n // 2))
    assert Fraction("0.00039988004798080767692922") < newton
    assert newton < Fraction("0.00039988004798080767692923")

    results = {
        "status": "all exact-rational assertions passed",
        "scope": (
            "Numerical consequences of the separately derived CUE variance "
            "formula. Operator theory, the formula's derivation, monotonicity, "
            "and the asymptotic expansion are not checked by this program."
        ),
        "arithmetic": "Python fractions.Fraction; no floating-point arithmetic",
        "pi_method": {
            "identity": "pi = 16 atan(1/5) - 4 atan(1/239)",
            "terms": {"atan_1_over_5": 50, "atan_1_over_239": 15},
            "certified_interval_width_less_than": "4.02e-72",
            "enclosure": decimal_enclosure(pi_lower, pi_upper, 70),
        },
        "variance_formula": (
            "V_N = N/4 - (2/pi^2) * sum((N-d)/d^2 for odd 1 <= d < N)"
        ),
        "variances": {
            str(n): decimal_enclosure(*interval)
            for n, interval in intervals.items()
        },
        "threshold_checks": {
            "V_1996_less_than_1": True,
            "V_1998_greater_than_1": True,
            "first_even_N_given_separately_proved_strict_monotonicity": 1998,
        },
        "N_10000_curvature_lower_bound_1_over_4V": decimal_enclosure(
            *curvature_interval
        ),
        "N_10000_Newton_bound": {
            "exact_numerator": newton.numerator,
            "exact_denominator": newton.denominator,
            "enclosure": decimal_enclosure(newton, newton),
        },
    }

    lines = [
        "PASS: all exact-rational assertions passed.",
        "Scope: formula and monotonicity were checked algebraically separately.",
        "pi interval width < 4.02e-72 (rational Machin/alternating-series bounds).",
    ]
    for n in intervals:
        enclosure = decimal_enclosure(*intervals[n], digits=30)
        lines.append(f"V_{n} in [{enclosure['lower']}, {enclosure['upper']}]")
    lines.append("Certified V_1996 < 1 < V_1998.")
    lines.append("Strict monotonicity therefore gives first even N = 1998.")
    for label, interval in (
        ("N=10000: 1/(4V)", curvature_interval),
        ("N=10000: Newton bound", (newton, newton)),
    ):
        enclosure = decimal_enclosure(*interval, digits=30)
        lines.append(f"{label} in [{enclosure['lower']}, {enclosure['upper']}]")
    log = "\n".join(lines) + "\n"
    script = Path(__file__).resolve()
    script.with_suffix(".json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8"
    )
    script.with_suffix(".log").write_text(log, encoding="utf-8")
    print(log, end="")


if __name__ == "__main__":
    main()
