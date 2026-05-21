"""Proof certificate for Route-2 on the pure spider lane S(2^a).

This is intentionally conservative.  It verifies the finite range exactly and
checks the elementary envelope used for all a >= 30:

    |Delta_a - Delta_a^0|
      <= 4 a^2 / 2^((2a-5)/3) + 4 a (3/4)^(a-1)
      < 1/6 <= Delta_a^0.

The exact finite check covers 2 <= a <= 29.
"""

from __future__ import annotations

from fractions import Fraction

from pure_spider_route2_asymptotic import exact_record


THRESHOLD = 30


def simple_error_bound(a: int) -> float:
    return 4 * a * a / (2 ** ((2 * a - 5) / 3)) + 4 * a * ((3 / 4) ** (a - 1))


def main() -> None:
    finite = [exact_record(a) for a in range(2, THRESHOLD)]
    min_finite = min(finite, key=lambda rec: rec["slack"])
    finite_failures = [rec for rec in finite if rec["slack"] <= 0]

    b0 = simple_error_bound(THRESHOLD)
    rational_b0_upper = Fraction(4 * THRESHOLD * THRESHOLD, 2**18) + Fraction(
        4 * THRESHOLD * 3**29, 4**29
    )
    term1_ratio = ((THRESHOLD + 1) / THRESHOLD) ** 2 * 2 ** (-Fraction(2, 3))
    term2_ratio = ((THRESHOLD + 1) / THRESHOLD) * Fraction(3, 4)

    print(f"finite exact range: a=2..{THRESHOLD - 1}")
    print(f"finite failures: {len(finite_failures)}")
    print(
        "minimum finite slack: "
        f"a={min_finite['a']} slack={float(min_finite['slack']):.12g}"
    )
    print()
    print(f"asymptotic threshold: a >= {THRESHOLD}")
    print(f"simple error bound at threshold: {b0:.12g}")
    print(
        "rational upper at threshold: "
        f"{rational_b0_upper} ~= {float(rational_b0_upper):.12g}"
    )
    print(f"reserve 1/6 - bound: {1 / 6 - b0:.12g}")
    print(
        "rational reserve: "
        f"{Fraction(1, 6) - rational_b0_upper} "
        f"~= {float(Fraction(1, 6) - rational_b0_upper):.12g}"
    )
    print(f"term1 ratio bound at threshold: {float(term1_ratio):.12g}")
    print(f"term2 ratio bound at threshold: {float(term2_ratio):.12g}")

    assert not finite_failures
    assert b0 < 1 / 6
    assert rational_b0_upper < Fraction(1, 6)
    assert term1_ratio < 1
    assert term2_ratio < 1


if __name__ == "__main__":
    main()
