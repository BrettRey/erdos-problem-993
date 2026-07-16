"""Exact certificate for the central-window T3 sequence.

No floating-point arithmetic or numerical root finding is used.
"""

from math import comb

from sympy import I, Poly, Rational, symbols


x = symbols("x")
Q_COEFFS = (1, 26, 15, 20, 15, 6, 1)
Q = Poly(sum(a * x**k for k, a in enumerate(Q_COEFFS)), x)

# Exact Gaussian-rational root counts for Q.
boxes = (
    (-3, -2),
    (-Rational(1, 25), -Rational(2, 51)),
    (-2 - 2 * I, -Rational(3, 2) - Rational(3, 2) * I),
    (-2 + Rational(3, 2) * I, -Rational(3, 2) + 2 * I),
    (Rational(1, 5) - Rational(5, 4) * I, Rational(1, 4) - I),
    (Rational(1, 5) + I, Rational(1, 4) + Rational(5, 4) * I),
)
root_counts = [Q.count_roots(a, b) for a, b in boxes]
assert root_counts == [1] * 6


def coefficients(s: int) -> list[int]:
    """Coefficients of Q(x)(1+x^2)^(2s), in increasing degree."""
    even_factor = [0] * (4 * s + 1)
    for t in range(2 * s + 1):
        even_factor[2 * t] = comb(2 * s, t)
    result = [0] * (4 * s + 7)
    for j, qj in enumerate(Q_COEFFS):
        for k, bk in enumerate(even_factor):
            result[j + k] += qj * bk
    return result


for s in range(3, 101):
    p = coefficients(s)
    alpha = 4 * s + 6
    window_end = (2 * alpha + 1) // 3
    assert len(p) == alpha + 1
    assert all(a > 0 for a in p)
    assert p[2 * s - 1] > p[2 * s] < p[2 * s + 1]
    assert 2 * s + 1 <= window_end

smallest = coefficients(3)
assert smallest == [
    1, 26, 21, 176, 120, 516, 336, 856, 546, 880,
    546, 576, 336, 236, 120, 56, 21, 6, 1,
]

print("Q root counts:", root_counts)
print("s=3 coefficients:", smallest)
print("exact central valleys verified for 3 <= s <= 100")
