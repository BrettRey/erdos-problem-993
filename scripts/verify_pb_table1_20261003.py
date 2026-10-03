#!/usr/bin/env python3
"""Numbers audit for Table 1 of the PB paper (balanced variance-one family).

For m >= 2, m success probabilities eps_m and m equal to 1 - eps_m, with
2 m eps_m (1 - eps_m) = 1, so V = 1 and D = m + 1. Computes, at D, in 40-digit
precision: the true deficit, the ULC bound (1.7), the Dumbgen-Wellner bound
1/(D+1), Johnson's odds-sum bound, and the bound from Pitman's (20) combined at
k-1 and k; prints each scaled by its claimed rate (2/m, 1/m, m^-2, sqrt2/m) and
the Skellam(1/2,1/2) limit of delta_D. W - m is the difference of two
independent Bin(m, eps_m) variables, which gives the pmf exactly.
"""
from mpmath import mp, mpf, sqrt, binomial, sinh, cosh, findroot, exp, besseli

mp.dps = 40


def row(m: int) -> dict:
    eps = (1 - sqrt(1 - mpf(2) / m)) / 2
    pa = [binomial(m, a) * eps**a * (1 - eps) ** (m - a) for a in range(m + 1)]

    def g(j: int):  # P(W - m = j)
        return sum(pa[a] * pa[a - j] for a in range(max(0, j), m + 1) if 0 <= a - j <= m)

    f0, f1, f2 = g(0), g(1), g(2)  # indices D-1 = m, D, D+1
    true_delta = 1 - f0 * f2 / f1**2
    ulc = mpf(2 * m + 1) / (m * (m + 2))
    dw = mpf(1) / (m + 2)
    odds = m * eps / (1 - eps) + m * (1 - eps) / eps
    johnson = (f2 / f1) / odds
    mu = lambda t: sinh(t) / (1 + (cosh(t) - 1) / m)  # tilted mean minus m
    lo = findroot(lambda t: mu(t) - (1 - mpf(1) / (m + 1)), mpf(0.88))
    hi = findroot(lambda t: mu(t) - (1 + mpf(1) / (m + 3)), mpf(0.88))
    pitman = 1 - exp(lo - hi)
    return dict(m=m, true_delta=true_delta, ulc_x_m_over2=ulc * m / 2, dw_x_m=dw * m,
                johnson_x_m2=johnson * m * m, pitman_x_m_over_sqrt2=pitman * m / sqrt(2))


if __name__ == "__main__":
    for m in (10, 100, 1000):
        r = row(m)
        print({k: (round(float(v), 5) if k != "m" else v) for k, v in r.items()})
    I0, I1, I2 = besseli(0, 1), besseli(1, 1), besseli(2, 1)
    print("Skellam limit of delta_D:", float(1 - I0 * I2 / I1**2))
    r = row(1000)
    assert abs(r["ulc_x_m_over2"] - 1) < 0.01 and abs(r["dw_x_m"] - 1) < 0.01
    assert abs(r["pitman_x_m_over_sqrt2"] - 1) < 0.01 and r["johnson_x_m2"] < 1
    assert abs(r["true_delta"] - (1 - I0 * I2 / I1**2)) < 0.01
    print("TABLE 1 CHECKS PASSED")
