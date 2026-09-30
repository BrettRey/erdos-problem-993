#!/usr/bin/env python3
"""Certified (Arb interval arithmetic) sharpening of Lemma 4.2 of arXiv:2609.20961.

Claim certified for each listed rational p and rational rho:
    for all lam in [1/4, 12] and all y >= 0, with r = lam e^{-y}, q = r/(1+r):
        y q^p / log(1+r) <= rho  (< 1).
Reduction (proved in REPORT.md): for fixed r the prefactor y = log(lam/r) is increasing
in lam, and r ranges over (0, lam], so the supremum over K x [0,inf) equals
sup_{y>=0} F(y) with lam = 12, F(y) = y q^p / log(1+r), r = 12 e^{-y}.
Tail: for y >= Y0, q <= log(1+r) and q <= r give F(y) <= y r^{p-1} = y 12^{p-1} e^{-(p-1)y},
which is decreasing for y >= 1/(p-1); we require Y0 >= 1/(p-1) and bound it at Y0.
On [0, Y0] we cover by intervals and evaluate F on Arb balls (y = 0 handled: F(0) = 0 and
the ball evaluation at y in [0, h] is fine since log(1+r) >= log(1 + 12 e^{-h}) > 0).
Run with the project venv: venv/bin/python certify_lemma42.py
"""
import json
from fractions import Fraction as Fr
from flint import arb, ctx

ctx.prec = 120


def F_ball(lo, hi, p):
    # lo, hi are dyadic (bisection of [0, Y0] with Y0 an integer), so mid and rad are exact;
    # the radius is still padded by a relative 1e-12 for safety.
    y = arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 1e-12))
    r = arb(12) * (-y).exp()
    q = r / (1 + r)
    return y * q ** p / (1 + r).log()


def certify(p_frac, rho_frac, Y0=None, N=None):
    p = arb(p_frac.numerator) / p_frac.denominator
    pm1 = float(p_frac) - 1
    if Y0 is None:
        Y0 = max(60.0, 1.0 / pm1 + 1)
    assert Y0 >= 1.0 / pm1
    # tail bound at Y0
    tailb = arb(Y0) * arb(12) ** (p - 1) * (-(p - 1) * arb(Y0)).exp()
    rho = arb(rho_frac.numerator) / rho_frac.denominator
    ok_tail = (tailb < rho)
    # adaptive cover of [0, Y0]
    stack = [(0.0, float(Y0))]
    worst = arb(0)
    nboxes = 0
    while stack:
        lo, hi = stack.pop()
        v = F_ball(lo, hi, p)
        nboxes += 1
        if v < rho:
            if (v.mid() + v.rad()) > (worst.mid() + worst.rad()):
                worst = v
            continue
        if hi - lo < 1e-9:
            return dict(ok=False, reason=f"cannot separate at [{lo},{hi}] value {v}")
        mid = (lo + hi) / 2
        stack.append((lo, mid)); stack.append((mid, hi))
    return dict(ok=bool(ok_tail), p=str(p_frac), rho_upper=str(rho_frac), Y0=Y0,
                tail_bound=str(tailb), boxes=nboxes, max_box_upper=str(worst.mid() + worst.rad()))


if __name__ == "__main__":
    targets = [(Fr(193, 100), Fr(98, 100)), (Fr(19, 10), Fr(995, 1000)), (Fr(195, 100), Fr(97, 100)),
               (Fr(2, 1), Fr(92, 100))]
    res = []
    for p, rho in targets:
        r = certify(p, rho)
        print(r)
        res.append(r)
    # also: the paper's own value b of eq. (4.5) vs the true sup at p = 2 (for the report)
    cert = [r for r in res if r.get("ok")]
    json.dump(dict(prec_bits=ctx.prec, certified=cert, all=res), open("data/lemma42_certificates.json", "w"), indent=1)
