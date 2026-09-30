#!/usr/bin/env python3
"""Independent kernel-constant check (lclt-check, resumed run).

(1) CERTIFIED (python-flint Arb): I3 := int |phi'''| = 2(phi(0) + 4 phi(sqrt3)), and
    c_* := phi(0)/I3 = 1/(2 + 8 e^{-3/2}); closed form plus an independent acb.integral
    of phi''' on the sign-definite pieces [0, sqrt3] and [sqrt3, 40] (tail beyond 40
    bounded by |phi''(40)|).
(2) EXACT (Fractions): TV(h_B) for h_B = -Delta^2 b_{B,p}, and the test
    TV(h_B)^2 v^3  >  U^2   with U = 151002/100000 > I3 (rational upper bound for I3),
    which would REFUTE the uniform bound TV(h_B) <= I3 / v^{3/2}.  Also records the exact
    rational TV^2 v^3 for a few (B, p) and the sign pattern of Delta^3 b.
Run with the project venv: venv/bin/python kernel_check_exact.py
"""
import json, sys
from fractions import Fraction as Fr
from math import comb

def tv_exact(B, p):
    q = 1 - p
    b = {j: comb(B, j) * p ** j * q ** (B - j) for j in range(B + 1)}
    g = lambda j: b.get(j, Fr(0))
    h = {j: 2 * g(j) - g(j - 1) - g(j + 1) for j in range(-2, B + 3)}
    H = lambda j: h.get(j, Fr(0))
    d = [H(j + 1) - H(j) for j in range(-3, B + 3)]
    TV = sum(abs(x) for x in d)
    signs = [1 if x > 0 else -1 for x in d if x != 0]
    ch = sum(1 for i in range(len(signs) - 1) if signs[i] != signs[i + 1])
    v = B * p * q
    return TV, v, ch

if __name__ == "__main__":
    out = {}
    try:
        from flint import arb, acb, ctx
        ctx.prec = 200
        two_pi = 2 * arb.pi()
        phi0 = 1 / two_pi.sqrt()
        s3 = arb(3).sqrt()
        phis3 = (arb(-3) / 2).exp() * phi0
        I3 = 2 * (phi0 + 4 * phis3)
        cstar = phi0 / I3
        cstar2 = 1 / (2 + 8 * (arb(-3) / 2).exp())
        # independent: acb.integral of phi'''(x) = (3x - x^3) phi(x)
        f = lambda z, analytic: (3 * z - z ** 3) * (-(z * z) / 2).exp() / acb(two_pi).sqrt()
        P1 = acb.integral(f, 0, acb(s3))                    # phi''' >= 0 here
        P2 = acb.integral(f, acb(s3), 40)                   # phi''' <= 0 here
        tail = (arb(40) ** 2 - 1) * (arb(-800)).exp() * phi0  # |phi''(40)| bounds int_40^inf |phi'''|
        I3_quad = 2 * (P1.real - P2.real) + 2 * tail * arb(0, 1)
        out.update(I3_closed=I3.str(25), I3_quad=I3_quad.str(25), c_star=cstar.str(25),
                   c_star_closed_form_1_over_2_plus_8e_minus_3half=cstar2.str(25),
                   I3_below_151002e_5=bool(I3 < arb(Fr(151002, 100000).numerator) / Fr(151002, 100000).denominator),
                   arb_note="python-flint %s, prec 200 bits" % __import__('flint').__version__)
    except ImportError as e:
        out["arb_error"] = repr(e)
    
    U = Fr(151002, 100000)
    
    rows = []
    cases = [(256, Fr(12, 13)), (128, Fr(12, 13)), (200, Fr(12, 13)), (300, Fr(12, 13)), (64, Fr(12, 13)),
             (256, Fr(1, 2)), (256, Fr(1, 4)), (64, Fr(1, 2)), (16, Fr(1, 2)), (1, Fr(1, 2)), (2, Fr(1, 2)), (4, Fr(1, 2))]
    for B, p in cases:
        TV, v, ch = tv_exact(B, p)
        val = TV * TV * v ** 3                       # exact rational = (TV v^{3/2})^2
        rows.append(dict(B=B, p=str(p), v=str(v), TVv32_float_diag=float(val) ** 0.5,
                         exceeds_U=bool(val > U * U), delta3_sign_changes=ch))
        print(rows[-1], flush=True)
    out["TV_exact_rows"] = rows
    json.dump(out, sys.stdout if len(sys.argv) < 2 else open(sys.argv[1], "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "TV_exact_rows"}, indent=1))
