#!/usr/bin/env python3
"""Kernel constants for claim C (DIAGNOSTIC except where marked certified).

1. int |phi'''| for phi the standard normal density:
   closed form 2(phi(0) + 4 phi(sqrt3)) (zeros of phi''' = (3x - x^3) phi at 0, +-sqrt3),
   checked (a) by mpmath quadrature on the monotone pieces, (b) by python-flint arb
   ball evaluation of the closed form (certified enclosure), (c) c_* = phi(0)/int|phi'''|.
2. TV(h_B) v^{3/2}, h_B = -Delta^2 b_{B,p}, v = B p (1-p), for B up to 10^4 and
   p in {1/4, 1/2, 3/4, 12/13}: raw sum of |differences| (mpmath, 60 digits), and the
   Polya-frequency formula TV = 2(max h + |min_left h| + |min_right h|) (valid when
   Delta^3 b has <= 3 sign changes, which is also counted).
"""
import json, sys
import mpmath as mp
mp.mp.dps = 60

phi = lambda x: mp.exp(-x * x / 2) / mp.sqrt(2 * mp.pi)
phi3 = lambda x: (3 * x - x ** 3) * phi(x)
s3 = mp.sqrt(3)
closed = 2 * (phi(0) + 4 * phi(s3))
quad = mp.quad(lambda x: abs(phi3(x)), [-mp.inf, -s3, 0, s3, mp.inf])
out = dict(int_abs_phi3_closed=mp.nstr(closed, 20), int_abs_phi3_quad=mp.nstr(quad, 20),
           c_star=mp.nstr(phi(0) / closed, 20))
try:
    from flint import arb
    a0 = 1 / arb.pi().__mul__(2).sqrt()
    a3 = (arb(-1.5)).exp() * a0
    I = 2 * (a0 + 4 * a3)
    out["int_abs_phi3_arb"] = I.str(20); out["c_star_arb"] = (a0 / I).str(20)
except Exception as e:
    out["arb_error"] = repr(e)

def tv_h(B, p):
    p = mp.mpf(p); q = 1 - p
    # pmf via recurrence in log space
    b = [mp.mpf(0)] * (B + 1)
    lb = [mp.loggamma(B + 1) - mp.loggamma(j + 1) - mp.loggamma(B - j + 1) + j * mp.log(p) + (B - j) * mp.log(q) for j in range(B + 1)]
    b = [mp.exp(x) for x in lb]
    ext = [mp.mpf(0)] * 3 + b + [mp.mpf(0)] * 3          # j = -3 .. B+3
    h = [2 * ext[i] - ext[i - 1] - ext[i + 1] for i in range(1, len(ext) - 1)]   # j = -2 .. B+2
    hh = [mp.mpf(0)] + h + [mp.mpf(0)]
    tv_raw = sum(abs(hh[i + 1] - hh[i]) for i in range(len(hh) - 1))
    d3 = [hh[i + 1] - hh[i] for i in range(len(hh) - 1)]
    tol = mp.mpf(10) ** (-50)
    signs = [1 if x > tol else -1 for x in d3 if abs(x) > tol]
    changes = sum(1 for i in range(len(signs) - 1) if signs[i] != signs[i + 1])
    imax = max(range(len(h)), key=lambda i: h[i])
    tv_pf = 2 * (h[imax] + abs(min(h[:imax + 1])) + abs(min(h[imax:])))
    v = B * p * q
    return dict(B=B, p=float(p), v=float(v), TVv32=float(tv_raw * v ** 1.5), TVv32_PF=float(tv_pf * v ** 1.5),
                dsign_changes=changes, maxh_v32=float(h[imax] * v ** 1.5))

rows = []
for p in [mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4, mp.mpf(12) / 13]:
    for B in [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 10000]:
        mp.mp.dps = 40
        r = tv_h(B, p); rows.append(r)
lim = float(closed)
for r in rows:
    r["excess"] = r["TVv32"] / lim - 1
    r["v_times_excess"] = r["excess"] * r["v"]
    r["maxh_excess"] = r["maxh_v32"] / float(phi(0)) - 1
out["TV_rows"] = rows
json.dump(out, sys.stdout, indent=1)
