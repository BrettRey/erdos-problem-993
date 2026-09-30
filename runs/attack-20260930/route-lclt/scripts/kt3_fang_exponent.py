#!/usr/bin/env python3
"""KT3: how small can Fang et al.'s exponent a be, and what n does the
dominant CLT error term n^{(a-1)/8} force?  Diagnostic, 80-digit mpmath.

Source (arXiv:2609.20961v1, Lemma 4.2 and its proof, pp. 5-6; rate display
end of Section 5, p. 10):
  * choose b in (0,1) with 12 e^{-1-3b/2} < b                    (proof of 4.2)
  * rho = (2 sqrt(12)/e)^{4-2p} * b^{2p-3} < 1,  p in (3/2, 2)    (proof of 4.2)
  * a = 2 - 2/p                                                  (proof of Prop 4.1)
  * CLT error  O( b/sqrt(n) + b^{(a-1)/2} + b^{a-1} + n^{a-1} ),
    with b = ceil(n^{1/4})                                       (p. 10)
Every implied constant is set to 1 here, so the n reported is a LOWER bound
for what that architecture needs (it ignores C, C0, C1, G, H, c, Raic's
Berry-Esseen constant and the Fourier step).
We also report the same computation for a general fugacity ceiling Lmax,
and for the best-case balance of b against n.
"""
import json, sys
from mpmath import mp, mpf, e, sqrt, log, findroot, exp, log10
mp.dps = 80

def exponent_for(Lmax):
    Lmax = mpf(Lmax)
    # b-constraint: Lmax * e^{-1-3b/2} < b.  (Paper uses 12, i.e. lambda <= 12,
    # through y r/(1+3r/2) <= b with r = lambda e^{-y}; maximum lambda e^{-1-3b/2}.)
    f = lambda b: Lmax*exp(-1-3*b/2) - b
    # find root in (0, 1] if exists
    if f(mpf(1)) > 0:
        return None  # exponent-2 inequality unavailable with b<1
    lo, hi = mpf('1e-30'), mpf(1)
    if f(lo) <= 0:
        bstar = lo
    else:
        for _ in range(300):
            mid = (lo+hi)/2
            if f(mid) > 0: lo = mid
            else: hi = mid
        bstar = hi
    A = log(2*sqrt(Lmax)/e)   # exponent-3/2 constant log
    lb = log(bstar)
    if A <= 0:
        # exponent 3/2 already gives rho<1: p can be taken at 3/2
        return dict(Lmax=float(Lmax), bstar=float(bstar), p_min=1.5, a_min=float(2-2/mpf(1.5)))
    # rho<1  <=>  (4-2p)A + (2p-3) lb < 0  <=> p > (4A - 3lb)/(2A - 2lb)
    pmin = (4*A - 3*lb)/(2*A - 2*lb)
    amin = 2 - 2/pmin
    return dict(Lmax=float(Lmax), bstar=float(bstar), A=float(A), p_min=float(pmin), a_min=float(amin))

out = {}
for L in [12, 8, 6, 4, 3, 2, mpf(e**2)/4]:
    r = exponent_for(L)
    if r is None:
        out[str(float(L))] = None; continue
    a = mpf(r['a_min'])
    # dominant term with b = n^{1/4}: n^{(a-1)/8} <= eps  <=> log10 n >= 8 log10(1/eps)/(1-a)
    for eps in ['0.1']:
        r['log10_n_for_term_le_0.1_b_eq_n^(1/4)'] = float(8*log10(1/mpf(eps))/(1-a))
    # best balance: choose b = n^beta to equalise b/sqrt(n) and b^{(a-1)/2}:
    # beta - 1/2 = beta (a-1)/2  => beta = 1/(3-a); rate n^{-(1-a)/(2(3-a))}
    rate = (1-a)/(2*(3-a))
    r['best_balance_rate_exponent'] = float(rate)
    r['log10_n_for_best_balance_le_0.1'] = float(1/rate)
    out[str(float(L))] = r

# Sensitivity: sharp exponent a for complete d-ary trees at the fixed point (heuristic,
# not a proof): delta ~ (d q*)^H, n ~ d^H  => q(1-q) delta^2 ~ n^{2 log(d q*)/log d}
sharp = {}
for d in range(2, 13):
    lo, hi = mpf(0), mpf(12)
    for _ in range(300):
        mid=(lo+hi)/2
        if mid*(1+mid)**d < 12: lo=mid
        else: hi=mid
    r = lo
    q = r/(1+r)
    sharp[d] = float(2*log(d*q)/log(d)) if d*q > 1 else 0.0
out['heuristic_fixed_point_exponent_complete_dary_lambda12'] = sharp
json.dump(out, sys.stdout, indent=1)
print()
