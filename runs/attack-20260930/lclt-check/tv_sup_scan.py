#!/usr/bin/env python3
"""FLOAT64 DIAGNOSTIC: sup of TV(h_B) * v^{3/2} over B = 1..Bmax and p = lam/(1+lam),
lam on a grid in [1/3, 12] (the fugacity range the window needs, Fang Prop 7.1).
h_B = -Delta^2 b_{B,p};  v = B p (1-p).  Limit value I3 = int|phi'''| = 1.5100130.
Purpose: an EFFECTIVE certificate needs a uniform constant C_TV >= sup TV v^{3/2};
c_eff = phi(0)/C_TV replaces c_* = phi(0)/I3.  The exact rational check in
kernel_check_exact.py already shows sup > I3 (B=256, p=12/13)."""
import json, sys
import numpy as np
from scipy.stats import binom
I3 = 1.5100130001304771
Bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
lams = sorted(set([1/3, 0.5, 1, 2, 3, 4, 5, 6, 8, 10, 11, 12] + list(np.linspace(1/3, 12, 60))))
best = []; overall = (0, None)
for lam in lams:
    p = lam / (1 + lam); q = 1 - p
    bl = (0, None)
    for B in range(1, Bmax + 1):
        b = binom.pmf(np.arange(B + 1), B, p)
        ext = np.concatenate([np.zeros(3), b, np.zeros(3)])
        h = 2 * ext[1:-1] - ext[:-2] - ext[2:]
        tv = np.abs(np.diff(np.concatenate([[0.0], h, [0.0]]))).sum()
        val = tv * (B * p * q) ** 1.5
        if val > bl[0]: bl = (val, B)
    best.append(dict(lam=lam, p=p, sup_TVv32=bl[0], argmax_B=bl[1]))
    if bl[0] > overall[0]: overall = (bl[0], (lam, bl[1]))
out = dict(Bmax=Bmax, I3=I3, overall_sup=overall[0], overall_arg=overall[1], per_lam=best)
json.dump(out, open("rerun/tv_sup_scan.json", "w"), indent=1)
print("overall sup TV v^1.5 =", overall, " ratio to I3 =", overall[0] / I3)
for r in best[::6]: print(r)
