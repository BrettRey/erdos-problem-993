#!/usr/bin/env python3
"""DIAGNOSTIC (floats): growth of the root quantity q(1-q) delta^2 (Prop 4.1) on complete
d-ary trees of depth h, at lam in {1/4, lam_c(d), 12}. Recursions (4.1)-(4.2):
 r_v = lam prod (1+r_i)^{-1}, q = r/(1+r), delta_v = 1 - sum q_i delta_i.
Levels are homogeneous, so each level is one number. n = (d^{h+1}-1)/(d-1).
Prop 4.1 claims q(1-q) delta^2 <= C n^a with a<1; here we just look at the observed growth."""
import math, json
def run(d, lam, H):
    r, delta = lam, 1.0         # leaf
    rows=[]
    for h in range(1, H+1):
        q = r/(1+r)
        r_new = lam/(1+r)**d
        delta_new = 1 - d*q*delta
        r, delta = r_new, delta_new
        qq = r/(1+r)
        n = (d**(h+1)-1)//(d-1) if d>1 else h+1
        rows.append((h, n, qq*(1-qq)*delta*delta, delta))
    return rows
out={}
for d in [2,3,4,5,6]:
    lamc = d**d/(d-1)**(d+1)
    for lam in [0.25, lamc, 12.0]:
        H = int(60/math.log2(d))
        rows = run(d, lam, H)
        mx = max(x[2] for x in rows)
        last = rows[-1]
        out[f"d={d},lam={lam:.4f}"] = dict(max_qq_delta2=mx, depth=last[0], log10_n=math.log10(last[1]), last=last[2])
        print(f"d={d} lam={lam:7.4f} depth={last[0]:3d} log10n={math.log10(last[1]):6.2f}  max q(1-q)delta^2 over depths={mx:10.4f}  at max depth={last[2]:.4f}")
json.dump(out, open("data/delta_diag.json","w"), indent=1)
