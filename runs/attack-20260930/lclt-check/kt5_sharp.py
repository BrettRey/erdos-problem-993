#!/usr/bin/env python3
"""FLOAT64 DIAGNOSTIC. Rerun route-lclt's KT5 census at one n with the sharp tail term.
Route-lclt certificate:  P(E0)[E h(k-Gd) - TV(h) dK] - 2 P(E0^c)
Sharp form (exact lower bound, verified in exact_AC.py):
                         P(E0)[E h(k-Gd) - TV(h) dK] - (a_{k-1} + a_{k+1}),
                         a_j = P(X = j, E0^c) = sum_{B<B0} sum_s P[s,B] b_{B,p}(j-s).
Reuses route-lclt's joint_JB (2-D float DP) and cert_one unchanged; only the tail is recomputed.
Selection per (tree, k): best over both sides and eta in (0.1, 0.2, 0.35, 0.5), separately
for each form (as route-lclt did for the old form)."""
import sys, json, math, subprocess
sys.path.insert(0, "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/runs/attack-20260930/route-lclt/scripts")
import numpy as np
from scipy.stats import binom
import kt4_mixing_clt as W
from lclt_lib import parent_line_to_adj, tree_polys, window

def sharp_tail(P, p, k, B0):
    tot = 0.0
    for B in range(0, min(B0, P.shape[1])):
        col = P[:, B]
        for s in np.nonzero(col > 0)[0]:
            tot += col[s] * (binom.pmf(k - 1 - s, B, p) + binom.pmf(k + 1 - s, B, p))
    return tot

def run(n, npts=5, etas=(0.1, 0.2, 0.35, 0.5)):
    out = subprocess.run(["gentreeg", "-p", "-q", str(n)], capture_output=True, text=True).stdout
    st = dict(n=n, trees=0, pairs=0, fail_old=0, fail_sharp=0, fail_sharp_trees=0, worst_sharp=None)
    for line in out.splitlines():
        par = list(map(int, line.split()))
        if len(par) != n: continue
        adj = parent_line_to_adj(par)
        Z, ZL, col = tree_polys(adj); alpha = len(Z) - 1
        lo, hi, q = window(n, alpha)
        if hi < lo: continue
        st["trees"] += 1
        ks = sorted(set(lo + round(i * (hi - lo) / (npts - 1)) for i in range(npts)))
        tree_fail = False
        for k in ks:
            lam = math.exp(0.5 * (math.log(Z[k - 1]) - math.log(Z[k + 1]))); p = lam / (1 + lam)
            best_old = -1e9; best_sharp = -1e9
            for even in (True, False):
                P = W.joint_JB(adj, even, lam)
                EB = (P.sum(axis=0) * np.arange(P.shape[1])).sum()
                for eta in etas:
                    B0 = max(1, int(math.floor((1 - eta) * EB)))
                    c = W.cert_one(P, p, k, B0)
                    if c is None or not c["EhG"] > 0: continue
                    main = c["PE0"] * (c["EhG"] - c["TV_h"] * c["dK_Mprime"])
                    old = (main - 2 * (1 - c["PE0"])) / c["EhG"]
                    sh = (main - sharp_tail(P, p, k, B0)) / c["EhG"]
                    best_old = max(best_old, old); best_sharp = max(best_sharp, sh)
            st["pairs"] += 1
            st["fail_old"] += int(best_old <= 0); st["fail_sharp"] += int(best_sharp <= 0)
            if best_sharp <= 0: tree_fail = True
            if st["worst_sharp"] is None or best_sharp < st["worst_sharp"][0]:
                st["worst_sharp"] = (float(best_sharp), line.strip(), int(k))
        st["fail_sharp_trees"] += int(tree_fail)
    return st

if __name__ == "__main__":
    print(json.dumps(run(int(sys.argv[1]))))
