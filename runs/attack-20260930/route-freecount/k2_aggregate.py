"""K2: per-vertex PV defects on LC-failing census trees and on the hub-star
counterexample family; aggregate diagnostics and one natural repair.

D_v(k) = k i_{k-1}(T) j^v_k - (k+1) i_k(T) j^v_{k-1},  j^v = I(T - N[v]).
Exact: sum_v D_v(k) = k (k+1) (i_{k-1} i_{k+1} - i_k^2).
PV_v(k) <=> D_v(k) <= 0.

Repair R1 (closed-neighbourhood averaging): L_u = sum_{v in N[u]} D_v/|N[v]|.
sum_u L_u = sum_v D_v, so "L_u <= 0 for all u" is a sufficient condition for LC.

Usage: python3 k2_aggregate.py census|hubstars > out.jsonl
"""

import glob
import json
import sys
from fractions import Fraction

from pv_lib import co, parse_parent_array, tree_data, window
from k3_families import hub_stars

ROOT = "/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993"


def analyse(adj, label, ks_extra=("am1", "am2")):
    n = len(adj)
    I, J, E, alpha = tree_data(adj)
    lo, q = window(n, alpha)
    lo = max(lo, 1)
    ks = list(range(lo, q + 1))
    extra = []
    if "am1" in ks_extra and alpha - 1 >= 1:
        extra.append(alpha - 1)
    if "am2" in ks_extra and alpha - 2 >= 1:
        extra.append(alpha - 2)
    out = {"label": label, "n": n, "alpha": alpha, "lo": lo, "q": q, "levels": []}
    for k in sorted(set(ks + extra)):
        D = [k * co(I, k - 1) * co(J[v], k) - (k + 1) * co(I, k) * co(J[v], k - 1)
             for v in range(n)]
        S = sum(D)
        assert S == k * (k + 1) * (co(I, k - 1) * co(I, k + 1) - co(I, k) ** 2)
        pos = sum(d for d in D if d > 0)
        neg = -sum(d for d in D if d < 0)
        fails = [v for v in range(n) if D[v] > 0]
        maxdeg = max(len(a) for a in adj)
        # R1
        L = []
        for u in range(n):
            s = Fraction(D[u], len(adj[u]) + 1)
            for v in adj[u]:
                s += Fraction(D[v], len(adj[v]) + 1)
            L.append(s)
        r1_fail = [u for u in range(n) if L[u] > 0]
        pv_weight = [Fraction(co(J[v], k - 1), co(I, k - 1)) if co(I, k - 1) else 0 for v in range(n)]
        out["levels"].append({
            "k": k, "in_window": lo <= k <= q, "k_minus_alpha": k - alpha,
            "lc_holds": S <= 0,
            "n_pv_fail": len(fails),
            "fail_degrees": sorted(len(adj[v]) for v in fails),
            "fail_is_maxdeg": [len(adj[v]) == maxdeg for v in fails],
            # diagnostics only (floats):
            "pos_over_neg": (float(Fraction(pos, neg)) if neg else None),
            "fail_weight_share": (float(sum(pv_weight[v] for v in fails) / sum(pv_weight)) if fails else 0.0),
            "n_r1_fail": len(r1_fail),
        })
    return out


def census_trees():
    for f in sorted(glob.glob(f"{ROOT}/results/lc_census_20260814/n*_s*of64.txt")):
        for line in open(f):
            if line.startswith("FAIL"):
                par = line.split("par=")[1].split()[0].split(",")
                yield f"census:{f.split('/')[-1]}", parse_parent_array(par)


def main():
    mode = sys.argv[1]
    if mode == "census":
        seen = 0
        stride = int(sys.argv[2]) if len(sys.argv) > 2 else 1
        for i, (lab, adj) in enumerate(census_trees()):
            if i % stride:
                continue
            print(json.dumps(analyse(adj, lab)), flush=True)
            seen += 1
    elif mode == "hubstars":
        for t in [2, 3]:
            for m in range(6, 26 if t == 2 else 20):
                print(json.dumps(analyse(hub_stars(m, t), f"hubstars_{m}_{t}")), flush=True)


if __name__ == "__main__":
    main()
