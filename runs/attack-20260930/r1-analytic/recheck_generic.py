"""Independent recheck of an R1 counterexample with the wave-1 generic tree DP (pv_lib),
on the explicit tree (parent array written to disk). Exact integers / Fractions only."""
import sys, json
from fractions import Fraction
sys.path.insert(0, "../route-freecount")
from pv_lib import tree_data, co, window
import fam

h, m, s = (int(a) for a in sys.argv[1:4])
out = sys.argv[4]
par = [-1]
for _ in range(h):
    hub = len(par); par.append(0)
    for _ in range(m):
        c = len(par); par.append(hub)
        for _ in range(s):
            par.append(c)
adj = fam._adj(par)
n = len(adj)
I, J, _E, alpha = tree_data(adj)
lo, q = window(n, alpha)
res = {"family": f"SH({h},{m},{s})", "n": n, "alpha": alpha, "window": [lo, q],
       "parent_array_1indexed_0root": [p + 1 for p in par],
       "I": [str(c) for c in I], "levels": []}
for k in range(lo, q + 1):
    D = [k * co(I, k - 1) * co(Jv, k) - (k + 1) * co(I, k) * co(Jv, k - 1) for Jv in J]
    lc = co(I, k - 1) * co(I, k + 1) <= co(I, k) ** 2
    bad = []
    for u in range(n):
        L = Fraction(D[u], len(adj[u]) + 1) + sum(Fraction(D[v], len(adj[v]) + 1) for v in adj[u])
        if L > 0:
            bad.append({"u": u, "deg": len(adj[u]), "L_num": str(L.numerator), "L_den": str(L.denominator)})
    res["levels"].append({"k": k, "lc_holds": lc, "lc_margin": str(co(I, k) ** 2 - co(I, k - 1) * co(I, k + 1)),
                          "n_r1_fail": len(bad), "r1_fail": bad[:5],
                          "D_centre": str(D[0]), "D_hub": str(D[1])})
json.dump(res, open(out, "w"), indent=1)
print(json.dumps([(l["k"], l["lc_holds"], l["n_r1_fail"], [b["u"] for b in l["r1_fail"]]) for l in res["levels"] if l["n_r1_fail"]]))
