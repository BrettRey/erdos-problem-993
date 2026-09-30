"""Independent recheck (wave-1 generic tree DP, pv_lib) of a heterogeneous star-of-hubs tree.
argv: JSON list of [count, m, s], output path. Full window scanned; exact integers/Fractions only."""
import sys, json
from fractions import Fraction
sys.path.insert(0, "../route-freecount")
from pv_lib import tree_data, co, window
import fam, mixed
types = [tuple(t) for t in json.loads(sys.argv[1])]; out = sys.argv[2]
par = [-1]
for c, m, s in types:
    for _ in range(c):
        hub = len(par); par.append(0)
        for _ in range(m):
            v = len(par); par.append(hub)
            for _ in range(s):
                par.append(v)
adj = fam._adj(par); n = len(adj)
I, J, _E, alpha = tree_data(adj)
lo, q = window(n, alpha)
F = mixed.mixed(types)
assert [int(F.I[i]) for i in range(F.I.degree() + 1)] == I and F.n == n
res = {"types_[count,m,s]": [list(t) for t in types], "n": n, "alpha": alpha, "window": [lo, q],
       "parent_array_1indexed_0root": [p + 1 for p in par],
       "edge_list_0indexed": [[i, p] for i, p in enumerate(par) if p >= 0],
       "I": [str(c) for c in I], "levels": []}
for k in range(lo, q + 1):
    D = [k * co(I, k - 1) * co(Jv, k) - (k + 1) * co(I, k) * co(Jv, k - 1) for Jv in J]
    lc = co(I, k - 1) * co(I, k + 1) <= co(I, k) ** 2
    bad = []
    for u in range(n):
        L = Fraction(D[u], len(adj[u]) + 1) + sum(Fraction(D[v], len(adj[v]) + 1) for v in adj[u])
        if L > 0:
            nb = sum(Fraction(abs(D[v]), len(adj[v]) + 1) for v in adj[u] + [u])
            bad.append({"u": u, "deg": len(adj[u]), "L": f"{L.numerator}/{L.denominator}",
                        "L_over_nbhd_abs_mass_float": float(L / nb)})
    Dc, Lc, *_ = F.level(k)
    Lg0 = Fraction(D[0], len(adj[0]) + 1) + sum(Fraction(D[v], len(adj[v]) + 1) for v in adj[0])
    assert Lg0 == Lc["z"], ("closed form vs generic mismatch", k)
    res["levels"].append({"k": k, "lc_holds": lc, "n_r1_fail": len(bad), "r1_fail": bad,
                          "n_pv_fail": sum(1 for d in D if d > 0)})
json.dump(res, open(out, "w"), indent=1)
print(n, alpha, [lo, q], "LC all window:", all(l["lc_holds"] for l in res["levels"]),
      [(l["k"], [(b["u"], b["deg"], round(b["L_over_nbhd_abs_mass_float"], 5)) for b in l["r1_fail"]]) for l in res["levels"] if l["n_r1_fail"]])
