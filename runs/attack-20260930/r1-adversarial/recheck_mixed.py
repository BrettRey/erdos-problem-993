"""Independent recheck for a star of hubs with non-uniform hubs: centre (vertex 0)
joined to hubs w_1..w_h, hub w_i carrying m_i cherries (vertex + 2 leaves).
Method 1 closed form (pure-Python integer polynomials):
  Ic = (1+x)^2 + x;  A_i = Ic^{m_i} + x (1+x)^{2 m_i};  B_i = Ic^{m_i}
  I(T) = prod A_i + x prod B_i;  J_centre = prod B_i;  J_{w_i} = (1+x)^{2 m_i} prod_{j != i} A_j
Method 2: route-freecount/pv_lib.indep_seq on the explicit tree.
Usage: python recheck_mixed.py m_1 m_2 ... m_h  (prints JSON)"""
import sys, json
from fractions import Fraction
sys.path.insert(0, '../route-freecount')
import pv_lib
from recheck_independent import pmul, padd, ppow, co  # noqa (module runs on import only under __main__ guard below)

def main():
    ms = [int(a) for a in sys.argv[1:]]
    h = len(ms)
    X = [0, 1]; onex = [1, 1]
    Ic = padd(ppow(onex, 2), X)
    A = [padd(ppow(Ic, m), pmul(X, ppow(onex, 2 * m))) for m in ms]
    Bs = [ppow(Ic, m) for m in ms]
    prodA = [1]
    for a in A: prodA = pmul(prodA, a)
    prodB = [1]
    for b in Bs: prodB = pmul(prodB, b)
    I = padd(prodA, pmul(X, prodB))
    Jc = prodB
    Jh = []
    for i, m in enumerate(ms):
        p = ppow(onex, 2 * m)
        for j, a in enumerate(A):
            if j != i: p = pmul(p, a)
        Jh.append(p)
    edges = []; nxt = 1; hubs = []
    for m in ms:
        w = nxt; nxt += 1; edges.append((0, w)); hubs.append(w)
        for _ in range(m):
            c = nxt; nxt += 1; edges.append((w, c))
            for _ in range(2):
                edges.append((c, nxt)); nxt += 1
    n = nxt
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    V = set(range(n))
    strip = lambda s: s[:max(i for i, x in enumerate(s) if x) + 1]
    assert strip(I) == strip(pv_lib.indep_seq(adj, V)), "I mismatch"
    assert strip(Jc) == strip(pv_lib.indep_seq(adj, V - {0} - set(adj[0]))), "Jc mismatch"
    for i, w in enumerate(hubs):
        assert strip(Jh[i]) == strip(pv_lib.indep_seq(adj, V - {w} - set(adj[w]))), ("Jh mismatch", i)
    alpha = len(strip(I)) - 1
    lo, q = -((-n) // 4), -((-(2 * alpha - 1)) // 3)
    res = {"hubs_m": ms, "n": n, "alpha": alpha, "window": [lo, q], "methods_agree": True, "levels": []}
    for k in range(lo, q + 1):
        D = lambda J: k * co(I, k - 1) * co(J, k) - (k + 1) * co(I, k) * co(J, k - 1)
        L = Fraction(D(Jc), h + 1) + sum(Fraction(D(Jh[i]), ms[i] + 2) for i in range(h))
        if L > 0:
            lcm = co(I, k) ** 2 - co(I, k - 1) * co(I, k + 1)
            res["levels"].append({"k": k, "L_centre": str(L), "D_centre": str(D(Jc)), "D_hubs": [str(D(J)) for J in Jh],
                                  "lc_margin_ik2_minus_prod": str(lcm), "lc_margin_rel": float(Fraction(lcm, co(I, k) ** 2)),
                                  "L_over_k_ik1_ik": float(L / (k * co(I, k - 1) * co(I, k)))})
    res["edges"] = edges
    print(json.dumps(res))

if __name__ == '__main__':
    main()
