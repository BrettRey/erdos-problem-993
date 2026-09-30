"""Independent recheck of an R1 violation on the star of hubs S(h,m,t).
Method 1 (closed form, pure-Python integer polynomials, no flint, no tree DP):
  cherry c = vertex with t leaves:   Ic = (1+x)^t + x
  hub-star A = I(H(m,t)) = Ic^m + x (1+x)^{tm};  B = I(H(m,t) - hub) = Ic^m
  I(T)        = A^h + x B^h                      (centre out / centre in)
  J_centre    = I(T - N[centre]) = Ic^{hm}
  J_hub       = I(T - N[hub])    = (1+x)^{tm} A^{h-1}
Method 2: route-freecount/pv_lib.indep_seq (fresh tree DP) on the explicit
edge list for T, T - N[centre] and T - N[hub_i] for every hub i.
Then L_centre(k) = D_centre/(h+1) + sum_i D_hub_i/(m+2) is evaluated with Fractions.
Usage: python recheck_independent.py h m t
"""
import sys, json
from fractions import Fraction
sys.path.insert(0, '../route-freecount')
import pv_lib

def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out
def padd(a, b):
    if len(a) < len(b): a, b = b, a
    out = list(a)
    for i, y in enumerate(b): out[i] += y
    return out
def ppow(a, e):
    r = [1]
    for _ in range(e): r = pmul(r, a)
    return r
def co(s, k): return s[k] if 0 <= k < len(s) else 0

def main():
    h, m, t = map(int, sys.argv[1:4])
    X = [0, 1]
    onex = [1, 1]
    Ic = padd(ppow(onex, t), X)
    A = padd(ppow(Ic, m), pmul(X, ppow(onex, t * m)))
    Bm = ppow(Ic, m)
    I = padd(ppow(A, h), pmul(X, ppow(Bm, h)))
    Jc = ppow(Ic, h * m)
    Jh = pmul(ppow(onex, t * m), ppow(A, h - 1))

    # explicit tree: vertex 0 = centre
    edges = []
    nxt = 1
    hubs = []
    for _ in range(h):
        hub = nxt; nxt += 1; edges.append((0, hub)); hubs.append(hub)
        for _ in range(m):
            c = nxt; nxt += 1; edges.append((hub, c))
            for _ in range(t):
                edges.append((c, nxt)); nxt += 1
    n = nxt
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    V = set(range(n))
    I2 = pv_lib.indep_seq(adj, V)
    Jc2 = pv_lib.indep_seq(adj, V - {0} - set(adj[0]))
    Jh2 = [pv_lib.indep_seq(adj, V - {w} - set(adj[w])) for w in hubs]
    strip = lambda s: s[:max(i for i, x in enumerate(s) if x) + 1]
    assert strip(I) == strip(I2), "I mismatch"
    assert strip(Jc) == strip(Jc2), "J_centre mismatch"
    assert all(strip(Jh) == strip(j) for j in Jh2), "J_hub mismatch"
    alpha = len(strip(I)) - 1
    lo, q = -((-n) // 4), -((-(2 * alpha - 1)) // 3)
    res = {"h": h, "m": m, "t": t, "n": n, "alpha": alpha, "window": [lo, q], "methods_agree": True, "levels": []}
    for k in range(lo, q + 1):
        D = lambda J: k * co(I, k - 1) * co(J, k) - (k + 1) * co(I, k) * co(J, k - 1)
        Dc, Dh = D(Jc), D(Jh)
        L = Fraction(Dc, h + 1) + h * Fraction(Dh, m + 2)
        lcm = co(I, k) ** 2 - co(I, k - 1) * co(I, k + 1)
        if L > 0:
            res["levels"].append({"k": k, "L_centre": str(L), "D_centre": str(Dc), "D_hub": str(Dh),
                                  "lc_margin_ik2_minus_prod": str(lcm), "lc_margin_rel": float(Fraction(lcm, co(I, k) ** 2)),
                                  "L_over_k_ik1_ik": float(L / (k * co(I, k - 1) * co(I, k)))})
    res["edges"] = edges
    print(json.dumps(res))


if __name__ == '__main__':
    main()
