#!/usr/bin/env python3
"""KT1 + KT2 on adversarial families at large n.

Polynomials Z, ZL are exact integers (lclt_lib.tree_polys).  The sign of
delta_k = 1 - i_{k-1}i_{k+1}/i_k^2 is decided in exact integer arithmetic.
Gamma_k = sigma^2 * delta_k and rho_k are DIAGNOSTICS evaluated with 60-digit
mpmath at lam in {i_{k-1}/i_k, i_k/i_{k+1}} (all sums have positive terms;
sigma^2 loses ~log10(n^2) digits to cancellation).  Cross-checked against the
exact Fraction path for n <= 120 in selfcheck().
k is sampled on the central window [ceil(n/4), min(q, alpha-1)] at 9 points.
"""
import sys, json, math, random
from fractions import Fraction
from mpmath import mp, mpf
from lclt_lib import tree_polys, window, moments_at
mp.dps = 60

# ---------- generators (adjacency lists, 0-indexed) ----------
def new(n): return [[] for _ in range(n)]
def edge(adj, u, v): adj[u].append(v); adj[v].append(u)

def star(n):
    a = new(n)
    for v in range(1, n): edge(a, 0, v)
    return a

def path(n):
    a = new(n)
    for v in range(1, n): edge(a, v - 1, v)
    return a

def spider(legs):
    n = 1 + sum(legs); a = new(n); nxt = 1
    for L in legs:
        prev = 0
        for _ in range(L):
            edge(a, prev, nxt); prev = nxt; nxt += 1
    return a

def complete_dary(d, h):
    n = sum(d ** i for i in range(h + 1)); a = new(n)
    for v in range(1, n): edge(a, (v - 1) // d, v)
    return a

def subdivided_binary(h, s):
    """complete binary tree of height h, every edge subdivided s times."""
    base = complete_dary(2, h); nb = len(base)
    edges = [(u, v) for u in range(nb) for v in base[u] if u < v]
    n = nb + s * len(edges); a = new(n); nxt = nb
    for u, v in edges:
        prev = u
        for _ in range(s):
            edge(a, prev, nxt); prev = nxt; nxt += 1
        edge(a, prev, v)
    return a

def bouquet(m, t):
    """T_{m,t}: root with m children, each child carries t pendant P_2 legs."""
    n = 1 + m * (1 + 2 * t); a = new(n); nxt = 1
    for _ in range(m):
        c = nxt; nxt += 1; edge(a, 0, c)
        for _ in range(t):
            x = nxt; y = nxt + 1; nxt += 2
            edge(a, c, x); edge(a, x, y)
    return a

def corona_path(m):
    n = 2 * m; a = new(n)
    for v in range(1, m): edge(a, v - 1, v)
    for v in range(m): edge(a, v, m + v)
    return a

def caterpillar(m, s):
    """spine P_m, each spine vertex with s pendant leaves."""
    n = m * (1 + s); a = new(n); nxt = m
    for v in range(1, m): edge(a, v - 1, v)
    for v in range(m):
        for _ in range(s): edge(a, v, nxt); nxt += 1
    return a

def broom(plen, s):
    n = plen + s; a = new(n)
    for v in range(1, plen): edge(a, v - 1, v)
    for j in range(s): edge(a, plen - 1, plen + j)
    return a

def double_star(s1, s2):
    n = 2 + s1 + s2; a = new(n); edge(a, 0, 1); nxt = 2
    for _ in range(s1): edge(a, 0, nxt); nxt += 1
    for _ in range(s2): edge(a, 1, nxt); nxt += 1
    return a

def star_of_stars(m, s):
    """root joined to m centres, each centre with s leaves (hub-of-hubs)."""
    n = 1 + m * (1 + s); a = new(n); nxt = 1
    for _ in range(m):
        c = nxt; nxt += 1; edge(a, 0, c)
        for _ in range(s): edge(a, c, nxt); nxt += 1
    return a

def random_tree(n, seed):
    rng = random.Random(seed)
    if n <= 2: return path(n)
    pr = [rng.randrange(n) for _ in range(n - 2)]
    deg = [1] * n
    for x in pr: deg[x] += 1
    a = new(n)
    import heapq
    leaves = [i for i in range(n) if deg[i] == 1]; heapq.heapify(leaves)
    for x in pr:
        l = heapq.heappop(leaves); edge(a, l, x); deg[x] -= 1
        if deg[x] == 1: heapq.heappush(leaves, x)
    u = heapq.heappop(leaves); v = heapq.heappop(leaves); edge(a, u, v)
    return a

# ---------- evaluation ----------
def mom_mp(Z, ZL, lam):
    S0 = S1 = S2 = SL = mpf(0); x = mpf(1)
    for j, c in enumerate(Z):
        t = c * x
        S0 += t; S1 += j * t; S2 += j * j * t
        if j < len(ZL): SL += ZL[j] * x
        x *= lam
    mu = S1 / S0; s2 = S2 / S0 - mu * mu; muL = SL / S0
    return mu, s2, muL

def analyse_family(adj, npts=9):
    n = len(adj)
    Z, ZL, _ = tree_polys(adj)
    alpha = len(Z) - 1
    lo, hi, q = window(n, alpha)
    if hi < lo: return None
    ks = sorted(set(lo + round(i * (hi - lo) / (npts - 1)) for i in range(npts)))
    out = []
    for k in ks:
        im, i0, ip = Z[k - 1], Z[k], Z[k + 1]
        dnum = i0 * i0 - im * ip       # exact sign
        delta = mpf(dnum) / (mpf(i0) * mpf(i0))
        g = []; rh = []
        for lam in (mpf(im) / mpf(i0), mpf(i0) / mpf(ip)):
            mu, s2, muL = mom_mp(Z, ZL, lam)
            g.append(s2 * delta)
            rh.append(max(muL, mu - muL) / ((1 + lam) * s2))
        out.append(dict(k=k, lc_exact=dnum > 0, Gamma_lo=float(min(g)), Gamma_hi=float(max(g)),
                        rho=float(min(rh)), lam_hi=float(mpf(i0) / mpf(ip))))
    return dict(n=n, alpha=alpha, window=[lo, hi, q], rows=out,
                min_Gamma=min(r["Gamma_lo"] for r in out),
                max_dev=max(max(1 - r["Gamma_lo"], r["Gamma_hi"] - 1) for r in out),
                min_rho=min(r["rho"] for r in out),
                max_lam=max(r["lam_hi"] for r in out),
                lc_all=all(r["lc_exact"] for r in out))

def selfcheck():
    from lclt_lib import analyse
    for adj in (complete_dary(3, 3), bouquet(4, 3), spider([2] * 12), random_tree(60, 3), subdivided_binary(3, 2)):
        Z, alpha, w, rows = analyse(adj)
        fam = analyse_family(adj, npts=50)
        ex = {r["k"]: r for r in rows}
        for r in fam["rows"]:
            e = ex[r["k"]]
            assert abs(float(e["Gamma_lo"]) - r["Gamma_lo"]) < 1e-12 and abs(float(e["rho"]) - r["rho"]) < 1e-12
    print("selfcheck ok: mpmath diagnostics match exact Fractions", file=sys.stderr)

FAMILIES = {
    "star": [(star, (n,)) for n in (50, 100, 200, 400, 800)],
    "path": [(path, (n,)) for n in (50, 100, 200, 400, 800)],
    "spider_2^m": [(spider, ([2] * m,)) for m in (25, 50, 100, 200, 400)],
    "spider_3^m": [(spider, ([3] * m,)) for m in (17, 33, 67, 133, 267)],
    "corona_path": [(corona_path, (m,)) for m in (25, 50, 100, 200, 400)],
    "caterpillar_s2": [(caterpillar, (m, 2)) for m in (17, 33, 67, 133, 267)],
    "broom_half": [(broom, (n // 2, n // 2)) for n in (50, 100, 200, 400, 800)],
    "double_star": [(double_star, (n // 2 - 1, n // 2 - 1)) for n in (50, 100, 200, 400, 800)],
    "binary_complete": [(complete_dary, (2, h)) for h in (5, 6, 7, 8, 9)],
    "ternary_complete": [(complete_dary, (3, h)) for h in (3, 4, 5, 6)],
    "4ary_complete": [(complete_dary, (4, h)) for h in (2, 3, 4, 5)],
    "5ary_complete": [(complete_dary, (5, h)) for h in (2, 3, 4)],
    "7ary_complete": [(complete_dary, (7, h)) for h in (2, 3)],
    "binary_subdiv2": [(subdivided_binary, (h, 2)) for h in (3, 4, 5, 6, 7)],
    "binary_subdiv4": [(subdivided_binary, (h, 4)) for h in (3, 4, 5, 6)],
    "bouquet_T_m,t_balanced": [(bouquet, ((3 ** t) // (2 ** t), t)) for t in (4, 5, 6, 7, 8)],
    "bouquet_T_m,4": [(bouquet, (m, 4)) for m in (5, 10, 20, 40, 80)],
    "star_of_stars_s=m": [(star_of_stars, (m, m)) for m in (7, 10, 14, 20, 28)],
    "random": [(random_tree, (n, s)) for n in (100, 200, 400, 800) for s in (1, 2, 3)],
}

if __name__ == "__main__":
    selfcheck()
    names = sys.argv[1:] or list(FAMILIES)
    for name in names:
        for gen, args in FAMILIES[name]:
            adj = gen(*args)
            if len(adj) > 1600: continue
            r = analyse_family(adj)
            if r is None: continue
            r["family"] = name; r["args"] = [a if not isinstance(a, list) else f"{a[0]}^{len(a)}" for a in args]
            print(json.dumps(r), flush=True)
