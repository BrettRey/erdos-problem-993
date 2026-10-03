"""Probe: Vatter's first-letter argument (arXiv:2608.22147) transplanted to trees.

Fix a vertex order and classify independent sets by their least vertex:
i_{k+1}(T) = sum_v i_k(G_v), where G_v is the subgraph induced on the later
non-neighbours of v. Then rho_{k+1}(T) is a weighted average of the rho_k(G_v),
rho_k = i_k / i_{k-1}, with weights i_{k-1}(G_v). So T is log-concave whenever
some order has rho_k(G_v) <= rho_k(T) for every v and k with positive weight
(a "Vatter order"). Vatter's Claim 2 gives this for words; for trees it cannot
hold in general, because it implies log-concavity, which fails at order 26.

has_chain searches for a Vatter order by DFS over the set of vertices still to
come, memoised on that set; every tail is compared with the full tree's
polynomial. Exact integer and Fraction arithmetic throughout.

Run (from the project root): venv/bin/python scripts/probe_vatter_order_certificate_20261002.py N [N_start]
Result 2026-10-02: every tree with 2..13 vertices has a Vatter order.
"""
import sys
from fractions import Fraction
from functools import lru_cache
sys.path.insert(0, '.')
from trees import trees

def indpoly_subset(adj, verts):
    verts = tuple(sorted(verts))
    @lru_cache(None)
    def f(vs):
        if not vs: return (1,)
        v, rest = vs[0], vs[1:]
        a = f(rest)
        b = f(tuple(u for u in rest if u not in adj[v]))
        n = max(len(a), len(b) + 1)
        return tuple((a[i] if i < len(a) else 0) + (b[i-1] if 0 <= i-1 < len(b) else 0) for i in range(n))
    return f(verts)

def ratio_ok(small, big):
    # rho_k(small) <= rho_k(big) for all k >= 1 where small has i_{k-1} > 0
    for k in range(1, len(small)):
        if small[k-1] == 0: continue
        if k >= len(big) or big[k-1] == 0:
            if small[k] > 0: return False
            continue
        if Fraction(small[k], small[k-1]) > Fraction(big[k], big[k-1]): return False
    return True

def has_chain(adj, verts, P, memo):
    """Is there an order of verts (the vertices still to come) such that every vertex's tail
    G_u (later non-neighbours) satisfies rho_k(G_u) <= rho_k(P) for all k? P is the full tree."""
    key = frozenset(verts)
    if key in memo: return memo[key]
    if not verts: memo[key] = True; return True
    for v in verts:
        gv = [u for u in verts if u != v and u not in adj[v]]
        rest = [u for u in verts if u != v]
        if ratio_ok(indpoly_subset(adj, gv), P) and has_chain(adj, rest, P, memo):
            memo[key] = True; return True
    memo[key] = False; return False

for n in range(int(sys.argv[2]) if len(sys.argv) > 2 else 2, int(sys.argv[1]) + 1):
    total = fail = 0
    for _n, A in trees(n):
        adj = [set(a) for a in A]
        total += 1
        if not has_chain(adj, list(range(n)), indpoly_subset(adj, range(n)), {}):
            fail += 1
    print(n, total, "no-order trees:", fail, flush=True)
