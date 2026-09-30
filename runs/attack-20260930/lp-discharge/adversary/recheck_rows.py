"""Independent row-level recheck for large trees (fresh forest DP per vertex, code
adapted from ../../r1_third_check_20260930.py; shares nothing with slib/orbit/r1lib).
Input: {"n", "edges", "rows": [[u, k], ...]}. For each row computes, in exact Fractions,
D_w(k) for w in N[u] (a fresh DP on T - N[w] for each), then
  L_u(sigma) = sigma D_u + (1-sigma) A_u,  A_u = sum_{v~u} D_v/deg v,
and the one-sided bound on sigma the row imposes. Two rows of one tree whose lower
bound exceeds the other's upper bound certify that NO constant sigma works on that tree."""
import json, sys
from fractions import Fraction
from math import ceil
sys.set_int_max_str_digits(0)
rec = json.load(open(sys.argv[1])); edges = rec['edges']; n = rec['n']
adj = [[] for _ in range(n)]
for a, b in edges: adj[a].append(b); adj[b].append(a)
assert len(edges) == n - 1
seen = {0}; st = [0]
while st:
    v = st.pop()
    for w in adj[v]:
        if w not in seen: seen.add(w); st.append(w)
assert len(seen) == n, 'not connected'
def pmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q): r[i + j] += x * y
    return r
def padd(p, q):
    m = max(len(p), len(q)); return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(m)]
def indpoly(alive):
    alive = set(alive); seen = set(); total = [1]
    for r in sorted(alive):
        if r in seen: continue
        order, parent, stack = [], {r: None}, [r]; seen.add(r)
        while stack:
            v = stack.pop(); order.append(v)
            for w in adj[v]:
                if w in alive and w not in seen: seen.add(w); parent[w] = v; stack.append(w)
        out, inn = {}, {}
        for v in reversed(order):
            ex, ic = [1], [0, 1]
            for w in adj[v]:
                if w in alive and parent.get(w) == v: ex = pmul(ex, padd(out[w], inn[w])); ic = pmul(ic, out[w])
            out[v], inn[v] = ex, ic
        total = pmul(total, padd(out[r], inn[r]))
    while len(total) > 1 and total[-1] == 0: total.pop()
    return total
V = set(range(n)); I = indpoly(V); alpha = len(I) - 1
q = ceil(Fraction(2 * alpha - 1, 3)); lo = ceil(Fraction(n, 4)); hi = min(q, alpha - 1)
co = lambda p, k: p[k] if 0 <= k < len(p) else 0
Jc = {}
def J(w):
    if w not in Jc: Jc[w] = indpoly(V - {w} - set(adj[w]))
    return Jc[w]
print(f'n={n} alpha={alpha} window=[{lo},{hi}]', flush=True)
res = []
for u, k in rec['rows']:
    assert lo <= k <= hi, 'row outside window'
    Dw = {w: k * I[k - 1] * co(J(w), k) - (k + 1) * I[k] * co(J(w), k - 1) for w in [u] + adj[u]}
    A = sum(Fraction(Dw[v], len(adj[v])) for v in adj[u]); a = Dw[u] - A
    side = 'upper' if a > 0 else ('lower' if a < 0 else 'none')
    bound = (-A / a) if a != 0 else None
    L23 = Fraction(2, 3) * Dw[u] + Fraction(1, 3) * A
    lc = I[k] ** 2 - I[k - 1] * I[k + 1]
    print(f'row u={u} deg={len(adj[u])} k={k}: sigma {side} bound {float(bound):.9f}; L_u(2/3) {"> 0 (S(2/3) VIOLATED)" if L23 > 0 else "<= 0"}; LC slack i_k^2 - i_(k-1)i_(k+1) > 0: {lc > 0}', flush=True)
    res.append(dict(u=u, k=k, deg=len(adj[u]), side=side, bound=[bound.numerator, bound.denominator] if bound is not None else None,
                    bound_f=float(bound) if bound is not None else None, L23_positive=L23 > 0))
lows = [Fraction(*r['bound']) for r in res if r['side'] == 'lower']; ups = [Fraction(*r['bound']) for r in res if r['side'] == 'upper']
verdict = bool(lows and ups and max(lows) > min(ups))
print('EMPTY INTERVAL CERTIFIED (max lower bound > min upper bound, exact):', verdict)
json.dump(dict(n=n, alpha=alpha, window=[lo, hi], rows=res, empty_certified=verdict), open(sys.argv[2], 'w'))
