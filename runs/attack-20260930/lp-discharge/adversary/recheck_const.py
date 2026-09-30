"""Independent recheck of the constant rule S(sigma) on an explicit tree.
Adapted from ../../r1_third_check_20260930.py (fresh forest DP on every T - N[v];
shares no code with slib, orbit, r1lib or pv_lib). Reads {"n":..,"edges":[[a,b],..]}.
Exact Fractions: reports every (u,k) with L_u(k) > 0 at sigma = 2/3 on the window
W = [ceil(n/4), min(q, alpha-1)], and the exact feasible sigma-interval."""
import json, sys
from fractions import Fraction
from math import ceil
rec = json.load(open(sys.argv[1])); edges = rec['edges']; n = rec['n']
adj = [[] for _ in range(n)]
for a, b in edges: adj[a].append(b); adj[b].append(a)
assert len(edges) == n - 1
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
    for r in alive:
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
V = range(n); I = indpoly(V); alpha = len(I) - 1
q = ceil(Fraction(2 * alpha - 1, 3)); lo = ceil(Fraction(n, 4)); hi = min(q, alpha - 1)
J = {v: indpoly(set(V) - {v} - set(adj[v])) for v in V}
co = lambda p, k: p[k] if 0 <= k < len(p) else 0
for k in range(1, alpha):
    assert sum(co(J[v], k - 1) for v in V) == k * I[k]
sig = Fraction(2, 3); viol = []; low = None; up = None; worst = None
for k in range(lo, hi + 1):
    D = {v: k * I[k - 1] * co(J[v], k) - (k + 1) * I[k] * co(J[v], k - 1) for v in V}
    assert sum(D.values()) == k * (k + 1) * (I[k - 1] * I[k + 1] - I[k] ** 2)
    for u in V:
        A = sum(Fraction(D[v], len(adj[v])) for v in adj[u])
        L = sig * D[u] + (1 - sig) * A
        nm = L / (k * I[k - 1] * I[k])
        if worst is None or nm > worst[0]: worst = (nm, u, k)
        if L > 0: viol.append((u, len(adj[u]), k, float(nm)))
        a = D[u] - A
        if a > 0:
            x = -A / a; up = (x, u, k) if up is None or x < up[0] else up
        elif a < 0:
            x = A / (-a); low = (x, u, k) if low is None or x > low[0] else low
print(f'n={n} alpha={alpha} window=[{lo},{hi}] I[:4]={I[:4]}')
print(f'S(2/3) violations: {len(viol)}'); [print('   u=%d deg=%d k=%d  L/(k i_{k-1} i_k)=%.3e' % x) for x in viol[:12]]
print('max normalised L at 2/3: %.6e at u=%d k=%d' % (float(worst[0]), worst[1], worst[2]))
print('interval: low=%s (%.6f) at u=%d k=%d ; up=%s (%.6f) at u=%d k=%d' % (str(low[0])[:60], float(low[0]), low[1], low[2], str(up[0])[:60], float(up[0]), up[1], up[2]))
json.dump(dict(n=n, alpha=alpha, window=[lo, hi], n_viol=len(viol), viol=viol[:50], low=[low[0].numerator, low[0].denominator, low[1], low[2]],
               up=[up[0].numerator, up[0].denominator, up[1], up[2]], I=I), open(sys.argv[2], 'w'))
