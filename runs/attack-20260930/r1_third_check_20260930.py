"""Third, independent check of the R1 counterexample MSH(4x9, 4x10), n = 237.
Written fresh on 2026-09-30 after the pause: shares no code with pv_lib, r1lib or
the analytic agent's closed forms. Reads only the edge list; exact integers/Fractions."""
import json, sys
from fractions import Fraction
from math import ceil

rec = json.load(open(sys.argv[1]))
edges = rec['edge_list_0indexed']
n = 1 + max(max(e) for e in edges)
adj = [[] for _ in range(n)]
for a, b in edges:
    adj[a].append(b); adj[b].append(a)
assert len(edges) == n - 1

def pmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                r[i + j] += x * y
    return r

def padd(p, q):
    m = max(len(p), len(q)); return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(m)]

def indpoly(alive):
    """Independence polynomial of the forest induced on vertex set `alive` (iterative tree DP)."""
    alive = set(alive); seen = set(); total = [1]
    for r in alive:
        if r in seen: continue
        order, parent, stack = [], {r: None}, [r]; seen.add(r)
        while stack:
            v = stack.pop(); order.append(v)
            for w in adj[v]:
                if w in alive and w not in seen:
                    seen.add(w); parent[w] = v; stack.append(w)
        out, inn = {}, {}
        for v in reversed(order):
            ex, ic = [1], [0, 1]
            for w in adj[v]:
                if w in alive and parent.get(w) == v:
                    ex = pmul(ex, padd(out[w], inn[w])); ic = pmul(ic, out[w])
            out[v], inn[v] = ex, ic
        total = pmul(total, padd(out[r], inn[r]))
    while len(total) > 1 and total[-1] == 0: total.pop()
    return total

V = range(n)
I = indpoly(V)
alpha = len(I) - 1
q = ceil((2 * alpha - 1) / 3); lo = ceil(n / 4)
assert I == [int(x) for x in rec['I']], 'polynomial disagrees with the recorded one'
J = {v: indpoly(set(V) - {v} - set(adj[v])) for v in V}
co = lambda p, k: p[k] if 0 <= k < len(p) else 0
# double-count identities as an internal consistency check
for k in range(1, alpha):
    assert sum(co(J[v], k) for v in V) == (k + 1) * I[k + 1]
    assert sum(co(J[v], k - 1) for v in V) == k * I[k]
viol, worst = [], None
for k in range(lo, q + 1):
    D = {v: k * I[k - 1] * co(J[v], k) - (k + 1) * I[k] * co(J[v], k - 1) for v in V}
    assert sum(D.values()) == k * (k + 1) * (I[k - 1] * I[k + 1] - I[k] ** 2)
    lc = I[k] ** 2 - I[k - 1] * I[k + 1]
    for u in V:
        L = sum(Fraction(D[v], len(adj[v]) + 1) for v in [u] + adj[u])
        norm = float(L / (k * I[k - 1] * I[k]))
        if worst is None or norm > worst[0]: worst = (norm, u, k)
        if L > 0: viol.append((u, len(adj[u]), k, str(L), lc > 0))
print(f'n={n} alpha={alpha} window=[{lo},{q}]  polynomial matches record')
print('R1 violations (u, deg u, k, L_u, LC holds at k):')
for x in viol: print('  ', x[0], x[1], x[2], x[3][:40] + '...', x[4])
print('max normalized L over window:', worst)
