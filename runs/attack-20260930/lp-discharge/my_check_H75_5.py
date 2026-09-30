"""Fresh parent-side check of the S(2/3) violation at H(75,5), n = 451, built from
the construction (not the adversary's edge list). Polynomial code copied from
../r1_third_check_20260930.py (written independently of the adversary's slib/orbit)."""
from fractions import Fraction
from math import ceil
m, s = 75, 5
n = 1 + m + m * s
adj = [[] for _ in range(n)]; c = 1; supports = []
for _ in range(m):
    sp = c; c += 1; supports.append(sp); adj[0].append(sp); adj[sp].append(0)
    for _ in range(s):
        adj[sp].append(c); adj[c].append(sp); c += 1
assert c == n
def pmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q): r[i + j] += x * y
    return r
def padd(p, q):
    L = max(len(p), len(q)); return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(L)]
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
V = set(range(n)); I = indpoly(V); a = len(I) - 1; q = ceil((2 * a - 1) / 3); lo = ceil(n / 4)
co = lambda p, k: p[k] if 0 <= k < len(p) else 0
Jh = indpoly(V - {0} - set(adj[0])); sp = supports[0]; Js = indpoly(V - {sp} - set(adj[sp]))
print(f'n={n} alpha={a} window=[{lo},{q}]')
for k in range(218, 227):
    D = lambda J: k * co(I, k - 1) * co(J, k) - (k + 1) * co(I, k) * co(J, k - 1)
    Dh, Ds = D(Jh), D(Js)
    L = Fraction(2, 3) * Dh + Fraction(1, 3) * m * Fraction(Ds, len(adj[sp]))   # hub row: 75 support neighbours, each of degree s+1
    lc = co(I, k) ** 2 - co(I, k - 1) * co(I, k + 1)
    print(f'k={k}: L_hub {"> 0 VIOLATION" if L > 0 else "<= 0"}  normalized {float(L / (k * co(I, k - 1) * co(I, k))):+.3e}  LC holds: {lc > 0}')
