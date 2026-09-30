import sys; sys.path.insert(0, '../lp-discharge')
from forestdp import indpoly
from layered import layered, coeffs
def build(b):
    adj = [[]]; level = [0]; frontier = [0]
    for j, bj in enumerate(b):
        nf = []
        for p in frontier:
            for _ in range(bj):
                c = len(adj); adj.append([p]); adj[p].append(c); level.append(j + 1); nf.append(c)
        frontier = nf
    return adj, level
ok = True
for b in ([3], [4, 2], [3, 2, 2], [2, 3, 1, 2], [3, 1, 2], [2, 2, 2, 2]):
    I, J, N, deg = layered(b); adj, level = build(b); n = len(adj)
    this = coeffs(I) == indpoly(adj, set(range(n))) and sum(N) == n
    for i in range(len(b) + 1):
        v = level.index(i)
        this &= coeffs(J[i]) == indpoly(adj, set(range(n)) - {v} - set(adj[v])) and len(adj[v]) == deg[i]
    ok &= this; print(b, 'match' if this else 'MISMATCH')
print('ALL OK' if ok else 'FAIL')
