"""Exact max-flow test on the killer trees: can every positive defect be cancelled
by negative defect within meeting distance 2 (i.e. a tree-specific radius-1
certificate)? Positive v -> negative w if dist(v,w) <= 2; capacities D_v and -D_w
(exact integers). Feasible iff max flow = total positive defect."""
import networkx as nx
from math import ceil
from layered import layered, coeffs
co = lambda p, k: p[k] if 0 <= k < len(p) else 0
def build(b):
    adj = [[]]; level = [0]; frontier = [0]
    for j, bj in enumerate(b):
        nf = []
        for p in frontier:
            for _ in range(bj):
                c = len(adj); adj.append([p]); adj[p].append(c); level.append(j + 1); nf.append(c)
        frontier = nf
    return adj, level
for name, b in [('H(9,2)', [9, 2]), ('H(75,5)', [75, 5]), ('MSH(8;10,2) [R1 killer family]', [8, 10, 2]),
                ('MSH(38;11,2)', [38, 11, 2]), ('MSH(96;11,2)', [96, 11, 2]), ('depth-6 worst', [2, 1, 4, 3, 6, 2])]:
    I, J, N, deg = layered(b); I = coeffs(I); Jc = [coeffs(x) for x in J]
    adj, level = build(b); n = len(adj); a = len(I) - 1; q = ceil((2 * a - 1) / 3); lo = ceil(n / 4)
    G = nx.Graph([(u, v) for u in range(n) for v in adj[u] if u < v])
    tested = infeasible = 0; minslack = None
    for k in range(lo, min(q, a - 1) + 1):
        D = [k * co(I, k - 1) * co(Jc[i], k) - (k + 1) * co(I, k) * co(Jc[i], k - 1) for i in range(len(b) + 1)]
        if max(D) <= 0: continue
        tested += 1
        F = nx.DiGraph(); pos_total = 0
        for v in range(n):
            d = D[level[v]]
            if d > 0:
                F.add_edge('s', v, capacity=d); pos_total += d
                for w, dist in nx.single_source_shortest_path_length(G, v, cutoff=2).items():
                    if D[level[w]] < 0: F.add_edge(v, ('n', w), capacity=10**400)
            elif d < 0:
                F.add_edge(('n', v), 't', capacity=-d)
        flow = nx.maximum_flow_value(F, 's', 't')
        if flow < pos_total: infeasible += 1
        sl = sum(-D[level[w]] for w in range(n) if D[level[w]] < 0) / pos_total
        minslack = sl if minslack is None else min(minslack, sl)
    print(f'{name}: n={n}, window levels with positive defect={tested}, radius-1 infeasible levels={infeasible}, min (total negative / total positive)={minslack:.1f}' if tested else f'{name}: n={n}, no positive defect in window', flush=True)
