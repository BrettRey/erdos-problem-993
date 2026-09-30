"""Repair depth for nested stars: centre -> c hubs -> each hub m sub-hubs -> each
sub-hub s leaves (the R1-killer shape, one level deeper than SoS). Least k with a
k-chain of proper GTS moves reaching central margin <= m(T); iso-dedup per level;
k capped at 5."""
import sys, time, networkx as nx
sys.argv = [sys.argv[0], '0']
exec(open('gts_margin_test.py').read().split("N = int(sys.argv[1])")[0])
def canon(T):
    def enc(v, p):
        return '(' + ''.join(sorted(enc(w, v) for w in T[v] if w != p)) + ')'
    return min(enc(c, None) for c in nx.center(T))
def nested(c, m, s):
    T = nx.Graph(); T.add_node(0); k = 1
    for _ in range(c):
        h = k; T.add_edge(0, h); k += 1
        for _ in range(m):
            g = k; T.add_edge(h, g); k += 1
            for _ in range(s): T.add_edge(g, k); k += 1
    return T
for (c, m, s) in [(2,2,2),(2,3,2),(3,2,2),(2,2,3),(3,3,2),(2,3,3),(3,2,3)]:
    t0 = time.time(); T = nested(c, m, s); n = T.number_of_nodes(); mT = margin(T, 4)
    dist = (n - 1) - sum(1 for v in T if T.degree(v) == 1)
    level = {canon(T): T}; kmin = None; sizes = []
    for k in range(1, min(dist, 5) + 1):
        nxt = {}
        for U in level.values():
            for _, V in proper_gts_images(U): nxt.setdefault(canon(V), V)
        sizes.append(len(nxt))
        if any(margin(V, 4) <= mT for V in nxt.values()): kmin = k; break
        level = nxt
    print(f'nested(c={c},m={m},s={s}) n={n} steps-to-star={dist} least k={kmin if kmin else ">"+str(min(dist,5))} level sizes {sizes} ({time.time()-t0:.0f}s)', flush=True)
