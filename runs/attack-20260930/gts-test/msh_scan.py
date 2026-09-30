"""Repair depth across the R1-killer family MSH(c; m, s): a centre joined to c
hub-stars H(m, s) (hub -> m middle vertices -> s leaves each). Least k (<= 7)."""
import sys, time, networkx as nx
sys.argv = [sys.argv[0], '0']
exec(open('gts_margin_test.py').read().split("N = int(sys.argv[1])")[0])
sys.setrecursionlimit(10000)
def canon(T):
    def enc(v, p):
        return '(' + ''.join(sorted(enc(w, v) for w in T[v] if w != p)) + ')'
    return min(enc(c, None) for c in nx.center(T))
def msh(c, m, s):
    T = nx.Graph(); T.add_node(0); k = 1
    for _ in range(c):
        h = k; T.add_edge(0, h); k += 1
        for _ in range(m):
            g = k; T.add_edge(h, g); k += 1
            for _ in range(s): T.add_edge(g, k); k += 1
    return T
for (c, m, s) in [(2,9,2),(4,9,2),(6,9,2),(8,9,2),(10,9,2),(8,5,2),(8,13,2)]:
    t0 = time.time(); T = msh(c, m, s); n = T.number_of_nodes(); mT = margin(T, 4)
    level = {canon(T): T}; kmin = None; sizes = []
    for k in range(1, 8):
        nxt = {}
        for U in level.values():
            for _, V in proper_gts_images(U): nxt.setdefault(canon(V), V)
        sizes.append(len(nxt))
        if any(margin(V, 4) <= mT for V in nxt.values()): kmin = k; break
        level = nxt
        if time.time() - t0 > 600: break
    print(f'MSH(c={c};m={m},s={s}) n={n} m/m_star={float(mT/margin(nx.star_graph(n-1),4)):.3f} least k={kmin if kmin else ">"+str(len(sizes))} level sizes {sizes} ({time.time()-t0:.0f}s)', flush=True)
