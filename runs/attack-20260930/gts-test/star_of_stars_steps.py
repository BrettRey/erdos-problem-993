"""Star of stars SoS(m,s): a centre joined to m hubs, each hub with s leaves.
For each, find the least k such that some chain of k proper GTS moves reaches a
tree with central margin <= m(SoS). The star is (number of leaves to gain)
steps away; if the least k equals that distance, the only repair is the whole
chain to the star (circular for a local lemma). Exact Fractions; isomorphism
classes deduplicated by an AHU canonical string rooted at the tree centre(s)."""
import sys, networkx as nx
sys.argv = [sys.argv[0], '0']
exec(open('gts_margin_test.py').read().split("N = int(sys.argv[1])")[0])

def canon(T):
    def enc(v, p):
        return '(' + ''.join(sorted(enc(w, v) for w in T[v] if w != p)) + ')'
    return min(enc(c, None) for c in nx.center(T))

def sos(m, s):
    T = nx.Graph(); T.add_node(0); c = 1
    for _ in range(m):
        h = c; T.add_edge(0, h); c += 1
        for _ in range(s): T.add_edge(h, c); c += 1
    return T

for m in range(3, 7):
    for s in range(2, 5):
        T = sos(m, s); n = T.number_of_nodes(); mT = margin(T, 4)
        dist = (n - 1) - sum(1 for v in T if T.degree(v) == 1)
        level = {canon(T): T}; kmin = None
        for k in range(1, min(dist, 4) + 1):
            nxt = {}
            for U in level.values():
                for _, V in proper_gts_images(U):
                    nxt.setdefault(canon(V), V)
            if any(margin(V, 4) <= mT for V in nxt.values()): kmin = k; break
            level = nxt
        print(f'SoS(m={m},s={s}) n={n} m/m_star={float(mT/margin(nx.star_graph(n-1),4)):.4f} steps to star={dist} least repairing k={kmin if kmin else "> " + str(min(dist,4))} (level size {len(level)})', flush=True)
