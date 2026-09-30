"""Repair depth (least k <= 3) for the R1 counterexample MSH(4x9, 4x10), n=237."""
import sys, time, json, networkx as nx
_K = sys.argv[1] if len(sys.argv) > 1 else "9"; sys.argv = [sys.argv[0], "0"]
exec(open('gts_margin_test.py').read().split("N = int(sys.argv[1])")[0])
def canon(T):
    def enc(v, p):
        return '(' + ''.join(sorted(enc(w, v) for w in T[v] if w != p)) + ')'
    return min(enc(c, None) for c in nx.center(T))
sys.setrecursionlimit(10000)
d = json.load(open('../r1-analytic/data/recheck_MSH_4x9_4x10_generic.json'))
T = nx.Graph(d['edge_list_0indexed']); n = T.number_of_nodes(); t0 = time.time()
mT = margin(T, 4); ms = margin(nx.star_graph(n - 1), 4)
print(f'n={n} m(T)/m(star)={float(mT/ms):.4f} ({time.time()-t0:.0f}s)', flush=True)
level = {canon(T): T}
for k in range(1, int(_K) + 1):
    nxt = {}
    for U in level.values():
        for _, V in proper_gts_images(U): nxt.setdefault(canon(V), V)
    ok = [V for V in nxt.values() if margin(V, 4) <= mT]
    print(f'k={k}: {len(nxt)} iso classes, {len(ok)} with margin <= m(T) ({time.time()-t0:.0f}s)', flush=True)
    if ok: break
    level = nxt
