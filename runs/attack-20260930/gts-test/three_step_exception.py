"""The one n=13 tree that two-step GTS moves cannot move without raising the
central margin: identify it and test three steps."""
import sys, networkx as nx
sys.argv = [sys.argv[0], '0']
exec(open('gts_margin_test.py').read().split("N = int(sys.argv[1])")[0])
T = nx.from_graph6_bytes(b'Ls`?GGCO??_A?C')
ms = margin(nx.star_graph(12), 4); mT = margin(T, 4)
degs = sorted((d for _, d in T.degree()), reverse=True)
hubs = [v for v in T if T.degree(v) >= 3]
print('degrees', degs, 'alpha', len(indpoly(T)) - 1, 'diameter', nx.diameter(T), 'hubs', [(h, T.degree(h)) for h in hubs], 'hub distance', nx.shortest_path_length(T, hubs[0], hubs[1]) if len(hubs) == 2 else None)
print('m(T)/m(star) =', float(mT / ms))
lvl1 = [U for _, U in proper_gts_images(T)]
print('one-step images: m/m_star =', sorted(round(float(margin(U, 4) / ms), 4) for U in lvl1))
found = None
for U in lvl1:
    for _, V in proper_gts_images(U):
        for _, W in proper_gts_images(V):
            if margin(W, 4) <= mT: found = W; break
        if found: break
    if found: break
print('three-step repair exists:', found is not None, ('-> degrees ' + str(sorted((d for _, d in found.degree()), reverse=True))) if found else '')
