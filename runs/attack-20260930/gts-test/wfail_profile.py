"""Profile the weak-form failures: trees whose every proper GTS image has a
larger minimum central margin. How close are they to the star's margin?"""
import sys, json, networkx as nx
from fractions import Fraction
sys.argv = [sys.argv[0], '0']
exec(open('gts_margin_test.py').read().split("N = int(sys.argv[1])")[0])
N = 16
for n in range(7, N + 1):
    trees = list(nx.nonisomorphic_trees(n)); star = nx.star_graph(n - 1); ms = margin(star, 4)
    fails = []; allmins = []
    for T in trees:
        if max(d for _, d in T.degree()) == n - 1: continue
        mT = margin(T, 4); allmins.append(mT / ms)
        imgs = [margin(U, 4) for _, U in proper_gts_images(T)]
        if all(mU > mT for mU in imgs):
            degs = sorted((d for _, d in T.degree()), reverse=True)
            a = len(indpoly(T)) - 1
            fails.append((float(mT / ms), degs[:4], a, nx.diameter(T)))
    fails.sort()
    print(f'n={n} W-fail={len(fails)} min m/m_star over W-fails={fails[0][0] if fails else None:.4} | min m/m_star over all non-stars={float(min(allmins)):.4} | closest W-fails: {fails[:3]}', flush=True)
