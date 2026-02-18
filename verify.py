import sys
sys.path.insert(0, '.')
from indpoly import independence_poly
import networkx as nx
import random

rng = random.Random(42)
violations = 0
bound1_pass = 0
bound2_pass = 0
total = 0

for _ in range(100):
    n = rng.randint(8, 15)
    G = nx.random_labeled_tree(n, seed=rng.randint(0, 10000))
    adj = [list(G.neighbors(v)) for v in range(n)]
    poly = independence_poly(n, adj)
    
    d = -1
    for i in range(1, len(poly)):
        if poly[i] < poly[i-1]:
            d = i
            break
    if d == -1:
        continue
    
    for u, v in G.edges():
        try:
            G2 = G.copy()
            w = n
            G2.remove_edge(u, v)
            G2.add_edge(u, w)
            G2.add_edge(w, v)
            adj2 = [list(G2.neighbors(v)) for v in range(n+1)]
            poly2 = independence_poly(n+1, adj2)
            
            max_len = max(len(poly), len(poly2))
            poly = poly + [0] * (max_len - len(poly))
            poly2 = poly2 + [0] * (max_len - len(poly2))
            A = [poly2[i] - poly[i] for i in range(max_len)]
            
            for k in range(d+1, min(len(A), len(poly))-1):
                total += 1
                if k > 0 and A[k] <= poly[k-1]:
                    bound1_pass += 1
                if k+1 < len(A) and A[k+1] <= A[k]:
                    bound2_pass += 1
                S = [poly[i] + A[i] for i in range(max_len)]
                for i in range(1, len(S)-1):
                    if S[i-1] > S[i] < S[i+1]:
                        violations += 1
        except:
            pass

print(f'Total: {total}')
print(f'Bound1: {bound1_pass}/{total} = {100*bound1_pass/total:.1f}%')
print(f'Bound2: {bound2_pass}/{total} = {100*bound2_pass/total:.1f}%')
print(f'Violations: {violations}')
