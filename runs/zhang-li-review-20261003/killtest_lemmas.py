"""Exact kill-tests of Zhang-Li (Zenodo 22999166) Lemma 3.1 (terminal segment via matching
ratio domination) and Theorem 6.1 (activity-two mean bounds (47), (48)). Trees only: all three
statements are additive or componentwise for forests."""
import sys
from fractions import Fraction as Fr
from math import comb, floor
sys.path.insert(0,'scripts'); sys.path.insert(0,'.')
import vatter_lib_20261003 as V
from trees import trees
import networkx as nx

def check(adj):
    n=len(adj); I=V.forest_poly(adj,range(n)); a=I.degree(); c=[V.coeff(I,k) for k in range(a+1)]
    v=n-a
    # matching number check: a + nu = n for forests
    G=nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from((i,j) for i in range(n) for j in adj[i] if i<j)
    nu=len(nx.max_weight_matching(G,maxcardinality=True)); assert nu==v,(nu,v)
    # Q = (1+2x)^v (1+x)^(a-v)
    q=[sum(comb(v,i)*2**i*comb(a-v,k-i) for i in range(0,min(k,v)+1)) for k in range(a+1)]
    ratio_ok=all(c[k+1]*q[k] <= c[k]*q[k+1] for k in range(a))
    U=floor((n+2*a)/6)+1
    tail_ok=all(c[k]>=c[k+1] for k in range(U,a))
    Z=sum(ck*2**k for k,ck in enumerate(c)); M=Fr(sum(k*ck*2**k for k,ck in enumerate(c)),Z)
    e47 = M - (Fr(n,3)+Fr(2,15))          # connected: c(G)=1
    e48 = 3*M - (a+Fr(29*n,60))
    return ratio_ok, tail_ok, e47, e48, a, U

N=int(sys.argv[1])
worst47=worst48=None; bad=[]
for n in range(2,N+1):
    for _n,adj in trees(n):
        r,t,e47,e48,a,U=check(adj)
        if not (r and t and e47>=0 and e48>=0): bad.append((n,r,t,float(e47),float(e48)))
        if worst48 is None or e48/n<worst48[0]: worst48=(e48/n,n,a)
        if worst47 is None or e47<worst47[0]: worst47=(e47,n)
    print(n,"bad so far",len(bad),"min e48/n",float(worst48[0]),"at n",worst48[1],"min e47",float(worst47[0]),flush=True)
print("families:")
for name,adj in [("path P_200",V.path(200)),("path P_201",V.path(201)),("spider S(50^4)",V.spider([50]*4)),("spider S(1^60)",V.spider([1]*60)),
                 ("caterpillar 60x1",V.caterpillar(60,1)),("caterpillar 60x2",V.caterpillar(60,2)),("H(9,2)",V.hubstar(9,2)),("H(75,5)",V.hubstar(75,5)),
                 ("MSH(4x9,4x10)",V.msh([9]*4+[10]*4,2)),("spider S(2^40)",V.spider([2]*40)),("spider S(3^30)",V.spider([3]*30))]:
    r,t,e47,e48,a,U=check(adj); n=len(adj)
    print(f"  {name}: n={n} ratio_dom={r} tail_from_U={t} e47={float(e47):.4f} e48/n={float(e48)/n:.5f}")
print("BAD:",bad[:10])
