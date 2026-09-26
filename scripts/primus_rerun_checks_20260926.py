"""Exact checks of the computational claims in the 2026-09-26 Primus rerun report (paper-2):
path margins P33..P38, the three order-8 residual-regime trees, Lemma 6 by brute force
(all trees n<=12), and the spider deficiency formula (Prop. 10). See
gpt_attack/primus_rerun_audit_2026-09-26.md."""
import itertools, networkx as nx
from math import comb
def indsets(G):
    V=list(G); out=[]
    for r in range(len(V)+1):
        for S in itertools.combinations(V,r):
            if all(not G.has_edge(a,b) for a,b in itertools.combinations(S,2)): out.append(frozenset(S))
    return out
def profile(G):
    I=indsets(G); a=max(len(S) for S in I); M=[S for S in I if len(S)==a]
    e=[0]*5; b=[0]*5
    for S in I:
        d=a-len(S)
        if d<=4:
            (e if any(S<=m for m in M) else b)[d]+=1
    return a,e,b,I
def D(e,b): return b[3]**2-b[2]*b[4]+2*e[3]*b[3]-e[2]*b[4]-e[4]*b[2]
# (a) path margins, depth 3, n=33..38 via i_k = C(n-k+1,k)
for n in range(33,39):
    i=[comb(n-k+1,k) for k in range(n+1)]; i=[x for x in i if x>0]; a=len(i)-1
    s=lambda d:i[a-d]
    print('P%d alpha=%d delta=%d margin=%d'%(n,a,2*a-n,s(3)**2-s(2)*s(4)))
# (b) order-8 trees in the residual regime (alpha=4)
hits=[]
for T in nx.nonisomorphic_trees(8):
    a,e,b,_=profile(T)
    if b[1]>0 and 3*(a-3)*b[4]>(a-7)*e[4] and D(e,b)<0:
        hits.append((a,2*a-8,tuple(b[1:]),tuple(e[1:]),D(e,b)))
print('order-8 R trees:',len(hits)); [print(' ',h) for h in hits]
# (c) Lemma 6 brute force: E(S)>=2d => extendable, all trees n<=12
bad=0; checked=0
for n in range(2,13):
    for T in nx.nonisomorphic_trees(n):
        a,_,_,I=profile(T); M=[S for S in I if len(S)==a]
        for S in I:
            d=a-len(S)
            if d>=1:
                N=set(S)|{u for v in S for u in T[v]}
                if n-len(N)>=2*d:
                    checked+=1
                    if not any(S<=m for m in M): bad+=1
print('Lemma 6 checks',checked,'violations',bad)
# (d) Prop 10 spider deficiency, all leg multisets with 3..5 legs, lengths 1..7
viol=0
for k in range(3,6):
    for legs in itertools.combinations_with_replacement(range(1,8),k):
        T=nx.Graph(); T.add_node(0); c=1
        for L in legs:
            prev=0
            for _ in range(L): T.add_edge(prev,c); prev=c; c+=1
        n=T.number_of_nodes(); a=len(nx.max_weight_matching(nx.complement(T)) ) if False else None
        # alpha of a tree = n - max matching size
        a=n-len(nx.max_weight_matching(T,maxcardinality=True))
        q=sum(L%2 for L in legs); pred=q-1 if q>=1 else 1
        if 2*a-n!=pred: viol+=1
print('Prop 10 violations',viol)
