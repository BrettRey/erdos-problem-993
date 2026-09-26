"""Census of every spider (>=3 legs) and path in the depth-three window (33<=n<=38,
alpha in {17,18,19}, 2 alpha - n <= 5), using the packet's exact analyzer, to test whether
the Primus rerun's spider exclusion removes any tree from the residual regime.
Run from gpt_attack/primus_depth3_2026-09-23/. Result 2026-09-26: 24,343 window spiders,
0 in the residual regime; the six paths P33..P38 also not in it."""
import sys, time, networkx as nx
sys.path.insert(0,'.')
from replay import analyze
def partitions(n, k_min=3, maxpart=None):
    # partitions of n into parts (non-increasing), at least k_min parts
    def rec(rem, mx):
        if rem==0: yield []; return
        for p in range(min(rem,mx),0,-1):
            for rest in rec(rem-p,p): yield [p]+rest
    for P in rec(n, maxpart or n):
        if len(P)>=k_min: yield P
def spider(legs):
    T=nx.Graph(); T.add_node(0); c=1
    for L in legs:
        prev=0
        for _ in range(L): T.add_edge(prev,c); prev=c; c+=1
    return T
import networkx as _nx
t0=time.time(); tot=0; inwin=0; b1pos=0; R=[]; negD=0
for n in range(33,39):
    for legs in partitions(n-1):
        q=sum(L%2 for L in legs); delta=q-1 if q>=1 else 1
        alpha=(n+delta)//2
        if not (17<=alpha<=19 and delta<=5 and 2*alpha-n==delta): continue
        r=analyze(spider(legs)); tot+=1
        assert r['alpha']==alpha and r['in_target_window']
        if r['b_0_to_4'][1]>0: b1pos+=1
        if r['combined_correction']<0: negD+=1
        if r['in_residual_regime_R']: R.append((n,tuple(legs),r['b_0_to_4'],r['combined_correction']))
print('window spiders:',tot,'b1>0:',b1pos,'D<0:',negD,'in residual regime:',len(R),'time %.0fs'%(time.time()-t0))
for x in R[:10]: print(' ',x)
for n in range(33,39):
    r=analyze(_nx.path_graph(n)); print("P%d"%n, "in residual regime" if r["in_residual_regime_R"] else "not in residual regime")
