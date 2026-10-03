import sys, time
sys.path.insert(0,'scripts')
import vatter_lib_20261003 as V
def centre(adj):
    n=len(adj); deg=[len(a) for a in adj]; layer=[v for v in range(n) if deg[v]<=1]; rem=n; gone=set()
    while rem>2:
        nxt=[]
        for v in layer:
            gone.add(v); rem-=1
            for w in adj[v]:
                if w not in gone:
                    deg[w]-=1
                    if deg[w]==1: nxt.append(w)
        layer=nxt
    return [v for v in range(n) if v not in gone][0]
def orders(adj):
    c=centre(adj); b=V.bfs_order(adj,c)
    deg=[len(a) for a in adj]
    return {"bfs-centre":b, "reverse-bfs":b[::-1],
            "deg-desc":sorted(range(len(adj)),key=lambda v:(-deg[v],v)),
            "deg-asc":sorted(range(len(adj)),key=lambda v:(deg[v],v))}
fams = [
 ("path P_40", V.path(40)), ("caterpillar 12x2", V.caterpillar(12,2)), ("spider S(3^10)", V.spider([3]*10)),
 ("H(9,2) n=28", V.hubstar(9,2)), ("H(20,3) n=81", V.hubstar(20,3)),
 ("MSH(4x9,4x10) n=237", V.msh([9]*4+[10]*4,2)),
]
for name, adj in fams:
    n=len(adj); a=V.alpha(adj); W=V.window(n,a); t=time.time()
    res=[]
    for oname,o in orders(adj).items():
        f=V.check_order(adj,o,W)
        bad=[k for k in W if f[k]]
        res.append(f"{oname}: {'PASS' if not bad else f'fails at {len(bad)}/{len(W)} k'}")
    print(f"{name} n={n} window=[{W.start},{W.stop-1}] | "+" | ".join(res)+f" ({time.time()-t:.1f}s)", flush=True)
