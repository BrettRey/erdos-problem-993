import sys, time, json
sys.path.insert(0,'scripts'); sys.path.insert(0,'.')
import vatter_lib_20261003 as V
from trees import trees
sys.path.insert(0,'/private/tmp/claude-502/-Users-brettreynolds-projects-LLM-CLI-projects-papers-queue-erdos-problem-993/eec3f615-7ddb-4430-8ddb-9f052413b82b/scratchpad')
centre = V.centre
N=int(sys.argv[1])
fails={"revbfs-window":[], "degasc-window":[], "revbfs-allk":[], "degasc-allk":[]}
tot=0
for n in range(int(sys.argv[2]) if len(sys.argv)>2 else 4, N+1):
    t=time.time(); cnt=0
    for _n, adj in trees(n):
        cnt+=1; tot+=1
        I=V.forest_poly(adj,range(n)); a=I.degree(); W=V.window(n,a); K=range(1,a)
        deg=[len(x) for x in adj]
        for oname,o in (("revbfs", V.bfs_order(adj,centre(adj))[::-1]), ("degasc", sorted(range(n),key=lambda v:(deg[v],v)))):
            f=V.check_order(adj,o,K, assert_identity=(n<=10))
            if any(f[k] for k in W): fails[oname+"-window"].append((n,[k for k in W if f[k]]))
            if any(f[k] for k in K): fails[oname+"-allk"].append((n,[k for k in K if f[k]]))
    print(n, cnt, {k:len(v) for k,v in fails.items()}, f"{time.time()-t:.1f}s", flush=True)
print(json.dumps({k:v[:5] for k,v in fails.items()}))
