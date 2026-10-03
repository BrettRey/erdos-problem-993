"""Run both leaves-first Vatter orders over the 30 Sep adversary's family zoo."""
import sys, time, json
P='/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993'
sys.path.insert(0,P+'/scripts'); sys.path.insert(0,P+'/runs/attack-20260930/lp-discharge/adversary')
import vatter_lib_20261003 as V
import families as F
nmax=int(sys.argv[1]); out=sys.argv[2]
S=F.specs(nmax)
t0=time.time(); stats={"specs":len(S),"checked":0,"fail_window":{"revbfs":0,"degasc":0},"fail_any":{"revbfs":0,"degasc":0}}
fails=[]
with open(out,'w') as fh:
    for i,(fam,par,(ct,args)) in enumerate(S):
        adj=F.CT[ct](*args); n=len(adj)
        I=V.forest_poly(adj,range(n)); a=I.degree(); W=V.window(n,a); K=range(1,a)
        rec={"family":fam,"params":par,"n":n,"alpha":a}
        for oname,fn in (("revbfs",V.order_revbfs),("degasc",V.order_degasc)):
            f=V.check_order(adj,fn(adj),K,assert_identity=(i%50==0))
            fw=[k for k in W if f[k]]; fa=[k for k in K if f[k]]
            rec[oname]={"fail_window":fw,"fail_any":fa}
            if fw: stats["fail_window"][oname]+=1
            if fa: stats["fail_any"][oname]+=1
        stats["checked"]+=1
        fh.write(json.dumps(rec)+"\n")
        if rec["revbfs"]["fail_window"] or rec["degasc"]["fail_window"]:
            fails.append(rec); print("WINDOW FAIL", rec, flush=True)
        if i%200==0: print(i,len(S),stats,f"{time.time()-t0:.0f}s",flush=True)
print("DONE",json.dumps(stats),f"{time.time()-t0:.0f}s")
