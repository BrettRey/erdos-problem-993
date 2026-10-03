"""Simulated annealing against the leaves-first Vatter certificate on the window.
Mutation operator and seed families reused from the 30 Sep adversary (sa_r1.mutate, sa_const seeds).
Objective (float, steering only): max over tails v and window k of rho_k(tail)/rho_k(T) - 1.
Positive means a violation; every reported violation is rechecked exactly.
Usage: sa_vatter.py ORDER NMIN NMAX SECONDS SEED OUT.json"""
import sys, os, json, math, random, time
from fractions import Fraction
P='/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993'
sys.path.insert(0,P+'/scripts'); sys.path.insert(0,P+'/runs/attack-20260930/r1-adversarial'); sys.path.insert(0,P+'/runs/attack-20260930/lp-discharge/adversary')
import vatter_lib_20261003 as V
from sa_r1 import mutate, edges_of
order_name, nmin, nmax, secs, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
ORDER = V.order_revbfs if order_name=='revbfs' else V.order_degasc
def score(adj):
    n=len(adj); I=V.forest_poly(adj,range(n)); a=I.degree(); W=V.window(n,a)
    tails=V.order_tails(adj,ORDER(adj)); best=-1e9; arg=None
    for k in W:
        b1,b0=V.coeff(I,k),V.coeff(I,k-1)
        for i,t in enumerate(tails):
            a1,a0=V.coeff(t,k),V.coeff(t,k-1)
            if a0==0: continue
            r=float(Fraction(a1*b0, b1*a0))-1.0
            if r>best: best,arg=r,(i,k)
    return best,arg,n,a
def hstar(m,s):
    adj=V.hubstar(m,s); return adj
SEEDS=[]
for s in (2,3,4,5):
    for m in range(8,60,4):
        if nmin<=1+m*(s+1)<=nmax: SEEDS.append(V.hubstar(m,s))
for s in (1,2,3):
    for m in (4,6,8,10,12):
        for h in (2,4,6,8,12):
            if nmin<=1+h*(1+m*(s+1))<=nmax: SEEDS.append(V.msh([m]*h,s))
rng=random.Random(seed); t0=time.time()
scored=sorted(((score(a)[0],i) for i,a in enumerate(SEEDS)),reverse=True)
top=[SEEDS[i] for _,i in scored[:8]]
best=None; viols=[]; evals=0; restarts=0
while time.time()-t0<secs:
    restarts+=1; adj=top[restarts%len(top)]; cur,arg,n,a=score(adj); T0=0.02*(abs(cur)+1e-4)
    for step in range(400):
        if time.time()-t0>=secs: break
        T=T0*(0.01)**(step/400)
        try: adj2=mutate(adj,rng,nmin,nmax)
        except (IndexError,ValueError): continue
        if not (nmin<=len(adj2)<=nmax): continue
        v2,arg2,n2,a2=score(adj2); evals+=1
        if v2>0:
            # exact recheck
            n=len(adj2); W=V.window(n,a2); f=V.check_order(adj2,ORDER(adj2),W)
            bad=[k for k in W if f[k]]
            if bad and len(viols)<10: viols.append(dict(n=n,alpha=a2,ks=bad[:5],edges=edges_of(adj2))); print("EXACT VIOLATION",n,a2,bad[:5],flush=True)
        if v2>=cur or rng.random()<math.exp((v2-cur)/T):
            adj,cur=adj2,v2
            if best is None or cur>best['val']:
                best=dict(val=cur,n=n2,alpha=a2,at=arg2,edges=edges_of(adj),degseq=sorted((len(x) for x in adj),reverse=True)[:6])
json.dump(dict(order=order_name,nmin=nmin,nmax=nmax,seed=seed,elapsed=time.time()-t0,evals=evals,restarts=restarts,seed_best=scored[0][0],best=best,violations=viols),open(out,'w'))
print("DONE",order_name,"evals",evals,"best slack",best and best['val'],"violations",len(viols))
