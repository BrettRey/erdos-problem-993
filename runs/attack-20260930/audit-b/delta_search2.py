#!/usr/bin/env python3
"""DIAGNOSTIC (floats): larger hill-climb for sup q(1-q)delta^2 over spherically symmetric trees
with root a centroid (root arity >= 2 forces this), lam in K.  Fits the local log-log slope,
a lower-bound proxy for the exponent a in Prop 4.1 (true sup over all trees is >= what is found)."""
import math, random, json
random.seed(99301)
def evalseq(ds, lam):
    r, delta = lam, 1.0
    for d in reversed(ds):
        q = r/(1+r)
        r, delta = lam/(1+r)**d, 1 - d*q*delta
    qq = r/(1+r)
    return qq*(1-qq)*delta*delta
def size(ds):
    n, level = 1, 1
    for d in ds:
        level *= d; n += level
    return n
res=[]
for e in range(3, 16):
    cap=10**e
    best=(0,None,None)
    for trial in range(200):
        lam = math.exp(random.uniform(math.log(0.25), math.log(12)))
        ds = [random.randint(2,12)] + [random.randint(1,12) for _ in range(random.randint(1,40))]
        while size(ds) > cap and len(ds)>1: ds.pop()
        if size(ds)>cap: continue
        v = evalseq(ds, lam)
        for it in range(1500):
            nd = list(ds); nl = lam
            m = random.random()
            if m < 0.4: nd[random.randrange(len(nd))] = random.randint(1,14)
            elif m < 0.6: nd.insert(random.randrange(1,len(nd)+1), random.randint(1,14))
            elif m < 0.75 and len(nd)>1: nd.pop(random.randrange(1,len(nd)))
            else: nl = min(12, max(0.25, lam*math.exp(random.gauss(0,0.1))))
            if nd[0] < 2 or size(nd) > cap: continue
            nv = evalseq(nd, nl)
            if nv >= v: ds, lam, v = nd, nl, nv
        if v > best[0]: best=(v, list(ds), lam)
    n=size(best[1]); res.append((e, n, best[0], best[2], best[1]))
    print(f"cap 1e{e:2d}: best={best[0]:.4g} n={n:.3e} lam={best[2]:.3f} arities={best[1]}")
slopes=[]
for i in range(1,len(res)):
    s=(math.log(res[i][2])-math.log(res[i-1][2]))/(math.log(res[i][1])-math.log(res[i-1][1]))
    slopes.append(s)
print("local slopes:", [round(s,3) for s in slopes])
x=[math.log(r[1]) for r in res[2:]]; y=[math.log(r[2]) for r in res[2:]]
mx=sum(x)/len(x); my=sum(y)/len(y)
fit=sum((a-mx)*(b-my) for a,b in zip(x,y))/sum((a-mx)**2 for a in x)
print("LS slope (caps 1e5..1e15):", round(fit,4))
json.dump(dict(rows=[dict(cap_exp=r[0],n=r[1],best=r[2],lam=r[3],arities=r[4]) for r in res],
               local_slopes=slopes, ls_slope=fit), open("data/delta_search2.json","w"), indent=1)
