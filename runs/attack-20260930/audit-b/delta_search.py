#!/usr/bin/env python3
"""DIAGNOSTIC (floats): adversarial hill-climb for large q(1-q) delta^2 over spherically
symmetric trees (level h has arity d_h, root level first), lam in K = [1/4,12].
Scalar per-level recursion, evaluated from the leaves up. Reports best value found per size cap.
Not a proof of anything; it tests whether the n^a growth allowed by Prop 4.1 shows up."""
import math, random, json
random.seed(993)
def evalseq(ds, lam):
    # ds[0] = arity of root, ds[-1] = arity of deepest internal level; leaves below
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
res={}
for cap in [10**3, 10**4, 10**5, 10**6, 10**8]:
    best=(0,None,None)
    for trial in range(300):
        lam = math.exp(random.uniform(math.log(0.25), math.log(12)))
        ds = [random.randint(1,4) for _ in range(random.randint(2,30))]
        while size(ds) > cap: ds.pop()
        if not ds: continue
        v = evalseq(ds, lam)
        for it in range(400):
            nd = list(ds); nl = lam
            m = random.random()
            if m < 0.4 and nd: nd[random.randrange(len(nd))] = random.randint(1,12)
            elif m < 0.6: nd.insert(random.randrange(len(nd)+1), random.randint(1,12))
            elif m < 0.75 and len(nd)>1: nd.pop(random.randrange(len(nd)))
            else: nl = min(12, max(0.25, lam*math.exp(random.gauss(0,0.2))))
            if size(nd) > cap: continue
            nv = evalseq(nd, nl)
            if nv >= v: ds, lam, v = nd, nl, nv
        if v > best[0]: best=(v, list(ds), lam)
    res[cap]=dict(best=best[0], arities=best[1], lam=best[2], n=size(best[1]))
    print(f"n<={cap:>10}: best q(1-q)delta^2 = {best[0]:.4f}  n={size(best[1])} lam={best[2]:.4f} arities={best[1][:40]}")
json.dump({str(k):v for k,v in res.items()}, open("data/delta_search.json","w"), indent=1)
