#!/usr/bin/env python3
"""Exact cross-check (Fractions, brute force over independent sets) that the level recursion in
delta_search.py computes delta = E(X|root in) - E(X|root out) and q for a small spherically
symmetric tree, then reprints the growth of the hill-climb optimum with a larger size range."""
from fractions import Fraction as Fr
import itertools
def build(ds):
    # vertices 0..n-1, root 0, arities per level
    edges=[]; level=[0]; nxt=1
    for d in ds:
        new=[]
        for v in level:
            for _ in range(d):
                edges.append((v,nxt)); new.append(nxt); nxt+=1
        level=new
    return nxt, edges
def brute(ds, lam):
    n, E = build(ds)
    adj=[set() for _ in range(n)]
    for a,b in E: adj[a].add(b); adj[b].add(a)
    Zin=Zout=Fr(0); Xin=Xout=Fr(0)
    for mask in range(1<<n):
        ok=True
        for a,b in E:
            if (mask>>a)&1 and (mask>>b)&1: ok=False; break
        if not ok: continue
        k=bin(mask).count("1"); w=lam**k
        if mask&1: Zin+=w; Xin+=w*k
        else: Zout+=w; Xout+=w*k
    q=Zin/(Zin+Zout); delta=Xin/Zin - Xout/Zout
    return q, delta
def rec(ds, lam):
    r, delta = lam, Fr(1)
    for d in reversed(ds):
        q = r/(1+r)
        r, delta = lam/(1+r)**d, 1 - d*q*delta
    return r/(1+r), delta
for ds, lam in [([3,2],Fr(12)), ([2,1,2],Fr(5,2)), ([4,1],Fr(1,4)), ([2,3,1],Fr(10))]:
    b=brute(ds,lam); r=rec(ds,lam)
    print(ds, lam, "brute==rec:", b==r, float(b[0]), float(b[1]))
