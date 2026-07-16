#!/usr/bin/env python3
"""Depth-3 recursive gadget valley sweep.

Architecture: root joined to A copies of a super-gadget plus B plain
spiders S(2^t). A super-gadget is a center carrying s sub-spiders
S(2^u) (its own centers carry legs). Three nested inclusion phases:
root / super-centers / sub-centers. Tests whether phase stacking can
beat the two-phase deficit barrier 1-V ~ C/n, C in [4,6].

Exact arithmetic (Kronecker-packed convolutions).
"""

import sys
import time

sys.path.insert(0, ".")

from scripts.valley_scaling_probe import kmul, kpow, kadd, shift, path_poly
from scripts.valley_search import valley_score


def spider_EI(u: int, cnt: int):
    """(E, I) for center with cnt legs of length u."""
    E = kpow(path_poly(u), cnt)
    I = shift(kpow(path_poly(u - 1), cnt))
    return E, I


def super_gadget_EI(s: int, u: int, own_legs: int = 0):
    """Center carrying s sub-spiders S(2^u) and own_legs legs of length 1."""
    subE, subI = spider_EI(u, 1)
    # s copies of sub-spider S(u^1)? No: sub-spider is S(2^u) = center with
    # u legs of length 2.
    raise NotImplementedError


def spider2_EI(u: int):
    """(E, I) for the sub-spider S(2^u): center with u legs of length 2."""
    E = kpow(path_poly(2), u)
    I = shift(kpow(path_poly(1), u))
    return E, I


def superg_EI(s: int, u: int, legs1: int = 0):
    """(E, I) for super-center with s sub-spiders S(2^u) + legs1 leaves."""
    subE, subI = spider2_EI(u)
    both = kadd(subE, subI)
    E = kpow(both, s)
    I = shift(kpow(subE, s))
    if legs1:
        E = kmul(E, kpow(path_poly(1), legs1))
        # leaves excluded when center included: factor 1 each
    return E, I


def tree_poly(A: int, s: int, u: int, B: int, t: int):
    """Root + A super-gadgets(s,u) + B spiders S(2^t)."""
    ex = [1]
    inc = [1]
    if A:
        E, I = superg_EI(s, u)
        ex = kmul(ex, kpow(kadd(E, I), A))
        inc = kmul(inc, kpow(E, A))
    if B:
        E, I = spider2_EI(t)
        ex = kmul(ex, kpow(kadd(E, I), B))
        inc = kmul(inc, kpow(E, B))
    return kadd(ex, shift(inc))


def size(A, s, u, B, t):
    return 1 + A * (1 + s * (1 + 2 * u)) + B * (1 + 2 * t)


def selftest():
    from indpoly import independence_poly

    def build_adj(A, s, u, B, t):
        adj = [[]]

        def add(parent):
            adj.append([])
            v = len(adj) - 1
            adj[parent].append(v)
            adj[v].append(parent)
            return v

        for _ in range(A):
            sc = add(0)
            for _ in range(s):
                c = add(sc)
                for _ in range(u):
                    x = add(c)
                    add(x)
        for _ in range(B):
            c = add(0)
            for _ in range(t):
                x = add(c)
                add(x)
        return len(adj), adj

    for args in ((1, 2, 2, 1, 3), (2, 3, 2, 0, 2), (0, 2, 2, 3, 4),
                 (2, 2, 4, 2, 2)):
        n, adj = build_adj(*args)
        assert n == size(*args), (args, n, size(*args))
        assert tree_poly(*args) == independence_poly(n, adj), args
    print("depth3 selftest passed", flush=True)


def main():
    selftest()
    best = []
    print(f"{'(A,s,u,B,t)':>22} {'n':>6} {'alpha':>6} {'V':>14} {'n(1-V)':>9} "
          f"{'b/alpha':>8} {'win':>5}")
    t0 = time.time()
    count = 0
    for u in (2, 3, 4, 6):
        for s in (2, 3, 4, 6, 9, 14):
            for A in (1, 2, 3, 5, 8, 13, 21, 34):
                for t in (2, 4, 8, 12):
                    for B in (0, 1, 2, 5, 13, 34, 89):
                        n = size(A, s, u, B, t)
                        if not (120 <= n <= 2600):
                            continue
                        poly = tree_poly(A, s, u, B, t)
                        vs = valley_score(poly)
                        count += 1
                        alpha = len(poly) - 1
                        rec = ((vs["window"], vs["ratio"]),
                               (A, s, u, B, t), n, alpha, vs)
                        best.append(rec)
                        if vs["witness"]:
                            print(f"*** WITNESS *** {(A,s,u,B,t)} n={n}")
                            print("poly=", poly, flush=True)
                            return
    best.sort(reverse=True)
    for key, args, n, alpha, vs in best[:25]:
        print(f"{str(args):>22} {n:>6} {alpha:>6} {vs['ratio']:>14.10f} "
              f"{n*(1-vs['ratio']):>9.4f} {vs['pos']/alpha:>8.4f} "
              f"{str(vs['window']):>5s}")
    print(f"\n{count} configs in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
