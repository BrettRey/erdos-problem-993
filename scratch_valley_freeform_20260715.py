#!/usr/bin/env python3
"""Free-form evolutionary valley search over arbitrary trees.

Escapes the spider-bouquet grammar: starts from the 2026-07-15 valley
campaign champions and applies unrestricted local tree mutations,
accepting improvements in (window, valley_ratio). Exact arithmetic.
"""

import random
import sys
import time

sys.path.insert(0, ".")

from indpoly import independence_poly
from scripts.valley_search import bouquet_adj, valley_score
from scripts.product_valley_search import dumbbell_adj


def mutate_tree(n, adj, rng):
    """Random SPR / leaf move / subdivide / contract on a copy."""
    adj2 = [list(x) for x in adj]
    op = rng.randrange(4)
    if op == 0:  # move a random leaf to a random vertex
        leaves = [v for v in range(n) if len(adj2[v]) == 1]
        v = rng.choice(leaves)
        p = adj2[v][0]
        w = rng.randrange(n)
        if w == v or w == p:
            return None
        adj2[p].remove(v)
        adj2[v] = [w]
        adj2[w].append(v)
        return n, adj2
    if op == 1:  # subtree prune and regraft (small subtree)
        u = rng.randrange(n)
        if len(adj2[u]) < 2:
            return None
        c = rng.choice(adj2[u])
        # collect subtree of c away from u; cap size
        comp = {c}
        stack = [c]
        while stack:
            x = stack.pop()
            for y in adj2[x]:
                if y != u and y not in comp:
                    comp.add(y)
                    stack.append(y)
            if len(comp) > n // 3:
                return None
        target = rng.randrange(n)
        if target in comp or target == u:
            return None
        adj2[u].remove(c)
        adj2[c].remove(u)
        adj2[c].append(target)
        adj2[target].append(c)
        return n, adj2
    if op == 2:  # subdivide a random edge (n grows by 1)
        u = rng.randrange(n)
        if not adj2[u]:
            return None
        v = rng.choice(adj2[u])
        w = n
        adj2.append([u, v])
        adj2[u][adj2[u].index(v)] = w
        adj2[v][adj2[v].index(u)] = w
        return n + 1, adj2
    if op == 3:  # contract a random degree-2 vertex (n shrinks by 1)
        deg2 = [v for v in range(n) if len(adj2[v]) == 2]
        if not deg2:
            return None
        v = rng.choice(deg2)
        a, b = adj2[v]
        if b in adj2[a]:
            return None
        adj2[a][adj2[a].index(v)] = b
        adj2[b][adj2[b].index(v)] = a
        # remove v by swapping with last vertex
        last = n - 1
        if v != last:
            for w2 in adj2[last]:
                adj2[w2][adj2[w2].index(last)] = v
            adj2[v] = adj2[last]
        adj2.pop()
        return n - 1, adj2
    return None


def seeds():
    out = []
    # champion bouquet: S(2^3) + S(2^77) + S(2^75 + 3^3)
    g = ((2,) * 3, (2,) * 77, tuple([2] * 75 + [3] * 3))
    out.append(("champ_bouquet", bouquet_adj(g)))
    # champion dumbbell: 12xS(3^6) --2-- 12xS(3^6)
    sA = (tuple([(3,) * 6] * 12), (), 0)
    out.append(("champ_dumbbell", dumbbell_adj(sA, sA, 2)))
    # hybrid: 34xS(2^2)+7xS(2^10)
    g = tuple([(2,) * 2] * 34 + [(2,) * 10] * 7)
    out.append(("champ_hybrid", bouquet_adj(g)))
    return out


def main():
    rng = random.Random(15)
    budget_s = float(sys.argv[1]) if len(sys.argv) > 1 else 900
    t0 = time.time()
    results = []
    for name, (n, adj) in seeds():
        poly = independence_poly(n, adj)
        vs = valley_score(poly)
        cur = (n, adj, (vs["window"], vs["ratio"]))
        best = cur
        evals = 0
        stall = 0
        share = budget_s / 3
        t1 = time.time()
        while time.time() - t1 < share:
            m = mutate_tree(cur[0], cur[1], rng)
            if m is None:
                continue
            n2, adj2 = m
            if n2 < 30 or n2 > 600:
                continue
            poly = independence_poly(n2, adj2)
            vs = valley_score(poly)
            evals += 1
            key = (vs["window"], vs["ratio"])
            if vs["witness"]:
                print(f"*** WITNESS *** seed={name} n={n2}")
                print("edges=", [(u, v) for u in range(n2)
                                 for v in adj2[u] if u < v])
                print("poly=", poly, flush=True)
                return
            if key > cur[2] or (stall > 400 and rng.random() < 0.02):
                cur = (n2, adj2, key)
                stall = 0
                if key > best[2]:
                    best = cur
            else:
                stall += 1
        print(f"{name}: evals={evals} best=(win={best[2][0]}, "
              f"V={best[2][1]:.10f}) n={best[0]}", flush=True)
        results.append((name, best))
    print("done", time.time() - t0, "s")


if __name__ == "__main__":
    main()
