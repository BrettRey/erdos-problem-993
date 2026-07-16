#!/usr/bin/env python3
"""Attack: map the dominance-angle Pareto frontier by evolution.

Instead of maximizing angle at one dominance cutoff (which the 611k-eval
run gamed by riding the boundary), maintain per-threshold champions
f_tau(T) = max certified angle among roots with modulus ratio <= tau,
for tau in {1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0}. A child enters the
archive if it improves any f_tau. Mutations: free tree moves plus
subtree crossover between archive members. Certified Arb roots in the
loop; exact witness + R_gap checks on every child.

Seeds: T_{6,6,1}, the 0.444 champion, failure trees, and the best
aperiodic-composition constructions if present.

Run with the scratchpad flint venv.
"""

import json
import math
import random
import sys
import time

sys.path.insert(0, ".")

from flint import fmpz_poly

from indpoly import independence_poly
from scripts.valley_search import bouquet_adj, valley_score
from scratch_valley_freeform_20260715 import mutate_tree
from scratch_rgap_rescore_20260715 import rgap

TAUS = (1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0)
OUT = "results/pareto_root_frontier_20260716.json"


def f_vector(poly):
    P = fmpz_poly(list(poly))
    roots = P.complex_roots()
    mods = [abs(complex(float(z.real), float(z.imag))) for z, m in roots]
    zmin = min(mods)
    pts = []
    for (z, m), mod in zip(roots, mods):
        if z.imag.rad() >= abs(z.imag.mid()):
            continue
        phi = math.pi - abs(math.atan2(float(z.imag), float(z.real)))
        pts.append((mod / zmin, phi))
    return [max((phi for ra, phi in pts if ra <= tau), default=0.0)
            for tau in TAUS]


# ---------------------------------------------------------------------------
# subtree surgery
# ---------------------------------------------------------------------------

def _component(adj, c, avoid):
    comp = {c}
    stack = [c]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y != avoid and y not in comp:
                comp.add(y)
                stack.append(y)
    return comp


def crossover(adjA, adjB, rng):
    nA, nB = len(adjA), len(adjB)
    for _ in range(20):
        # cut a small subtree out of A
        p = rng.randrange(nA)
        if not adjA[p]:
            continue
        c = rng.choice(adjA[p])
        compA = _component(adjA, c, p)
        if len(compA) > nA // 3:
            continue
        # extract a small subtree from B
        p2 = rng.randrange(nB)
        if not adjB[p2]:
            continue
        c2 = rng.choice(adjB[p2])
        compB = _component(adjB, c2, p2)
        if len(compB) > nB // 3:
            continue
        keep = [v for v in range(nA) if v not in compA]
        remapA = {v: i for i, v in enumerate(keep)}
        newadj = [[remapA[w] for w in adjA[v] if w not in compA]
                  for v in keep]
        offs = len(newadj)
        listB = sorted(compB)
        remapB = {v: offs + i for i, v in enumerate(listB)}
        for v in listB:
            newadj.append([remapB[w] for w in adjB[v] if w in compB])
        attach_at = rng.randrange(offs)
        newadj[attach_at].append(remapB[c2])
        newadj[remapB[c2]].append(attach_at)
        return len(newadj), newadj
    return None


def seeds():
    out = []

    def bq(g, p, label):
        spec = (tuple(tuple(sorted(x)) for x in g), tuple(p), 0)
        n, adj = bouquet_adj(*spec)
        return [label, n, adj]

    out.append(bq([(2,) * 6] * 6, [], "T_{6,6,1}"))
    out.append(bq([(3,) * 6] * 4, [], "4xS(3^6)"))
    try:
        d = json.load(open("results/certified_root_evolution_20260716.json"))
        n = d["n"]
        adj = [[] for _ in range(n)]
        for u, v in d["edges"]:
            adj[u].append(v)
            adj[v].append(u)
        out.append(["champ0.444", n, adj])
    except Exception:
        pass
    try:
        recs = json.load(open("results/aperiodic_composition_20260716.json"))
        bestper = {}
        for r in recs:
            key = r["word"]
            top = max((phi for ra, phi in r["points"] if ra <= 1.5),
                      default=0.0)
            if top > bestper.get(key, (0.0, None))[0]:
                bestper[key] = (top, r)
        for key, (top, r) in bestper.items():
            n = r["n"]
            adj = [[] for _ in range(n)]
            for u, v in r["edges"]:
                adj[u].append(v)
                adj[v].append(u)
            out.append([f"aper-{key}", n, adj])
    except Exception as e:
        print("no aperiodic seeds:", e)
    return out


def main():
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 2700
    rng = random.Random(1618)
    champs = {}     # tau index -> [f, n, adj, label]
    def consider(fv, n, adj, label):
        improved = False
        for i, v in enumerate(fv):
            if v > champs.get(i, (0.0,))[0]:
                champs[i] = (v, n, adj, label)
                improved = True
        return improved

    for label, n, adj in seeds():
        fv = f_vector(independence_poly(n, adj))
        consider(fv, n, adj, label)
        print(f"seed {label:>16}: " +
              " ".join(f"{v:.3f}" for v in fv), flush=True)

    t0 = time.time()
    evals = 0
    last_report = 0
    while time.time() - t0 < budget:
        pool = list(champs.values())
        a = pool[rng.randrange(len(pool))]
        if len(pool) > 1 and rng.random() < 0.35:
            b = pool[rng.randrange(len(pool))]
            m = crossover(a[2], b[2], rng)
        else:
            m = mutate_tree(a[1], a[2], rng)
        if m is None:
            continue
        n2, adj2 = m
        if n2 < 30 or n2 > 200:
            continue
        poly = independence_poly(n2, adj2)
        vs = valley_score(poly)
        if vs["witness"]:
            print("*** WITNESS ***")
            print("poly=", poly)
            json.dump({"witness": True, "n": n2,
                       "edges": [(u, v) for u in range(n2)
                                 for v in adj2[u] if u < v]},
                      open(OUT, "w"))
            return
        try:
            fv = f_vector(poly)
        except Exception:
            continue
        evals += 1
        if consider(fv, n2, adj2, "evolved"):
            now = time.time() - t0
            if now - last_report > 30:
                last_report = now
                r = rgap(poly, thetas=(0.001,))[0.001]
                rr = (f"R={r['R']:.4f} cb={r['c']-r['b']}"
                      if r else "nodesc")
                print(f"[{now:6.0f}s] frontier: " +
                      " ".join(f"{champs.get(i,(0,))[0]:.3f}"
                               for i in range(len(TAUS))) +
                      f"  (n={n2}, {rr}, evals={evals})", flush=True)

    print(f"\ndone: evals={evals} in {time.time()-t0:.0f}s")
    print(f"{'tau':>6} {'phi*':>8} {'n':>5} {'label':>12}")
    final = {}
    for i, tau in enumerate(TAUS):
        if i in champs:
            v, n, adj, label = champs[i]
            print(f"{tau:>6} {v:>8.5f} {n:>5} {label:>12}")
            final[str(tau)] = {"phi": v, "n": n, "label": label,
                               "edges": [(u, w) for u in range(n)
                                         for w in adj[u] if u < w]}
    json.dump(final, open(OUT, "w"))
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
