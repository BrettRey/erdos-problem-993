#!/usr/bin/env python3
"""Attack 2: root herding — search in root-geometry space.

A non-unimodal coefficient sequence with rebound distance c-b >> 1
requires a conjugate root pair contributing a long-wavelength
oscillation with near-dominant amplitude: small angular deviation from
... precisely, oscillation wavelength in k is 2*pi/theta where theta is
the root's angle seen from the negative real axis, and amplitude decays
like (|z_min|/|z|)^k. Known families' roots hug the negative axis
(phi = pi - |arg z| near 0) or Bautista's circle. This search maximizes
phi among near-dominant roots directly.

Phase 1: empirical frontier over known objects (failure trees, brooms,
paths, bouquets, random trees).
Phase 2: evolutionary maximization of phi with unrestricted mutations.
Roots via numpy on normalized coefficients (degree <= ~75), spot-checked
with mpmath on champions.
"""

import json
import random
import sys
import time

import networkx as nx
import numpy as np

sys.path.insert(0, ".")

from indpoly import independence_poly
from scripts.valley_search import bouquet_adj, valley_score
from scratch_valley_freeform_20260715 import mutate_tree
from scratch_rgap_rescore_20260715 import rgap


def root_score(poly):
    """(phi_dom, phi_max, modulus_ratio) among non-real roots.

    phi = pi - |arg z|: 0 on the negative real axis, pi/2 at +-i,
    -> pi approaching the positive real axis.
    phi_dom restricts to near-dominant roots (|z| <= 2 |z|_min).
    """
    c = np.array(poly, dtype=float)
    c = c / c.max()
    roots = np.roots(c[::-1])
    if len(roots) == 0:
        return 0.0, 0.0, 0.0
    zmin = np.min(np.abs(roots))
    phi = np.pi - np.abs(np.angle(roots))
    nonreal = np.abs(roots.imag) > 1e-9 * np.abs(roots)
    if not nonreal.any():
        return 0.0, 0.0, zmin
    phi_max = float(phi[nonreal].max())
    dom = nonreal & (np.abs(roots) <= 2.0 * zmin)
    phi_dom = float(phi[dom].max()) if dom.any() else 0.0
    return phi_dom, phi_max, float(zmin)


def load_seeds():
    seeds = []
    d = json.load(open("results/analysis_n28_modal_lc_nm.json"))
    for i, f in enumerate(d["top_lc_failures"][:6]):
        G = nx.from_graph6_bytes(f["graph6"].encode())
        adj = [sorted(G.neighbors(v)) for v in range(G.number_of_nodes())]
        seeds.append((f"F28_{i}", len(adj), adj))
    def bq(g, p, label):
        spec = (tuple(tuple(sorted(x)) for x in g), tuple(p), 0)
        n, adj = bouquet_adj(*spec)
        return (label, n, adj)
    seeds.append(bq([], [1] * 30 + [12], "broom(12,30)"))
    seeds.append(bq([], [100], "P101"))
    seeds.append(bq([(2,) * 6] * 6, [], "T_{6,6,1}"))
    seeds.append(bq([(2,) * 3, (2,) * 10, (2,) * 11], [], "T_{3,10,11}"))
    seeds.append(bq([(3,) * 6] * 4, [], "4xS(3^6)"))
    seeds.append(bq([], [2] * 25 + [3] * 8, "substar(25x2,8x3)"))
    return seeds


def main():
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 900
    rng = random.Random(2026)
    seeds = load_seeds()

    print("=== phase 1: empirical phi frontier ===")
    print(f"{'object':>24} {'n':>5} {'phi_dom':>9} {'phi_max':>9} "
          f"{'|z|min':>8}")
    pool = []
    for label, n, adj in seeds:
        poly = independence_poly(n, adj)
        pd, pm, zm = root_score(poly)
        print(f"{label:>24} {n:>5} {pd:>9.5f} {pm:>9.5f} {zm:>8.4f}",
              flush=True)
        pool.append([pd, pm, n, adj, label])
    for i in range(150):
        n = rng.randint(30, 120)
        G = nx.random_labeled_tree(n, seed=rng.randint(0, 10**9))
        adj = [sorted(G.neighbors(v)) for v in range(n)]
        poly = independence_poly(n, adj)
        pd, pm, zm = root_score(poly)
        pool.append([pd, pm, n, adj, f"rand{i}"])
    pool.sort(key=lambda r: -r[0])
    print(f"random-tree phi_dom range: "
          f"{min(r[0] for r in pool):.5f}..{max(r[0] for r in pool):.5f}")

    print("\n=== phase 2: evolutionary phi maximization ===", flush=True)
    population = pool[:12]
    best = population[0][0]
    t0 = time.time()
    evals = 0
    while time.time() - t0 < budget:
        parent = population[rng.randrange(len(population))]
        m = mutate_tree(parent[2], parent[3], rng)
        if m is None:
            continue
        n2, adj2 = m
        if n2 < 30 or n2 > 150:
            continue
        poly = independence_poly(n2, adj2)
        vs = valley_score(poly)
        evals += 1
        if vs["witness"]:
            print("*** WITNESS ***", poly)
            return
        pd, pm, _ = root_score(poly)
        if pd > population[-1][0]:
            population.append([pd, pm, n2, adj2, "evolved"])
            population.sort(key=lambda r: -r[0])
            population = population[:12]
            if pd > best:
                best = pd
    print(f"evals={evals}  best phi_dom={best:.5f}", flush=True)
    print("\ntop evolved/seed objects, with corrected rebound metric:")
    for pd, pm, n, adj, label in population[:8]:
        poly = independence_poly(n, adj)
        r = rgap(poly, thetas=(0.001,))[0.001]
        rr = f"R={r['R']:.6f} c-b={r['c']-r['b']}" if r else "no descent"
        print(f"  phi_dom={pd:.5f} phi_max={pm:.5f} n={n:>4} {label:>10} "
              f"{rr}", flush=True)
    # dump best adjacency for reproducibility
    bestrec = population[0]
    json.dump({"phi_dom": bestrec[0], "n": bestrec[2],
               "edges": [(u, v) for u in range(bestrec[2])
                         for v in bestrec[3][u] if u < v]},
              open("results/root_herding_best_20260716.json", "w"))
    print("best tree saved to results/root_herding_best_20260716.json")


if __name__ == "__main__":
    main()
