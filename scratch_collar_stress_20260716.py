#!/usr/bin/env python3
"""Kill-tests for the real-collar conjecture (Conjecture A).

Phase 1: mass certified scan over random, structured, and mutated
trees: record the minimum modulus ratio of any certified non-real zero
(collar margin) and the minimum |arg z| over non-real zeros
(positive-axis sector margin, relevant to the mid-band bridge).
Phase 2: adversarial evolution MINIMIZING the non-real ratio.

Any certified non-real zero at ratio <= 1.1 kills Conjecture A as
stated. Run with the scratchpad flint venv.
"""

import json
import math
import random
import sys
import time

sys.path.insert(0, ".")

import networkx as nx
from flint import fmpz_poly

from indpoly import independence_poly
from scripts.valley_search import bouquet_adj
from scratch_valley_freeform_20260715 import mutate_tree

OUT = "results/collar_stress_20260716.json"


def spectrum_margins(poly):
    """(min nonreal modulus ratio, min |arg| over nonreal zeros)."""
    P = fmpz_poly(list(poly))
    roots = P.complex_roots()
    mods = [abs(complex(float(z.real), float(z.imag))) for z, m in roots]
    zmin = min(mods)
    ratio_min = float("inf")
    arg_min = float("inf")
    for (z, m), mod in zip(roots, mods):
        if z.imag.rad() >= abs(z.imag.mid()):
            continue
        ratio_min = min(ratio_min, mod / zmin)
        arg_min = min(arg_min, abs(math.atan2(float(z.imag),
                                              float(z.real))))
    return ratio_min, arg_min


def random_tree(n, rng):
    G = nx.random_labeled_tree(n, seed=rng.randint(0, 10 ** 9))
    return [sorted(G.neighbors(v)) for v in range(n)]


def structured_sample():
    out = []

    def bq(g, p):
        spec = (tuple(tuple(sorted(x)) for x in g), tuple(p), 0)
        return bouquet_adj(*spec)

    for t in range(2, 9):
        for m in range(2, 9):
            out.append(bq([(2,) * t] * m, []))
    for a in range(3, 30, 3):
        for b in range(0, 20, 4):
            out.append(bq([], [2] * a + [3] * b))
    for w in range(5, 40, 5):
        for L in (5, 10, 20):
            out.append(bq([], [1] * w + [L]))
    return out


def main():
    budget2 = float(sys.argv[1]) if len(sys.argv) > 1 else 900
    rng = random.Random(1101)
    worst = {"ratio": float("inf"), "arg": float("inf")}
    worst_trees = {}
    scanned = 0
    t0 = time.time()

    def check(n, adj, label):
        nonlocal scanned
        poly = independence_poly(n, adj)
        r, a = spectrum_margins(poly)
        scanned += 1
        for key, val in (("ratio", r), ("arg", a)):
            if val < worst[key]:
                worst[key] = val
                worst_trees[key] = {
                    "label": label, "n": n, key: val,
                    "edges": [(u, v) for u in range(n)
                              for v in adj[u] if u < v]}
                if key == "ratio" and val <= 1.1:
                    print(f"*** COLLAR VIOLATION *** ratio={val:.5f} "
                          f"{label} n={n}", flush=True)
        return r

    print("phase 1: mass scan", flush=True)
    for n, adj in structured_sample():
        check(n, adj, "structured")
    for i in range(12000):
        n = rng.randint(30, 140)
        check(n, random_tree(n, rng), f"rand{i}")
        if i % 2000 == 0:
            print(f"  {scanned} scanned, min ratio={worst['ratio']:.4f}, "
                  f"min arg={worst['arg']:.4f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    d28 = json.load(open("results/analysis_n28_modal_lc_nm.json"))
    fail_seeds = []
    for f in d28["top_lc_failures"][:5]:
        G = nx.from_graph6_bytes(f["graph6"].encode())
        fail_seeds.append([sorted(G.neighbors(v))
                           for v in range(G.number_of_nodes())])
    for i in range(2000):
        adj = fail_seeds[i % len(fail_seeds)]
        m = mutate_tree(len(adj), adj, rng)
        if m is None:
            continue
        n2, adj2 = m
        if 30 <= n2 <= 140:
            check(n2, adj2, f"mut{i}")
    print(f"phase 1 done: {scanned} trees, min nonreal ratio "
          f"{worst['ratio']:.5f}, min |arg| {worst['arg']:.5f}",
          flush=True)

    print("\nphase 2: adversarial minimization of nonreal ratio",
          flush=True)
    seed = worst_trees["ratio"]
    n = seed["n"]
    adj = [[] for _ in range(n)]
    for u, v in seed["edges"]:
        adj[u].append(v)
        adj[v].append(u)
    cur = (n, adj, worst["ratio"])
    best = cur[2]
    t1 = time.time()
    evals = 0
    stall = 0
    while time.time() - t1 < budget2:
        m = mutate_tree(cur[0], cur[1], rng)
        if m is None:
            continue
        n2, adj2 = m
        if not (30 <= n2 <= 160):
            continue
        r = check(n2, adj2, "adv")
        evals += 1
        if r < cur[2] or (stall > 300 and rng.random() < 0.03):
            cur = (n2, adj2, r)
            stall = 0
            if r < best:
                best = r
                print(f"  [{time.time()-t1:5.0f}s] adv min ratio="
                      f"{r:.5f} n={n2} (evals={evals})", flush=True)
        else:
            stall += 1
    print(f"\nphase 2 done: evals={evals}, adversarial min ratio "
          f"{best:.5f}", flush=True)
    print(f"global minima: nonreal ratio {worst['ratio']:.5f}, "
          f"|arg| {worst['arg']:.5f} over {scanned} certified trees")
    json.dump({"scanned": scanned, "min_nonreal_ratio": worst["ratio"],
               "min_arg": worst["arg"], "extremal": worst_trees},
              open(OUT, "w"))
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
