#!/usr/bin/env python3
"""Attack 1: interference experiments with the exact LC-failure trees.

Blocks: the 19 exact n=28 LC-failure trees and the 2 n=26 ones (the only
known trees whose coefficient sequences have a genuine convexity bump,
worst ratio 1.503 at k=14).

Experiments:
  A. All pairwise joins of failure trees at ALL vertex pairs, connector
     lengths 1-3. (Bump lands at ~0.93 of composite alpha: outside the
     window, but the subtractive correction acts mid-band; exact check.)
  B. Failure tree x large smooth partner (paths, brooms, subdivided
     stars) sized so the translated bump ridge lands mid-band of the
     composite; all failure-tree vertices x partner reps.
Scores: exact witness, R_gap(0.001), rise distance c-b.
"""

import json
import sys
import time

import networkx as nx

sys.path.insert(0, ".")

from indpoly import independence_poly
from scripts.valley_search import bouquet_adj, valley_score
from scratch_rgap_rescore_20260715 import rgap


def load_failure_trees():
    out = []
    d = json.load(open("results/analysis_n28_modal_lc_nm.json"))
    for i, f in enumerate(d["top_lc_failures"]):
        G = nx.from_graph6_bytes(f["graph6"].encode())
        adj = [sorted(G.neighbors(v)) for v in range(G.number_of_nodes())]
        out.append((f"F28_{i}(lc={f['lc_ratio']:.3f})", len(adj), adj))
    d26 = json.load(open("results/analysis_n26.json"))
    for i, f in enumerate(d26["lc_failures"]):
        g6 = f["graph6"] if isinstance(f, dict) else f
        G = nx.from_graph6_bytes(g6.encode())
        adj = [sorted(G.neighbors(v)) for v in range(G.number_of_nodes())]
        out.append((f"F26_{i}", len(adj), adj))
    return out


def join(nA, adjA, u, nB, adjB, v, ell):
    inter = max(0, ell - 1)
    n = nA + nB + inter
    adj = ([list(x) for x in adjA]
           + [[] for _ in range(inter)]
           + [[x + nA + inter for x in row] for row in adjB])

    def edge(a, b):
        adj[a].append(b)
        adj[b].append(a)

    prev = u
    for i in range(inter):
        edge(prev, nA + i)
        prev = nA + i
    edge(prev, v + nA + inter)
    return n, adj


def score(n, adj):
    poly = independence_poly(n, adj)
    vs = valley_score(poly)
    if vs["witness"]:
        return ("WITNESS", poly, vs)
    r = rgap(poly, thetas=(0.001,))[0.001]
    return ("ok", poly, r)


def main():
    trees = load_failure_trees()
    print(f"{len(trees)} failure-tree blocks loaded", flush=True)

    results = []
    witnesses = []

    # --- Experiment A: all pairwise joins, all vertex pairs, ell 1..3
    t0 = time.time()
    count = 0
    for i in range(len(trees)):
        lA, nA, adjA = trees[i]
        for j in range(i, len(trees)):
            lB, nB, adjB = trees[j]
            for u in range(nA):
                for v in range(nB):
                    if j == i and v < u:
                        continue
                    for ell in (1, 2, 3):
                        n, adj = join(nA, adjA, u, nB, adjB, v, ell)
                        tag, poly, r = score(n, adj)
                        count += 1
                        if tag == "WITNESS":
                            witnesses.append((lA, u, lB, v, ell, poly))
                            print(f"*** WITNESS *** {lA}[{u}]--{ell}--"
                                  f"{lB}[{v}]")
                            print("poly=", poly, flush=True)
                        elif r:
                            results.append((r["R"], r["c"] - r["b"], n,
                                            f"{lA}[{u}]--{ell}--{lB}[{v}]"))
        if i % 5 == 0:
            print(f"  expA block {i}: {count} joins, "
                  f"{time.time()-t0:.0f}s", flush=True)
    print(f"expA: {count} joins in {time.time()-t0:.0f}s, "
          f"witnesses={len(witnesses)}", flush=True)

    # --- Experiment B: failure tree x smooth partner, bump mid-band
    partners = []
    for L in (30, 60, 100, 160):
        partners.append((f"P{L+1}", *bouquet_adj(((), ), (L,), 0)[0:2])) \
            if False else None
    # build partners via bouquet_adj spec: (gadgets, paths, leaves)
    def partner(spec_g, spec_p, label):
        spec = (tuple(tuple(sorted(g)) for g in spec_g), tuple(spec_p), 0)
        n, adj = bouquet_adj(*spec)
        return (label, n, adj)

    partners = [
        partner([], [40], "P41"),
        partner([], [80], "P81"),
        partner([], [140], "P141"),
        partner([], [1] * 30 + [10], "broom(10,30)"),
        partner([], [1] * 60 + [15], "broom(15,60)"),
        partner([], [2] * 20 + [3] * 8, "substar(20x2,8x3)"),
        partner([], [2] * 40 + [3] * 15, "substar(40x2,15x3)"),
        partner([(2,) * 6] * 6, [], "T_{6,6,1}"),
    ]
    t1 = time.time()
    countB = 0
    for lA, nA, adjA in trees:
        for lB, nB, adjB in partners:
            # partner reps: root(0), a mid vertex, an end vertex
            reps = [0, nB // 2, nB - 1]
            for u in range(nA):
                for v in reps:
                    for ell in (1, 2):
                        n, adj = join(nA, adjA, u, nB, adjB, v, ell)
                        tag, poly, r = score(n, adj)
                        countB += 1
                        if tag == "WITNESS":
                            witnesses.append((lA, u, lB, v, ell, poly))
                            print(f"*** WITNESS *** {lA}[{u}]--{ell}--"
                                  f"{lB}[{v}]")
                            print("poly=", poly, flush=True)
                        elif r:
                            results.append((r["R"], r["c"] - r["b"], n,
                                            f"{lA}[{u}]--{ell}--{lB}[{v}]"))
    print(f"expB: {countB} joins in {time.time()-t1:.0f}s, "
          f"witnesses={len(witnesses)}", flush=True)

    results.sort(key=lambda t: (-(t[1] > 1), -t[0]))
    multi = [r for r in results if r[1] > 1]
    print(f"\njoins with rise distance c-b > 1: {len(multi)}")
    for R, gap, n, label in (multi[:10] + results[:12]):
        print(f"  R={R:.8f} n(1-R)={n*(1-R):8.2f} c-b={gap} n={n}  {label}")


if __name__ == "__main__":
    main()
