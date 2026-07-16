#!/usr/bin/env python3
"""Forest/product valley search + tree realization via path joins.

Lane 5 from README ("Forest/product search route"), run with a direct
valley score instead of a log-concavity defect, and extended by a tree
realization step: for the best product valleys, build the actual TREE
consisting of the two component bouquets joined root-to-root by a path
of length ell, and score it exactly.

Rationale: I(T1 ⊔ T2) = I(T1) I(T2) (forest). Joining by a path of
length ell perturbs the product by lower-order terms that decay with
ell, so a robust product valley is a candidate valley for a genuine
tree. All arithmetic exact.
"""

import argparse
import heapq
import json
import os
import sys
import time

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

from indpoly import independence_poly  # noqa: E402
from scripts.valley_search import (  # noqa: E402
    bouquet_adj,
    bouquet_poly,
    bouquet_size,
    spec_label,
    valley_score,
    _polymul,
)


# ---------------------------------------------------------------------------
# Component library: (label, gadgets, extra) with diverse density profiles
# ---------------------------------------------------------------------------

def build_library(max_comp_n: int):
    lib = []

    def add(gadgets, paths=(), leaves=0):
        spec = (tuple(tuple(sorted(g)) for g in gadgets), tuple(paths), leaves)
        if bouquet_size(spec[0], spec[1], spec[2]) <= max_comp_n:
            lib.append(spec)

    # Galvin bouquets T_{m,t,1}, including LC-failing shapes (t <= m)
    for t in range(2, 11):
        for m in (t, t + 1, t + 2, 2 * t, 3 * t, 4 * t):
            add([(2,) * t] * m)
    # Kadrawi-Levit diagonal
    for k in range(3, 12):
        add([(2,) * 3, (2,) * k, (2,) * (k + 1)])
        add([(2, 2, 4), (2,) * k, (2,) * k])
    # subdivided stars: dense-low profile
    for a2 in (6, 12, 20, 30, 45):
        for a3 in (0, 6, 12, 20):
            if a2 + a3 >= 3:
                add([], [2] * a2 + [3] * a3)
    # stars / near-stars: mode at k=1..2, sharply decreasing (low profile)
    for w in (10, 20, 40, 80):
        add([], [1] * w)
        add([(1,) * w])
    # paths: flat middle profile
    for L in (20, 40, 80):
        add([], [L])
    # brooms: mode ~ 0.5 alpha
    for L in (10, 20, 30):
        for w in (20, 40, 80):
            add([], [1] * w + [L])
    # deep double-layer: gadget legs long
    for t in (4, 6, 8):
        for m in (4, 8, 12):
            add([(3,) * t] * m)
            add([(4,) * t] * m)

    # dedupe
    seen = set()
    out = []
    for spec in lib:
        if spec not in seen:
            seen.add(spec)
            out.append(spec)
    return out


# ---------------------------------------------------------------------------
# Tree realization: two bouquets joined root-to-root by a path of length ell
# ---------------------------------------------------------------------------

def dumbbell_adj(specA, specB, ell: int):
    nA, adjA = bouquet_adj(*specA)
    nB, adjB = bouquet_adj(*specB)
    # path of ell edges means ell-1 intermediate vertices when ell >= 1
    inter = max(0, ell - 1)
    n = nA + nB + inter
    adj = [list(x) for x in adjA] + [[] for _ in range(inter)] \
        + [[x + nA + inter for x in row] for row in adjB]

    def edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    rootA, rootB = 0, nA + inter
    if ell == 0:
        raise ValueError("ell >= 1 required")
    prev = rootA
    for i in range(inter):
        edge(prev, nA + i)
        prev = nA + i
    edge(prev, rootB)
    return n, adj


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-comp-n", type=int, default=260)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--ells", default="1,2,3,4,6,9")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    lib = build_library(args.max_comp_n)
    print(f"library: {len(lib)} components (n <= {args.max_comp_n})",
          flush=True)
    polys = []
    for spec in lib:
        p = bouquet_poly(list(spec[0]), spec[1], spec[2])
        polys.append(p)

    t0 = time.time()
    heap = []
    counter = 0
    product_witnesses = []
    for i in range(len(lib)):
        for j in range(i, len(lib)):
            prod = _polymul(polys[i], polys[j])
            vs = valley_score(prod)
            counter += 1
            k = (vs["window"], vs["ratio"])
            item = (k, counter, i, j, vs)
            if vs["witness"]:
                product_witnesses.append(item)
                print(f"*** PRODUCT VALLEY *** {spec_label(lib[i])} x "
                      f"{spec_label(lib[j])} ratio={vs['ratio']:.6f}",
                      flush=True)
            if len(heap) < args.top:
                heapq.heappush(heap, item)
            elif k > heap[0][0]:
                heapq.heapreplace(heap, item)
    print(f"scored {counter} pairwise products in {time.time()-t0:.1f}s; "
          f"product-level witnesses: {len(product_witnesses)}", flush=True)

    top = sorted(heap, key=lambda it: it[0], reverse=True)
    ells = [int(x) for x in args.ells.split(",")]
    tree_results = []
    tree_witnesses = []
    for k, _, i, j, vs in (product_witnesses + top):
        for ell in ells:
            n, adj = dumbbell_adj(lib[i], lib[j], ell)
            poly = independence_poly(n, adj)
            tvs = valley_score(poly)
            rec = {
                "A": spec_label(lib[i]), "B": spec_label(lib[j]),
                "specA": [[list(g) for g in lib[i][0]], list(lib[i][1]), lib[i][2]],
                "specB": [[list(g) for g in lib[j][0]], list(lib[j][1]), lib[j][2]],
                "ell": ell, "n": n, "alpha": len(poly) - 1,
                "product_ratio": vs["ratio"], **tvs,
            }
            if tvs["witness"]:
                rec["poly"] = poly
                tree_witnesses.append(rec)
                print(f"\n*** TREE VALLEY WITNESS *** {rec['A']} --{ell}-- "
                      f"{rec['B']} n={n}")
                print(f"    poly={poly}\n", flush=True)
            tree_results.append(rec)

    tree_results.sort(key=lambda r: (r["window"], r["ratio"]), reverse=True)
    print(f"\ntop tree (dumbbell) valley ratios:")
    for r in tree_results[:20]:
        print(f"  V={r['ratio']:.10f} win={str(r['window']):5s} "
              f"b={r['pos']:>3} rise={r['rise_pos']:>3} n={r['n']:>4} "
              f"ell={r['ell']} prodV={r['product_ratio']:.6f}")
        print(f"      A={r['A']}")
        print(f"      B={r['B']}")

    if args.out:
        with open(args.out, "w") as f:
            json.dump({
                "library_size": len(lib),
                "pair_products": counter,
                "product_witnesses": len(product_witnesses),
                "tree_witnesses": tree_witnesses,
                "top_trees": [
                    {k2: v for k2, v in r.items() if k2 != "poly"}
                    for r in tree_results[:60]
                ],
            }, f, indent=2)
        print(f"saved {args.out}")


if __name__ == "__main__":
    main()
