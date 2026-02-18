#!/usr/bin/env python3
"""Check tail ratio monotonicity for the perturbation lemma.

Given a tree T with a leaf v whose neighbor u has degree 2 (other neighbor w),
define:

  T'' = T - {u, v}
  f = I(T'')
  g = I(T'' - w)

We test whether r_k = g_k / f_k is nonincreasing for k >= d(f),
where d(f) is the first descent index of f.
"""

from __future__ import annotations

import argparse
import json
from typing import Dict, List, Tuple

from indpoly import independence_poly
from trees import trees


def encode_graph6_small(adj: List[List[int]]) -> str:
    n = len(adj)
    if n >= 63:
        raise ValueError("n too large for small graph6 encoder")
    aset = [set(nei) for nei in adj]
    bits: List[int] = []
    for j in range(1, n):
        sj = aset[j]
        for i in range(j):
            bits.append(1 if i in sj else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k : k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)


def first_descent(seq: List[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return -1


def mode_index(seq: List[int]) -> int:
    maxv = max(seq) if seq else 0
    idx = -1
    for i, v in enumerate(seq):
        if v == maxv:
            idx = i
    return idx


def induced_subgraph(adj: List[List[int]], removed: List[bool]) -> List[List[int]]:
    n = len(adj)
    mapping = [-1] * n
    keep: List[int] = []
    for i in range(n):
        if not removed[i]:
            mapping[i] = len(keep)
            keep.append(i)
    m = len(keep)
    out: List[List[int]] = [[] for _ in range(m)]
    for old in keep:
        i = mapping[old]
        for nb in adj[old]:
            j = mapping[nb]
            if j != -1:
                out[i].append(j)
    return out


def tail_ratio_ok(f: List[int], g: List[int], start: int) -> Tuple[bool, int]:
    # check g_{k+1}/f_{k+1} <= g_k/f_k for k >= start
    # return (ok, first_bad_k)
    if start < 0:
        return True, -1
    max_k = min(len(g), len(f)) - 2
    if max_k < start:
        return True, -1
    for k in range(start, max_k + 1):
        # compare g_{k+1} * f_k <= g_k * f_{k+1}
        if g[k + 1] * f[k] > g[k] * f[k + 1]:
            return False, k
    return True, -1


def check_tree(
    n: int, adj: List[List[int]], max_witnesses: int, witnesses: List[Dict]
) -> Tuple[int, int, int, int]:
    """Return (leaf_cases, ratio_failures, mode_failures, tail_h_failures, pre_h_failures)."""
    leaf_cases = 0
    ratio_failures = 0
    mode_failures = 0
    tail_h_failures = 0
    pre_h_failures = 0
    for v in range(n):
        if len(adj[v]) != 1:
            continue
        u = adj[v][0]
        if len(adj[u]) != 2:
            continue
        # u has exactly two neighbors: v and w
        w = adj[u][0] if adj[u][1] == v else adj[u][1]
        leaf_cases += 1

        removed_f = [False] * n
        removed_f[u] = True
        removed_f[v] = True
        adj_f = induced_subgraph(adj, removed_f)
        f = independence_poly(len(adj_f), adj_f)

        removed_g = [False] * n
        removed_g[u] = True
        removed_g[v] = True
        removed_g[w] = True
        adj_g = induced_subgraph(adj, removed_g)
        g = independence_poly(len(adj_g), adj_g)

        d = first_descent(f)
        ok, bad_k = tail_ratio_ok(f, g, d if d != -1 else len(f))
        if not ok:
            ratio_failures += 1
            if len(witnesses) < max_witnesses:
                witnesses.append(
                    {
                        "n": n,
                        "graph6": encode_graph6_small(adj),
                        "leaf": v,
                        "u": u,
                        "w": w,
                        "d_f": d,
                        "bad_k": bad_k,
                        "f": f,
                        "g": g,
                    }
                )

        # mode(g) <= d(f) check (tail nonincreasing for g past d(f))
        mode_g = mode_index(g)
        if d != -1 and mode_g > d:
            mode_failures += 1

        # tail check on h = (1+x)f + xg for k >= d+1
        if d != -1:
            deg_h = max(len(f), len(g) + 1) - 1
            h = [0] * (deg_h + 1)
            for k, fk in enumerate(f):
                h[k] += fk
                if k + 1 < len(h):
                    h[k + 1] += fk
            for k, gk in enumerate(g):
                if k + 1 < len(h):
                    h[k + 1] += gk
            for k in range(d + 1, len(h) - 1):
                if h[k + 1] > h[k]:
                    tail_h_failures += 1
                    break
            for k in range(0, d):
                if h[k + 1] < h[k]:
                    pre_h_failures += 1
                    break
    return leaf_cases, ratio_failures, mode_failures, tail_h_failures, pre_h_failures


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--backend", default="auto", choices=["auto", "geng", "networkx"])
    ap.add_argument("--out", default="")
    ap.add_argument("--max-witnesses", type=int, default=5)
    args = ap.parse_args()

    stats = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "backend": args.backend,
        "trees": 0,
        "leaf_cases": 0,
        "ratio_failures": 0,
        "mode_failures": 0,
        "tail_h_failures": 0,
        "pre_h_failures": 0,
        "by_n": [],
        "witnesses": [],
    }

    for n in range(args.min_n, args.max_n + 1):
        tree_count = 0
        leaf_cases = 0
        ratio_failures = 0
        mode_failures = 0
        tail_h_failures = 0
        pre_h_failures = 0
        for _, adj in trees(n, backend=args.backend):
            tree_count += 1
            lc, rf, mf, thf, phf = check_tree(n, adj, args.max_witnesses, stats["witnesses"])
            leaf_cases += lc
            ratio_failures += rf
            mode_failures += mf
            tail_h_failures += thf
            pre_h_failures += phf
        stats["trees"] += tree_count
        stats["leaf_cases"] += leaf_cases
        stats["ratio_failures"] += ratio_failures
        stats["mode_failures"] += mode_failures
        stats["tail_h_failures"] += tail_h_failures
        stats["pre_h_failures"] += pre_h_failures
        stats["by_n"].append(
            {
                "n": n,
                "trees": tree_count,
                "leaf_cases": leaf_cases,
                "ratio_failures": ratio_failures,
                "mode_failures": mode_failures,
                "tail_h_failures": tail_h_failures,
                "pre_h_failures": pre_h_failures,
            }
        )

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2)
    else:
        print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
