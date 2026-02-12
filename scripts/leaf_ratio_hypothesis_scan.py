#!/usr/bin/env python3
"""Creative hypothesis scan for leaf-step ratio monotonicity.

Hypothesis H_R (leaf case):
  For f = I(T), g = I(T-w) with w a leaf, the ratio sequence
    R_k := g_k / f_k
  is nonincreasing for k = 0..d(g)-2.

Equivalent check without division:
  g_{k+1} * f_k <= g_k * f_{k+1}.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import _polymul, independence_poly
from trees import trees


def first_descent(seq: list[int]) -> int:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return len(seq)


def encode_graph6_small(adj: list[list[int]]) -> str:
    n = len(adj)
    if n >= 63:
        return "<n>=63"
    aset = [set(nei) for nei in adj]
    bits: list[int] = []
    for j in range(1, n):
        sj = aset[j]
        for i in range(j):
            bits.append(1 if i in sj else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for t in range(0, len(bits), 6):
        v = 0
        for b in bits[t : t + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)


def forest_poly_removed_set(adj: list[list[int]], removed: set[int]) -> list[int]:
    n = len(adj)
    remaining = [i for i in range(n) if i not in removed]
    if not remaining:
        return [1]
    rem_set = set(remaining)
    seen: set[int] = set()
    out = [1]
    for start in remaining:
        if start in seen:
            continue
        comp: list[int] = []
        stack = [start]
        seen.add(start)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v in rem_set and v not in seen:
                    seen.add(v)
                    stack.append(v)
        mapping = {old: i for i, old in enumerate(comp)}
        cadj: list[list[int]] = [[] for _ in range(len(comp))]
        for old in comp:
            i = mapping[old]
            for v in adj[old]:
                j = mapping.get(v)
                if j is not None:
                    cadj[i].append(j)
        out = _polymul(out, independence_poly(len(comp), cadj))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=2)
    ap.add_argument("--max-n", type=int, default=18)
    ap.add_argument("--backend", default="auto", choices=["auto", "geng", "networkx"])
    ap.add_argument("--progress-every", type=int, default=10000)
    ap.add_argument("--max-witnesses", type=int, default=8)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    stats: dict[str, Any] = {
        "hypothesis": "For leaf w, g_k/f_k is nonincreasing for k<=d(g)-2",
        "min_n": args.min_n,
        "max_n": args.max_n,
        "backend": args.backend,
        "trees": 0,
        "leaf_cases": 0,
        "checked_pairs": 0,
        "violations": 0,
        "violations_full_prefix": 0,
        "gap_counts": {},
        "witnesses": [],
        "by_n": [],
    }

    for n in range(args.min_n, args.max_n + 1):
        n_trees = 0
        n_leaf_cases = 0
        n_pairs = 0
        n_viol = 0
        n_viol_full = 0
        n_gap_counts: dict[str, int] = {}

        for _, adj in trees(n, backend=args.backend):
            n_trees += 1
            f = independence_poly(n, adj)
            leaves = [w for w in range(n) if len(adj[w]) == 1]
            if not leaves:
                continue

            for w in leaves:
                n_leaf_cases += 1
                g = forest_poly_removed_set(adj, {w})
                d_g = first_descent(g)
                d_f = first_descent(f)
                gap = d_f - d_g
                n_gap_counts[str(gap)] = n_gap_counts.get(str(gap), 0) + 1
                stats["gap_counts"][str(gap)] = stats["gap_counts"].get(str(gap), 0) + 1

                # target prefix: k <= d(g)-2
                kmax = min(d_g - 2, min(len(f), len(g)) - 2)
                bad_k = None
                if kmax >= 0:
                    for k in range(kmax + 1):
                        n_pairs += 1
                        if g[k + 1] * f[k] > g[k] * f[k + 1]:
                            bad_k = k
                            n_viol += 1
                            break

                # full common prefix (diagnostic)
                full_kmax = min(len(f), len(g)) - 2
                if full_kmax >= 0:
                    for k in range(full_kmax + 1):
                        if g[k + 1] * f[k] > g[k] * f[k + 1]:
                            n_viol_full += 1
                            break

                if bad_k is not None and len(stats["witnesses"]) < args.max_witnesses:
                    stats["witnesses"].append(
                        {
                            "n": n,
                            "graph6": encode_graph6_small(adj),
                            "leaf": w,
                            "bad_k": bad_k,
                            "d_f": d_f,
                            "d_g": d_g,
                            "f_k": f[bad_k],
                            "f_k1": f[bad_k + 1],
                            "g_k": g[bad_k],
                            "g_k1": g[bad_k + 1],
                        }
                    )

            if args.progress_every > 0 and n_trees % args.progress_every == 0:
                print(
                    f"n={n} trees={n_trees:,} leaf_cases={n_leaf_cases:,} "
                    f"pairs={n_pairs:,} viol={n_viol:,} full_viol={n_viol_full:,}"
                )

        stats["trees"] += n_trees
        stats["leaf_cases"] += n_leaf_cases
        stats["checked_pairs"] += n_pairs
        stats["violations"] += n_viol
        stats["violations_full_prefix"] += n_viol_full
        stats["by_n"].append(
            {
                "n": n,
                "trees": n_trees,
                "leaf_cases": n_leaf_cases,
                "checked_pairs": n_pairs,
                "violations": n_viol,
                "violations_full_prefix": n_viol_full,
                "gap_counts": n_gap_counts,
            }
        )

        print(
            f"done n={n}: trees={n_trees:,} leaf_cases={n_leaf_cases:,} "
            f"pairs={n_pairs:,} viol={n_viol:,} full_viol={n_viol_full:,}"
        )

    status = "PASS" if stats["violations"] == 0 else "FAIL"
    print(
        f"{status}: trees={stats['trees']:,} leaf_cases={stats['leaf_cases']:,} "
        f"pairs={stats['checked_pairs']:,} viol={stats['violations']:,} "
        f"full_viol={stats['violations_full_prefix']:,}"
    )

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fobj:
            json.dump(stats, fobj, indent=2)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
