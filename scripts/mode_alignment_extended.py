#!/usr/bin/env python3
"""Extended mode-alignment check: mode(I(T-w)) <= d(I(T)) for n=19..21.

For every tree T on n vertices and every vertex w, we verify that:
    mode(I(T-w)) <= d(I(T))
where:
    mode(P) = last index achieving max coefficient of polynomial P
    d(I(T)) = first descent index of the independence polynomial of T

Uses multiprocessing with geng partitions for speed.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool
from typing import Any

# Add parent dir to path so we can import project modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import independence_poly, _polymul, _polyadd
from trees import trees_geng_raw


def mode_index(seq: list[int]) -> int:
    """Index of the last occurrence of the maximum value."""
    if not seq:
        return -1
    maxv = max(seq)
    idx = -1
    for i, v in enumerate(seq):
        if v == maxv:
            idx = i
    return idx


def first_descent(seq: list[int]) -> int:
    """First index where seq[i] < seq[i-1]. Returns len(seq) if no descent."""
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return len(seq)


def encode_graph6_small(adj: list[list[int]]) -> str:
    """Encode adjacency list to graph6 string for small graphs (n < 63)."""
    n = len(adj)
    if n >= 63:
        raise ValueError("n too large for small graph6 encoder")
    aset = [set(nei) for nei in adj]
    bits: list[int] = []
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


def forest_poly_removed(adj: list[list[int]], w: int) -> list[int]:
    """Compute I(T - w) by removing vertex w and multiplying component polys."""
    n = len(adj)
    remaining = [i for i in range(n) if i != w]
    if not remaining:
        return [1]
    rem_set = set(remaining)
    seen = set()
    out = [1]
    for start in remaining:
        if start in seen:
            continue
        # BFS component
        comp = []
        stack = [start]
        seen.add(start)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v in rem_set and v not in seen:
                    seen.add(v)
                    stack.append(v)
        # Build adjacency for component
        mapping = {old: i for i, old in enumerate(comp)}
        cadj: list[list[int]] = [[] for _ in range(len(comp))]
        for old in comp:
            i = mapping[old]
            for v in adj[old]:
                j = mapping.get(v)
                if j is not None:
                    cadj[i].append(j)
        poly_c = independence_poly(len(comp), cadj)
        out = _polymul(out, poly_c)
    return out


def check_tree(n: int, adj: list[list[int]], g6: str) -> dict[str, Any]:
    """Check mode(I(T-w)) <= d(I(T)) for all vertices w in one tree.

    Returns a dict with per-tree statistics.
    """
    f = independence_poly(n, adj)
    mode_f = mode_index(f)
    d_f = first_descent(f)

    vertex_cases = 0
    failures = 0
    equality_count = 0
    max_mode_gap = -(10**9)   # max(mode(I(T-w)) - d(I(T)))
    max_abs_mode_diff = 0     # max |mode(I(T-w)) - mode(I(T))|
    tight_examples: list[dict] = []

    for w in range(n):
        vertex_cases += 1
        g = forest_poly_removed(adj, w)
        mode_g = mode_index(g)

        gap = mode_g - d_f
        if gap > max_mode_gap:
            max_mode_gap = gap

        abs_diff = abs(mode_g - mode_f)
        if abs_diff > max_abs_mode_diff:
            max_abs_mode_diff = abs_diff

        if gap > 0:
            failures += 1

        if gap == 0:
            equality_count += 1
            if len(tight_examples) < 3:
                tight_examples.append({
                    "graph6": g6,
                    "vertex": w,
                    "mode_Tw": mode_g,
                    "d_IT": d_f,
                    "mode_IT": mode_f,
                })

    return {
        "vertex_cases": vertex_cases,
        "failures": failures,
        "equality_count": equality_count,
        "max_mode_gap": max_mode_gap,
        "max_abs_mode_diff": max_abs_mode_diff,
        "tight_examples": tight_examples,
    }


def worker_partition(args: tuple[int, int, int]) -> dict[str, Any]:
    """Worker: enumerate and check one geng partition for a given n."""
    n, res, mod = args
    tree_count = 0
    vertex_cases = 0
    failures = 0
    equality_count = 0
    max_mode_gap = -(10**9)
    max_abs_mode_diff = 0
    tight_examples: list[dict] = []

    for tree_n, adj, raw in trees_geng_raw(n, res=res, mod=mod):
        tree_count += 1
        g6 = raw.decode("ascii", errors="replace")
        result = check_tree(tree_n, adj, g6)

        vertex_cases += result["vertex_cases"]
        failures += result["failures"]
        equality_count += result["equality_count"]
        if result["max_mode_gap"] > max_mode_gap:
            max_mode_gap = result["max_mode_gap"]
        if result["max_abs_mode_diff"] > max_abs_mode_diff:
            max_abs_mode_diff = result["max_abs_mode_diff"]
        # Keep up to 10 tight examples total across all workers
        if len(tight_examples) < 10:
            tight_examples.extend(result["tight_examples"][:10 - len(tight_examples)])

    return {
        "tree_count": tree_count,
        "vertex_cases": vertex_cases,
        "failures": failures,
        "equality_count": equality_count,
        "max_mode_gap": max_mode_gap,
        "max_abs_mode_diff": max_abs_mode_diff,
        "tight_examples": tight_examples,
    }


# Expected tree counts (OEIS A000055)
A000055 = {
    19: 317955,
    20: 823065,
    21: 2144505,
    22: 5623756,
}


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Extended mode-alignment check: mode(I(T-w)) <= d(I(T))"
    )
    ap.add_argument("--min-n", type=int, default=19)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument(
        "--out",
        default="",
        help="Output JSON path (default: results/mode_alignment_n{max_n}.json)",
    )
    args = ap.parse_args()

    if not args.out:
        args.out = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "results",
            f"mode_alignment_n{args.max_n}.json",
        )

    print(f"Mode-alignment check: mode(I(T-w)) <= d(I(T))")
    print(f"Range: n = {args.min_n} to {args.max_n}, workers = {args.workers}")
    print()

    grand_start = time.time()
    total_trees = 0
    total_vertex_cases = 0
    total_failures = 0
    total_equality = 0
    global_max_mode_gap = -(10**9)
    global_max_abs_mode_diff = 0
    all_tight_examples: list[dict] = []
    per_n_counts: dict[str, Any] = {}

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        expected_str = f"{A000055[n]:,}" if n in A000055 else "?"
        print(f"n={n}: enumerating {expected_str} trees with {args.workers} workers...")

        tasks = [(n, res, args.workers) for res in range(args.workers)]

        n_trees = 0
        n_vertex_cases = 0
        n_failures = 0
        n_equality = 0
        n_max_mode_gap = -(10**9)
        n_max_abs_mode_diff = 0
        n_tight: list[dict] = []

        with Pool(args.workers) as pool:
            for result in pool.imap_unordered(worker_partition, tasks):
                n_trees += result["tree_count"]
                n_vertex_cases += result["vertex_cases"]
                n_failures += result["failures"]
                n_equality += result["equality_count"]
                if result["max_mode_gap"] > n_max_mode_gap:
                    n_max_mode_gap = result["max_mode_gap"]
                if result["max_abs_mode_diff"] > n_max_abs_mode_diff:
                    n_max_abs_mode_diff = result["max_abs_mode_diff"]
                if len(n_tight) < 10:
                    n_tight.extend(result["tight_examples"][:10 - len(n_tight)])

        elapsed = time.time() - t0

        # Verify tree count
        expected = A000055.get(n)
        if expected is not None and n_trees != expected:
            print(f"  WARNING: tree count mismatch! got {n_trees:,}, expected {expected:,}")

        status = "PASS" if n_failures == 0 else f"FAIL ({n_failures} failures)"
        print(
            f"  {n_trees:>10,} trees, {n_vertex_cases:>12,} vertex-cases, "
            f"max_gap={n_max_mode_gap}, equality={n_equality:,}, "
            f"max|diff|={n_max_abs_mode_diff} [{status}] ({elapsed:.1f}s)"
        )

        total_trees += n_trees
        total_vertex_cases += n_vertex_cases
        total_failures += n_failures
        total_equality += n_equality
        if n_max_mode_gap > global_max_mode_gap:
            global_max_mode_gap = n_max_mode_gap
        if n_max_abs_mode_diff > global_max_abs_mode_diff:
            global_max_abs_mode_diff = n_max_abs_mode_diff
        if len(all_tight_examples) < 10:
            all_tight_examples.extend(n_tight[:10 - len(all_tight_examples)])

        per_n_counts[str(n)] = {
            "trees": n_trees,
            "vertex_cases": n_vertex_cases,
            "failures": n_failures,
            "equality_count": n_equality,
            "max_mode_gap": n_max_mode_gap,
            "max_abs_mode_diff": n_max_abs_mode_diff,
        }

    grand_elapsed = time.time() - grand_start

    output = {
        "claim": "mode(I(T-w)) <= d(I(T)) for all trees T and vertices w",
        "max_n": args.max_n,
        "min_n": args.min_n,
        "total_trees": total_trees,
        "total_vertex_cases": total_vertex_cases,
        "max_mode_gap": global_max_mode_gap,
        "max_abs_mode_diff": global_max_abs_mode_diff,
        "equality_count": total_equality,
        "failures": total_failures,
        "verified": total_failures == 0,
        "tightest_examples": all_tight_examples[:10],
        "per_n_counts": per_n_counts,
        "wall_time_seconds": round(grand_elapsed, 1),
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {args.out}")

    print(f"\nSummary:")
    print(f"  Trees checked: {total_trees:,}")
    print(f"  Vertex-cases: {total_vertex_cases:,}")
    print(f"  Max mode gap (mode(I(T-w)) - d(I(T))): {global_max_mode_gap}")
    print(f"  Max |mode(I(T-w)) - mode(I(T))|: {global_max_abs_mode_diff}")
    print(f"  Equality cases: {total_equality:,}")
    print(f"  Failures: {total_failures}")
    print(f"  Wall time: {grand_elapsed:.1f}s")

    if total_failures == 0:
        print(f"\n  Claim VERIFIED for n <= {args.max_n}.")
    else:
        print(f"\n  CLAIM FAILS for n <= {args.max_n}! See failures above.")


if __name__ == "__main__":
    main()
