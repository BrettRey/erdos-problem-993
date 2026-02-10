#!/usr/bin/env python3
"""Search for a counterexample to Erdős Problem #993.

Enumerates all non-isomorphic trees up to a given vertex count, computes
their independence polynomials, and checks unimodality. Any non-unimodal
tree is a counterexample to the conjecture.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool

from graph6 import parse_graph6
from indpoly import independence_poly, is_unimodal
from trees import trees, trees_geng


def check_tree(args: tuple[int, list[list[int]]]) -> tuple[bool, list[int]] | None:
    """Check a single tree. Returns (False, poly) if non-unimodal, else None."""
    n, adj = args
    poly = independence_poly(n, adj)
    if not is_unimodal(poly):
        return (False, poly)
    return None


def check_trees_single(n: int, backend: str) -> tuple[int, list | None]:
    """Check all trees on n vertices, single-process.

    Returns (count, counterexample_or_None).
    """
    count = 0
    for tree_n, adj in trees(n, backend=backend):
        poly = independence_poly(tree_n, adj)
        if not is_unimodal(poly):
            return count + 1, {
                "n": tree_n,
                "adj": adj,
                "poly": poly,
            }
        count += 1
    return count, None


def _worker_geng_partition(args: tuple[int, int, int]) -> tuple[int, dict | None]:
    """Worker: enumerate and check one geng partition."""
    n, res, mod = args
    count = 0
    for tree_n, adj in trees_geng(n, res=res, mod=mod):
        poly = independence_poly(tree_n, adj)
        if not is_unimodal(poly):
            return count + 1, {
                "n": tree_n,
                "adj": adj,
                "poly": poly,
            }
        count += 1
    return count, None


def check_trees_parallel(n: int, workers: int) -> tuple[int, dict | None]:
    """Check all trees on n vertices using parallel geng partitions.

    Requires geng on PATH. Falls back to single-process networkx if unavailable.
    """
    if not shutil.which("geng"):
        print(f"  (geng not found, falling back to single-process networkx)")
        return check_trees_single(n, backend="networkx")

    tasks = [(n, res, workers) for res in range(workers)]
    total = 0
    with Pool(workers) as pool:
        for count, counterexample in pool.imap_unordered(_worker_geng_partition, tasks):
            total += count
            if counterexample is not None:
                pool.terminate()
                return total, counterexample
    return total, None


def save_counterexample(result: dict, directory: str = "results") -> str:
    """Save a counterexample to a JSON file. Returns the file path."""
    os.makedirs(directory, exist_ok=True)
    n = result["n"]
    ts = int(time.time())
    path = os.path.join(directory, f"counterexample_n{n}_{ts}.json")
    with open(path, "w") as f:
        json.dump(result, f, indent=2)
    return path


def main():
    parser = argparse.ArgumentParser(
        description="Search for counterexamples to Erdős Problem #993"
    )
    parser.add_argument(
        "--min-n", type=int, default=1, help="Minimum vertex count (default: 1)"
    )
    parser.add_argument(
        "--max-n", type=int, default=20, help="Maximum vertex count (default: 20)"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1; >1 requires geng)",
    )
    parser.add_argument(
        "--backend",
        choices=["auto", "geng", "networkx"],
        default="auto",
        help="Tree enumeration backend (default: auto)",
    )
    args = parser.parse_args()

    print(f"Erdős Problem #993: searching for non-unimodal tree independence sequences")
    print(f"Range: n = {args.min_n} to {args.max_n}, workers = {args.workers}")
    print(f"Backend: {args.backend}")
    print()

    grand_total = 0
    grand_start = time.time()

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()

        if args.workers > 1:
            count, counterexample = check_trees_parallel(n, args.workers)
        else:
            count, counterexample = check_trees_single(n, backend=args.backend)

        elapsed = time.time() - t0
        grand_total += count

        if counterexample is not None:
            print(f"\n*** COUNTEREXAMPLE FOUND at n = {n} ***")
            print(f"Independence polynomial: {counterexample['poly']}")
            print(f"Adjacency list: {counterexample['adj']}")
            path = save_counterexample(counterexample)
            print(f"Saved to: {path}")
            print(f"\nThis disproves the conjecture!")
            sys.exit(0)

        print(
            f"n={n:>3}: {count:>12,} trees checked in {elapsed:>8.2f}s "
            f"— all unimodal"
        )

    grand_elapsed = time.time() - grand_start
    print(f"\nDone. {grand_total:,} trees checked in {grand_elapsed:.2f}s. "
          f"All unimodal (conjecture holds for n ≤ {args.max_n}).")


if __name__ == "__main__":
    main()
