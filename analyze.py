#!/usr/bin/env python3
"""Analyze trees for log-concavity failures and unimodality near-misses.

Re-scans all trees on n vertices, collecting:
  1. Every tree whose independence polynomial fails log-concavity
  2. The top-K trees closest to violating unimodality (by near-miss ratio)

Results saved to results/analysis_n{N}.json.
"""

import argparse
import heapq
import json
import os
import time
from multiprocessing import Pool

from indpoly import (
    independence_poly,
    is_log_concave,
    is_unimodal,
    log_concavity_ratio,
    near_miss_ratio,
)
from trees import trees_geng_raw


def _worker(args: tuple[int, int, int, int]) -> dict:
    """Analyze one geng partition.

    Returns a dict with:
      - count: total trees processed
      - lc_failures: list of log-concavity failures
      - top_near_misses: list of (neg_ratio, data) for heap merging
      - non_unimodal: list of any unimodality failures (the prize)
      - stats: summary statistics
    """
    n, res, mod, top_k = args
    count = 0
    lc_failures = []
    non_unimodal = []
    # Min-heap of (-ratio, data) to keep top-K largest ratios
    near_miss_heap: list[tuple[float, dict]] = []

    for tree_n, adj, g6 in trees_geng_raw(n, res=res, mod=mod):
        poly = independence_poly(tree_n, adj)
        count += 1

        # Check unimodality
        uni = is_unimodal(poly)
        if not uni:
            non_unimodal.append({
                "graph6": g6.decode(),
                "poly": poly,
                "n": tree_n,
            })
            continue  # also a near-miss and lc failure, but save separately

        # Check log-concavity
        lc = is_log_concave(poly)
        if not lc:
            lc_ratio, lc_pos = log_concavity_ratio(poly)
            lc_failures.append({
                "graph6": g6.decode(),
                "poly": poly,
                "n": tree_n,
                "lc_ratio": lc_ratio,
                "lc_pos": lc_pos,
            })

        # Near-miss ratio
        nm_ratio, nm_pos = near_miss_ratio(poly)
        if nm_ratio > 0:
            entry = {
                "graph6": g6.decode(),
                "poly": poly,
                "n": tree_n,
                "nm_ratio": nm_ratio,
                "nm_pos": nm_pos,
            }
            if len(near_miss_heap) < top_k:
                heapq.heappush(near_miss_heap, (nm_ratio, count, entry))
            elif nm_ratio > near_miss_heap[0][0]:
                heapq.heapreplace(near_miss_heap, (nm_ratio, count, entry))

    top_near_misses = sorted(near_miss_heap, key=lambda x: -x[0])
    return {
        "count": count,
        "lc_failures": lc_failures,
        "non_unimodal": non_unimodal,
        "top_near_misses": [(r, d) for r, _, d in top_near_misses],
    }


def analyze_n(n: int, workers: int, top_k: int = 100) -> dict:
    """Run full analysis on all trees of order n."""
    tasks = [(n, res, workers, top_k) for res in range(workers)]

    total_count = 0
    all_lc_failures = []
    all_non_unimodal = []
    # Merge heaps across workers
    merged_heap: list[tuple[float, int, dict]] = []
    merge_counter = 0

    t0 = time.time()
    with Pool(workers) as pool:
        for result in pool.imap_unordered(_worker, tasks):
            total_count += result["count"]
            all_lc_failures.extend(result["lc_failures"])
            all_non_unimodal.extend(result["non_unimodal"])

            for ratio, data in result["top_near_misses"]:
                merge_counter += 1
                if len(merged_heap) < top_k:
                    heapq.heappush(merged_heap, (ratio, merge_counter, data))
                elif ratio > merged_heap[0][0]:
                    heapq.heapreplace(merged_heap, (ratio, merge_counter, data))

    elapsed = time.time() - t0

    # Sort results
    all_lc_failures.sort(key=lambda x: -x["lc_ratio"])
    top_near_misses = sorted(merged_heap, key=lambda x: -x[0])
    top_near_misses = [{"nm_ratio": r, **d} for r, _, d in top_near_misses]

    return {
        "n": n,
        "total_trees": total_count,
        "elapsed_s": round(elapsed, 2),
        "non_unimodal_count": len(all_non_unimodal),
        "non_unimodal": all_non_unimodal,
        "lc_failure_count": len(all_lc_failures),
        "lc_failures": all_lc_failures,
        "top_near_misses": top_near_misses,
    }


def print_summary(result: dict):
    n = result["n"]
    total = result["total_trees"]
    elapsed = result["elapsed_s"]

    print(f"\n{'='*70}")
    print(f"Analysis of n = {n}: {total:,} trees in {elapsed:.1f}s")
    print(f"{'='*70}")

    # Unimodality violations
    if result["non_unimodal"]:
        print(f"\n*** {result['non_unimodal_count']} UNIMODALITY VIOLATION(S) ***")
        for item in result["non_unimodal"]:
            print(f"  graph6: {item['graph6']}")
            print(f"  poly:   {item['poly']}")
    else:
        print(f"\nUnimodality: ALL {total:,} trees are unimodal.")

    # Log-concavity failures
    lc_count = result["lc_failure_count"]
    print(f"\nLog-concavity failures: {lc_count}")
    if lc_count > 0:
        print(f"  Worst violation ratio: {result['lc_failures'][0]['lc_ratio']:.10f}")
        print(f"  Showing top 10:")
        for item in result["lc_failures"][:10]:
            alpha = len(item["poly"]) - 1
            print(
                f"    ratio={item['lc_ratio']:.10f} at k={item['lc_pos']} "
                f"(alpha={alpha})  graph6={item['graph6']}"
            )

    # Near-misses
    print(f"\nTop 10 unimodality near-misses (ratio a_{{j+1}}/a_j in descending tail):")
    for item in result["top_near_misses"][:10]:
        alpha = len(item["poly"]) - 1
        print(
            f"  ratio={item['nm_ratio']:.10f} at j={item['nm_pos']} "
            f"(alpha={alpha})  graph6={item['graph6']}"
        )
        # Show the relevant part of the polynomial
        j = item["nm_pos"]
        lo = max(0, j - 2)
        hi = min(len(item["poly"]), j + 3)
        excerpt = item["poly"][lo:hi]
        labels = list(range(lo, hi))
        print(f"    poly[{lo}..{hi-1}]: {dict(zip(labels, excerpt))}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze trees for log-concavity and near-miss metrics"
    )
    parser.add_argument("n", type=int, help="Number of vertices")
    parser.add_argument(
        "--workers", type=int, default=8, help="Parallel workers (default: 8)"
    )
    parser.add_argument(
        "--top-k", type=int, default=100, help="Keep top K near-misses (default: 100)"
    )
    args = parser.parse_args()

    print(f"Analyzing all trees on n = {args.n} vertices...")
    print(f"Workers: {args.workers}, Top-K: {args.top_k}")

    result = analyze_n(args.n, args.workers, args.top_k)

    # Save
    os.makedirs("results", exist_ok=True)
    path = f"results/analysis_n{args.n}.json"
    with open(path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nFull results saved to {path}")

    print_summary(result)


if __name__ == "__main__":
    main()
