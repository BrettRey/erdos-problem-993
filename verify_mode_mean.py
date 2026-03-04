#!/usr/bin/env python3
"""Verify Conjecture mode(I(T)) ≤ ⌈μ(T)⌉ for all trees through n=max_n.

Uses the same geng/indpoly infrastructure as search.py.
For each tree, computes the full IS polynomial and checks:
  mode ≤ ceil(μ)   where μ = I'(1)/I(1)
"""

import argparse
import math
import shutil
import time
from multiprocessing import Pool

from graph6 import parse_graph6
from indpoly import independence_poly
from trees import trees, trees_geng

# OEIS A000055 counts for verification
A000055 = {
    1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 11, 8: 23, 9: 47, 10: 106,
    11: 235, 12: 551, 13: 1301, 14: 3159, 15: 7741, 16: 19320, 17: 48629,
    18: 123867, 19: 317955, 20: 823065, 21: 2144505, 22: 5623756,
    23: 14828074, 24: 39299897, 25: 104636890, 26: 279793450, 27: 751065460,
}


def check_mode_mean(poly: list[int]) -> tuple[bool, int, float]:
    """Check mode ≤ ceil(μ) for a polynomial coefficient list.

    Returns (passes, mode, mu).
    """
    # mode = argmax
    mode = max(range(len(poly)), key=lambda k: poly[k])

    # μ = I'(1)/I(1) = sum(k * i_k) / sum(i_k)
    total = sum(poly)
    if total == 0:
        return True, mode, 0.0
    weighted = sum(k * c for k, c in enumerate(poly))
    mu = weighted / total

    return mode <= math.ceil(mu), mode, mu


def worker_partition(args):
    """Worker: check one geng partition for mode-mean."""
    n, res, mod = args
    count = 0
    violations = []
    worst_margin = float('inf')  # smallest ceil(mu) - mode
    for tree_n, adj in trees_geng(n, res=res, mod=mod):
        poly = independence_poly(tree_n, adj)
        passes, mode, mu = check_mode_mean(poly)
        margin = math.ceil(mu) - mode
        if margin < worst_margin:
            worst_margin = margin
        if not passes:
            violations.append({
                'n': tree_n, 'adj': adj, 'poly': poly,
                'mode': mode, 'mu': mu
            })
        count += 1
    return count, violations, worst_margin


def main():
    ap = argparse.ArgumentParser(description="Verify mode ≤ ⌈μ⌉ for all trees")
    ap.add_argument("--min-n", type=int, default=1)
    ap.add_argument("--max-n", type=int, default=27)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    print(f"Verifying Conjecture mode(I(T)) ≤ ⌈μ(T)⌉")
    print(f"Range: n = {args.min_n} to {args.max_n}, workers = {args.workers}")
    print()

    grand_total = 0
    grand_violations = 0
    grand_start = time.time()

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()

        if n <= 2 or not shutil.which("geng"):
            # Single-process for tiny n or no geng
            count = 0
            violations = []
            worst = float('inf')
            for tree_n, adj in trees(n, backend="auto"):
                poly = independence_poly(tree_n, adj)
                passes, mode, mu = check_mode_mean(poly)
                margin = math.ceil(mu) - mode
                if margin < worst:
                    worst = margin
                if not passes:
                    violations.append({
                        'n': tree_n, 'mode': mode, 'mu': mu
                    })
                count += 1
        else:
            tasks = [(n, res, args.workers) for res in range(args.workers)]
            count = 0
            violations = []
            worst = float('inf')
            with Pool(args.workers) as pool:
                for c, v, w in pool.imap_unordered(worker_partition, tasks):
                    count += c
                    violations.extend(v)
                    if w < worst:
                        worst = w

        elapsed = time.time() - t0
        grand_total += count
        grand_violations += len(violations)

        expected = A000055.get(n)
        if expected is not None and count != expected:
            print(f"  WARNING: count mismatch at n={n}: {count} vs {expected}")

        status = "OK" if not violations else f"*** {len(violations)} VIOLATIONS ***"
        margin_str = f"min margin={worst}" if worst != float('inf') else ""
        print(f"n={n:>3}: {count:>12,} trees, {elapsed:>8.2f}s  {status}  {margin_str}")

        if violations:
            for v in violations[:5]:
                print(f"  VIOLATION: mode={v['mode']}, μ={v['mu']:.6f}, "
                      f"⌈μ⌉={math.ceil(v['mu'])}")

    grand_elapsed = time.time() - grand_start
    print(f"\nDone. {grand_total:,} trees checked in {grand_elapsed:.1f}s.")
    if grand_violations == 0:
        print(f"Conjecture mode ≤ ⌈μ⌉ holds for all trees on n ≤ {args.max_n}.")
    else:
        print(f"*** {grand_violations} total violations found! ***")


if __name__ == "__main__":
    main()
