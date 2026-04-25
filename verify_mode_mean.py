#!/usr/bin/env python3
"""Verify Conjecture mode(I(T)) ≤ ⌈μ(T)⌉ for all trees through n=max_n.

Uses the same geng/indpoly infrastructure as search.py.
For each tree, computes the full IS polynomial and checks:
  mode ≤ ceil(μ)   where μ = I'(1)/I(1)
"""

import argparse
import shutil
import time
from fractions import Fraction
from multiprocessing import Pool

from indpoly import independence_poly
from search import A000055_COUNTS
from trees import trees, trees_geng


def check_mode_mean(poly: list[int]) -> tuple[bool, int, int, int]:
    """Check mode ≤ ceil(μ) for a polynomial coefficient list.

    Returns (passes, mode, weighted, total), so μ=weighted/total exactly.
    """
    # mode = argmax
    mode = max(range(len(poly)), key=lambda k: poly[k])

    # μ = I'(1)/I(1) = sum(k * i_k) / sum(i_k)
    total = sum(poly)
    if total == 0:
        return True, mode, 0, 1
    weighted = sum(k * c for k, c in enumerate(poly))
    ceil_mu = (weighted + total - 1) // total

    return mode <= ceil_mu, mode, weighted, total


def worker_partition(args):
    """Worker: check one geng partition for mode-mean."""
    n, res, mod = args
    count = 0
    violations = []
    worst_margin = float('inf')  # smallest ceil(mu) - mode
    for tree_n, adj in trees_geng(n, res=res, mod=mod):
        poly = independence_poly(tree_n, adj)
        passes, mode, weighted, total = check_mode_mean(poly)
        ceil_mu = (weighted + total - 1) // total
        margin = ceil_mu - mode
        if margin < worst_margin:
            worst_margin = margin
        if not passes:
            violations.append({
                'n': tree_n, 'adj': adj, 'poly': poly,
                'mode': mode, 'mu': [weighted, total]
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
                passes, mode, weighted, total = check_mode_mean(poly)
                ceil_mu = (weighted + total - 1) // total
                margin = ceil_mu - mode
                if margin < worst:
                    worst = margin
                if not passes:
                    violations.append({
                        'n': tree_n, 'mode': mode,
                        'mu': [weighted, total],
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

        expected = A000055_COUNTS.get(n)
        if expected is not None and count != expected:
            print(f"  WARNING: count mismatch at n={n}: {count} vs {expected}")

        status = "OK" if not violations else f"*** {len(violations)} VIOLATIONS ***"
        margin_str = f"min margin={worst}" if worst != float('inf') else ""
        print(f"n={n:>3}: {count:>12,} trees, {elapsed:>8.2f}s  {status}  {margin_str}")

        if violations:
            for v in violations[:5]:
                mu = Fraction(v["mu"][0], v["mu"][1])
                ceil_mu = (mu.numerator + mu.denominator - 1) // mu.denominator
                print(f"  VIOLATION: mode={v['mode']}, μ={mu} "
                      f"(≈{float(mu):.6f}), ⌈μ⌉={ceil_mu}")

    grand_elapsed = time.time() - grand_start
    print(f"\nDone. {grand_total:,} trees checked in {grand_elapsed:.1f}s.")
    if grand_violations == 0:
        print(f"Conjecture mode ≤ ⌈μ⌉ holds for all trees on n ≤ {args.max_n}.")
    else:
        print(f"*** {grand_violations} total violations found! ***")


if __name__ == "__main__":
    main()
