#!/usr/bin/env python3
"""Count homeomorphically irreducible (series-reduced) trees.

A tree is HI if it has no vertex of degree exactly 2.
If the subdivision lemma holds, any minimal counterexample must be HI.

Combined with the PNP framework (Conjecture A reduction):
- A minimal counterexample must also satisfy d_leaf(v) <= 1 for all v.
- An HI tree with d_leaf <= 1 has every internal vertex of degree >= 3
  with at most 1 leaf child.

This script counts HI trees and HI+d_leaf<=1 trees through various n,
and verifies unimodality/LC for all of them.
"""
import subprocess
import sys
import time
from collections import Counter

from indpoly import independence_poly, is_log_concave, is_unimodal

GENG = "/opt/homebrew/bin/geng"


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode("ascii")]
    n = data[0]
    idx = 1
    adj = [[] for _ in range(n)]
    bit = 5
    word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5
                idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj


def is_hi(n, adj):
    """Check if tree is homeomorphically irreducible (no degree-2 vertex)."""
    for v in range(n):
        if len(adj[v]) == 2:
            return False
    return True


def d_leaf_max(n, adj):
    """Max d_leaf(v) = number of leaf-children of any vertex."""
    mx = 0
    for v in range(n):
        count = sum(1 for u in adj[v] if len(adj[u]) == 1)
        if count > mx:
            mx = count
    return mx


def first_descent(seq):
    for k in range(len(seq) - 1):
        if seq[k] > seq[k + 1]:
            return k
    return len(seq) - 1


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    print(f"Counting HI trees and HI+d_leaf<=1 trees, n up to {max_n}", flush=True)
    print("=" * 90, flush=True)
    print(f"{'n':>3} {'trees':>12} {'HI':>10} {'HI%':>7} "
          f"{'HI_dl1':>10} {'HI_dl1%':>8} {'uni_fail':>8} {'lc_fail':>8} {'mode':>6}",
          flush=True)
    print("-" * 90, flush=True)

    t0 = time.time()
    total_trees = 0
    total_hi = 0
    total_hi_dl1 = 0
    total_uni_fail = 0
    total_lc_fail = 0

    for n in range(2, max_n + 1):
        tn = time.time()
        n_trees = 0
        n_hi = 0
        n_hi_dl1 = 0
        n_uni_fail = 0
        n_lc_fail = 0
        mode_dist = Counter()

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            n_trees += 1

            if not is_hi(nn, adj):
                continue

            n_hi += 1
            dl = d_leaf_max(nn, adj)

            if dl <= 1:
                n_hi_dl1 += 1
                # Verify unimodality/LC for these constrained trees
                I = independence_poly(nn, adj)
                d = first_descent(I)
                mode_dist[d] += 1

                if not is_unimodal(I):
                    n_uni_fail += 1
                if not is_log_concave(I):
                    n_lc_fail += 1

        proc.wait()
        elapsed = time.time() - tn
        total_trees += n_trees
        total_hi += n_hi
        total_hi_dl1 += n_hi_dl1
        total_uni_fail += n_uni_fail
        total_lc_fail += n_lc_fail

        hi_pct = 100 * n_hi / n_trees if n_trees > 0 else 0
        dl1_pct = 100 * n_hi_dl1 / n_trees if n_trees > 0 else 0
        mode_str = f"{max(mode_dist.keys()) if mode_dist else 0}"

        print(f"{n:3d} {n_trees:12,} {n_hi:10,} {hi_pct:6.2f}% "
              f"{n_hi_dl1:10,} {dl1_pct:7.3f}% {n_uni_fail:8d} {n_lc_fail:8d} "
              f"≤{mode_str:>4}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 90, flush=True)
    print(f"Total trees: {total_trees:,}", flush=True)
    print(f"Total HI: {total_hi:,} ({100*total_hi/max(total_trees,1):.2f}%)", flush=True)
    print(f"Total HI+d_leaf<=1: {total_hi_dl1:,} ({100*total_hi_dl1/max(total_trees,1):.3f}%)",
          flush=True)
    print(f"Unimodality failures in HI+d_leaf<=1: {total_uni_fail}", flush=True)
    print(f"Log-concavity failures in HI+d_leaf<=1: {total_lc_fail}", flush=True)
    print(f"Time: {total_time:.1f}s", flush=True)

    if total_uni_fail == 0:
        print(flush=True)
        print("ALL HI + d_leaf<=1 trees are unimodal!", flush=True)
        if total_lc_fail == 0:
            print("ALL HI + d_leaf<=1 trees are also log-concave!", flush=True)


if __name__ == "__main__":
    main()
