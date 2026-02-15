#!/usr/bin/env python3
"""Analyze the structure of high-mode trees to understand why they exclude low-priv IS.

A "high-mode" tree has mode > floor(n/3) + 1. For these trees, the 1-Private bound
k >= (n+1)/3 is insufficient (mode exceeds the bound). We need to prove that no
maximal IS with all priv <= 1 exists in such trees.

This script characterizes high-mode trees:
- What do they look like? (degree sequences, diameter, etc.)
- Why don't they admit low-priv IS?
- Is there a structural invariant we can prove?
"""

import subprocess
import sys
import time
from collections import Counter

sys.path.insert(0, ".")
from indpoly import independence_poly


def parse_graph6(s):
    s = s.strip()
    if not s:
        return 0, []
    idx = 0
    if ord(s[0]) - 63 < 63:
        n = ord(s[0]) - 63
        idx = 1
    else:
        idx = 1
        n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63)
            idx += 1
    adj = [[] for _ in range(n)]
    bits = []
    for c in s[idx:]:
        val = ord(c) - 63
        for b in range(5, -1, -1):
            bits.append((val >> b) & 1)
    bit_idx = 0
    for j in range(n):
        for i in range(j):
            if bit_idx < len(bits) and bits[bit_idx]:
                adj[i].append(j)
                adj[j].append(i)
            bit_idx += 1
    return n, adj


def tree_diameter(n, adj):
    """BFS-based diameter computation."""
    def bfs_farthest(start):
        dist = [-1] * n
        dist[start] = 0
        queue = [start]
        farthest = start
        for v in queue:
            for w in adj[v]:
                if dist[w] == -1:
                    dist[w] = dist[v] + 1
                    queue.append(w)
                    if dist[w] > dist[farthest]:
                        farthest = w
        return farthest, dist[farthest]
    far1, _ = bfs_farthest(0)
    far2, diam = bfs_farthest(far1)
    return diam


def degree_sequence(n, adj):
    return sorted([len(adj[v]) for v in range(n)], reverse=True)


def count_leaves(n, adj):
    return sum(1 for v in range(n) if len(adj[v]) == 1)


def max_degree(n, adj):
    return max(len(adj[v]) for v in range(n))


def main():
    print("HIGH-MODE TREE STRUCTURAL ANALYSIS")
    print("=" * 70)
    print()

    # Collect statistics on high-mode vs non-high-mode trees
    for n in range(8, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 2  # mode > floor(n/3) + 1

        n_trees = len(lines)
        n_high = 0

        # Stats for high-mode trees
        high_diameters = []
        high_max_degrees = []
        high_leaf_counts = []
        high_modes = []
        high_degree_seqs = []

        # Stats for all trees
        all_diameters = []
        all_modes = []

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            all_modes.append(mode)

            diam = tree_diameter(tn, adj_data)
            all_diameters.append(diam)

            if mode >= threshold:
                n_high += 1
                high_diameters.append(diam)
                high_max_degrees.append(max_degree(tn, adj_data))
                high_leaf_counts.append(count_leaves(tn, adj_data))
                high_modes.append(mode)
                if n_high <= 3 and n <= 12:
                    ds = degree_sequence(tn, adj_data)
                    high_degree_seqs.append(ds)

        elapsed = time.time() - t0

        avg_mode = sum(all_modes) / len(all_modes) if all_modes else 0
        print(f"n={n:2d}: {n_trees:6d} trees, threshold={threshold}, "
              f"high_mode={n_high:5d} ({100*n_high/max(1,n_trees):.1f}%), "
              f"avg_mode={avg_mode:.2f} ({elapsed:.1f}s)")

        if high_diameters:
            avg_diam = sum(high_diameters) / len(high_diameters)
            min_diam = min(high_diameters)
            max_diam_val = max(high_diameters)
            avg_maxdeg = sum(high_max_degrees) / len(high_max_degrees)
            min_maxdeg = min(high_max_degrees)
            avg_leaves = sum(high_leaf_counts) / len(high_leaf_counts)
            min_leaves = min(high_leaf_counts)
            max_leaves = max(high_leaf_counts)
            print(f"  HIGH-MODE: diam=[{min_diam},{max_diam_val}] avg={avg_diam:.1f}, "
                  f"maxdeg=[{min_maxdeg},...] avg={avg_maxdeg:.1f}, "
                  f"leaves=[{min_leaves},{max_leaves}] avg={avg_leaves:.1f}, "
                  f"modes={Counter(high_modes).most_common(5)}")

            # Key question: do high-mode trees have high max degree?
            # (If so, the IS vertex at the high-degree hub has many privs)
            pct_maxdeg_ge3 = sum(1 for d in high_max_degrees if d >= 3) / len(high_max_degrees)
            pct_maxdeg_ge4 = sum(1 for d in high_max_degrees if d >= 4) / len(high_max_degrees)
            print(f"  maxdeg>=3: {100*pct_maxdeg_ge3:.1f}%, maxdeg>=4: {100*pct_maxdeg_ge4:.1f}%")

            # Key question: do high-mode trees have short diameter?
            pct_diam_le_n_half = sum(1 for d in high_diameters if d <= n // 2) / len(high_diameters)
            print(f"  diam <= n/2: {100*pct_diam_le_n_half:.1f}%")

            if high_degree_seqs:
                for i, ds in enumerate(high_degree_seqs[:3]):
                    print(f"  example degree seq: {ds}")

    print()
    print("=" * 70)
    print("KEY QUESTION: What structural property distinguishes high-mode trees?")
    print("If high-mode trees are 'star-like' (high max degree, many leaves),")
    print("then any maximal IS must include the hub, which has many private")
    print("neighbors, making all priv <= 1 impossible.")
    print("=" * 70)


if __name__ == "__main__":
    main()
