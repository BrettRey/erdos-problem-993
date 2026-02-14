#!/usr/bin/env python3
"""
Check Hall's degree condition for the augmented bipartite graph.

Hall's degree condition: If every left vertex has degree ≥ d and every
right vertex has degree ≤ d, then Hall's condition holds.

In the containment-only graph:
  Right degree = exactly k+1 (remove each vertex of T to get S ⊂ T)
  Left degree = number of free vertices (0 for maximal IS)

In the augmented graph:
  Right degree = k+1 + (swap neighbors from right side)
  Left degree = containment + swap degree

Check: does left min-degree ≥ right max-degree in the augmented graph?
If so, Hall's degree condition gives the matching for free.

Also check: left min-degree ≥ k+1? (since right has degree ≤ k+1 + swaps)
"""

import subprocess
import sys
import time
from collections import defaultdict

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


def enumerate_independent_sets(n, adj):
    nbr = [set(adj[v]) for v in range(n)]
    levels = defaultdict(list)
    def backtrack(v, current, forbidden):
        if v == n:
            levels[len(current)].append(frozenset(current))
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])
    backtrack(0, [], set())
    return levels


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def analyze_degrees(levels, k, n, adj):
    """Compute left and right degree distributions in augmented graph."""
    left = levels.get(k, [])
    right = levels.get(k + 1, [])
    if not left or not right:
        return None

    nbr = [set(adj[v]) for v in range(n)]
    left_idx = {s: i for i, s in enumerate(left)}
    right_idx = {s: j for j, s in enumerate(right)}
    n_left = len(left)
    n_right = len(right)

    left_deg = [0] * n_left
    right_deg = [0] * n_right

    # Containment edges
    for j, t in enumerate(right):
        for v in t:
            s = t - {v}
            if s in left_idx:
                i = left_idx[s]
                left_deg[i] += 1
                right_deg[j] += 1

    # Swap edges
    for i, s in enumerate(left):
        for u in s:
            s_minus_u = s - {u}
            forbidden_by_remaining = set()
            for x in s_minus_u:
                forbidden_by_remaining |= nbr[x]
            candidates = [v for v in range(n)
                          if v not in s and v not in forbidden_by_remaining]
            for ci, v in enumerate(candidates):
                for w in candidates[ci + 1:]:
                    if w not in nbr[v]:
                        t = s_minus_u | {v, w}
                        if t in right_idx:
                            j = right_idx[t]
                            left_deg[i] += 1
                            right_deg[j] += 1

    return {
        "left_min": min(left_deg),
        "left_max": max(left_deg),
        "left_avg": sum(left_deg) / n_left,
        "right_min": min(right_deg),
        "right_max": max(right_deg),
        "right_avg": sum(right_deg) / n_right,
        "n_left": n_left,
        "n_right": n_right,
        "k": k,
    }


def main():
    max_n = 14  # Keep moderate — degree analysis is O(swap edges)

    print("=" * 70)
    print("HALL'S DEGREE CONDITION CHECK")
    print("For augmented bipartite graph: left_min ≥ right_max?")
    print("=" * 70)

    total_levels = 0
    hall_degree_holds = 0
    hall_degree_fails = 0

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_levels = 0
        n_holds = 0
        n_fails = 0
        worst_ratio = float("inf")  # left_min / right_max
        worst_example = None

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)
            levels = enumerate_independent_sets(n, adj_data)

            for k in range(mode):
                lk = levels.get(k, [])
                if not lk:
                    continue

                result = analyze_degrees(levels, k, n, adj_data)
                if result is None:
                    continue

                n_levels += 1
                total_levels += 1

                ratio = result["left_min"] / result["right_max"] if result["right_max"] > 0 else float("inf")

                if result["left_min"] >= result["right_max"]:
                    n_holds += 1
                    hall_degree_holds += 1
                else:
                    n_fails += 1
                    hall_degree_fails += 1

                if ratio < worst_ratio:
                    worst_ratio = ratio
                    worst_example = {
                        "n": n, "tree": tidx, **result,
                    }

        elapsed = time.time() - t0
        print(f"\nn={n:2d} ({len(lines)} trees, {elapsed:.1f}s):", flush=True)
        print(f"  Levels: {n_levels}  Hall degree holds: {n_holds}  fails: {n_fails}",
              flush=True)
        if worst_example:
            ex = worst_example
            print(f"  Worst ratio (left_min/right_max): {worst_ratio:.3f}")
            print(f"    left_min={ex['left_min']} left_max={ex['left_max']} "
                  f"left_avg={ex['left_avg']:.1f}")
            print(f"    right_min={ex['right_min']} right_max={ex['right_max']} "
                  f"right_avg={ex['right_avg']:.1f}")
            print(f"    k={ex['k']} |left|={ex['n_left']} |right|={ex['n_right']}")

    print(f"\n\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total levels: {total_levels}")
    print(f"Hall degree condition holds: {hall_degree_holds}/{total_levels} "
          f"({100*hall_degree_holds/total_levels:.1f}%)")
    print(f"Hall degree condition fails: {hall_degree_fails}/{total_levels} "
          f"({100*hall_degree_fails/total_levels:.1f}%)")

    if hall_degree_fails > 0:
        print("\nHall's degree condition is NOT sufficient.")
        print("The augmented matching works anyway (0 failures in Hopcroft-Karp),")
        print("so the graph has expansion beyond what the degree condition captures.")


if __name__ == "__main__":
    main()
