#!/usr/bin/env python3
"""
Investigate WHY the augmented matching always works.

Key questions:
1. What is the minimum augmented degree for maximal IS? (Lower bound?)
2. What do the "hardest" cases look like structurally?
3. Can we prove the key lemma: every maximal IS below mode has ≥1 swap neighbor?
4. What is the relationship between IS size, tree structure, and swap degree?

Also: explore the exchange structure more carefully.
For a maximal IS S of size k and any IS T* of size k+1:
  |S ∩ T*| = k - m, |S \ T*| = m, |T* \ S| = m + 1
  S is a swap neighbor of T* when m = 1.
  How often is m = 1 vs m > 1?
"""

import subprocess
import sys
import time
from collections import defaultdict, Counter

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


def analyze_maximal_swap(n, adj, levels, mode):
    """For each maximal IS below mode, compute detailed swap structure."""
    nbr = [set(adj[v]) for v in range(n)]
    deg = [len(adj[v]) for v in range(n)]

    results = []

    for k in range(mode):
        left = levels.get(k, [])
        right_set = set(levels.get(k + 1, []))

        for s in left:
            # Check if maximal
            can_extend = any(v not in s and not (nbr[v] & s) for v in range(n))
            if can_extend:
                continue

            # This IS is maximal at size k < mode
            # Compute swap neighbors
            swap_neighbors = set()
            swap_detail = []  # (u_removed, v_added, w_added) tuples

            for u in s:
                s_minus_u = s - {u}
                # Vertices not adjacent to any remaining vertex
                forbidden = set()
                for x in s_minus_u:
                    forbidden |= nbr[x]
                candidates = [v for v in range(n)
                              if v not in s and v not in forbidden]
                # All pairs of non-adjacent candidates
                for ci, v in enumerate(candidates):
                    for w in candidates[ci + 1:]:
                        if w not in nbr[v]:
                            t = s_minus_u | {v, w}
                            if t in right_set:
                                swap_neighbors.add(t)
                                swap_detail.append((u, v, w))

            # Analyze the IS structure
            s_degs = sorted([deg[v] for v in s], reverse=True)
            max_deg_in_s = max(deg[v] for v in s)
            sum_deg_in_s = sum(deg[v] for v in s)

            results.append({
                "k": k,
                "s": sorted(s),
                "swap_degree": len(swap_neighbors),
                "swap_triples": len(swap_detail),
                "max_deg_in_s": max_deg_in_s,
                "sum_deg_in_s": sum_deg_in_s,
                "deg_seq": s_degs,
                "i_k": len(left),
                "i_k1": len(levels.get(k + 1, [])),
            })

    return results


def analyze_exchange_distance(n, adj, levels, mode):
    """For each maximal IS S below mode and each IS T* of size k+1,
    compute the exchange distance |S \\ T*|.
    Returns distribution of m values."""
    m_dist = Counter()

    for k in range(mode):
        left = levels.get(k, [])
        right = levels.get(k + 1, [])
        nbr = [set(adj[v]) for v in range(n)]

        for s in left:
            can_extend = any(v not in s and not (nbr[v] & s) for v in range(n))
            if can_extend:
                continue

            # Maximal IS at level k
            for t in right:
                m = len(s - t)  # |S \ T*|
                m_dist[m] += 1

    return m_dist


def main():
    max_n = 14

    print("=" * 70)
    print("ANALYSIS: Why does the augmented matching always work?")
    print("=" * 70)

    # Track minimum swap degree across all trees
    global_min_swap_deg = float("inf")
    global_min_swap_example = None
    swap_deg_hist = Counter()
    total_maximal = 0
    m_dist_global = Counter()

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_maximal = 0
        min_swap_deg_n = float("inf")
        min_swap_example_n = None

        for tidx, line in enumerate(lines):
            tn, adj = parse_graph6(line)
            poly = independence_poly(n, adj)
            mode = compute_mode(poly)
            levels = enumerate_independent_sets(n, adj)

            # Detailed swap analysis for maximal IS
            max_results = analyze_maximal_swap(n, adj, levels, mode)
            for r in max_results:
                n_maximal += 1
                total_maximal += 1
                sd = r["swap_degree"]
                swap_deg_hist[sd] += 1
                if sd < min_swap_deg_n:
                    min_swap_deg_n = sd
                    min_swap_example_n = {
                        "n": n, "tree": tidx, "g6": line.strip(),
                        **r,
                    }
                if sd < global_min_swap_deg:
                    global_min_swap_deg = sd
                    global_min_swap_example = {
                        "n": n, "tree": tidx, "g6": line.strip(),
                        **r,
                    }

            # Exchange distance analysis (only for small n, expensive)
            if n <= 10:
                m_dist = analyze_exchange_distance(n, adj, levels, mode)
                for m, cnt in m_dist.items():
                    m_dist_global[m] += cnt

        elapsed = time.time() - t0
        print(f"\nn={n} ({len(lines)} trees, {elapsed:.1f}s):", flush=True)
        print(f"  Maximal IS below mode: {n_maximal}", flush=True)
        if n_maximal > 0:
            print(f"  Min swap degree at n={n}: {min_swap_deg_n}", flush=True)
            if min_swap_example_n:
                ex = min_swap_example_n
                print(f"    Example: k={ex['k']} IS={ex['s']} "
                      f"swap_deg={ex['swap_degree']} "
                      f"deg_seq={ex['deg_seq']} "
                      f"i_k={ex['i_k']} i_k+1={ex['i_k1']}", flush=True)

    print(f"\n\n{'='*70}")
    print("GLOBAL SUMMARY")
    print(f"{'='*70}")

    print(f"\nTotal maximal IS below mode analyzed: {total_maximal}")
    print(f"Global minimum swap degree: {global_min_swap_deg}")
    if global_min_swap_example:
        ex = global_min_swap_example
        print(f"  Achieved at: n={ex['n']} tree={ex['tree']} k={ex['k']}")
        print(f"  IS = {ex['s']}, deg_seq = {ex['deg_seq']}")
        print(f"  swap_degree = {ex['swap_degree']}, i_k = {ex['i_k']}, i_k+1 = {ex['i_k1']}")

    print(f"\nSwap degree distribution:")
    for sd in sorted(swap_deg_hist.keys()):
        cnt = swap_deg_hist[sd]
        pct = 100 * cnt / total_maximal
        bar = "#" * min(50, int(pct))
        print(f"  deg={sd:4d}: {cnt:6d} ({pct:5.1f}%) {bar}")

    if m_dist_global:
        print(f"\nExchange distance distribution (|S \\ T*| for maximal S, n≤10):")
        total_pairs = sum(m_dist_global.values())
        for m in sorted(m_dist_global.keys()):
            cnt = m_dist_global[m]
            pct = 100 * cnt / total_pairs
            print(f"  m={m}: {cnt:6d} ({pct:5.1f}%)")
        # m=1 means direct swap neighbor
        if 1 in m_dist_global:
            print(f"\n  Fraction of (maximal S, T*) pairs at distance 1 (direct swap): "
                  f"{100*m_dist_global[1]/total_pairs:.1f}%")


if __name__ == "__main__":
    main()
