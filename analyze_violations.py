#!/usr/bin/env python3
"""
Detailed analysis of the 20 cases where n < 3k but pigeonhole still works.
These are maximal IS of size k < mode where the bound P >= n-2k+1 is too weak
for pigeonhole, yet the actual private neighbor distribution still gives max_priv >= 2.

Why? What structural property makes these ISes have more private neighbors
than the bound predicts?
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


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def find_maximal_is_below_mode(n, adj, mode):
    nbr = [set(adj[v]) for v in range(n)]
    results = []

    def backtrack(v, current, forbidden):
        if v == n:
            k = len(current)
            if k >= mode:
                return
            s = frozenset(current)
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                results.append((s, k))
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return results


def analyze_is_detail(n, adj, s):
    """Detailed analysis of an independent set in a tree."""
    nbr = [set(adj[v]) for v in range(n)]
    deg = [len(adj[v]) for v in range(n)]
    k = len(s)

    # Private neighbor analysis
    priv_count = {}  # u -> private neighbor count
    priv_vertices = {}  # u -> set of private neighbors
    for u in s:
        privs = set()
        for v in nbr[u]:
            if v not in s:
                s_nbrs = nbr[v] & s
                if len(s_nbrs) == 1:  # v's only S-neighbor is u
                    privs.add(v)
        priv_count[u] = len(privs)
        priv_vertices[u] = privs

    total_P = sum(priv_count.values())
    bound_P = n - 2 * k + 1

    # Non-S vertices with exactly 1 S-neighbor (private) vs >=2 (shared)
    non_s = set(range(n)) - s
    private_non_s = set()
    shared_non_s = set()
    for v in non_s:
        s_nbrs = nbr[v] & s
        if len(s_nbrs) == 1:
            private_non_s.add(v)
        elif len(s_nbrs) >= 2:
            shared_non_s.add(v)
        # len == 0 shouldn't happen (S is dominating)

    # Edge count from S to V\S
    edges_s_to_non_s = sum(deg[u] for u in s)  # All edges from S go to V\S (S is independent)

    # Degree sequence of S-vertices
    s_degs = {u: deg[u] for u in s}

    # Check swap availability for vertices with priv >= 2
    swap_available = {}
    for u in s:
        if priv_count[u] >= 2:
            privs = sorted(priv_vertices[u])
            # Check all pairs: in a tree, two private nbrs of same vertex are always non-adjacent
            # (because u-v1 and u-v2 are edges, and v1-v2 edge would create a triangle)
            # But let's verify
            non_adj_pairs = 0
            for i, v1 in enumerate(privs):
                for v2 in privs[i+1:]:
                    if v2 not in nbr[v1]:
                        non_adj_pairs += 1
            swap_available[u] = non_adj_pairs

    return {
        "k": k,
        "s": sorted(s),
        "s_degs": {v: deg[v] for v in sorted(s)},
        "total_P": total_P,
        "bound_P": bound_P,
        "slack": total_P - bound_P,
        "priv_count": {v: priv_count[v] for v in sorted(s)},
        "priv_vertices": {v: sorted(priv_vertices[v]) for v in sorted(s)},
        "n_private_non_s": len(private_non_s),
        "n_shared_non_s": len(shared_non_s),
        "edges_s_to_non_s": edges_s_to_non_s,
        "max_priv": max(priv_count.values()),
        "swap_available": swap_available,
    }


def print_tree_structure(n, adj, s):
    """Print tree structure with IS membership."""
    deg = [len(adj[v]) for v in range(n)]
    print(f"    Vertices: {n}, Edges: {n-1}")
    for v in range(n):
        marker = " *" if v in s else "  "
        nbrs = sorted(adj[v])
        print(f"      v={v:2d} deg={deg[v]}{marker} -> {nbrs}")


def main():
    print("=" * 70)
    print("DETAILED ANALYSIS OF n < 3k VIOLATION CASES")
    print("=" * 70)

    for n in range(14, 17):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            max_is_list = find_maximal_is_below_mode(n, adj_data, mode)
            for s, k in max_is_list:
                excess = n - 3 * k
                if excess < 0:  # Violation!
                    print(f"\n{'='*60}")
                    print(f"VIOLATION: n={n}, tree={tidx}, k={k}, mode={mode}, n-3k={excess}")
                    print(f"g6: {line.strip()}")
                    print(f"Poly: {poly}")
                    print(f"{'='*60}")

                    detail = analyze_is_detail(n, adj_data, s)
                    print(f"\n  IS = {detail['s']}")
                    print(f"  IS vertex degrees: {detail['s_degs']}")
                    print(f"  Actual P = {detail['total_P']}, Bound P >= {detail['bound_P']}, Slack = {detail['slack']}")
                    print(f"  Per-vertex private neighbors:")
                    for v in detail['s']:
                        print(f"    u={v}: priv={detail['priv_count'][v]} "
                              f"privs={detail['priv_vertices'][v]}")
                    print(f"  Max priv: {detail['max_priv']}")
                    print(f"  Private non-S: {detail['n_private_non_s']}, "
                          f"Shared non-S: {detail['n_shared_non_s']}")
                    print(f"  Total edges S -> V\\S: {detail['edges_s_to_non_s']}")
                    print(f"  Swap pairs available: {detail['swap_available']}")

                    print(f"\n  Tree structure (* = in IS):")
                    print_tree_structure(n, adj_data, s)

    print(f"\n\n{'='*70}")
    print("ANALYSIS")
    print(f"{'='*70}")
    print("""
In the violation cases (n < 3k), the bound P >= n-2k+1 gives small P
(e.g., P >= 5 for n=16, k=6), but the ACTUAL P is much larger.

Why? The bound uses: edges from S to V\\S <= n-1 (tree edge count).
But when IS vertices have high degree, they "waste" edges on shared
non-S vertices, leaving more private neighbors than expected.

Actually the OPPOSITE: high-degree IS vertices mean more edges from S,
which could mean more sharing. But the tree structure constrains sharing.

The key structural property may be: in a tree, the number of shared
non-S vertices is limited by the tree structure. Specifically, each
shared vertex v with s >= 2 S-neighbors creates s-1 "extra" edges
beyond the domination minimum. The total extra edges equal
sum(s_nbrs(v) - 1) for shared v = edges_from_S - (n-k).
Since edges_from_S <= n-1, extra edges <= k-1.
So #shared_non_S * 1 <= k-1 (each shared saves at least 1 vs private).
Thus P = n-k - #shared >= n-k - (k-1) = n-2k+1. [This IS the bound!]

The actual P exceeds the bound when the average sharing multiplicity > 2
for shared vertices. In trees, this is common because shared vertices
tend to lie on paths between IS vertices and have exactly 2 S-neighbors.
""")


if __name__ == "__main__":
    main()
