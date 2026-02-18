#!/usr/bin/env python3
"""
Investigate the n=18 no-MLSV tree with a below-mode maximal IS.

g6 = Q??????_A?C?C?a?G_@C?CO?N_?
n=18, i(T)=5, mode=6

Questions:
1. What does the tree look like?
2. What are the maximal IS of size < 6?
3. Do they have any vertex with >= 2 leaf-neighbors?
   (If no MLSV, this seems impossible... unless leaf-neighbors
   are leaves adjacent to a vertex, not leaf-CHILDREN.)

WAIT: I need to distinguish:
- MLSV: vertex with >= 2 leaf-CHILDREN (in a rooting)
- Leaf-neighbor: vertex adjacent to a leaf (in unrooted tree)

Actually, in an unrooted tree, "leaf-children" of v = leaves adjacent to v.
A leaf has degree 1. If v has >= 2 neighbors with degree 1, v is an MLSV.

So if the tree has no MLSV, no vertex is adjacent to >= 2 leaves.
Then no vertex in any IS can have >= 2 leaf-neighbors.

THIS WOULD BREAK THE LNP.

Let me check very carefully.
"""

import subprocess
import sys
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


def find_all_maximal_is(n, adj):
    nbr = [set(adj[v]) for v in range(n)]
    results = []

    def backtrack(v, current, forbidden):
        if v == n:
            s = frozenset(current)
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                results.append(s)
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return results


g6 = "Q??????_A?C?C?a?G_@C?CO?N_?"
n, adj = parse_graph6(g6)

print(f"n = {n}")
print(f"Adjacency list:")
for v in range(n):
    print(f"  {v}: {adj[v]}")

deg = [len(adj[v]) for v in range(n)]
print(f"\nDegree sequence: {deg}")
print(f"Sorted: {sorted(deg, reverse=True)}")

leaves = [v for v in range(n) if deg[v] == 1]
print(f"\nLeaves: {leaves}")

# Check MLSV
for v in range(n):
    lc = sum(1 for w in adj[v] if deg[w] == 1)
    if lc >= 2:
        print(f"  MLSV: vertex {v} with {lc} leaf-children")
    elif lc == 1:
        leaf_nbr = [w for w in adj[v] if deg[w] == 1]
        print(f"  Vertex {v}: 1 leaf-child ({leaf_nbr[0]})")

poly = independence_poly(n, adj)
mode = compute_mode(poly)
print(f"\nIndependence polynomial: {poly}")
print(f"Mode: {mode} (i_{mode} = {poly[mode]})")

# Find all maximal IS
all_mis = find_all_maximal_is(n, adj)
print(f"\nTotal maximal IS: {len(all_mis)}")

# Below mode
below = [s for s in all_mis if len(s) < mode]
print(f"Maximal IS below mode ({mode}): {len(below)}")

for s in below:
    print(f"\n  S = {sorted(s)}  (size {len(s)})")
    # Check leaf-neighbors
    nbr = [set(adj[v]) for v in range(n)]
    for u in sorted(s):
        leaf_nbrs = [w for w in adj[u] if deg[w] == 1 and w not in s]
        print(f"    u={u} (deg={deg[u]}): "
              f"leaf-neighbors outside S = {leaf_nbrs} "
              f"(count={len(leaf_nbrs)})")

    max_leaf_nbrs = max(
        sum(1 for w in adj[u] if deg[w] == 1 and w not in s)
        for u in s
    )
    print(f"    MAX leaf-neighbors: {max_leaf_nbrs}")
    if max_leaf_nbrs < 2:
        print(f"    *** LNP FAILURE! ***")
    else:
        print(f"    LNP holds")
