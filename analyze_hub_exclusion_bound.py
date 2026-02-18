#!/usr/bin/env python3
"""For high-mode trees: compute the hub-exclusion bound.

When hub v is excluded from a maximal IS S:
  - All pendant leaves of v must be in S: contributes ell = leaf_deg(v)
  - Each non-leaf subtree T_w (hanging off non-leaf neighbor w) has
    S ∩ T_w being a maximal IS of T_w: contributes >= gamma(T_w)
  - Total: k >= ell + sum(gamma(T_w))

Key question: is ell + sum(gamma(T_w)) >= mode for all high-mode trees?
If so, excluding the hub forces k >= mode, which combined with
"including hub gives priv >= 2" proves PNP for high-mode trees.
"""

import subprocess
import sys
import time
from collections import deque

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


def find_domination_number(n, adj):
    """Find gamma(T) = minimum maximal IS size."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    nbr = [set(adj[v]) for v in range(n)]
    best = [n]  # upper bound

    def backtrack(v, current, forbidden):
        if len(current) >= best[0]:
            return
        if v == n:
            s = frozenset(current)
            for w in range(n):
                if w not in s and not (nbr[w] & s):
                    return
            best[0] = len(current)
            return
        # Try without v
        backtrack(v + 1, current, forbidden)
        # Try with v
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return best[0]


def get_subtree(n, adj, root, exclude):
    """Get the connected component containing root in T - exclude."""
    visited = set()
    queue = deque([root])
    visited.add(root)
    while queue:
        v = queue.popleft()
        for w in adj[v]:
            if w != exclude and w not in visited:
                visited.add(w)
                queue.append(w)
    return visited


def subtree_gamma(n, adj, vertices):
    """Compute gamma for the induced subtree on given vertices."""
    vlist = sorted(vertices)
    if not vlist:
        return 0
    remap = {v: i for i, v in enumerate(vlist)}
    m = len(vlist)
    sub_adj = [[] for _ in range(m)]
    for v in vlist:
        for w in adj[v]:
            if w in remap:
                sub_adj[remap[v]].append(remap[w])
    return find_domination_number(m, sub_adj)


def main():
    print("HUB EXCLUSION BOUND FOR HIGH-MODE TREES")
    print("=" * 70)
    print()
    print("bound = ell + sum(gamma(T_w)), where ell = pendant leaf-deg of hub")
    print()
    print(f"{'n':>3} {'high':>5} {'all_ok':>7} {'min_gap':>8} {'min_ell':>8} "
          f"{'min_bound':>10} {'min_mode':>9} {'time':>6}")
    print("-" * 60)

    for n in range(8, 21):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 2

        n_high = 0
        all_ok = True
        min_gap = float('inf')
        min_ell = float('inf')
        min_bound = float('inf')
        min_mode = float('inf')
        violation_examples = []

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode < threshold:
                continue

            n_high += 1

            # Find hub (max degree vertex)
            hub = max(range(tn), key=lambda v: len(adj_data[v]))
            hub_deg = len(adj_data[hub])

            # Count pendant leaves of hub
            ell = sum(1 for w in adj_data[hub] if len(adj_data[w]) == 1)

            # Find non-leaf neighbors and their subtrees
            non_leaf_nbrs = [w for w in adj_data[hub] if len(adj_data[w]) > 1]
            total_gamma = 0
            for w in non_leaf_nbrs:
                subtree_verts = get_subtree(tn, adj_data, w, hub)
                g = subtree_gamma(tn, adj_data, subtree_verts)
                total_gamma += g

            bound = ell + total_gamma
            gap = bound - mode

            min_gap = min(min_gap, gap)
            min_ell = min(min_ell, ell)
            min_bound = min(min_bound, bound)
            min_mode = min(min_mode, mode)

            if gap < 0:
                all_ok = False
                if len(violation_examples) < 3:
                    violation_examples.append(
                        (tn, mode, ell, total_gamma, bound, gap, line.strip()))

        elapsed = time.time() - t0
        if n_high == 0:
            print(f"{n:3d} {n_high:5d}       -        -        -          -         - {elapsed:6.1f}s")
            continue

        ok_str = "YES" if all_ok else "NO"
        print(f"{n:3d} {n_high:5d} {ok_str:>7} {min_gap:8d} {min_ell:8d} "
              f"{min_bound:10d} {min_mode:9d} {elapsed:6.1f}s")

        for ex in violation_examples:
            print(f"  VIOLATION: n={ex[0]}, mode={ex[1]}, ell={ex[2]}, "
                  f"sum_gamma={ex[3]}, bound={ex[4]}, gap={ex[5]}")

    print()
    print("=" * 70)
    print("If all_ok = YES for all n:")
    print("  ell + sum(gamma(T_w)) >= mode for every high-mode tree.")
    print("  This proves: excluding the hub forces k >= mode.")
    print("  Combined with 'including hub gives priv >= 2': PNP proved!")
    print("=" * 70)


if __name__ == "__main__":
    main()
