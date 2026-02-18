#!/usr/bin/env python3
"""Test: mode(I(T)) ≤ F + ceil((n - F - h + c) / 3) for trees with d_leaf ≥ 2.

F = total forced leaves (sum of d_leaf(v) for v with d_leaf ≥ 2)
h = number of excluded hubs (vertices with d_leaf ≥ 2)
c = components of residual forest T' = T - hubs - forced leaves
    (c ≥ 1 always; often c ≥ h)

If this bound holds, Case B of PNP follows from Part 1 applied to the residual.
The bound simplifies to: mode ≤ F + ceil((n - F) / 3) since c ≥ h ≥ 1.
"""

import subprocess
import sys
import time
import math

sys.path.insert(0, ".")
from indpoly import independence_poly


def parse_graph6(s):
    s = s.strip()
    if not s:
        return 0, []
    idx = 0
    if ord(s[0]) - 63 < 63:
        n = ord(s[0]) - 63; idx = 1
    else:
        idx = 1; n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63); idx += 1
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
                adj[i].append(j); adj[j].append(i)
            bit_idx += 1
    return n, adj


def compute_mode(poly):
    """Return the mode index (largest k with maximum i_k)."""
    mx = max(poly)
    for k in range(len(poly) - 1, -1, -1):
        if poly[k] == mx:
            return k
    return 0


def compute_dleaf_info(n, adj):
    """Compute d_leaf for each vertex and return forced leaf info."""
    leaves = set()
    for v in range(n):
        if len(adj[v]) == 1:
            leaves.add(v)

    dleaf = {}
    for v in range(n):
        if v in leaves:
            dleaf[v] = 0
        else:
            dleaf[v] = sum(1 for u in adj[v] if u in leaves)

    hubs = [v for v in range(n) if dleaf[v] >= 2]
    forced_leaves = set()
    for v in hubs:
        for u in adj[v]:
            if u in leaves:
                forced_leaves.add(u)

    F = len(forced_leaves)
    h = len(hubs)

    # Components of residual
    removed = set(hubs) | forced_leaves
    remaining = [v for v in range(n) if v not in removed]
    if not remaining:
        c = 0
    else:
        # BFS to count components
        visited = set()
        c = 0
        for start in remaining:
            if start in visited:
                continue
            c += 1
            queue = [start]
            visited.add(start)
            while queue:
                v = queue.pop()
                for u in adj[v]:
                    if u not in visited and u not in removed:
                        visited.add(u)
                        queue.append(u)

    return F, h, c, dleaf


def main():
    print("FORCED INCLUSION BOUND TEST", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    print("For trees with d_leaf >= 2 at some vertex:", flush=True)
    print("  k >= F + ceil((n - F - h + c) / 3)", flush=True)
    print("Test: mode(I(T)) <= F + ceil((n - F - h + c) / 3)?", flush=True)
    print(flush=True)

    total_tested = 0
    total_violations = 0
    min_surplus_all = float('inf')

    for n in range(5, 23):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_tested = 0
        n_violations = 0
        min_surplus = float('inf')
        tightest_info = None

        for line in lines:
            tn, adj_data = parse_graph6(line)
            F, h, c, dleaf = compute_dleaf_info(tn, adj_data)

            if F == 0:
                continue  # Only test trees with d_leaf >= 2

            n_tested += 1
            poly = independence_poly(tn, adj_data)
            mode = compute_mode(poly)

            # The bound
            n_residual = tn - F - h
            bound = F + math.ceil((n_residual + c) / 3) if n_residual + c > 0 else F

            surplus = bound - mode
            if surplus < min_surplus:
                min_surplus = surplus
                deg_seq = sorted([len(adj_data[v]) for v in range(tn)],
                                 reverse=True)
                dleaf_vals = sorted([dleaf[v] for v in range(tn)
                                     if dleaf[v] >= 2], reverse=True)
                tightest_info = {
                    "deg": deg_seq[:6],
                    "dleaf": dleaf_vals,
                    "F": F, "h": h, "c": c,
                    "mode": mode, "bound": bound,
                }

            if surplus < 0:
                n_violations += 1

        elapsed = time.time() - t0
        total_tested += n_tested
        total_violations += n_violations
        if min_surplus < min_surplus_all:
            min_surplus_all = min_surplus

        if n_tested > 0:
            status = "ALL OK" if n_violations == 0 else f"VIOLATIONS: {n_violations}"
            tight_str = ""
            if tightest_info:
                t = tightest_info
                tight_str = (f" tight: F={t['F']},h={t['h']},c={t['c']},"
                            f"mode={t['mode']},bound={t['bound']},"
                            f"dleaf={t['dleaf']}")
            print(f"n={n:2d}: {n_tested:6d} trees, min_surplus={min_surplus:2d} | "
                  f"{status}{tight_str} ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print(f"TOTAL: {total_tested} trees with d_leaf >= 2", flush=True)
    print(f"Violations: {total_violations}", flush=True)
    print(f"Minimum surplus: {min_surplus_all}", flush=True)

    if total_violations == 0:
        print(flush=True)
        print("BOUND HOLDS for all tested trees!", flush=True)
        print("=> PNP Case B follows from Part 1 + Hub Exclusion + transfer.", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
