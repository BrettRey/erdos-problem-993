#!/usr/bin/env python3
"""Verify the forced inclusion bound: for trees with d_leaf >= 2,
any 1-Private maximal IS has k >= F + gamma(T'), and F + gamma(T') >= mode.

F = total forced leaves (leaf-children of vertices with d_leaf >= 2)
gamma(T') = minimum maximal IS of residual graph (T minus multi-leaf
            support vertices and their leaf-children)

Also check: does the "d/2 surplus" principle hold?
  mode <= F/2 + mode(I(T'))  (forced leaves add F to k but F/2 to mode)
"""

import subprocess
import sys
import time
from itertools import combinations

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


def get_residual(n, adj):
    """Compute forced leaves and residual graph for 1-Private IS.

    Returns (F, n_res, adj_res, vertex_map) where:
    - F = number of forced leaves (sum of d_leaf for vertices with d_leaf >= 2)
    - n_res, adj_res = residual graph after removing multi-leaf supports + their leaves
    - vertex_map = mapping from residual indices to original indices
    """
    leaves = set()
    for v in range(n):
        if len(adj[v]) == 1:
            leaves.add(v)

    # Find multi-leaf support vertices
    multi_support = set()
    for v in range(n):
        if v in leaves:
            continue
        d_leaf = sum(1 for u in adj[v] if u in leaves)
        if d_leaf >= 2:
            multi_support.add(v)

    # Forced leaves: leaf-children of multi-support vertices
    forced_leaves = set()
    for v in multi_support:
        for u in adj[v]:
            if u in leaves:
                forced_leaves.add(u)

    F = len(forced_leaves)

    # Residual: remove multi_support and forced_leaves
    removed = multi_support | forced_leaves
    remaining = [v for v in range(n) if v not in removed]

    if not remaining:
        return F, 0, [], {}

    # Build residual adjacency
    old_to_new = {v: i for i, v in enumerate(remaining)}
    n_res = len(remaining)
    adj_res = [[] for _ in range(n_res)]
    for v in remaining:
        for u in adj[v]:
            if u in old_to_new:
                adj_res[old_to_new[v]].append(old_to_new[u])

    return F, n_res, adj_res, old_to_new


def min_maximal_is(n, adj):
    """Find minimum maximal IS size in a graph (small n only)."""
    if n == 0:
        return 0

    nbr = [set(adj[v]) for v in range(n)]
    best = [n + 1]

    def backtrack(v, current, forbidden):
        if v == n:
            s = set(current)
            # Check maximality
            for w in range(n):
                if w not in s and not (nbr[w] & s):
                    return  # not maximal
            if len(current) < best[0]:
                best[0] = len(current)
            return

        if len(current) >= best[0]:
            return  # prune

        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return best[0]


def main():
    print("FORCED INCLUSION BOUND VERIFICATION", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    print("For trees with d_leaf >= 2:", flush=True)
    print("  Check: F + gamma(T') >= mode(I(T))", flush=True)
    print("  Check: mode(I(T)) <= F/2 + mode(I(T'))", flush=True)
    print(flush=True)

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 1
        n_trees = 0
        n_highmode = 0
        n_has_multi = 0  # high-mode trees with d_leaf >= 2
        n_bound_ok = 0   # F + gamma >= mode
        n_surplus_ok = 0 # mode <= F/2 + mode(T')
        n_viol_bound = 0
        n_viol_surplus = 0
        min_surplus = None  # min (F + gamma - mode)

        for line in lines:
            tn, adj_data = parse_graph6(line)
            n_trees += 1

            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode <= threshold:
                continue
            n_highmode += 1

            # Check for multi-leaf support
            leaves = set()
            for v in range(tn):
                if len(adj_data[v]) == 1:
                    leaves.add(v)

            has_multi = any(
                sum(1 for u in adj_data[v] if u in leaves) >= 2
                for v in range(tn) if v not in leaves
            )

            if not has_multi:
                continue  # Case A (d_leaf <= 1, should be low-mode)
            n_has_multi += 1

            # Compute forced inclusion bound
            F, n_res, adj_res, _ = get_residual(tn, adj_data)

            if n_res > 0:
                gamma = min_maximal_is(n_res, adj_res)
                poly_res = independence_poly(n_res, adj_res)
                mode_res = max(range(len(poly_res)),
                              key=lambda i: poly_res[i])
            else:
                gamma = 0
                mode_res = 0

            bound = F + gamma
            surplus_bound = F // 2 + mode_res  # F/2 + mode(T')

            if bound >= mode:
                n_bound_ok += 1
            else:
                n_viol_bound += 1
                print(f"  BOUND FAIL n={tn}: F={F}, gamma={gamma}, "
                      f"F+gamma={bound} < mode={mode}", flush=True)

            if surplus_bound >= mode:
                n_surplus_ok += 1
            else:
                n_viol_surplus += 1

            surplus = bound - mode
            if min_surplus is None or surplus < min_surplus:
                min_surplus = surplus

        elapsed = time.time() - t0
        ms = f"min_surplus={min_surplus}" if min_surplus is not None else "N/A"
        print(f"n={n:2d}: {n_highmode:4d} high-mode, {n_has_multi:4d} with d_leaf>=2, "
              f"bound_ok={n_bound_ok}, surplus_ok={n_surplus_ok}, "
              f"viol={n_viol_bound}/{n_viol_surplus}, {ms} ({elapsed:.1f}s)",
              flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
