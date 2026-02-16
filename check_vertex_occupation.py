#!/usr/bin/env python3
"""Check per-vertex occupation probability P(v in S) for d_leaf <= 1 trees.

If P(v in S) <= 1/3 for ALL vertices in ALL d_leaf <= 1 trees,
then μ = sum P(v in S) <= n/3, and mode <= ceil(μ) <= floor(n/3)+1.

Compute P(v in S) = (# IS containing v) / (# total IS) for each vertex.
"""

import subprocess
import sys
import time

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


def vertex_is_counts(n, adj):
    """For each vertex v, count IS containing v.

    Uses DP: I(T;x) with v forced in = x * I(T - N[v]; x).
    Evaluate at x=1 to get count.
    """
    # Total IS count
    total_poly = independence_poly(n, adj)
    total = sum(total_poly)

    counts = []
    nbr = [set(adj[v]) for v in range(n)]

    for v in range(n):
        # IS containing v: must exclude all neighbors of v
        # Remaining vertices: V \ N[v]
        remaining = [u for u in range(n) if u != v and u not in nbr[v]]
        if not remaining:
            # Only v can be in the IS
            counts.append(1)
            continue

        # Build subgraph on remaining vertices
        old_to_new = {u: i for i, u in enumerate(remaining)}
        m = len(remaining)
        sub_adj = [[] for _ in range(m)]
        for u in remaining:
            for w in adj[u]:
                if w in old_to_new:
                    sub_adj[old_to_new[u]].append(old_to_new[w])

        sub_poly = independence_poly(m, sub_adj)
        counts.append(sum(sub_poly))  # 1 * I(T-N[v]; 1)

    return counts, total


def main():
    print("PER-VERTEX OCCUPATION PROBABILITY FOR d_leaf <= 1 TREES",
          flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    print("Checking: P(v in S) <= 1/3 for every vertex v.", flush=True)
    print("If true: μ = Σ P(v in S) <= n/3, so mode <= floor(n/3)+1.",
          flush=True)
    print(flush=True)

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_dleaf1 = 0
        n_violations = 0
        max_prob = 0.0
        max_prob_info = None

        for line in lines:
            tn, adj_data = parse_graph6(line)

            # Check d_leaf condition
            leaves = set()
            for v in range(tn):
                if len(adj_data[v]) == 1:
                    leaves.add(v)

            ok = True
            for v in range(tn):
                if v in leaves:
                    continue
                d_leaf_v = sum(1 for u in adj_data[v] if u in leaves)
                if d_leaf_v >= 2:
                    ok = False
                    break
            if not ok:
                continue

            n_dleaf1 += 1
            counts, total = vertex_is_counts(tn, adj_data)

            for v in range(tn):
                p = counts[v] / total
                if p > max_prob:
                    max_prob = p
                    deg_v = len(adj_data[v])
                    is_leaf = v in leaves
                    ds = sorted([len(adj_data[u]) for u in range(tn)],
                                reverse=True)
                    max_prob_info = (tn, v, deg_v, is_leaf, ds[:6])

                if p > 1/3 + 1e-12:
                    n_violations += 1
                    if n_violations <= 5:
                        print(f"  VIOLATION n={tn}: v={v}, P={p:.6f}, "
                              f"deg={len(adj_data[v])}, "
                              f"leaf={'Y' if v in leaves else 'N'}",
                              flush=True)

        elapsed = time.time() - t0
        info_str = ""
        if max_prob_info:
            tn, v, dv, il, ds = max_prob_info
            info_str = (f" max_vertex: deg={dv}, "
                        f"{'leaf' if il else 'internal'}, tree_deg={ds}")
        print(f"n={n:2d}: {n_dleaf1:5d} trees, max_P={max_prob:.6f}, "
              f"viol={n_violations} ({elapsed:.1f}s){info_str}", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
