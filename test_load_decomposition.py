#!/usr/bin/env python3
"""Test the load decomposition approach for proving μ < n/3.

In the hard-core model at λ=1, μ = Σ P(v ∈ S) where P is the
hard-core measure. Define load L(v) = P(v) - 1/3. Then μ < n/3
iff Σ L(v) < 0.

For d_leaf ≤ 1 trees, pair each leaf w with its support vertex v:
  P(w) = (1 - P(v)) / 2
  L(w) + L(v) = P(v)/2 - 1/6

This is ≤ 0 iff P(v) ≤ 1/3 (i.e., R(v) ≤ 1/2 in the tree recursion).

For support vertices with ≥ 1 non-leaf child: R(v) = 1/(2·∏(1+R(c)))
where product is over non-leaf children. Since ∏(1+R(c)) > 1 when there
are non-leaf children, R(v) < 1/2, so P(v) < 1/3, and the pair load < 0.

For pendant edges (support vertex has ONLY the leaf child, degree 2):
R(v) = 1/2, P(v) = 1/3, pair load = 0 (tight but ok).

After pairing: remaining "core" vertices have d_leaf = 0 (no leaf children).
Need: Σ_core L(v) < 0, i.e., Σ_core P(v) < |core|/3.

This script verifies the decomposition computationally.
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


def check_dleaf(n, adj):
    leaves = set()
    for v in range(n):
        if len(adj[v]) == 1:
            leaves.add(v)
    for v in range(n):
        if v in leaves:
            continue
        d_leaf = sum(1 for u in adj[v] if u in leaves)
        if d_leaf >= 2:
            return False
    return True


def compute_vertex_probs(n, adj):
    """Compute P(v ∈ S) for each vertex using the IS polynomial."""
    # For each vertex v: P(v ∈ S) = I_{T-N[v]}(1) / I_T(1)
    # where I_{T-N[v]} counts IS in the graph with v and its neighbors removed,
    # and the IS includes v (so we count IS containing v).
    # Actually: I(T;x) = I(T-v; x) + x * I(T-N[v]; x)
    # At x=1: I(1) = I(T-v;1) + I(T-N[v];1)
    # IS containing v = I(T-N[v];1), so P(v) = I(T-N[v];1) / I(T;1)

    poly = independence_poly(n, adj)
    total = sum(poly)

    probs = []
    for v in range(n):
        # Compute I(T-N[v]; 1) = number of IS containing v
        # Build subgraph T - N[v]
        nv = set(adj[v]) | {v}
        remaining = [u for u in range(n) if u not in nv]
        if not remaining:
            # v is connected to everyone; only IS containing v is {v}
            probs.append(1.0 / total)
            continue

        # Build adjacency for remaining vertices
        remap = {u: i for i, u in enumerate(remaining)}
        rn = len(remaining)
        radj = [[] for _ in range(rn)]
        for u in remaining:
            for w in adj[u]:
                if w in remap:
                    radj[remap[u]].append(remap[w])

        rpoly = independence_poly(rn, radj)
        is_with_v = sum(rpoly)
        probs.append(is_with_v / total)

    return probs


def main():
    print("LOAD DECOMPOSITION ANALYSIS", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    print("For d_leaf ≤ 1 trees, decompose μ - n/3 into:", flush=True)
    print("  (1) Leaf-support pair loads: L(w)+L(v) = P(v)/2 - 1/6", flush=True)
    print("  (2) Core vertex loads: L(v) = P(v) - 1/3", flush=True)
    print(flush=True)

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_tested = 0
        max_core_load = -float('inf')  # max Σ_core L(v) / |core|
        max_pair_load = -float('inf')  # max pair load
        max_total_load = -float('inf')  # max total load (should be < 0)
        tightest = None

        for line in lines:
            tn, adj_data = parse_graph6(line)
            if not check_dleaf(tn, adj_data):
                continue

            n_tested += 1
            probs = compute_vertex_probs(tn, adj_data)

            # Classify vertices
            leaves = set()
            for v in range(tn):
                if len(adj_data[v]) == 1:
                    leaves.add(v)

            support = set()
            for v in range(tn):
                if v in leaves:
                    continue
                for u in adj_data[v]:
                    if u in leaves:
                        support.add(v)
                        break

            core = set(range(tn)) - leaves - support

            # Pair loads
            pair_load_total = 0.0
            for v in support:
                pair_load_total += probs[v] / 2 - 1.0 / 6

            # Core loads
            core_load_total = sum(probs[v] - 1.0/3 for v in core)
            core_size = len(core)

            total_load = pair_load_total + core_load_total

            if core_size > 0:
                avg_core_load = core_load_total / core_size
                if avg_core_load > max_core_load:
                    max_core_load = avg_core_load

            if pair_load_total > max_pair_load:
                max_pair_load = pair_load_total

            if total_load > max_total_load:
                max_total_load = total_load
                deg_seq = sorted([len(adj_data[v]) for v in range(tn)],
                                 reverse=True)
                mu = sum(probs)
                tightest = {
                    "deg": deg_seq[:6],
                    "ell": len(leaves),
                    "core": core_size,
                    "mu": mu,
                    "ratio": mu / (tn / 3),
                    "pair_load": pair_load_total,
                    "core_load": core_load_total,
                    "total_load": total_load,
                }

        elapsed = time.time() - t0

        if n_tested > 0 and tightest:
            t = tightest
            status = "OK" if t["total_load"] < -1e-12 else (
                "TIGHT" if abs(t["total_load"]) < 1e-6 else "VIOLATION!")
            print(f"n={n:2d}: {n_tested:5d} trees | "
                  f"max_total_load={max_total_load:.6f} "
                  f"(pair={max_pair_load:.4f}, "
                  f"core_avg={max_core_load:.4f}) | "
                  f"{status} ({elapsed:.1f}s)", flush=True)
            if n <= 12 or n % 3 == 0:
                print(f"      tightest: deg={t['deg']}, ℓ={t['ell']}, "
                      f"core={t['core']}, μ/(n/3)={t['ratio']:.4f}",
                      flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
