#!/usr/bin/env python3
"""Investigate WHY d_leaf ≤ 1 trees always have low mode.

Key idea: tree IS polynomials are real-rooted (Chudnovsky-Seymour 2007).
For real-rooted polynomials with positive coefficients, |mode - mean| ≤ 1.
So proving mean ≤ n/3 + 1 would imply mode ≤ n/3 + 2 ≈ threshold.

Compute:
1. Mean μ = I'(1)/I(1) for all d_leaf ≤ 1 trees
2. Compare μ to n/3 and threshold
3. Identify tightest cases
4. Analyze whether path P_n is always the extremal case
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


def poly_mean(poly):
    """Compute mean μ = I'(1)/I(1) of the IS size distribution."""
    total = sum(poly)  # I(1)
    weighted = sum(k * poly[k] for k in range(len(poly)))  # I'(1)
    if total == 0:
        return 0.0
    return weighted / total


def classify_tree(n, adj):
    """Classify a d_leaf ≤ 1 tree by structure."""
    leaves = set()
    for v in range(n):
        if len(adj[v]) == 1:
            leaves.add(v)

    deg_seq = sorted([len(adj[v]) for v in range(n)], reverse=True)
    max_deg = deg_seq[0] if deg_seq else 0

    # Count branch points (degree >= 3)
    branch_pts = sum(1 for v in range(n) if len(adj[v]) >= 3)

    # Is it a path?
    is_path = all(len(adj[v]) <= 2 for v in range(n))

    # Is it a caterpillar? (all vertices within distance 1 of spine)
    # Simplified check: remove leaves, check if remaining is a path
    inner = [v for v in range(n) if v not in leaves]
    is_caterpillar = False
    if inner:
        inner_set = set(inner)
        inner_degs = {v: sum(1 for u in adj[v] if u in inner_set)
                      for v in inner}
        is_caterpillar = all(d <= 2 for d in inner_degs.values())

    if is_path:
        return "path"
    elif is_caterpillar:
        return f"caterpillar(b={branch_pts})"
    else:
        return f"other(b={branch_pts},Δ={max_deg})"


def main():
    print("INVESTIGATING d_leaf ≤ 1 TREES: MEAN vs MODE vs THRESHOLD",
          flush=True)
    print("=" * 70, flush=True)
    print(flush=True)
    print("Real-rooted ⟹ |mode - mean| ≤ 1.", flush=True)
    print("So if mean ≤ n/3, then mode ≤ n/3 + 1 = threshold.", flush=True)
    print(flush=True)

    # Also compute for paths specifically
    print("PATH BASELINE:", flush=True)
    for n in range(5, 30):
        adj = [[] for _ in range(n)]
        for i in range(n - 1):
            adj[i].append(i + 1)
            adj[i + 1].append(i)
        poly = independence_poly(n, adj)
        mode = max(range(len(poly)), key=lambda i: poly[i])
        mu = poly_mean(poly)
        thr = n // 3 + 1
        print(f"  P_{n:2d}: mode={mode:2d}, μ={mu:.4f}, n/3={n/3:.2f}, "
              f"thr={thr}, μ/(n/3)={mu/(n/3):.4f}", flush=True)

    print(flush=True)
    print("EXHAUSTIVE SCAN:", flush=True)
    print("-" * 70, flush=True)

    all_tight = []  # Collect tightest cases for later analysis

    for n in range(5, 21):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 1
        n_dleaf1 = 0
        max_mu = 0.0
        max_mode = 0
        max_mu_ratio = 0.0  # μ / (n/3)
        tightest = None  # (μ, mode, tree_type, deg_seq)
        path_mu = None

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
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            mu = poly_mean(poly)
            mu_ratio = mu / (tn / 3) if tn > 0 else 0

            tree_type = classify_tree(tn, adj_data)
            if tree_type == "path":
                path_mu = mu

            if mu_ratio > max_mu_ratio:
                max_mu_ratio = mu_ratio
                max_mu = mu
                max_mode = mode
                ds = sorted([len(adj_data[v]) for v in range(tn)],
                            reverse=True)
                tightest = (mu, mode, tree_type, ds[:8])

            if mu > max_mu:
                max_mu = mu

        elapsed = time.time() - t0

        thr_third = n / 3
        tight_desc = ""
        if tightest:
            mu_t, mode_t, ttype, ds = tightest
            tight_desc = (f"  tightest: μ={mu_t:.4f}, mode={mode_t}, "
                          f"type={ttype}, deg={ds}")
            if mu_t > thr_third:
                tight_desc += " [μ > n/3!]"
            all_tight.append((n, mu_t, mode_t, ttype, ds))

        pm_str = f"path_μ={path_mu:.4f}" if path_mu else ""
        print(f"n={n:2d}: {n_dleaf1:6d} trees, max_μ={max_mu:.4f}, "
              f"max_μ/(n/3)={max_mu_ratio:.4f}, max_mode={max_mode}, "
              f"thr={threshold} ({elapsed:.1f}s) {pm_str}", flush=True)
        if tight_desc:
            print(tight_desc, flush=True)

    # Summary
    print(flush=True)
    print("TIGHTEST CASES SUMMARY:", flush=True)
    print("-" * 70, flush=True)
    for n, mu, mode, ttype, ds in all_tight:
        thr = n / 3
        print(f"n={n:2d}: μ={mu:.4f} (n/3={thr:.2f}), mode={mode}, "
              f"thr={n//3+1}, type={ttype}", flush=True)

    # Key question: is the path always the tightest?
    print(flush=True)
    print("KEY QUESTION: Is P_n always the tightest d_leaf ≤ 1 tree?",
          flush=True)
    print("(If yes, proving μ(P_n) ≤ n/3 suffices for Conjecture A.)",
          flush=True)
    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
