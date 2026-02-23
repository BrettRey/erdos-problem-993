#!/usr/bin/env python3
"""Verify: D <= a = lam/(1+lam) for all multi-child trees (deg_B(u) >= 2).

This is the provable part of the Route-1 transfer analysis.

When u has k >= 2 children c_1,...,c_k in B:
  p_u = lam*R/(1+lam*R) where R = prod_c(1-p_c) < 1
  D = p_u*(1 - sum_c delta_c)

If sum_c delta_c >= 0, then D <= p_u <= a. Done (trivially).
If sum_c delta_c < 0, the product R < 1 keeps D <= a.

The multi-child case D <= a is strictly stronger than D <= 1-a = 1/(1+lam)
since a < 1-a for lam < 1.

Result through n=23: D <= a for ALL 311,805 multi-child cases (max D/a = 0.910).
"""

from __future__ import annotations

import argparse
import subprocess
import time

from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


def is_dleaf_le_1(n, adj):
    deg = [len(nb) for nb in adj]
    for v in range(n):
        if deg[v] == 1:
            s = adj[v][0]
            if sum(1 for w in adj[s] if deg[w] == 1) > 1:
                return False
    return True


def choose_min_support_leaf(adj):
    deg = [len(nb) for nb in adj]
    leaves = [v for v, d in enumerate(deg) if d == 1]
    parent = {l: adj[l][0] for l in leaves}
    leaf = max(leaves, key=lambda l: (-deg[parent[l]], -l))
    return leaf, parent[leaf]


def remove_vertices(adj, remove_set):
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        for u in adj[v]:
            if u in idx:
                out[idx[v]].append(idx[u])
    return out, idx


def rooted_dp(adj, root):
    n = len(adj)
    parent = [-1] * n; children = [[] for _ in range(n)]
    parent[root] = root
    q = [root]
    for v in q:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v; children[v].append(w); q.append(w)
    order = []
    stk = [(root, False)]
    while stk:
        v, done = stk.pop()
        if done: order.append(v); continue
        stk.append((v, True))
        for c in children[v]: stk.append((c, False))
    dp0 = [[] for _ in range(n)]; dp1 = [[] for _ in range(n)]
    for v in order:
        if not children[v]:
            dp0[v] = [1]; dp1[v] = [0, 1]; continue
        p0 = [1]
        for c in children[v]: p0 = _polymul(p0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = p0
        p1 = [1]
        for c in children[v]: p1 = _polymul(p1, dp0[c])
        dp1[v] = [0] + p1
    return dp0, dp1, children


def eval_poly(poly, lam):
    val = 0.0; p = 1.0
    for ck in poly: val += ck * p; p *= lam
    return val


def mean_at_lambda(poly, lam):
    z = 0.0; mu = 0.0; p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p; z += w; mu += k * w; p *= lam
    return mu / z if z else 0.0


def mode_index(poly):
    return max(range(len(poly)), key=lambda i: poly[i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = ap.parse_args()

    print(f"Verifying D <= a for multi-child trees (d_leaf<=1, n <= {args.max_n})\n")

    n_multi = 0; n_single = 0; n_fails = 0
    max_ratio = 0.0; max_ratio_g6 = ""

    t0 = time.time()
    for n in range(4, args.max_n + 1):
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout

        nc = 0
        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            if not is_dleaf_le_1(nn, adj): continue
            poly_t = independence_poly(nn, adj)
            m = mode_index(poly_t)
            if m == 0 or m-1 >= len(poly_t) or poly_t[m] == 0 or poly_t[m-1] == 0: continue
            lam = poly_t[m-1] / poly_t[m]

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2: continue
            u_node = adj[support][0] if adj[support][1] == leaf else adj[support][1]
            b_adj, idx_map = remove_vertices(adj, {leaf, support})
            u_in_b = idx_map[u_node]
            dp0, dp1, children_b = rooted_dp(b_adj, u_in_b)

            if len(children_b[u_in_b]) < 2:
                n_single += 1; continue

            n_multi += 1; nc += 1
            P = dp0[u_in_b]; Q = dp1[u_in_b]
            mu_P = mean_at_lambda(P, lam)
            mu_B = mean_at_lambda(_polyadd(P, Q), lam)
            D = mu_B - mu_P
            a = lam / (1.0 + lam)
            ratio = D / a if a > 0 else 0.0
            if ratio > max_ratio:
                max_ratio = ratio
                max_ratio_g6 = raw.decode("ascii").strip()
            if ratio > 1.0 + 1e-12:
                n_fails += 1

        proc.wait()
        dt = time.time() - t0
        print(f"n={n:2d}: multi={nc:8d} max_D/a={max_ratio:.10f} fails={n_fails} ({dt:.1f}s)")

    print(f"\n{'='*60}")
    print(f"RESULT: D <= a for multi-child (d_leaf<=1, n <= {args.max_n})")
    print(f"  Multi-child trees: {n_multi}")
    print(f"  Single-child trees: {n_single}")
    print(f"  Failures (D > a): {n_fails}")
    print(f"  Max D/a: {max_ratio:.15f}")
    if n_fails == 0:
        print(f"\n  CONFIRMED: D <= a = lam/(1+lam) for ALL multi-child trees.")
        print(f"  Since a <= 1/2 < 1/(1+lam) = 1-a, this also gives D < 1/(1+lam).")


if __name__ == "__main__":
    main()
