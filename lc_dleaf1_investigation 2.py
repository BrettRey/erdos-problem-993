#!/usr/bin/env python3
"""
Investigate log-concavity of d_leaf<=1 IS polynomials.

Three checks:
1. ULC (ultra-log-concave) for d_leaf<=1 trees: (i_k/C(α,k))^2 >= i_{k-1}/C(α,k-1) * i_{k+1}/C(α,k+1)
2. Leaf-removal cross-term: h_k = f_k + g_{k-1}. LC of h requires 2f_k g_{k-1} >= f_{k-1}g_k + g_{k-2}f_{k+1}
3. Whether T-l (leaf removed) is always LC (even if not d_leaf<=1)

Also checks: does T-{l,s} stay d_leaf<=1? (should be yes)
"""
from __future__ import annotations
import subprocess
import sys
from math import comb
from indpoly import independence_poly as indpoly

def parse_graph6(line: bytes) -> tuple[int, list[list[int]]]:
    data = line.strip()
    if data[0] == ord('~'):
        n = (data[1] - 63) << 12 | (data[2] - 63) << 6 | (data[3] - 63)
        bits_start = 4
    else:
        n = data[0] - 63
        bits_start = 1
    adj: list[list[int]] = [[] for _ in range(n)]
    bit_idx = 0
    for b in data[bits_start:]:
        val = b - 63
        for i in range(5, -1, -1):
            if bit_idx >= n * (n - 1) // 2:
                break
            if (val >> i) & 1:
                col = int((1 + (1 + 8 * bit_idx) ** 0.5) / 2)
                row = bit_idx - col * (col - 1) // 2
                adj[row].append(col)
                adj[col].append(row)
            bit_idx += 1
    return n, adj

def is_dleaf_le1(n: int, adj: list[list[int]]) -> bool:
    for v in range(n):
        leaf_count = sum(1 for u in adj[v] if len(adj[u]) == 1)
        if leaf_count > 1:
            return False
    return True


def is_lc(poly: list[int | float]) -> bool:
    for k in range(1, len(poly) - 1):
        if poly[k] ** 2 < poly[k-1] * poly[k+1] - 1e-9:
            return False
    return True

def is_ulc(poly: list[int | float]) -> bool:
    """Ultra-log-concave: (i_k/C(d,k))^2 >= (i_{k-1}/C(d,k-1))*(i_{k+1}/C(d,k+1))"""
    d = len(poly) - 1
    for k in range(1, d):
        if poly[k-1] == 0 or poly[k+1] == 0:
            continue
        lhs = (poly[k] / comb(d, k)) ** 2
        rhs = (poly[k-1] / comb(d, k-1)) * (poly[k+1] / comb(d, k+1))
        if lhs < rhs - 1e-9:
            return False
    return True

def subgraph_induced(n: int, adj: list[list[int]], keep: set[int]) -> tuple[int, list[list[int]]]:
    """Induced subgraph on vertices in 'keep', re-indexed."""
    vmap = {v: i for i, v in enumerate(sorted(keep))}
    m = len(keep)
    new_adj: list[list[int]] = [[] for _ in range(m)]
    for v in keep:
        for u in adj[v]:
            if u in keep:
                new_adj[vmap[v]].append(vmap[u])
    return m, new_adj

def cross_term_lc(f: list, g: list, tol: float = 1e-9) -> tuple[bool, float]:
    """Check cross-term condition for h = f + x*g being LC.
    h_k = f_k + g_{k-1}.
    LC(h): h_k^2 >= h_{k-1}*h_{k+1} for all k.
    = (f_k + g_{k-1})^2 >= (f_{k-1}+g_{k-2})*(f_{k+1}+g_k)
    This expands to: f_k^2 - f_{k-1}f_{k+1} + g_{k-1}^2 - g_{k-2}g_k + 2f_k g_{k-1} - f_{k-1}g_k - g_{k-2}f_{k+1} >= 0
    Define cross = 2f_k*g_{k-1} - f_{k-1}*g_k - g_{k-2}*f_{k+1}.
    LC(f) and LC(g) give the first two terms. We only need cross >= -(lc_f_margin + lc_g_margin).
    But more directly, just check h is LC.
    Returns (h_is_lc, min_lc_ratio_of_h).
    """
    # Build h
    max_len = max(len(f), len(g) + 1)
    h = [0.0] * max_len
    for k in range(len(f)):
        h[k] += f[k]
    for k in range(len(g)):
        h[k+1] += g[k]

    # Check LC of h
    min_ratio = float('inf')
    for k in range(1, len(h) - 1):
        if h[k-1] > 0 and h[k+1] > 0:
            ratio = h[k]**2 / (h[k-1] * h[k+1])
            min_ratio = min(min_ratio, ratio)
        elif h[k-1] == 0 or h[k+1] == 0:
            pass  # vacuously ok
    return min_ratio >= 1 - tol, min_ratio

def main():
    max_n = 22  # adjust as needed
    geng = "/opt/homebrew/bin/geng"

    # Counters
    total = 0
    dleaf1_total = 0
    lc_fail = 0
    ulc_fail = 0
    ulc_fail_first = None

    # Leaf-removal analysis
    lr_total = 0       # total leaf-removal instances checked
    lr_Tl_lc_fail = 0  # T-l not LC
    lr_Tls_lc_fail = 0  # T-{l,s} not LC
    lr_h_lc_fail = 0   # h = I(T-l) + x*I(T-{l,s}) not LC (sanity: should match I(T) LC)
    lr_crossterm_fail = 0  # cross-term < 0 (but h might still be LC due to self-terms)
    lr_Tls_dleaf1_fail = 0  # T-{l,s} not d_leaf<=1
    lr_Tl_dleaf1_fail = 0   # T-l not d_leaf<=1 (expected to fail sometimes)

    worst_min_ratio = float('inf')
    worst_ulc_ratio = float('inf')

    print(f"Checking d_leaf<=1 trees n=3..{max_n}")
    print(f"{'n':>4} {'trees':>8} {'dleaf1':>8} {'lc_fail':>8} {'ulc_fail':>8} {'lr_Tl_lc':>9} {'lr_Tls_lc':>10} {'lr_cross_fail':>14}")
    print("-" * 90)

    for n in range(3, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout

        n_total = 0
        n_dleaf1 = 0
        n_lc_fail = 0
        n_ulc_fail = 0
        n_lr_Tl_lc = 0
        n_lr_Tls_lc = 0
        n_lr_cross = 0
        n_lr_Tls_d1_fail = 0
        n_lr_Tl_d1_fail = 0

        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_total += 1

            if not is_dleaf_le1(n0, adj):
                continue
            n_dleaf1 += 1

            poly = indpoly(n0, adj)

            # LC check
            if not is_lc(poly):
                n_lc_fail += 1

            # ULC check
            if not is_ulc(poly):
                n_ulc_fail += 1
                if ulc_fail_first is None:
                    ulc_fail_first = (n0, list(poly))

            # Leaf-removal analysis
            leaves = [v for v in range(n0) if len(adj[v]) == 1]
            for l in leaves:
                s = adj[l][0]  # unique support

                # Build T-l (remove vertex l)
                keep_Tl = set(range(n0)) - {l}
                m_Tl, adj_Tl = subgraph_induced(n0, adj, keep_Tl)
                f = indpoly(m_Tl, adj_Tl)

                # Build T-{l,s} (remove both l and s)
                keep_Tls = set(range(n0)) - {l, s}
                if len(keep_Tls) == 0:
                    continue  # trivial (n=2)
                m_Tls, adj_Tls = subgraph_induced(n0, adj, keep_Tls)
                g = indpoly(m_Tls, adj_Tls)

                lr_total += 1

                # Check d_leaf<=1 for T-l
                if not is_dleaf_le1(m_Tl, adj_Tl):
                    n_lr_Tl_d1_fail += 1

                # Check d_leaf<=1 for T-{l,s}
                if not is_dleaf_le1(m_Tls, adj_Tls):
                    n_lr_Tls_d1_fail += 1

                # LC check for T-l
                if not is_lc(f):
                    n_lr_Tl_lc += 1

                # LC check for T-{l,s}
                if not is_lc(g):
                    n_lr_Tls_lc += 1

                # Cross-term check: is h = f + x*g LC? (should be = I(T))
                h_lc, min_ratio = cross_term_lc(f, g)
                if not h_lc:
                    n_lr_cross += 1

                # Also check pure cross-term sign (2f_k g_{k-1} - f_{k-1}g_k - g_{k-2}f_{k+1})
                # for each k
                # (separate from LC of h, which includes the self-LC terms of f and g)

        proc.wait()

        total += n_total
        dleaf1_total += n_dleaf1
        lc_fail += n_lc_fail
        ulc_fail += n_ulc_fail
        lr_Tl_lc_fail += n_lr_Tl_lc
        lr_Tls_lc_fail += n_lr_Tls_lc
        lr_h_lc_fail += n_lr_cross
        lr_Tls_dleaf1_fail += n_lr_Tls_d1_fail
        lr_Tl_dleaf1_fail += n_lr_Tl_d1_fail

        print(f"{n:>4} {n_total:>8} {n_dleaf1:>8} {n_lc_fail:>8} {n_ulc_fail:>8} "
              f"{n_lr_Tl_lc:>9} {n_lr_Tls_lc:>10} {n_lr_cross:>14}", flush=True)

    print()
    print("=" * 90)
    print(f"TOTALS: trees={total:,} d_leaf<=1={dleaf1_total:,}")
    print(f"  LC failures: {lc_fail}")
    print(f"  ULC failures: {ulc_fail}")
    if ulc_fail_first:
        print(f"    First ULC failure: n={ulc_fail_first[0]}, poly={ulc_fail_first[1]}")
    print(f"  Leaf-removal ({lr_total:,} pairs):")
    print(f"    T-l NOT d_leaf<=1: {lr_Tl_dleaf1_fail} (expected sometimes)")
    print(f"    T-{{l,s}} NOT d_leaf<=1: {lr_Tls_dleaf1_fail} (should be 0)")
    print(f"    T-l NOT LC: {lr_Tl_lc_fail}")
    print(f"    T-{{l,s}} NOT LC: {lr_Tls_lc_fail}")
    print(f"    h = I(T-l)+x*I(T-{{l,s}}) NOT LC: {lr_h_lc_fail} (should match lc_fail)")

if __name__ == "__main__":
    main()
