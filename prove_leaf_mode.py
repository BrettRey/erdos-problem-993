#!/usr/bin/env python3
"""Investigate Leaf-Mode Inequality proof strategies.

Strategy: for a leaf v with support u, I(T;x) = I(T-v;x) + x*I(T-N[v];x).
Track how mode relates to ℓ through this recursion.

Also check: for the tightest cases, what's the exact polynomial structure?
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


def remove_vertex(n, adj, v):
    """Return (new_n, new_adj) with vertex v removed."""
    mapping = {}
    new_idx = 0
    for i in range(n):
        if i != v:
            mapping[i] = new_idx
            new_idx += 1
    new_n = n - 1
    new_adj = [[] for _ in range(new_n)]
    for i in range(n):
        if i == v:
            continue
        for j in adj[i]:
            if j != v:
                new_adj[mapping[i]].append(mapping[j])
    return new_n, new_adj


def remove_closed_nbhd(n, adj, v):
    """Return (new_n, new_adj) with N[v] removed. May be disconnected."""
    to_remove = set([v]) | set(adj[v])
    mapping = {}
    new_idx = 0
    for i in range(n):
        if i not in to_remove:
            mapping[i] = new_idx
            new_idx += 1
    new_n = n - len(to_remove)
    if new_n == 0:
        return 0, []
    new_adj = [[] for _ in range(new_n)]
    for i in range(n):
        if i in to_remove:
            continue
        for j in adj[i]:
            if j not in to_remove:
                new_adj[mapping[i]].append(mapping[j])
    return new_n, new_adj


def poly_add(p, q):
    """Add two polynomials (as lists)."""
    n = max(len(p), len(q))
    result = [0] * n
    for i, c in enumerate(p):
        result[i] += c
    for i, c in enumerate(q):
        result[i] += c
    return result


def poly_shift(p):
    """Multiply by x (shift right)."""
    return [0] + list(p)


def main():
    print("LEAF-MODE INEQUALITY: RECURSION ANALYSIS", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)

    # For each high-mode tree, check the recursion I(T) = I(T-v) + x*I(T-N[v])
    # at the tightest leaf v. Track how mode/ℓ relate.

    for n in range(8, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 1
        n_checked = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode <= threshold:
                continue

            leaves_count = sum(1 for v in range(tn) if len(adj_data[v]) == 1)
            gap = leaves_count - mode

            # Only look at tightest cases (gap <= 4)
            if gap > 4:
                continue

            n_checked += 1
            if n_checked > 3:
                continue

            # Find a leaf and analyze recursion
            leaf = next(v for v in range(tn) if len(adj_data[v]) == 1)
            support = adj_data[leaf][0]

            # T - v
            n1, adj1 = remove_vertex(tn, adj_data, leaf)
            if n1 > 0:
                poly1 = independence_poly(n1, adj1)
                mode1 = max(range(len(poly1)), key=lambda i: poly1[i])
                l1 = sum(1 for v in range(n1) if len(adj1[v]) == 1)
            else:
                poly1, mode1, l1 = [1], 0, 0

            # T - N[v] = T - {v, support}
            n2, adj2 = remove_closed_nbhd(tn, adj_data, leaf)
            if n2 > 0:
                poly2 = independence_poly(n2, adj2)
                mode2 = max(range(len(poly2)), key=lambda i: poly2[i])
                l2 = sum(1 for v in range(n2) if len(adj2[v]) == 1)
            else:
                poly2, mode2, l2 = [1], 0, 0

            # Verify: poly = poly1 + x * poly2
            reconstructed = poly_add(poly1, poly_shift(poly2))
            assert reconstructed == poly, f"Reconstruction failed!"

            deg_seq = sorted([len(adj_data[v]) for v in range(tn)], reverse=True)
            print(f"n={tn}, ℓ={leaves_count}, mode={mode}, gap={gap}", flush=True)
            print(f"  deg_seq={deg_seq[:8]}", flush=True)
            print(f"  T-v: n={n1}, ℓ={l1}, mode={mode1}", flush=True)
            print(f"  T-N[v]: n={n2}, ℓ={l2}, mode={mode2}", flush=True)
            print(f"  I(T) at mode: {poly[mode]}", flush=True)
            print(f"  I(T-v) at mode: {poly1[mode] if mode < len(poly1) else 0}", flush=True)
            print(f"  x*I(T-N[v]) at mode: {poly2[mode-1] if mode-1 < len(poly2) and mode >= 1 else 0}", flush=True)
            print(flush=True)

        elapsed = time.time() - t0

    # Now check: for trees where mode > ℓ (which shouldn't happen for high-mode),
    # what's the i_ℓ value? This tells us about the growth beyond ℓ.
    print("=" * 60, flush=True)
    print("CHECKING i_ℓ vs i_{ℓ+1} FOR ALL TREES", flush=True)
    print("=" * 60, flush=True)

    for n in range(8, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_above = 0  # trees with mode > ℓ
        max_mode_minus_l = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            leaves_count = sum(1 for v in range(tn) if len(adj_data[v]) == 1)

            if mode > leaves_count:
                n_above += 1
                diff = mode - leaves_count
                max_mode_minus_l = max(max_mode_minus_l, diff)

                # What is the ratio i_{ℓ+1}/i_ℓ?
                if leaves_count < len(poly) and leaves_count + 1 < len(poly):
                    ratio = poly[leaves_count + 1] / poly[leaves_count] if poly[leaves_count] > 0 else float('inf')
                else:
                    ratio = 0

        elapsed = time.time() - t0
        threshold = n // 3 + 1
        print(f"n={n:2d}: {n_above} trees with mode > ℓ, "
              f"max(mode-ℓ)={max_mode_minus_l}, "
              f"threshold={threshold} ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
