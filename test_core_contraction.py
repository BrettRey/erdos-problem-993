#!/usr/bin/env python3
"""Test: does contracting a core edge always increase μ?

A "core edge" is an edge {u,v} where neither u nor v is a leaf.
Contracting it: merge u and v into one vertex, preserving all other edges.
The result has n-1 vertices.

If μ(T/e) > μ(T) for every core edge e, then repeatedly contracting
core edges transforms any tree into a spider while monotonically
increasing μ. This would prove spider extremality.

But we compare μ/n (ratio) since n changes.

Actually, a better comparison: we compare μ directly, since both
trees have different n. We want to show that the RATIO μ/n increases:
μ(T/e)/(n-1) >= μ(T)/n.
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
    total = sum(poly)
    weighted = sum(k * poly[k] for k in range(len(poly)))
    return weighted / total if total else 0.0


def contract_edge(n, adj, u, v):
    """Contract edge {u,v}: merge v into u, remove v."""
    # New vertex set: all except v
    # u inherits v's neighbors (except u itself)
    new_n = n - 1
    # Map old vertices to new: v maps to u, others keep index (shifted)
    old_to_new = {}
    for w in range(n):
        if w == v:
            old_to_new[w] = old_to_new.get(u, u)  # v -> u
        elif w < v:
            old_to_new[w] = w
        else:
            old_to_new[w] = w - 1

    # Actually, let's be more careful
    new_adj = [[] for _ in range(new_n)]
    # Remap: vertex v becomes vertex u, vertices > v shift down by 1
    mapping = {}
    for w in range(n):
        if w == v:
            mapping[w] = mapping.get(u, u if u < v else u)
        elif w < v:
            mapping[w] = w
        else:
            mapping[w] = w - 1

    # Reindex u properly
    u_new = u if u < v else u - 1

    seen_edges = set()
    for w in range(n):
        for x in adj[w]:
            if w == v or x == v:
                # Edge involving v: redirect to u
                a = mapping[w]
                b = mapping[x]
            else:
                a = mapping[w]
                b = mapping[x]

            if a == b:
                continue  # self-loop from contraction
            edge = (min(a, b), max(a, b))
            if edge not in seen_edges:
                seen_edges.add(edge)
                new_adj[a].append(b)
                new_adj[b].append(a)

    return new_n, new_adj


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


def main():
    print("CORE EDGE CONTRACTION: DOES μ/n INCREASE?", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)

    total_contractions = 0
    ratio_increases = 0
    ratio_decreases = 0
    worst_decrease = 0

    for n in range(6, 17):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_tested = 0
        n_increase = 0
        n_decrease = 0
        worst_n = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            if not check_dleaf(tn, adj_data):
                continue

            leaves = set()
            for v in range(tn):
                if len(adj_data[v]) == 1:
                    leaves.add(v)

            poly = independence_poly(tn, adj_data)
            mu = poly_mean(poly)
            ratio = mu / tn

            # Find all core edges (neither endpoint is a leaf)
            core_edges = set()
            for v in range(tn):
                if v in leaves:
                    continue
                for u in adj_data[v]:
                    if u in leaves:
                        continue
                    if u > v:
                        core_edges.add((v, u))

            for u, v in core_edges:
                cn, cadj = contract_edge(tn, adj_data, u, v)
                cpoly = independence_poly(cn, cadj)
                cmu = poly_mean(cpoly)
                cratio = cmu / cn

                n_tested += 1
                total_contractions += 1

                if cratio > ratio + 1e-12:
                    n_increase += 1
                    ratio_increases += 1
                elif cratio < ratio - 1e-12:
                    n_decrease += 1
                    ratio_decreases += 1
                    decrease = ratio - cratio
                    if decrease > worst_n:
                        worst_n = decrease
                    if decrease > worst_decrease:
                        worst_decrease = decrease
                else:
                    n_increase += 1  # count ties as OK
                    ratio_increases += 1

        elapsed = time.time() - t0
        status = "ALL OK" if n_decrease == 0 else f"DECREASES: {n_decrease}"
        wd = f", worst={worst_n:.6f}" if n_decrease > 0 else ""
        print(f"n={n:2d}: {n_tested:5d} contractions, "
              f"inc={n_increase}, dec={n_decrease} | "
              f"{status}{wd} ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print(f"TOTAL: {total_contractions} contractions, "
          f"{ratio_increases} increases, {ratio_decreases} decreases",
          flush=True)
    if ratio_decreases > 0:
        print(f"Worst decrease: {worst_decrease:.8f}", flush=True)
        print("CORE CONTRACTION DOES NOT ALWAYS INCREASE μ/n!", flush=True)
    else:
        print("Core contraction ALWAYS increases μ/n!", flush=True)
        print("=> Spider extremality follows by induction.", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
