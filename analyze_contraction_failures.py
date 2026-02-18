#!/usr/bin/env python3
"""Analyze the core edge contractions that DECREASE μ/n.

Questions:
1. What do the failing trees look like?
2. Is the decrease always small?
3. Does contraction still move μ/n TOWARD the spider bound?
4. Is there a modified contraction that always increases μ/n?
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
    new_n = n - 1
    mapping = {}
    for w in range(n):
        if w == v:
            mapping[w] = u if u < v else u - 1
        elif w < v:
            mapping[w] = w
        else:
            mapping[w] = w - 1

    u_new = u if u < v else u - 1
    seen_edges = set()
    new_adj = [[] for _ in range(new_n)]
    for w in range(n):
        for x in adj[w]:
            a = mapping[w]
            b = mapping[x]
            if a == b:
                continue
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


def get_tree_structure(n, adj):
    """Describe the tree structure."""
    deg_seq = sorted([len(adj[v]) for v in range(n)], reverse=True)
    leaves = [v for v in range(n) if len(adj[v]) == 1]
    high_deg = [v for v in range(n) if len(adj[v]) >= 3]
    return {
        "n": n,
        "leaves": len(leaves),
        "high_deg_count": len(high_deg),
        "max_deg": deg_seq[0],
        "deg_seq": deg_seq[:8],
        "is_spider": len(high_deg) <= 1,
        "is_caterpillar": all(
            any(len(adj[u]) <= 2 for u in adj[v]) or len(adj[v]) <= 2
            for v in range(n)
        ),
    }


def main():
    print("CONTRACTION FAILURE ANALYSIS", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)

    # Collect all failures at n=11-14 for detailed analysis
    all_failures = []

    for n in range(9, 15):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        failures = []

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

            # Find all core edges
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

                if cratio < ratio - 1e-12:
                    decrease = ratio - cratio
                    tree_info = get_tree_structure(tn, adj_data)
                    contracted_info = get_tree_structure(cn, cadj)
                    failures.append({
                        "tree": tree_info,
                        "contracted": contracted_info,
                        "edge": (u, v),
                        "u_deg": len(adj_data[u]),
                        "v_deg": len(adj_data[v]),
                        "mu_before": mu,
                        "mu_after": cmu,
                        "ratio_before": ratio,
                        "ratio_after": cratio,
                        "decrease": decrease,
                    })

        all_failures.extend(failures)

        if failures:
            print(f"n={n}: {len(failures)} failures", flush=True)
            # Show worst failure
            worst = max(failures, key=lambda f: f["decrease"])
            print(f"  Worst: decrease={worst['decrease']:.6f}", flush=True)
            print(f"    Tree: deg={worst['tree']['deg_seq']}, "
                  f"spider={worst['tree']['is_spider']}", flush=True)
            print(f"    Edge degs: ({worst['u_deg']}, {worst['v_deg']})",
                  flush=True)
            print(f"    Contracted: deg={worst['contracted']['deg_seq']}, "
                  f"spider={worst['contracted']['is_spider']}", flush=True)

    print(flush=True)
    print("=" * 70, flush=True)
    print("PATTERN ANALYSIS", flush=True)
    print("=" * 70, flush=True)

    # Q1: What degree pairs cause failures?
    print(flush=True)
    print("Edge degree pairs causing failures:", flush=True)
    from collections import Counter
    deg_pairs = Counter()
    for f in all_failures:
        pair = tuple(sorted([f["u_deg"], f["v_deg"]]))
        deg_pairs[pair] += 1
    for pair, count in deg_pairs.most_common(10):
        print(f"  ({pair[0]}, {pair[1]}): {count}", flush=True)

    # Q2: Are failures always spider→spider?
    print(flush=True)
    print("Spider status:", flush=True)
    spider_to_spider = sum(1 for f in all_failures
                          if f["tree"]["is_spider"] and f["contracted"]["is_spider"])
    spider_to_non = sum(1 for f in all_failures
                       if f["tree"]["is_spider"] and not f["contracted"]["is_spider"])
    non_to_spider = sum(1 for f in all_failures
                       if not f["tree"]["is_spider"] and f["contracted"]["is_spider"])
    non_to_non = sum(1 for f in all_failures
                    if not f["tree"]["is_spider"] and not f["contracted"]["is_spider"])
    print(f"  Spider→Spider: {spider_to_spider}", flush=True)
    print(f"  Spider→Non: {spider_to_non}", flush=True)
    print(f"  Non→Spider: {non_to_spider}", flush=True)
    print(f"  Non→Non: {non_to_non}", flush=True)

    # Q3: Skip d_leaf recheck (need stored adjacency)

    # Q4: How large are the decreases relative to the gap to n/3?
    print(flush=True)
    print("Decrease vs gap to n/3:", flush=True)
    for f in sorted(all_failures, key=lambda f: -f["decrease"])[:10]:
        gap_before = f["tree"]["n"] / 3 - f["mu_before"]
        gap_after = f["contracted"]["n"] / 3 - f["mu_after"]
        print(f"  n={f['tree']['n']}: decrease={f['decrease']:.6f}, "
              f"gap_before={gap_before:.4f}, gap_after={gap_after:.4f}, "
              f"edge=({f['u_deg']},{f['v_deg']}), "
              f"deg={f['tree']['deg_seq'][:5]}", flush=True)

    # Q5: Does μ itself (not μ/n) always increase?
    print(flush=True)
    print("Does RAW μ (not μ/n) always increase?", flush=True)
    mu_increases = sum(1 for f in all_failures if f["mu_after"] > f["mu_before"])
    mu_decreases = sum(1 for f in all_failures if f["mu_after"] < f["mu_before"])
    print(f"  Among {len(all_failures)} ratio-decrease cases:", flush=True)
    print(f"    Raw μ increases: {mu_increases}", flush=True)
    print(f"    Raw μ decreases: {mu_decreases}", flush=True)

    # Q6: Does the gap n/3 - μ always decrease (good)?
    print(flush=True)
    print("Does gap (n/3 - μ) always decrease (approach 0)?", flush=True)
    gap_decreases = 0
    gap_increases = 0
    for f in all_failures:
        gap_before = f["tree"]["n"] / 3 - f["mu_before"]
        gap_after = f["contracted"]["n"] / 3 - f["mu_after"]
        if gap_after < gap_before - 1e-12:
            gap_decreases += 1
        elif gap_after > gap_before + 1e-12:
            gap_increases += 1
    print(f"  Gap decreases (good): {gap_decreases}", flush=True)
    print(f"  Gap increases (bad): {gap_increases}", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
