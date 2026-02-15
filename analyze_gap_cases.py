#!/usr/bin/env python3
"""Analyze the 43 gap cases at n=17 where neither pigeonhole nor degree bound covers PNP.

These are maximal IS of size k=6, mode=7, max_deg_in_S=6, n=17 < 3k=18.
What structural property gives priv >= 2?

For each gap case: identify the tree structure, what vertex provides max_priv,
and what mechanism gives it.
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
        n = ord(s[0]) - 63
        idx = 1
    else:
        idx = 1
        n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63)
            idx += 1
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
                adj[i].append(j)
                adj[j].append(i)
            bit_idx += 1
    return n, adj


def find_all_maximal_is(n, adj):
    nbr = [set(adj[v]) for v in range(n)]
    result = []
    def backtrack(v, current, forbidden):
        if v == n:
            s = frozenset(current)
            for w in range(n):
                if w not in s and not (nbr[w] & s):
                    return
            result.append(s)
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])
    backtrack(0, [], set())
    return result


def compute_priv(u, s, nbr):
    count = 0
    for v in nbr[u]:
        if v not in s:
            if len(nbr[v] & s) == 1:
                count += 1
    return count


def main():
    print("ANALYSIS OF GAP CASES (neither pigeonhole nor degree bound)")
    print("=" * 70)
    print()

    n = 17
    t0 = time.time()
    cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    lines = [l for l in proc.stdout.strip().split("\n") if l]

    gap_cases = []  # (g6, tree_info, IS_info)

    for line in lines:
        tn, adj_data = parse_graph6(line)
        poly = independence_poly(tn, adj_data)
        mode = max(range(len(poly)), key=lambda i: poly[i])
        nbr = [set(adj_data[v]) for v in range(tn)]

        all_mis = find_all_maximal_is(tn, adj_data)
        for s in all_mis:
            k = len(s)
            if k >= mode:
                continue

            ph = (tn >= 3 * k)
            max_deg = max(len(adj_data[u]) for u in s)
            deg = (max_deg >= k + 1)

            if not ph and not deg:
                # Gap case!
                privs = {u: compute_priv(u, s, nbr) for u in s}
                max_priv = max(privs.values())
                max_priv_v = max(privs, key=privs.get)
                deg_seq = sorted([len(adj_data[v]) for v in range(tn)], reverse=True)

                # How many vertices have deg >= 2 in S?
                high_deg_in_s = sum(1 for u in s if len(adj_data[u]) >= 2)

                # What's the sum of degrees in S?
                deg_sum_s = sum(len(adj_data[u]) for u in s)

                # Number of edges between S and V\S
                edges_sv = sum(len(nbr[u] - s) for u in s)

                # Number of edges within V\S
                vs = set(range(tn)) - s
                edges_vv = sum(1 for u in vs for v in nbr[u] if v in vs) // 2

                gap_cases.append({
                    'g6': line,
                    'k': k,
                    'mode': mode,
                    'max_deg_in_s': max_deg,
                    'max_priv': max_priv,
                    'max_priv_v': max_priv_v,
                    'deg_v': len(adj_data[max_priv_v]),
                    'privs': privs,
                    'deg_seq': deg_seq[:8],
                    'high_deg_in_s': high_deg_in_s,
                    'deg_sum_s': deg_sum_s,
                    'edges_sv': edges_sv,
                    'edges_vv': edges_vv,
                })

    elapsed = time.time() - t0
    print(f"Found {len(gap_cases)} gap cases at n={n} ({elapsed:.1f}s)\n")

    # Summarize
    print("Distribution of max_priv:")
    from collections import Counter
    mp_dist = Counter(c['max_priv'] for c in gap_cases)
    for mp in sorted(mp_dist):
        print(f"  max_priv={mp}: {mp_dist[mp]} cases")

    print(f"\nDistribution of deg(max_priv vertex):")
    dv_dist = Counter(c['deg_v'] for c in gap_cases)
    for dv in sorted(dv_dist):
        print(f"  deg={dv}: {dv_dist[dv]} cases")

    print(f"\nDistribution of edges_sv (edges between S and V\\S):")
    esv_dist = Counter(c['edges_sv'] for c in gap_cases)
    for esv in sorted(esv_dist):
        print(f"  edges_sv={esv}: {esv_dist[esv]} cases")

    print(f"\nDistribution of edges_vv (edges within V\\S):")
    evv_dist = Counter(c['edges_vv'] for c in gap_cases)
    for evv in sorted(evv_dist):
        print(f"  edges_vv={evv}: {evv_dist[evv]} cases")

    # Detailed look at first few
    print(f"\nDetailed first 5 cases:")
    for i, c in enumerate(gap_cases[:5]):
        print(f"\n  Case {i+1}: k={c['k']}, mode={c['mode']}, max_deg_in_S={c['max_deg_in_s']}")
        print(f"    max_priv={c['max_priv']} at vertex {c['max_priv_v']} (deg={c['deg_v']})")
        print(f"    deg_seq={c['deg_seq']}")
        print(f"    edges: S-VS={c['edges_sv']}, VS-VS={c['edges_vv']}")
        print(f"    privs: {dict(c['privs'])}")

    # Check: does every gap case have a vertex u in S where
    # deg(u) + priv(u) > some threshold?
    print(f"\n\nKey question: what structural property gives priv >= 2?")

    # Check: number of leaves adjacent to S
    n_leaves_adj_s = []
    for c in gap_cases:
        tn, adj_data = parse_graph6(c['g6'])
        nbr = [set(adj_data[v]) for v in range(tn)]
        # Reconstruct s from privs
        s = frozenset(c['privs'].keys())
        leaves_adj_s = 0
        for v in range(tn):
            if v not in s and len(adj_data[v]) == 1:
                if nbr[v] & s:
                    leaves_adj_s += 1
        n_leaves_adj_s.append(leaves_adj_s)

    print(f"\nLeaves adjacent to S: min={min(n_leaves_adj_s)}, max={max(n_leaves_adj_s)}")
    las_dist = Counter(n_leaves_adj_s)
    for las in sorted(las_dist):
        print(f"  leaves_adj_S={las}: {las_dist[las]} cases")

    # Check: does every gap case have total leaves > 2k?
    print(f"\nTotal leaves in tree:")
    n_leaves = []
    for c in gap_cases:
        tn, adj_data = parse_graph6(c['g6'])
        leaves = sum(1 for v in range(tn) if len(adj_data[v]) == 1)
        n_leaves.append(leaves)
    print(f"  min={min(n_leaves)}, max={max(n_leaves)}")
    nl_dist = Counter(n_leaves)
    for nl in sorted(nl_dist):
        print(f"  total_leaves={nl}: {nl_dist[nl]} cases")

    # Check: for the max_priv vertex, is it ALWAYS adjacent to leaves outside S?
    print(f"\nFor max_priv vertex: how many of its private neighbors are leaves?")
    leaf_priv_counts = []
    for c in gap_cases:
        tn, adj_data = parse_graph6(c['g6'])
        nbr = [set(adj_data[v]) for v in range(tn)]
        s = frozenset(c['privs'].keys())
        u = c['max_priv_v']
        leaf_priv = 0
        for v in nbr[u]:
            if v not in s and len(nbr[v] & s) == 1 and len(adj_data[v]) == 1:
                leaf_priv += 1
        leaf_priv_counts.append(leaf_priv)
    print(f"  min={min(leaf_priv_counts)}, max={max(leaf_priv_counts)}")
    lpc_dist = Counter(leaf_priv_counts)
    for lpc in sorted(lpc_dist):
        print(f"  leaf_privates={lpc}: {lpc_dist[lpc]} cases")


if __name__ == "__main__":
    main()
