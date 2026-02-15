#!/usr/bin/env python3
"""
For the n=18 LNP failure tree, check:
1. Does max_priv >= 2? (private neighbors, not just leaf-neighbors)
2. Can we still do a swap (remove u, add v, w)?
3. Are there other LNP failures at n=17..20?
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


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def find_maximal_is_below_mode(n, adj, mode):
    nbr = [set(adj[v]) for v in range(n)]
    results = []

    def backtrack(v, current, forbidden):
        if v == n:
            k = len(current)
            if k >= mode:
                return
            s = frozenset(current)
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                results.append((s, k))
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return results


# First: analyze the known failure tree
print("=" * 75)
print("ANALYZING n=18 LNP FAILURE TREE")
print("=" * 75)

g6 = "Q??????_A?C?C?a?G_@C?CO?N_?"
n, adj = parse_graph6(g6)
deg = [len(adj[v]) for v in range(n)]
nbr = [set(adj[v]) for v in range(n)]

S = frozenset([0, 1, 2, 3, 17])
k = len(S)

print(f"S = {sorted(S)}, k = {k}")
print(f"Tree: spider with center 17 (deg 5), 4 arms of length 2, leaf 8")

# Compute private and shared neighbors for each u in S
for u in sorted(S):
    privates = []
    shared_nbrs = []
    for v in nbr[u]:
        if v not in S:
            s_nbrs = nbr[v] & S
            if len(s_nbrs) == 1:
                privates.append(v)
            else:
                shared_nbrs.append(v)
    print(f"  u={u} (deg={deg[u]}): "
          f"private={privates} ({len(privates)}) "
          f"shared={shared_nbrs} ({len(shared_nbrs)})")

# Check max_priv
max_priv = 0
for u in S:
    p = sum(1 for v in nbr[u] if v not in S and len(nbr[v] & S) == 1)
    max_priv = max(max_priv, p)
print(f"\nmax_priv = {max_priv}")

if max_priv >= 2:
    print("max_priv >= 2: swap exists (not via leaf-neighbors, but via private neighbors)")
else:
    print("max_priv < 2: NO SWAP via private neighbors!")
    # Check: can we do a non-private swap?
    # A swap: remove u, add v, w where v, w not in S, v, w not adjacent,
    # and S' = (S \ {u}) ∪ {v, w} is independent
    print("\nChecking all possible swaps:")
    for u in sorted(S):
        freed = nbr[u] - S  # vertices freed by removing u
        candidates = set()
        for v in range(n):
            if v not in S and v != u:
                # Check if v can be added to S \ {u}
                s_minus_u = S - {u}
                if not (nbr[v] & s_minus_u):
                    candidates.add(v)
        # candidates: vertices that can be in S \ {u}
        # Need to find two non-adjacent candidates
        cand_list = sorted(candidates)
        pairs = []
        for i in range(len(cand_list)):
            for j in range(i + 1, len(cand_list)):
                v, w = cand_list[i], cand_list[j]
                if w not in nbr[v]:
                    pairs.append((v, w))
        if pairs:
            print(f"  Remove u={u}: {len(pairs)} valid (v,w) pairs")
            for v, w in pairs[:5]:
                print(f"    ({v}, {w})")
        else:
            print(f"  Remove u={u}: NO valid (v,w) pairs!")

# Second: scan for ALL LNP failures at n=17..18
print(f"\n{'='*75}")
print("SCANNING FOR ALL LNP FAILURES at n=17..18")
print(f"{'='*75}")

for n_scan in range(17, 19):
    t0 = time.time()
    cmd = f"/opt/homebrew/bin/geng {n_scan} {n_scan-1}:{n_scan-1} -c -q"
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    lines = [l for l in proc.stdout.strip().split("\n") if l]

    total_cases = 0
    lnp_fails = 0
    priv_fails = 0  # max_priv < 2
    swap_fails = 0  # no valid swap at all

    for tidx, line in enumerate(lines):
        tn, adj_data = parse_graph6(line)
        deg_t = [len(adj_data[v]) for v in range(tn)]
        leaves = {v for v in range(tn) if deg_t[v] == 1}
        nbr_t = [set(adj_data[v]) for v in range(tn)]

        poly = independence_poly(tn, adj_data)
        mode_t = compute_mode(poly)

        mis_below = find_maximal_is_below_mode(tn, adj_data, mode_t)

        for s, k in mis_below:
            total_cases += 1

            # LNP check
            max_leaf = max(
                sum(1 for w in nbr_t[u] if w in leaves and w not in s)
                for u in s
            )
            if max_leaf < 2:
                lnp_fails += 1

                # Private neighbor check
                max_p = 0
                for u in s:
                    p = sum(1 for v in nbr_t[u]
                            if v not in s and len(nbr_t[v] & s) == 1)
                    max_p = max(max_p, p)

                if max_p < 2:
                    priv_fails += 1

                    # Swap check
                    has_swap = False
                    for u in s:
                        s_minus = s - {u}
                        cands = [v for v in range(tn)
                                 if v not in s and v != u
                                 and not (nbr_t[v] & s_minus)]
                        for i in range(len(cands)):
                            for j in range(i + 1, len(cands)):
                                if cands[j] not in nbr_t[cands[i]]:
                                    has_swap = True
                                    break
                            if has_swap:
                                break
                        if has_swap:
                            break

                    if not has_swap:
                        swap_fails += 1
                        print(f"  *** NO SWAP: n={tn} tree={tidx} "
                              f"k={k} mode={mode_t} S={sorted(s)} "
                              f"g6={line.strip()}")

    elapsed = time.time() - t0
    print(f"n={n_scan}: {total_cases} cases, "
          f"LNP fails: {lnp_fails}, "
          f"priv<2: {priv_fails}, "
          f"no swap: {swap_fails} "
          f"({elapsed:.1f}s)")
