#!/usr/bin/env python3
"""
Verify the Steiner-tree peeling approach for the gap-formula WHNC.

n/3 - mu(T) = sum_{C\A}(1/3-P) + (1/2)*sum_A(1/3-P)

Define for S subset H_core:
  F_gap(S) = sum_{u in N_C(S) cap C\A} (1/3-P(u))
           + (1/2)*sum_{u in N_C(S) cap A} (1/3-P(u))
           - sum_{h in S} (P(h)-1/3)

Since n/3-mu >= F_gap(H_core) (extra supply from outside N_C(H_core) >= 0),
proving F_gap(S) >= 0 for all S subset H_core proves Conjecture A.

Steiner peeling: take h = leaf of Steiner tree of S.
  - h has degree >= 2 in T (h in H_core = C\A, no leaf children)
  - h has private C-neighbors (degree in Steiner tree = 1, degree in T >= 2)
  - Check: M_gap(h, S) = F_gap(S) - F_gap(S\{h}) >= 0?

Also checks:
  1. P(v) < 1/3 for all v in A (support vertices)
  2. P(v) = (1 - P(support))/2 for leaves (decimation identity)
"""
from __future__ import annotations
import subprocess
import sys
from fractions import Fraction

TOL = 1e-12

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
        leaf_children = sum(1 for u in adj[v] if len(adj[u]) == 1)
        if leaf_children > 1:
            return False
    return True

def hard_core_probs(n: int, adj: list[list[int]]) -> list[float]:
    """Belief propagation at fugacity 1."""
    msg = [[1.0] * n for _ in range(n)]
    for _ in range(200):
        new_msg = [[1.0] * n for _ in range(n)]
        for v in range(n):
            for u in adj[v]:
                prod = 1.0
                for w in adj[u]:
                    if w != v:
                        prod *= 1.0 / (1.0 + msg[w][u])
                new_msg[u][v] = prod
        msg = new_msg
    probs = []
    for v in range(n):
        m_in = 1.0
        for u in adj[v]:
            m_in *= 1.0 / (1.0 + msg[u][v])
        probs.append(m_in / (1.0 + m_in))
    return probs

def steiner_tree_leaves(S: list[int], adj: list[list[int]], n: int) -> list[int]:
    """
    Find leaves of the Steiner tree of S in the tree.
    A vertex h in S is a Steiner leaf iff it has exactly one Steiner-tree
    neighbor (i.e., the path to all other S-vertices goes through one direction).
    Equivalently: remove h from T; the number of components containing
    S\{h} members is 1.
    """
    if len(S) <= 1:
        return list(S)
    S_set = set(S)
    leaves = []
    for h in S:
        # BFS from h, counting S-components excluding h
        visited = [False] * n
        visited[h] = True
        # Count how many different subtrees of h contain S\{h} members
        components_with_S = 0
        for start_neighbor in adj[h]:
            if visited[start_neighbor]:
                continue
            # BFS component
            queue = [start_neighbor]
            visited[start_neighbor] = True
            has_S = start_neighbor in S_set
            qi = 1
            while qi <= len(queue) - 1:
                qi += 1
                v = queue[qi - 1]
                for w in adj[v]:
                    if not visited[w]:
                        visited[w] = True
                        queue.append(w)
                        if w in S_set:
                            has_S = True
            if has_S:
                components_with_S += 1
        if components_with_S <= 1:
            leaves.append(h)
    return leaves

def private_C_neighbors(h: int, S: list[int], adj: list[list[int]],
                         leaves: set[int]) -> list[int]:
    """Private C-neighbors of h w.r.t. S: in C (not a leaf), not in N(S\{h})."""
    S_minus_h = set(S) - {h}
    N_S_minus_h = set()
    for s in S_minus_h:
        N_S_minus_h.update(adj[s])
    return [u for u in adj[h] if u not in leaves and u not in N_S_minus_h]

def marginal_gap(h: int, S: list[int], adj: list[list[int]],
                  leaves: set[int], A: set[int], probs: list[float]) -> float:
    """
    M_gap(h, S) = F_gap(S) - F_gap(S\{h})
               = sum_{u private C-neighbor of h} s_gap(u) - demand(h)
    where s_gap(u) = 1/3-P(u) for u in C\A, (1/2)(1/3-P(u)) for u in A.
    """
    priv = private_C_neighbors(h, S, adj, leaves)
    supply = 0.0
    for u in priv:
        if u in A:
            supply += 0.5 * (1/3 - probs[u])
        else:
            supply += (1/3 - probs[u])
    demand = probs[h] - 1/3
    return supply - demand

def f_gap(S: list[int], adj: list[list[int]], leaves: set[int],
           A: set[int], probs: list[float]) -> float:
    """F_gap(S) = sum_{N_C(S) supply} - sum_{S demand}."""
    S_set = set(S)
    # N_C(S): C-neighbors of S not in S (= not leaves, not in S)
    N_C_S = set()
    for h in S:
        for u in adj[h]:
            if u not in leaves and u not in S_set:
                N_C_S.add(u)
    supply = sum(
        (0.5 if u in A else 1.0) * (1/3 - probs[u])
        for u in N_C_S
    )
    demand = sum(probs[h] - 1/3 for h in S)
    return supply - demand

def main():
    min_n, max_n = 3, 18
    geng = "/opt/homebrew/bin/geng"

    # Counters
    total_trees = 0
    total_dleaf = 0
    A_heavy_fails = 0          # P(v) >= 1/3 for v in A
    leaf_identity_fails = 0    # P(l) != (1-P(s))/2
    steiner_M_negative = 0     # Steiner leaf has M_gap < 0
    steiner_M_zero = 0         # M_gap = 0 (tight)
    F_gap_negative = 0         # F_gap(H_core) < 0
    n_div3_negative = 0        # n/3 - mu < 0

    worst_M = 0.0
    worst_M_tree = None

    for n in range(min_n, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout
        n_trees = 0
        n_dleaf = 0
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_trees += 1
            if not is_dleaf_le1(n0, adj):
                continue
            n_dleaf += 1
            probs = hard_core_probs(n0, adj)
            mu = sum(probs)

            # Partition vertices
            leaves_T = {v for v in range(n0) if len(adj[v]) == 1}
            A = {v for v in range(n0) if v not in leaves_T
                 and any(u in leaves_T for u in adj[v])}
            C_A = {v for v in range(n0) if v not in leaves_T and v not in A}

            # Check 1: P(v) < 1/3 for all v in A
            for v in A:
                if probs[v] >= 1/3 - TOL:
                    A_heavy_fails += 1

            # Check 2: P(l) = (1-P(s))/2 for leaves
            for v in leaves_T:
                s = adj[v][0]  # unique support
                expected = (1 - probs[s]) / 2
                if abs(probs[v] - expected) > 1e-9:
                    leaf_identity_fails += 1

            # Check 3: n/3 - mu >= 0
            gap_formula = n0/3 - mu
            if gap_formula < -TOL:
                n_div3_negative += 1

            # H_core: heavy non-support non-leaf vertices
            H_core = [v for v in C_A if probs[v] > 1/3 + TOL]
            if not H_core:
                continue

            # Check F_gap(H_core) >= 0
            fg = f_gap(H_core, adj, leaves_T, A, probs)
            if fg < -TOL:
                F_gap_negative += 1

            # Check Steiner peeling for all non-singleton subsets
            # (only up to |H_core| <= 8 for performance)
            if len(H_core) > 8:
                continue
            for mask in range(3, 1 << len(H_core)):
                S = [H_core[i] for i in range(len(H_core)) if (mask >> i) & 1]
                if len(S) < 2:
                    continue
                # Find Steiner leaves
                s_leaves = steiner_tree_leaves(S, adj, n0)
                if not s_leaves:
                    continue
                # Check if ANY Steiner leaf has M_gap >= 0
                best_M = None
                for h in s_leaves:
                    M = marginal_gap(h, S, adj, leaves_T, A, probs)
                    if best_M is None or M > best_M:
                        best_M = M
                if best_M is not None and best_M < -TOL:
                    steiner_M_negative += 1
                    if best_M < worst_M:
                        worst_M = best_M
                        worst_M_tree = (n0, line.decode().strip(), S)
                elif best_M is not None and abs(best_M) < TOL:
                    steiner_M_zero += 1

        proc.wait()
        total_trees += n_trees
        total_dleaf += n_dleaf
        print(f"n={n:2d}: trees={n_trees:8d} d_leaf<=1={n_dleaf:7d} "
              f"A_heavy={A_heavy_fails} leaf_id_fail={leaf_identity_fails} "
              f"n/3-mu<0={n_div3_negative} Fgap<0={F_gap_negative} "
              f"steiner_M<0={steiner_M_negative} M=0={steiner_M_zero}",
              flush=True)

    print()
    print(f"TOTAL trees={total_trees:,} d_leaf<=1={total_dleaf:,}")
    print(f"Checks:")
    print(f"  A heavy (P>=1/3): {A_heavy_fails} failures (expect 0)")
    print(f"  Leaf identity P(l)=(1-P(s))/2: {leaf_identity_fails} failures (expect 0)")
    print(f"  n/3-mu < 0: {n_div3_negative} failures (expect 0)")
    print(f"  F_gap(H_core) < 0: {F_gap_negative} failures (expect 0)")
    print(f"  Steiner leaf M_gap < 0: {steiner_M_negative} (key claim)")
    print(f"  Steiner leaf M_gap = 0: {steiner_M_zero} (tight cases)")
    if worst_M_tree:
        n0, g6, S = worst_M_tree
        print(f"  Worst M_gap={worst_M:.6f} at n={n0} g6={g6} S={S}")

if __name__ == "__main__":
    main()
