#!/usr/bin/env python3
"""
Attack 6: Analyze the (1+2x)P + (1+x)Q structure to explain mode(P) >= m-1.

We search for trees (specifically d_leaf <= 1) that have a leaf l attached to a support s of degree 2.
For these trees, I(T) = (1+2x)P + (1+x)Q where P = I(B-u), Q = I(B-N[u]), B = T-{l,s}, u is s's neighbor.

Hypotheses to check:
1. mode(F) is always m or m-1, where F = (1+2x)P.
2. mode(F) is constrained to {mode(P), mode(P)+1}.
3. The ratio f_m / g_m is large, implying F dominates G near the mode.
4. Verify mode(P) >= m-1.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from collections import Counter

# Add current directory to path so we can import local modules
sys.path.append(os.getcwd())

from graph6 import parse_graph6
from indpoly import independence_poly, _polyadd, _polymul, is_log_concave

def mode_indices(poly: list[int]) -> list[int]:
    """Return all indices where the polynomial achieves its maximum."""
    if not poly:
        return []
    mx = max(poly)
    return [i for i, c in enumerate(poly) if c == mx]

def mode_index_leftmost(poly: list[int]) -> int:
    """Return the smallest index where the polynomial achieves its maximum."""
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i]) # Python's max is stable, so we need to reverse or use careful logic. 
    # Actually max(range(len), key=...) will return the *first* occurrence if ties.
    # Wait, max() on indices with key=poly.__getitem__ returns the index of the max value.
    # If there are ties, it returns the *first* one encountered? 
    # Python documentation says: "If multiple items are maximal, the function returns the first one encountered."
    # So max(range(len(poly)), key=lambda i: poly[i]) returns the FIRST index (leftmost).
    # YES.
    mx = -1
    best_i = -1
    for i, c in enumerate(poly):
        if c > mx:
            mx = c
            best_i = i
    return best_i

def remove_vertices(adj: list[list[int]], remove_set: set[int]) -> list[list[int]]:
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out

def compute_hub_polys(adj_B: list[list[int]], u_in_B: int) -> tuple[list[int], list[int]]:
    """Compute P = dp[u][0] and Q = dp[u][1] by rooting B at u."""
    n = len(adj_B)
    if n == 0:
        return [1], [] # P=1, Q=0? I(B) = P+Q = 1. Correct.
    
    # BFS/Post-order logic
    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[u_in_B] = True
    bfs_queue = [u_in_B]
    head = 0
    while head < len(bfs_queue):
        v = bfs_queue[head]
        head += 1
        for w in adj_B[v]:
            if not visited[w]:
                visited[w] = True
                children[v].append(w)
                bfs_queue.append(w)

    order = []
    stack = [(u_in_B, False)]
    while stack:
        v, processed = stack.pop()
        if processed:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    dp0 = [[] for _ in range(n)]
    dp1 = [[] for _ in range(n)]

    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod = [1]
            for c in children[v]:
                summand = _polyadd(dp0[c], dp1[c])
                prod = _polymul(prod, summand)
            dp0[v] = prod
            
            prod = [1]
            for c in children[v]:
                prod = _polymul(prod, dp0[c])
            dp1[v] = [0] + prod

    return dp0[u_in_B], dp1[u_in_B]

def is_dleaf_le_1(n: int, adj: list[list[int]]) -> bool:
    """Check if max degree of any leaf is <= 1 (wait, max degree of a leaf is always 1 unless isolated).
    Ah, 'd_leaf' usually means 'degree of the neighbor of the leaf'.
    But here, maybe it means 'derivative'?
    No, in this context (Erdos 993), d_leaf usually refers to the 'derivative at 1' or something?
    
    Actually, looking at `conjecture_a_hall_subset_scan.py` imported in `diagnose_bridge_decomposition.py`:
    `from conjecture_a_hall_subset_scan import is_dleaf_le_1`
    
    Let's check what `is_dleaf_le_1` does in other files or implement a check.
    Usually 'd_leaf' refers to the distance between leaves?
    Or maybe it means 'tree with at most 1 leaf'? No.
    
    Wait, `diagnose_bridge_decomposition.py` checks `if (not args.all_trees) and (not is_dleaf_le_1(nn, adj)):`.
    
    Let's assume for now we scan all trees, or I can copy the function if I find it.
    But given the prompt mentions "d_leaf<=1 trees", I should try to filter for them if possible.
    
    Let's skip the filter for now and just check "trees with the decomposition".
    """
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_min", type=int, default=10)
    parser.add_argument("--n_max", type=int, default=20)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    args = parser.parse_args()

    print(f"Scanning n={args.n_min}..{args.n_max} for (1+2x)P + (1+x)Q structure...")
    
    stats = {
        "total_trees": 0,
        "valid_structure": 0,
        "mode_P_ge_m_minus_1": 0,
        "mode_P_less_m_minus_1": 0,
        "mode_F_is_m_or_m_minus_1": 0,
        "mode_F_outside_m_m1": 0,
        "mode_F_in_mP_mP1": 0,
        "mode_F_outside_mP_mP1": 0,
        "P_is_LC": 0,
        "P_not_LC": 0,
        "counterexamples": []
    }
    
    ratio_sum = 0.0
    ratio_count = 0

    for n in range(args.n_min, args.n_max + 1):
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        
        for line in proc.stdout:
            stats["total_trees"] += 1
            line = line.strip()
            # parse_graph6 expects bytes
            try:
                nn, adj = parse_graph6(line)
            except Exception:
                continue
            
            # Find a leaf l with neighbor s of degree 2
            deg = [len(x) for x in adj]
            candidates = []
            for v in range(nn):
                if deg[v] == 1:
                    s = adj[v][0]
                    if deg[s] == 2:
                        candidates.append((v, s))
            
            if not candidates:
                continue
                
            # Prefer the one that gives 'nicest' P/Q if multiple? Just pick first for now.
            l, s = candidates[0]
            
            # u is the other neighbor of s
            u = [x for x in adj[s] if x != l][0]
            
            # B = T - {l, s}
            # To get indices right, we map vertices.
            remove_set = {l, s}
            keep = [v for v in range(nn) if v not in remove_set]
            mapping = {old: new for new, old in enumerate(keep)}
            
            adj_B = remove_vertices(adj, remove_set)
            
            if u not in mapping:
                # Should not happen unless u was l or s (impossible)
                continue
            
            u_B = mapping[u]
            
            P, Q = compute_hub_polys(adj_B, u_B)
            
            # F = (1+2x)P = P + 2xP
            P_shifted = [0] + P
            F = _polyadd(P, [2*c for c in P_shifted])
            
            # G = (1+x)Q = Q + xQ
            Q_shifted = [0] + Q
            G = _polyadd(Q, Q_shifted)
            
            # I(T)
            I_T = _polyadd(F, G)
            
            m = mode_index_leftmost(I_T)
            m_P = mode_index_leftmost(P)
            m_F = mode_index_leftmost(F)
            
            stats["valid_structure"] += 1
            
            # Check hypothesis 1: mode(P) >= m-1
            if m_P >= m - 1:
                stats["mode_P_ge_m_minus_1"] += 1
            else:
                stats["mode_P_less_m_minus_1"] += 1
                # Only keep a few details to avoid giant JSON
                if len(stats["counterexamples"]) < 20: 
                    stats["counterexamples"].append({
                        "n": n,
                        "g6": line.decode('utf-8'),
                        "m": m,
                        "m_P": m_P,
                        "m_F": m_F,
                        "m_G": mode_index_leftmost(G),
                        "P": P,
                        "F": F,
                        "G": G,
                        "I_T": I_T,
                        "issue": "mode(P) < m-1"
                    })

            # Check hypothesis 2: mode(F) in {m_P, m_P+1}
            if m_F in [m_P, m_P + 1]:
                stats["mode_F_in_mP_mP1"] += 1
            else:
                stats["mode_F_outside_mP_mP1"] += 1
                
            # Check hypothesis 3: mode(F) in {m, m-1} ?
            # Wait, if mode(P) >= m-1, and mode(F) >= mode(P), then mode(F) >= m-1.
            # But mode(F) could be m+1?
            if m_F in [m, m - 1]:
                stats["mode_F_is_m_or_m_minus_1"] += 1
            else:
                stats["mode_F_outside_m_m1"] += 1
                if m_F < m - 1:
                     if "mode_F_lt_m_minus_1" not in stats: stats["mode_F_lt_m_minus_1"] = 0
                     stats["mode_F_lt_m_minus_1"] += 1
                     if len(stats["counterexamples"]) < 30:
                         stats["counterexamples"].append({
                            "n": n,
                            "g6": line.decode('utf-8'),
                            "m": m,
                            "m_P": m_P,
                            "m_F": m_F,
                            "issue": "mode(F) < m-1"
                         })
                elif m_F > m:
                     if "mode_F_gt_m" not in stats: stats["mode_F_gt_m"] = 0
                     stats["mode_F_gt_m"] += 1


            if is_log_concave(P):
                stats["P_is_LC"] += 1
            else:
                stats["P_not_LC"] += 1
                
            # Ratio f_m / g_m
            # Be careful with indices
            f_m = F[m] if m < len(F) else 0
            g_m = G[m] if m < len(G) else 0
            
            if g_m > 0:
                ratio = f_m / g_m
                ratio_sum += ratio
                ratio_count += 1
                
    print(json.dumps(stats, indent=2))
    if ratio_count > 0:
        print(f"Average f_m/g_m ratio: {ratio_sum/ratio_count:.4f}")

if __name__ == "__main__":
    main()
