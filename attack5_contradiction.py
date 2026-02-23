#!/usr/bin/env python3
"""
Attack 5: Verify the contradiction hypothesis for mode(P) >= m-1.

We assume mode(P) <= m-2, where m = mode(I(T)).
This implies P is decreasing at m-1 and m.
p_{m-2} > p_{m-1} > p_m.

The mode condition a_m >= a_{m-1} requires:
(p_m - p_{m-2}) + (p_{m-1} - p_{m-2}) + (q_m - q_{m-2}) >= 0.

Since P terms are negative, we need q_m - q_{m-2} to be large enough.
q_k is coeff of Q = x * G. So q_m = g_{m-1}, q_{m-2} = g_{m-3}.
We need g_{m-1} - g_{m-3} to be large.

This script scans d_leaf<=1 trees, finds a degree-2 support leaf,
computes P, Q, m=mode(I(T)), and checks:
1. Is mode(P) <= m-2 ever true?
2. What are the values of Drop_P = p_{m-2} - p_{m-1} and Gap_Q = q_m - q_{m-2}?
3. If mode(P) >= m-1 (which we expect), how close is it to failing?
   i.e. is p_{m-1} close to p_{m-2}?
"""

import sys
import subprocess
from indpoly import independence_poly, is_log_concave
from graph6 import parse_graph6

def get_mode_indices(poly):
    """Return all indices where poly achieves its maximum."""
    if not poly:
        return []
    max_val = max(poly)
    return [i for i, x in enumerate(poly) if x == max_val]

def remove_vertex(n, adj, v):
    """Return adj of graph with v removed, reindexed."""
    mapping = []
    new_adj = []
    curr = 0
    old_to_new = {}
    for i in range(n):
        if i != v:
            mapping.append(i)
            old_to_new[i] = curr
            curr += 1
            new_adj.append([])
    
    for i in mapping:
        u = old_to_new[i]
        for neighbor in adj[i]:
            if neighbor != v:
                v_new = old_to_new[neighbor]
                new_adj[u].append(v_new)
    return len(new_adj), new_adj

def remove_vertices(n, adj, vertices):
    """Return adj of graph with vertices removed."""
    current_n = n
    current_adj = adj
    # Remove one by one (inefficient but safe for small graphs)
    keep = [i for i in range(n) if i not in vertices]
    old_to_new = {old: new for new, old in enumerate(keep)}
    new_n = len(keep)
    new_adj = [[] for _ in range(new_n)]
    for old in keep:
        new = old_to_new[old]
        for neighbor in current_adj[old]:
            if neighbor in keep:
                new_adj[new].append(old_to_new[neighbor])
    return new_n, new_adj

def get_P_Q(n, adj, leaf, support):
    """
    Decompose based on leaf l and support s (deg(s)=2).
    Let u be the other neighbor of s.
    B = T - {l, s}.
    P = I(B - u).
    Q = x * I(B - N_B[u]).
    
    Returns P, Q (as coefficient lists).
    """
    # Identify u
    neighbors_s = adj[support]
    u = [x for x in neighbors_s if x != leaf][0]
    
    # B = T - {l, s}
    n_B, adj_B = remove_vertices(n, adj, {leaf, support})
    
    # Map u to B's indexing
    # We need to know which index u became in B.
    # remove_vertices preserves order of kept vertices.
    # count how many removed vertices are < u
    removed_count = sum(1 for v in {leaf, support} if v < u)
    u_in_B = u - removed_count
    
    # P = I(B - u)
    n_P, adj_P = remove_vertex(n_B, adj_B, u_in_B)
    P = independence_poly(n_P, adj_P)
    
    # Q = x * I(B - N_B[u])
    # neighbors of u in B
    neighbors_u_B = adj_B[u_in_B]
    n_Q_inner, adj_Q_inner = remove_vertices(n_B, adj_B, set(neighbors_u_B) | {u_in_B})
    Q_inner = independence_poly(n_Q_inner, adj_Q_inner)
    Q = [0] + Q_inner
    
    return P, Q

def solve():
    max_n = 20
    print(f"Checking trees with d_leaf <= 1 up to n={max_n}...")
    
    geng_cmd = ["/opt/homebrew/bin/geng", "-c", "-q"]
    
    counters = {
        "total": 0,
        "d_leaf_le_1": 0,
        "deg2_support": 0,
        "mode_P_le_m_minus_2": 0,
        "mode_P_eq_m_minus_1": 0,
        "mode_P_ge_m": 0,
        "mode_G_gt_mode_P": 0,
        "mode_G_gt_mode_P_plus_1": 0,
    }
    
    min_gap_Q_needed = float('inf')
    max_gap_Q_actual = float('-inf')
    max_gap_case = None

    for n in range(4, max_n + 1):
        cmd = geng_cmd + [str(n), f"{n-1}:{n-1}"] # Trees
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) # text=True for line iteration
        
        # Note: graph6 might expect bytes, but parse_graph6 handles bytes.
        # If we iterate text lines, we must encode back to bytes for parse_graph6.
        
        for line_str in process.stdout:
            counters["total"] += 1
            # Encode to bytes for graph6 parser
            line_bytes = line_str.encode('ascii')
            n_nodes, adj = parse_graph6(line_bytes)
            
            # Check d_leaf <= 1
            # d_leaf(v) is number of leaf neighbors of v
            is_valid = True
            leaves = [v for v in range(n_nodes) if len(adj[v]) == 1]
            leaf_set = set(leaves)
            
            for v in range(n_nodes):
                if v in leaf_set:
                    continue
                d_leaf = sum(1 for u in adj[v] if u in leaf_set)
                if d_leaf > 1:
                    is_valid = False
                    break
            
            if not is_valid:
                continue
                
            counters["d_leaf_le_1"] += 1
            
            # Find a leaf with support of degree 2
            target_leaf = -1
            target_support = -1
            
            for l in leaves:
                s = adj[l][0]
                if len(adj[s]) == 2:
                    target_leaf = l
                    target_support = s
                    break
            
            if target_leaf != -1:
                counters["deg2_support"] += 1
                
                # Compute I(T), P, Q
                I_T = independence_poly(n_nodes, adj)
                P, Q = get_P_Q(n_nodes, adj, target_leaf, target_support)
                
                modes_T = get_mode_indices(I_T)
                m = modes_T[0] 
                
                modes_P = get_mode_indices(P)
                mode_P = modes_P[-1] 
                
                # Q = x * G. So G is Q shifted by -1 (indices).
                G = Q[1:]
                modes_G = get_mode_indices(G)
                mode_G = modes_G[-1] if modes_G else 0

                if mode_G > mode_P:
                    counters["mode_G_gt_mode_P"] += 1
                    if mode_G > mode_P + 1:
                        counters["mode_G_gt_mode_P_plus_1"] += 1

                # Check condition
                if mode_P <= m - 2:
                    counters["mode_P_le_m_minus_2"] += 1
                    print(f"COUNTEREXAMPLE FOUND! n={n_nodes} mode(T)={m}, mode(P)={mode_P}")
                elif mode_P == m - 1:
                    counters["mode_P_eq_m_minus_1"] += 1
                else:
                    counters["mode_P_ge_m"] += 1
                    
                # Analyze Gap_Q
                val_p_m1 = P[m-1] if m-1 < len(P) else 0
                val_q_m = Q[m] if m < len(Q) else 0
                val_q_m2 = Q[m-2] if m-2 < len(Q) else 0
                
                term_Q = val_q_m - val_q_m2
                
                if term_Q > max_gap_Q_actual:
                    max_gap_Q_actual = term_Q
                    max_gap_case = {
                        "n": n_nodes,
                        "m": m,
                        "P": P,
                        "Q": Q,
                        "Gap_Q": term_Q,
                        "p_m1": val_p_m1,
                        "ratio": term_Q / val_p_m1 if val_p_m1 > 0 else 0
                    }

    print(f"Finished. Counters: {counters}")
    print(f"Max observed Gap_Q (q_m - q_{{m-2}}): {max_gap_Q_actual}")
    if max_gap_case:
        print("Details of max Gap_Q case:")
        print(f"  n={max_gap_case['n']}, m={max_gap_case['m']}")
        print(f"  Gap_Q={max_gap_case['Gap_Q']}, p_{{m-1}}={max_gap_case['p_m1']}")
        print(f"  Ratio Gap_Q/p_{{m-1}} = {max_gap_case['ratio']:.4f}")

if __name__ == "__main__":
    solve()
