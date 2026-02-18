#!/usr/bin/env python3
"""Verify Matching Cover Hypothesis for d_leaf <= 1 trees.

Hypothesis: There exists a matching M that covers all 'Heavy' vertices (P(v) > 1/3).
If true, then Unmatched vertices are all 'Light' (P <= 1/3).
Combined with Edge Bound P(u)+P(v) <= 2/3, this proves mu < n/3.
"""

import subprocess
import sys

# Use local modules
sys.path.insert(0, ".")
from indpoly import independence_poly
from graph6 import parse_graph6

def get_occupation_probs(n, adj):
    poly_T = independence_poly(n, adj)
    total_IS = sum(poly_T)
    probs = []
    for v in range(n):
        neighbors = set(adj[v])
        neighbors.add(v) 
        remaining = [u for u in range(n) if u not in neighbors]
        if not remaining:
            count_inc = 1 
        else:
            mapping = {old: new for new, old in enumerate(remaining)}
            new_adj = [[] for _ in range(len(remaining))]
            for old_u in remaining:
                for old_w in adj[old_u]:
                    if old_w in mapping: 
                         new_adj[mapping[old_u]].append(mapping[old_w])
            poly_sub = independence_poly(len(remaining), new_adj)
            count_inc = sum(poly_sub)
        probs.append(count_inc / total_IS)
    return probs

def check_heavy_cover(n, adj, probs):
    # Identify Heavy vertices
    heavy = [v for v in range(n) if probs[v] > 1/3 + 1e-9]
    if not heavy:
        return True, "No heavy vertices"

    # We need a matching that covers all 'heavy' vertices.
    # We can model this as a Maximum Matching problem.
    # To prioritize covering heavy vertices, we can find a matching in G
    # but we need to ensure specific nodes are covered.
    
    # Simple valid check:
    # Does G contain a matching covering S = heavy?
    # Hall's condition? No, that's bipartite.
    # Tutte's theorem?
    # Constructive approach:
    # Greedily match heavy vertices?
    # Better: Max Weight Matching where heavy edges have high weight?
    # Edges incident to heavy get weight 100. Others weight 1.
    # But an edge between two heavy nodes gets weight ??
    # We want to maximize number of covered heavy nodes.
    # Weight w(u,v):
    # - If u,v both heavy: 1000
    # - If one heavy: 100
    # - If neither: 1
    # Maximize weight.
    
    # Implement Max Weight Matching for Tree (Dynamic Programming)
    # State DP[u][0] = Max weight of matching in subtree u, u is unmatched
    # State DP[u][1] = Max weight of matching in subtree u, u is matched
    
    weights = {}
    for u in range(n):
        for v in adj[u]:
            if u > v: continue
            w = 0
            if u in heavy: w += 100
            if v in heavy: w += 100
            weights[(u,v)] = w
            weights[(v,u)] = w
            
    # Tree DP
    # Root at 0
    parent = [-1] * n
    children = [[] for _ in range(n)]
    order = []
    queue = [0]
    visited = {0}
    while queue:
        u = queue.pop(0)
        order.append(u)
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                parent[v] = u
                children[u].append(v)
                queue.append(v)
                
    order.reverse() # Process leaves first
    
    dp_unmatched = [0] * n
    dp_matched = [0] * n
    
    for u in order:
        # DP[u][0]: u is unmatched. Sum of max(DP[v][0], DP[v][1]) for children
        sum_children = 0
        for v in children[u]:
            sum_children += max(dp_unmatched[v], dp_matched[v])
        dp_unmatched[u] = sum_children
        
        # DP[u][1]: u is matched to ONE child v.
        # Weight = w(u,v) + DP[v][0] (v matched to u) + Sum_{others} max(DP[k][0], DP[k][1])
        # = w(u,v) + dp_unmatched[v] + (sum_children - max(dp_unmatched[v], dp_matched[v]))
        
        max_matched = -1 # Cannot be matched if no children... wait.
        # Actually DP[u][1] means u is matched DOWN.
        # But u could be matched UP.
        # Standard Independent Set DP logic but for Matching.
        # Let's clarify states.
        # DP[u][0]: Max weight in subtree rooted at u, u is FREE (available for parent).
        # DP[u][1]: Max weight in subtree rooted at u, u is COVERED (by child).
        
        # Valid matching in subtree:
        # If u free: all children must be covered or free.
        # DP[u][0] = Sum_{v in children} max(DP[v][0], DP[v][1])
        # Same as before.
        
        # If u covered by child v:
        # Edge (u,v) used. v must be free in its subtree (to match u).
        # Other children k can be free or covered.
        
        best_match_val = -float('inf')
        
        for v in children[u]:
            # Try matching u-v
            # Gain: w(u,v)
            # v contributes dp_unmatched[v] (since v matched to u, it was 'free' from below)
            # Others contribute max(...)
            
            val = weights[(u,v)] + dp_unmatched[v] + (sum_children - max(dp_unmatched[v], dp_matched[v]))
            if val > best_match_val:
                best_match_val = val
                
        dp_matched[u] = best_match_val
        
        # Special case: Leaf
        if not children[u]:
            dp_unmatched[u] = 0
            dp_matched[u] = -float('inf')
            
    # Final result at root
    total_weight = max(dp_unmatched[0], dp_matched[0])
    
    # Check if all heavy nodes covered
    # Each heavy node contributes 100 to any edge covering it.
    # Total possible weight from heavy nodes = 100 * |Heavy|.
    # Wait, if u,v both heavy, edge weight is 200. Covered weight is 200 (100 for u, 100 for v).
    # If u matched to non-heavy v, weight 100.
    # So if all heavy nodes are covered, the 'heavy portion' of the score should be 100 * |Heavy|.
    # Is it strictly additive?
    # Yes. W = Sum_(edges) (100*IsHeavy(u) + 100*IsHeavy(v)).
    #        = 100 * Sum_(v in Heavy) IsCovered(v).
    # So Total Weight / 100 (integer div) should equal |Heavy|.
    
    covered_heavy_count = int(total_weight) // 100
    
    if covered_heavy_count == len(heavy):
        return True, f"Covered all {len(heavy)}"
    else:
        return False, f"Covered {covered_heavy_count}/{len(heavy)}"

def main():
    print("ANALYSIS: Matching Cover for Heavy Vertices (P > 1/3)", flush=True)
    
    for n in range(3, 19): 
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True) 
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_trees = 0
        n_heavy = 0
        n_fail = 0
        
        for line in lines:
            tn, adj = parse_graph6(line)
            
            # Filter d_leaf <= 1
            leaves = [v for v in range(tn) if len(adj[v]) == 1]
            is_valid = True
            for v in range(tn):
                if v in leaves: continue
                leaf_neighbors = sum(1 for u in adj[v] if u in leaves)
                if leaf_neighbors > 1:
                    is_valid = False; break
            if not is_valid: continue
            
            n_trees += 1
            probs = get_occupation_probs(tn, adj)
            
            # Optimization: only run matching if there are heavy nodes
            if any(p > 1/3 + 1e-9 for p in probs):
                n_heavy += 1
                ok, msg = check_heavy_cover(tn, adj, probs)
                if not ok:
                    print(f"FAIL n={tn}: {msg}", flush=True)
                    print(f"Adj: {adj}", flush=True)
                    print(f"Probs: {['%.3f'%p for p in probs]}", flush=True)
                    n_fail += 1
                    
        print(f"n={n:2d}: {n_trees} valid trees. {n_heavy} with heavy nodes. {n_fail} failures.", flush=True)

if __name__ == "__main__":
    main()
