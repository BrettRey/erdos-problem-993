#!/usr/bin/env python3
"""Check Core Edge Contraction Hypothesis.

Hypothesis: Contracting an edge uv between two core vertices increases mu.
mu(G/uv) >= mu(G).
"""

import subprocess
import sys
from collections import defaultdict

# Use local modules
sys.path.insert(0, ".")
from indpoly import independence_poly
from graph6 import parse_graph6

def get_mu(n, adj):
    poly = independence_poly(n, adj)
    total = sum(poly)
    # mu = sum(k * a_k) / sum(a_k)
    weighted = sum(k * c for k, c in enumerate(poly))
    return weighted / total

def contract_edge(n, adj, u, v):
    # Contract v into u
    # New graph has n-1 vertices.
    # Map old indices to new.
    # We remove v.
    # Vertices 0..n-1.
    # Limit: keep u as the merged node.
    # Map v -> u.
    # Map w > v -> w-1.
    
    new_n = n - 1
    new_adj = [set() for _ in range(new_n)]
    
    # Mapping
    def map_node(x):
        if x == v: return map_node(u) # Map v to u's new index
        if x > v: return x - 1
        return x
        
    target_u = map_node(u)
    
    for old_i in range(n):
        if old_i == v: continue # Handled via u
        
        new_i = map_node(old_i)
        
        # Add edges
        for old_j in adj[old_i]:
            if old_j == v: new_j = map_node(u)
            else: new_j = map_node(old_j)
            
            if new_i != new_j: # Avoid self-loops
                new_adj[new_i].add(new_j)
                new_adj[new_j].add(new_i)
                
    # Also handle v's original edges (add to u)
    for old_j in adj[v]:
        if old_j == u: continue
        new_j = map_node(old_j)
        if target_u != new_j:
            new_adj[target_u].add(new_j)
            new_adj[new_j].add(target_u)

    return new_n, [list(s) for s in new_adj]

def check_contraction(n, adj):
    leaves = [x for x in range(n) if len(adj[x]) == 1]
    support_map = {} 
    for l in leaves:
        supp = adj[l][0]
        support_map[l] = supp
        
    used = set()
    for l in leaves:
        used.add(l)
        used.add(support_map[l])
        
    core = [x for x in range(n) if x not in used]
    
    if not core: return None, 0, 0
    
    # Identify core edges
    core_edges = []
    seen_edges = set()
    for u in core:
        for w in adj[u]:
            if w in core and w > u:
                edge = tuple(sorted((u, w)))
                if edge not in seen_edges:
                    seen_edges.add(edge)
                    core_edges.append(edge)
                    
    if not core_edges: return None, 0, 0
    
    mu_orig = get_mu(n, adj)
    
    violations = []
    
    for u, v in core_edges:
        # Contract uv
        new_n, new_adj = contract_edge(n, adj, u, v)
        mu_new = get_mu(new_n, new_adj)
        
        # User Hypothesis: mu_new >= mu_orig
        # Violation if mu_new < mu_orig
        
        if mu_new < mu_orig:
            violations.append((u, v, mu_new, mu_orig))
            
    return violations, mu_orig, len(core_edges)

def main():
    print("ANALYSIS: Core Edge Contraction (Hypothesis: mu increases on contraction)", flush=True)
    
    for n in range(3, 17): 
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True) 
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        total_viol = 0
        total_tested = 0
        
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
            
            viols, mu, n_edges = check_contraction(tn, adj)
            
            if viols is not None:
                total_tested += n_edges
                total_viol += len(viols)
                if len(viols) > 0:
                     # Print first few violations
                     u, v, m_new, m_old = viols[0]
                     # print(f"n={n} Violation: Contract {u}-{v} -> mu decreases ({m_old:.3f} -> {m_new:.3f})")
                     pass
                     
        if total_viol > 0:
            print(f"n={n:2d}: {total_viol}/{total_tested} contractions decreased mu (VIOLATION).", flush=True)
        else:
            print(f"n={n:2d}: All {total_tested} contractions increased mu (CONFIRMED).", flush=True)

if __name__ == "__main__":
    main()
