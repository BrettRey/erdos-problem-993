#!/usr/bin/env python3
"""Analyze load distribution on path segments in d_leaf <= 1 trees.

Tests specific structures:
1. Dumbbell: Leaf-Support --(Path)-- Support-Leaf
2. Star: Central Node --(Path)-- Support-Leaf (x3 arms)
3. H-Graph: (Leaf-Support-Path) -- Center -- (Path-Support-Leaf) (x4 arms)

Calculates GAP = n/3 - mu.
Checks if GAP is positive (hypothesis holds).
"""

import sys
sys.path.insert(0, ".")
from indpoly import independence_poly

def analyze_structure(name, n, adj):
    poly = independence_poly(n, adj)
    mu = 0
    total_w = sum(poly)
    # Compute mu = I'(1)/I(1)
    # I(x) = sum a_k x^k. I'(1) = sum k a_k.
    avg_size = sum(k * c for k, c in enumerate(poly)) / total_w
    
    gap = n/3 - avg_size
    return gap, avg_size

def make_dumbbell(k):
    """
    Leaf(0) - Support(1) - Path(k vertices) - Support(k+2) - Leaf(k+3)
    Internal path vertices: 2..k+1
    Total n = k + 4
    """
    n = k + 4
    adj = [[] for _ in range(n)]
    
    # Left end
    adj[0].append(1); adj[1].append(0)
    
    # Path from 1 to 2..k+1 to k+2
    current = 1
    for i in range(k + 1): # k+1 edges for k internal nodes
        # Edge current <-> current+1
        # Path nodes are 2, 3, ..., k+1
        # Support is 1. Next is 2.
        # Last internal is k+1. Next is k+2 (Support).
        # Vertices: 1, 2, ..., k+1, k+2.
        u, v = current, current + 1
        adj[u].append(v)
        adj[v].append(u)
        current += 1
        
    # Right end (Leaf) attached to k+2
    # Leaf is k+3.
    adj[k+2].append(k+3); adj[k+3].append(k+2)
    
    return n, adj

def make_star_3arm(k):
    """
    Center (0).
    3 Arms. Each arm has k internal path nodes, then Support, then Leaf.
    Arm structure: Center - P1 - ... - Pk - Support - Leaf.
    Length from Center to Support = k+1 edges.
    Nodes per arm: k (path) + 1 (supp) + 1 (leaf) = k+2.
    Total n = 1 + 3*(k+2).
    """
    n = 1 + 3*(k+2)
    adj = [[] for _ in range(n)]
    
    center = 0
    current_idx = 1
    
    for arm in range(3):
        # Center connects to first node of arm
        # If k=0, Connect Center - Support - Leaf
        # If k>0, Center - P1 ...
        
        last = center
        
        # Path nodes
        for _ in range(k):
            p = current_idx
            adj[last].append(p); adj[p].append(last)
            last = p
            current_idx += 1
            
        # Support
        supp = current_idx
        adj[last].append(supp); adj[supp].append(last)
        current_idx += 1
        
        # Leaf
        leaf = current_idx
        adj[supp].append(leaf); adj[leaf].append(supp)
        current_idx += 1
        
    return n, adj

def main():
    print("ANALYSIS: Path Segment Load (Gap = n/3 - mu)", flush=True)
    
    print("\n--- Dumbbell (Leaf-Supp ... Supp-Leaf) ---", flush=True)
    print("K = Internal Path Length. N = K+4.")
    for k in range(0, 21):
        n, adj = make_dumbbell(k)
        gap, mu = analyze_structure(f"Dumbbell-{k}", n, adj)
        print(f"K={k:2d}, N={n:2d}: Gap={gap:.5f} (Mu={mu:.2f})", flush=True)
        if gap < 0:
            print("  VIOLATION!", flush=True)

    print("\n--- Star-3Arm (Center ... Supp-Leaf) ---", flush=True)
    print("K = Path Nodes per Arm. N = 1 + 3(K+2).")
    for k in range(0, 11):
        n, adj = make_star_3arm(k)
        gap, mu = analyze_structure(f"Star3-{k}", n, adj)
        print(f"K={k:2d}, N={n:2d}: Gap={gap:.5f} (Mu={mu:.2f})", flush=True)
        if gap < 0:
            print("  VIOLATION!", flush=True)

if __name__ == "__main__":
    main()
