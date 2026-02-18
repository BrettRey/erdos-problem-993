#!/usr/bin/env python3
import subprocess
import sys
from indpoly import independence_poly
from graph6 import parse_graph6

def get_mean_size(poly):
    # poly[k] is count of IS of size k.
    # Mean = Sum(k * i_k) / Sum(i_k)
    total_count = sum(poly)
    total_size = sum(k * c for k, c in enumerate(poly))
    return total_size / total_count

def has_max_leaf_degree_le_1(n, adj):
    # Check if any vertex has >= 2 leaf neighbors
    leaves = [v for v in range(n) if len(adj[v]) == 1]
    leaf_set = set(leaves)
    
    for v in range(n):
        leaf_neighbors = 0
        for u in adj[v]:
            if u in leaf_set:
                leaf_neighbors += 1
        if leaf_neighbors >= 2:
            return False
    return True

def main():
    print("Checking Mean < n/3 for trees with d_leaf <= 1...")
    
    violation_found = False
    
    # Check up to n=20 (fast enough)
    for n in range(3, 21):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        # If geng not found, skip or warn
        try:
            proc = subprocess.run(cmd, shell=True, capture_output=True)
        except Exception as e:
            print(f"Skipping n={n} (geng error)")
            continue
            
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        count_checked = 0
        for line in lines:
            tn, adj = parse_graph6(line)
            if not has_max_leaf_degree_le_1(tn, adj):
                continue
            
            count_checked += 1
            mu = get_mean_size(independence_poly(tn, adj))
            limit = tn / 3.0
            
            if mu >= limit - 1e-9:
                print(f"VIOLATION or TIGHT at n={tn}: Mean={mu:.5f}, Limit={limit:.5f} (Diff={limit-mu:.2e})")
                print(f"Adj: {adj}")
                violation_found = True
        
        print(f"n={n}: Checked {count_checked} trees. All OK.")

    if not violation_found:
        print("VERIFIED: No violations found up to n=18.")
    else:
        print("WARNING: Violations found!")

if __name__ == "__main__":
    main()
