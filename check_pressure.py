#!/usr/bin/env python3
"""Check Tree Recursion Pressure Bound.

User Hypothesis: p = (1/n) * log(I(T)) satisfies p < log(phi)/2 approx 0.24.
Note: Path graph P_n has I(P_n) approx phi^n.
So p(P_n) approx log(phi) approx 0.48.
This suggests the user's bound might be too tight or interpreted differently.
"""

import subprocess
import sys
import math

# Use local modules
sys.path.insert(0, ".")
from indpoly import independence_poly
from graph6 import parse_graph6

PHI = (1 + math.sqrt(5)) / 2
LOG_PHI = math.log(PHI)
TARGET = LOG_PHI / 2

def check_pressure(n, adj):
    poly = independence_poly(n, adj)
    total_is = sum(poly)
    entropy = math.log(total_is)
    pressure = entropy / n
    return pressure

def main():
    print(f"ANALYSIS: Pressure Check (p = log(I(T))/n)", flush=True)
    print(f"Target Bound: < log(phi)/2 = {TARGET:.5f}", flush=True)
    print(f"Path Reference: log(phi) = {LOG_PHI:.5f}", flush=True)
    
    max_p = 0.0
    
    for n in range(3, 17): 
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True) 
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_viol = 0
        
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
            
            p = check_pressure(tn, adj)
            
            if p > max_p:
                max_p = p
                
            if p > TARGET:
                n_viol += 1
                
        print(f"n={n:2d}: Max Pressure = {max_p:.5f}. Violations over {TARGET:.3f}: {n_viol}/{len(lines)}", flush=True)

if __name__ == "__main__":
    main()
