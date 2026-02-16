#!/usr/bin/env python3
"""Check: does the path P_n have excess → 2/3 as n → ∞?

The path has a well-known closed-form independence polynomial related 
to Fibonacci numbers. The mean IS size for P_n should be approximately 
n * golden_ratio^{-1} ≈ n/phi. 

If mode ≈ floor(n/3) and mean ≈ n/phi, then:
  excess = mode - mean ≈ n(1/3 - 1/phi) = n(1/3 - 0.618) < 0

Wait, that's wrong. Let me just compute directly.
"""

import math
from indpoly import independence_poly

def build_path(n):
    adj = [[] for _ in range(n)]
    for i in range(n-1):
        adj[i].append(i+1)
        adj[i+1].append(i)
    return n, adj

def main():
    print("=" * 80)
    print("  PATH ANALYSIS: Mode - Mean Excess for P_n")
    print("=" * 80)
    print()
    print(f"{'n':>5} {'mode':>5} {'mean':>10} {'excess':>10} {'mode/n':>8} {'mean/n':>8} {'alpha':>5}")
    print("-" * 60)
    
    for n in [10, 20, 30, 40, 50, 60, 75, 100, 125, 150, 175, 200, 250, 300, 400, 500]:
        _, adj = build_path(n)
        poly = independence_poly(n, adj)
        total = sum(poly)
        mean = sum(k * poly[k] for k in range(len(poly))) / total
        mode = max(range(len(poly)), key=lambda k: poly[k])
        excess = mode - mean
        alpha = len(poly) - 1
        
        print(f"{n:5d} {mode:5d} {mean:10.4f} {excess:+10.6f} {mode/n:8.5f} {mean/n:8.5f} {alpha:5d}")
    
    print()
    
    # What's the asymptotic mode/n and mean/n?
    # Independence number of P_n = ceil(n/2)
    # Mean IS size: by the hard-core model at fugacity 1, 
    # the occupation probability on a path is p = 1 - 1/phi where phi = golden ratio.
    # So mean/n → 1 - 1/phi ≈ 0.38197
    
    phi = (1 + math.sqrt(5)) / 2
    p_occ = 1 - 1/phi  # ≈ 0.38197
    print(f"Asymptotic mean/n → 1 - 1/φ ≈ {p_occ:.5f}")
    print(f"Observed mode/n → mode/n column above")
    print()
    print("If mode/n → some limit L, then excess/n → L - (1 - 1/φ)")
    print("From the data, mode ≈ floor(n·p_occ) or ceil(n·p_occ), so excess → O(1)")
    print()
    
    # Detailed: what is mode/n converging to?
    print("Detailed mode_frac = mode/n:")
    for n in range(50, 501, 50):
        _, adj = build_path(n)
        poly = independence_poly(n, adj)
        total = sum(poly)
        mean = sum(k * poly[k] for k in range(len(poly))) / total
        mode = max(range(len(poly)), key=lambda k: poly[k])
        print(f"  P_{n}: mode/n = {mode/n:.6f}, mean/n = {mean/n:.6f}, excess = {mode-mean:+.4f}")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
