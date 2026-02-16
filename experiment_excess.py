#!/usr/bin/env python3
"""Focused investigation: WHY is mode - mean < 1 for all trees?

Key observation: 44% of trees have mode > mean, but the max excess
is only 0.619. If mode - mean < 1 always, then mode ≤ ceil(mean).

This script:
1. Tracks the trees with largest mode - mean at each n
2. Checks if the excess approaches 1 as n grows
3. Looks for structural patterns in the worst-case trees
"""

import subprocess
import math
from indpoly import independence_poly
from graph6 import parse_graph6

def leaves_and_degree(adj, n):
    """Return (num leaves, max degree) of a tree."""
    degs = [len(adj[v]) for v in range(n)]
    num_leaves = sum(1 for d in degs if d == 1)
    return num_leaves, max(degs)

def main():
    print("=" * 70)
    print("  WHY IS mode - mean < 1 FOR ALL TREES?")
    print("=" * 70)
    print()
    
    print(f"{'n':>3} {'max_excess':>10} {'mode':>5} {'mean':>8} {'alpha':>5} "
          f"{'leaves':>6} {'max_deg':>7} {'skew':>8} {'poly'}")
    print("-" * 100)
    
    worst_by_n = []
    
    for n in range(3, 23):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        max_excess = -float('inf')
        best_info = None
        
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            total = sum(poly)
            mean = sum(k * poly[k] for k in range(len(poly))) / total
            mode = max(range(len(poly)), key=lambda k: poly[k])
            excess = mode - mean
            
            if excess > max_excess:
                max_excess = excess
                var = sum(k**2 * poly[k] for k in range(len(poly))) / total - mean**2
                sd = math.sqrt(var)
                skew = sum((k - mean)**3 * poly[k] for k in range(len(poly))) / (total * sd**3) if sd > 0 else 0
                nl, md = leaves_and_degree(adj, tn)
                best_info = {
                    'excess': excess, 'mode': mode, 'mean': mean,
                    'alpha': len(poly)-1, 'leaves': nl, 'max_deg': md,
                    'skew': skew, 'poly': poly, 'g6': line,
                }
        
        worst_by_n.append((n, max_excess, best_info))
        
        info = best_info
        poly_str = str(info['poly']) if len(info['poly']) <= 12 else str(info['poly'][:12]) + '...'
        print(f"{n:3d} {info['excess']:10.6f} {info['mode']:5d} {info['mean']:8.4f} "
              f"{info['alpha']:5d} {info['leaves']:6d} {info['max_deg']:7d} "
              f"{info['skew']:8.4f}  {poly_str}")
    
    print()
    
    # Check the trend
    print("TREND: Does max_excess approach 1 as n grows?")
    print()
    for n, excess, _ in worst_by_n:
        bar = '█' * int(excess * 50)
        print(f"  n={n:2d}: {excess:.4f} {bar}│1.0")
    
    print()
    
    # Group by n parity
    odd_excesses = [e for n, e, _ in worst_by_n if n % 2 == 1]
    even_excesses = [e for n, e, _ in worst_by_n if n % 2 == 0]
    print(f"  Odd n:  max excess range [{min(odd_excesses):.4f}, {max(odd_excesses):.4f}]")
    print(f"  Even n: max excess range [{min(even_excesses):.4f}, {max(even_excesses):.4f}]")
    print()
    
    # Try to identify what the worst trees look like
    print("STRUCTURE: What do the worst-case trees look like?")
    print()
    for n, excess, info in worst_by_n:
        if info:
            print(f"  n={n}: leaves={info['leaves']}, max_deg={info['max_deg']}, "
                  f"alpha={info['alpha']}, mode/alpha={info['mode']/info['alpha']:.3f}")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
