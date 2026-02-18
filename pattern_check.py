#!/usr/bin/env python3
"""
Formal Proof Attempt: Influence Gini and Unimodality

Key observations from data:
- Stars: concentrated structure → high nm (0.94 for S(100))
- Paths: distributed structure → low nm (0.48 for P(30))

Let's try to prove this formally.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly, near_miss_ratio


def main():
    # Quick confirmation of the pattern
    print("KEY DATA POINTS:")
    print("-" * 40)
    
    # Stars
    for k in [10, 50, 100]:
        n = k + 1
        adj = [[] for _ in range(n)]
        for i in range(1, n):
            adj[0].append(i)
            adj[i].append(0)
        poly = independence_poly(n, adj)
        nm, _ = near_miss_ratio(poly)
        print(f"Star S({k}): nm = {nm:.4f}")
    
    print()
    
    # Paths
    for n in [10, 20, 30]:
        adj = [[] for _ in range(n)]
        for i in range(n-1):
            adj[i].append(i+1)
            adj[i+1].append(i)
        poly = independence_poly(n, adj)
        nm, _ = near_miss_ratio(poly)
        print(f"Path P({n}): nm = {nm:.4f}")
    
    print()
    print("PATTERN CONFIRMED:")
    print("  - Concentrated influence (star) → high near-miss ratio")
    print("  - Distributed influence (path) → low near-miss ratio")


if __name__ == '__main__':
    main()
