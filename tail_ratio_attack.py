#!/usr/bin/env python3
"""
Attack on the Subdivision Lemma: Tail Ratio Monotonicity

Goal: Prove that for forest independence polynomials, the coefficient ratios
are monotonic in the tail region.

Key definitions from subdivision lemma:
- p = P_u P_v (product of two "exclude-root" polynomials)
- j = q + r = Q_u Q_v + P_u Q_v + Q_u P_v (the "at least one root" polynomial)

We need to show tail ratios are nonincreasing for these polynomials.
"""

import sys
sys.path.insert(0, '.')

from indpoly import independence_poly
from collections import defaultdict
import random


def compute_ratio monotonicity(poly):
    """Compute coefficient ratios and check if they're decreasing in the tail."""
    # Find first descent
    d = -1
    for i in range(1, len(poly)):
        if poly[i] < poly[i-1]:
            d = i
            break
    
    if d == -1:
        return True, {}  # Sequence is nondecreasing
    
    # Check tail ratios from d onwards
    tail_ratios = []
    for i in range(d, len(poly)-1):
        if poly[i] > 0:
            tail_ratios.append(poly[i+1] / poly[i])
    
    # Check if ratios are nonincreasing
    is_monotonic = all(tail_ratios[i] >= tail_ratios[i+1] for i in range(len(tail_ratios)-1))
    
    return is_monotonic, {'d': d, 'ratios': tail_ratios}


def analyze_forest_polynomials():
    """Analyze ratio monotonicity for various forest structures."""
    print("=" * 70)
    print("TAIL RATIO MONOTONICITY FOR FOREST POLYNOMIALS")
    print("=" * 70)
    
    import networkx as nx
    
    # Test 1: Products of paths
    print("\n1. Products of paths (P_n * P_m):")
    for n in [5, 10, 15]:
        for m in [5, 10]:
            # Create two paths and compute their polynomials
            # Actually, we need forest polynomials
            pass
    
    # Test 2: Random forest polynomials
    print("\n2. Random forest polynomials:")
    for _ in range(100):
        # Create random forest (collection of trees)
        pass


def explore_ratio_patterns():
    """Explore what makes tail ratios monotonic or not."""
    print("=" * 70)
    print("EXPLORING RATIO PATTERNS")
    print("=" * 70)
    
    import networkx as nx
    from targeted import make_broom, make_spider
    
    # Key insight: The polynomials P_u and Q_u for rooted trees
    # might have special properties that guarantee tail monotonicity
    
    # Let me explore the actual coefficient behavior
    print("\nExploring broom coefficient ratios:")
    
    for s in [10, 50, 100]:
        n, adj = make_broom(5, s)
        poly = independence_poly(n, adj)
        
        # Find peak
        peak = max(poly)
        peak_idx = poly.index(peak)
        
        # Compute ratios around peak
        print(f"\nBroom(5,{s}): n={n}, peak at k={peak_idx}")
        print("  Ratios around peak:")
        for i in range(max(0, peak_idx-2), min(len(poly)-1, peak_idx+3)):
            if poly[i] > 0:
                ratio = poly[i+1] / poly[i]
                direction = "↓" if i >= peak_idx and ratio < 1 else "↑" if ratio > 1 else "="
                print(f"    r_{i} = {ratio:.3f} {direction}")


def test_specific_cases():
    """Test specific cases that might break monotonicity."""
    print("\n" + "=" * 70)
    print("TESTING SPECIFIC CASES")
    print("=" * 70)
    
    import networkx as nx
    
    # Test cases that might have non-monotonic ratios
    # The n=26 LC failures
    
    with open('results/analysis_n26.json') as f:
        import json
        data = json.load(f)
    
    # We can't easily reconstruct the trees from graph6, but let's try random
    print("\nRandom trees with various structures:")
    
    for n in [15, 20, 25]:
        monotonic_count = 0
        for _ in range(50):
            G = nx.random_labeled_tree(n, seed=random.randint(0, 10000))
            adj = [list(G.neighbors(v)) for v in range(n)]
            poly = independence_poly(n, adj)
            is_monotonic, _ = compute_ratio_monotonicity(poly)
            if is_monotonic:
                monotonic_count += 1
        
        print(f"  n={n}: {monotonic_count}/50 have monotonic tail ratios")


if __name__ == '__main__':
    explore_ratio_patterns()
    test_specific_cases()
