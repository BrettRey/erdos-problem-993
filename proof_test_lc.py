#!/usr/bin/env python3
"""Is mean ∈ [M-1, M] at tie points a property of ALL LC sequences?

Test with random LC sequences and see if the property holds.

A log-concave sequence a_0, a_1, ..., a_d satisfies a_k^2 ≥ a_{k-1} a_{k+1}.

We can generate an LC sequence by taking the coefficients of a product 
of linear factors (1 + t_i x) with t_i > 0. Such products are always LC.
"""

import random
import math

def make_lc_poly(factors):
    """Multiply (1 + t_1 x)(1 + t_2 x)...(1 + t_k x)."""
    poly = [1]
    for t in factors:
        new_poly = [0] * (len(poly) + 1)
        for i, c in enumerate(poly):
            new_poly[i] += c
            new_poly[i+1] += c * t
        poly = new_poly
    return poly

def check_tie_property(poly):
    """Check if mean ∈ [M-1, M] at every tie point."""
    d = len(poly) - 1
    violations = 0
    worst_excess = -100
    
    for M in range(1, d + 1):
        if poly[M] > 0 and poly[M-1] > 0:
            lam0 = poly[M-1] / poly[M]
            weights = [poly[k] * (lam0 ** k) for k in range(len(poly))]
            total_w = sum(weights)
            if total_w > 0:
                mean_val = sum(k * weights[k] for k in range(len(weights))) / total_w
                excess = M - mean_val
                if excess > worst_excess:
                    worst_excess = excess
                if mean_val < M - 1 - 1e-8:
                    violations += 1
    
    return violations, worst_excess

def main():
    print("=" * 80)
    print("  TEST: Is mean ∈ [M-1, M] at tie points for ALL LC polynomials?")
    print("=" * 80)
    print()
    
    random.seed(42)
    
    # Test 1: Products of (1 + t_i x) — these are always LC
    print("  Test 1: Products of (1 + t_i x), t_i ~ Uniform(0.01, 100)")
    total_polys = 0
    total_violations = 0
    max_excess = -100
    
    for trial in range(10000):
        k = random.randint(2, 20)
        factors = [random.uniform(0.01, 100) for _ in range(k)]
        poly = make_lc_poly(factors)
        v, e = check_tie_property(poly)
        total_polys += 1
        total_violations += v
        if e > max_excess:
            max_excess = e
    
    print(f"    Polys tested: {total_polys}")
    print(f"    Violations: {total_violations}")
    print(f"    Max excess at tie: {max_excess:.6f}")
    print()
    
    # Test 2: Random unimodal but NOT necessarily LC sequences
    print("  Test 2: Random UNIMODAL (but not LC) sequences")
    total_polys = 0
    total_violations = 0
    max_excess = -100
    worst_poly = None
    
    for trial in range(10000):
        d = random.randint(3, 15)
        mode = random.randint(1, d-1)
        
        # Generate unimodal sequence: increase to mode, decrease after
        poly = [0] * (d + 1)
        poly[mode] = random.uniform(50, 100)
        
        for k in range(mode - 1, -1, -1):
            poly[k] = poly[k+1] * random.uniform(0.3, 0.99)
        
        for k in range(mode + 1, d + 1):
            poly[k] = poly[k-1] * random.uniform(0.3, 0.99)
        
        v, e = check_tie_property(poly)
        total_polys += 1
        total_violations += v
        if e > max_excess:
            max_excess = e
            worst_poly = poly[:]
    
    print(f"    Polys tested: {total_polys}")
    print(f"    Violations: {total_violations}")
    print(f"    Max excess at tie: {max_excess:.6f}")
    
    if total_violations > 0 and worst_poly:
        print(f"    Worst poly: {[round(x, 2) for x in worst_poly]}")
    print()
    
    # Test 3: Specifically try to BREAK the property with pathological sequences
    print("  Test 3: Pathological sequences — try to get excess ≥ 1")
    
    # A sequence with a very sharp peak and heavy left tail
    for sharpness in [10, 50, 100, 500, 1000]:
        d = 20
        mode = 15  # mode far to the right
        poly = [0] * (d + 1)
        poly[mode] = sharpness
        
        # Heavy left tail
        for k in range(mode - 1, -1, -1):
            poly[k] = poly[k+1] * 0.95  # slow decay left
        
        # Sharp right drop
        for k in range(mode + 1, d + 1):
            poly[k] = poly[k-1] * 0.1  # fast decay right
        
        v, e = check_tie_property(poly)
        print(f"    sharpness={sharpness}: excess={e:.6f}, violations={v}")
    
    print()
    
    # Test 4: Extreme left-heavy distribution
    print("  Test 4: Extreme left-heavy — constant plateau then sharp peak")
    
    for plateau in [5, 10, 15]:
        for d in [20]:
            mode = d  # peak at the end
            poly = [1.0] * (plateau + 1)  # constant left part
            for k in range(plateau + 1, d):
                poly.append(poly[-1] * 1.1)  # slow rise
            poly.append(poly[-1] * 2.0)  # sharp jump at mode
            
            v, e = check_tie_property(poly)
            print(f"    plateau={plateau}, d={d}: excess={e:.6f}, violations={v}")
    
    print()
    
    # Test 5: The CRITICAL test — can we construct a unimodal sequence
    # where the tie-point excess exceeds 1?
    print("  Test 5: Systematic search for excess ≥ 1")
    
    found = False
    for trial in range(100000):
        d = random.randint(4, 25)
        mode = random.randint(2, d-1)
        
        poly = [0] * (d + 1)
        poly[mode] = 100
        
        # Left: slow decay (heavy left tail)
        for k in range(mode - 1, -1, -1):
            r = random.uniform(0.8, 0.999)
            poly[k] = poly[k+1] * r
        
        # Right: fast decay (light right tail)
        for k in range(mode + 1, d + 1):
            r = random.uniform(0.01, 0.5)
            poly[k] = poly[k-1] * r
        
        v, e = check_tie_property(poly)
        if v > 0:
            found = True
            print(f"    VIOLATION at trial {trial}: d={d}, mode={mode}, excess={e:.6f}")
            print(f"    poly = {[round(x, 4) for x in poly]}")
            break
    
    if not found:
        print(f"    No violations in 100K random unimodal sequences!")
    
    print()
    print("=" * 80)

if __name__ == "__main__":
    main()
