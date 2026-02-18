#!/usr/bin/env python3
"""Experiment 1: Check Lorentzian / ultra-log-concave properties of tree
independence polynomials.

Huh's idea: is there a positivity criterion (Hessian, Lorentzian, ULC)
that holds for trees? If it fails, WHERE does it fail? That's interesting.

A polynomial p(x) = sum a_k x^k is:
  - Log-concave:       a_k^2 >= a_{k-1} a_{k+1} for all k
  - Ultra-log-concave: a_k^2 / C(n,k)^2 >= (a_{k-1}/C(n,k-1)) * (a_{k+1}/C(n,k+1))
  - "Hessian-positive": The matrix H_{ij} = d^2/dx_i dx_j log(p) is negative
                        semidefinite (for homogeneous Lorentzian polynomials)

We check ULC since it's stronger than LC and is the signature of Lorentzian polys.
"""

import json
import math
from indpoly import independence_poly
from graph6 import parse_graph6

def ulc_ratios(poly):
    """Compute ultra-log-concavity ratios.
    
    ULC: (a_k / C(n,k))^2 >= (a_{k-1}/C(n,k-1)) * (a_{k+1}/C(n,k+1))
    
    Returns list of ratios R_k = LHS/RHS for k=1..len-2.
    ULC holds iff all R_k >= 1.
    """
    n = len(poly) - 1  # degree
    ratios = []
    for k in range(1, n):
        if poly[k-1] == 0 or poly[k+1] == 0:
            ratios.append(float('inf'))
            continue
        # Normalized coefficients: a_k / C(n, k)
        # Use log to avoid overflow
        log_ck = math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1)
        log_ck_minus = math.lgamma(n+1) - math.lgamma(k) - math.lgamma(n-k+2)
        log_ck_plus = math.lgamma(n+1) - math.lgamma(k+2) - math.lgamma(n-k)
        
        log_lhs = 2 * (math.log(poly[k]) - log_ck)
        log_rhs = (math.log(poly[k-1]) - log_ck_minus) + (math.log(poly[k+1]) - log_ck_plus)
        
        ratios.append(math.exp(log_lhs - log_rhs))
    return ratios

def lc_ratios(poly):
    """Standard log-concavity ratios: a_k^2 / (a_{k-1} * a_{k+1})."""
    ratios = []
    for k in range(1, len(poly) - 1):
        if poly[k-1] == 0 or poly[k+1] == 0:
            ratios.append(float('inf'))
            continue
        ratios.append(poly[k]**2 / (poly[k-1] * poly[k+1]))
    return ratios

def main():
    print("=" * 70)
    print("EXPERIMENT: Ultra-Log-Concavity of Tree Independence Polynomials")
    print("=" * 70)
    print()
    
    # Test on small trees first
    print("--- Small tree examples ---")
    
    # Path P_5
    adj_p5 = [[1], [0,2], [1,3], [2,4], [3]]
    poly = independence_poly(5, adj_p5)
    print(f"P_5: poly = {poly}")
    lr = lc_ratios(poly)
    ur = ulc_ratios(poly)
    print(f"  LC ratios:  {[f'{r:.4f}' for r in lr]}")
    print(f"  ULC ratios: {[f'{r:.4f}' for r in ur]}")
    print(f"  LC holds:  {all(r >= 1 - 1e-9 for r in lr)}")
    print(f"  ULC holds: {all(r >= 1 - 1e-9 for r in ur)}")
    print()
    
    # Star K_{1,4}
    adj_star = [[1,2,3,4], [0], [0], [0], [0]]
    poly = independence_poly(5, adj_star)
    print(f"K_1,4: poly = {poly}")
    lr = lc_ratios(poly)
    ur = ulc_ratios(poly)
    print(f"  LC ratios:  {[f'{r:.4f}' for r in lr]}")
    print(f"  ULC ratios: {[f'{r:.4f}' for r in ur]}")
    print(f"  LC holds:  {all(r >= 1 - 1e-9 for r in lr)}")
    print(f"  ULC holds: {all(r >= 1 - 1e-9 for r in ur)}")
    print()

    # Now scan the LC-failing trees at n=26
    print("--- n=26 LC-failing trees ---")
    with open("results/analysis_n26.json") as f:
        data = json.load(f)
    
    for i, tree in enumerate(data.get("lc_failures", [])):
        poly = tree["poly"]
        lr = lc_ratios(poly)
        ur = ulc_ratios(poly)
        min_lr = min(lr)
        min_ur = min(ur)
        min_lr_pos = lr.index(min_lr) + 1
        min_ur_pos = ur.index(min_ur) + 1
        print(f"  LC-failure #{i+1}:")
        print(f"    min LC ratio = {min_lr:.6f} at k={min_lr_pos}")
        print(f"    min ULC ratio = {min_ur:.6f} at k={min_ur_pos}")
        print(f"    ULC holds: {all(r >= 1 - 1e-9 for r in ur)}")
        print()

    # Scan top near-misses  
    print("--- Top 10 near-miss trees at n=26 ---")
    nm_trees = data.get("top_near_misses", [])[:10]
    ulc_failures = 0
    for i, tree in enumerate(nm_trees):
        poly = tree["poly"]
        ur = ulc_ratios(poly)
        min_ur = min(ur)
        min_ur_pos = ur.index(min_ur) + 1
        fails = not all(r >= 1 - 1e-9 for r in ur)
        if fails:
            ulc_failures += 1
        print(f"  NM #{i+1} (nm={tree['nm_ratio']:.4f}): min ULC = {min_ur:.4f} at k={min_ur_pos} {'FAIL' if fails else 'ok'}")
    
    print(f"\n  ULC failures among top 10 near-misses: {ulc_failures}/10")
    print()

    # Broader scan: check ULC for ALL trees at small n
    print("--- Exhaustive ULC check for n <= 14 ---")
    import subprocess
    for n in range(3, 15):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        ulc_fail_count = 0
        worst_ratio = float('inf')
        total = 0
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            ur = ulc_ratios(poly)
            total += 1
            if not all(r >= 1 - 1e-9 for r in ur):
                ulc_fail_count += 1
                worst_ratio = min(worst_ratio, min(ur))
        
        if ulc_fail_count > 0:
            print(f"  n={n}: {ulc_fail_count}/{total} ULC failures (worst ratio: {worst_ratio:.6f})")
        else:
            print(f"  n={n}: {total} trees, ALL ultra-log-concave ✓")

    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
