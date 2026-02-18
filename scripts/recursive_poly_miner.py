#!/usr/bin/env python3
"""
Rooted Forest Beam Search:
Iteratively constructs trees by merging subtrees at a new root.
State: (A, B) where A = I(T), B = I(T-root).
Transition: Pick k subtrees (A_i, B_i).
  New T has root v connected to roots of T_i.
  A_new = prod(A_i) + x * prod(B_i)
  B_new = prod(A_i)
Target: Maximize Log-Concavity violation of A_new.
"""

import argparse
import json
import random
import sys
import os
import math
from typing import List, Tuple, Dict

# Ensure we can import from parent directory
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import is_unimodal, log_concavity_ratio, _polymul, _polyadd

def prod_polys(polys: List[List[int]]) -> List[int]:
    if not polys:
        return [1]
    res = polys[0]
    for p in polys[1:]:
        res = _polymul(res, p)
    return res

def get_stats(poly: List[int]) -> Tuple[float, bool]:
    lc_ratio, _ = log_concavity_ratio(poly)
    is_uni = is_unimodal(poly)
    return lc_ratio, is_uni

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generations", type=int, default=50, help="Number of growth generations")
    parser.add_argument("--beam-width", type=int, default=2000, help="Max signatures to keep")
    parser.add_argument("--max-degree", type=int, default=5, help="Max subtrees to combine (degree of new root)")
    parser.add_argument("--out", type=str, default="results/forest_beam_search.json")
    args = parser.parse_args()

    # Pool: List of tuples (A, B, n)
    # A = I(T), B = I(T-r), n = size
    # Initial: Single vertex. I(K1)=[1,1], I(K1-r)=[1].
    pool: List[Tuple[List[int], List[int], int]] = [([1, 1], [1], 1)]
    
    best_lc_ratio = 0.0
    best_lc_poly = []
    best_n = 0
    
    print(f"Starting Rooted Forest Beam Search (Gens={args.generations}, Beam={args.beam_width})...")
    
    for gen in range(1, args.generations + 1):
        new_items = []
        
        # We want to form new trees.
        # Strategy: Randomly sample sets of k subtrees from the current pool.
        # Since exact enumeration is hard, we'll produce (beam_width * 5) candidates and prune.
        
        num_candidates = args.beam_width * 10
        
        for _ in range(num_candidates):
            # Pick degree k
            # Bias towards small k to allow path-like growth
            # Weights for k=1..5
            k = random.choices([1, 2, 3, 4, 5], weights=[0.4, 0.3, 0.1, 0.1, 0.1], k=1)[0]
            
            # Pick k components
            # Bias towards higher LC ratio (pool is sorted? we will sort it)
            # Tournament selection? Or just random from elite pool.
            
            components = random.choices(pool, k=k)
            
            As = [c[0] for c in components]
            Bs = [c[1] for c in components]
            ns = [c[2] for c in components]
            
            prod_A = prod_polys(As)
            prod_B = prod_polys(Bs)
            
            # New A = prod_A + x * prod_B
            # x * prod_B is shifted right by 1
            term2 = [0] + prod_B
            A_new = _polyadd(prod_A, term2)
            
            # New B = prod_A
            B_new = prod_A
            
            n_new = sum(ns) + 1
            
            # Check stats
            lc, uni = get_stats(A_new)
            
            if not uni:
                print(f"FOUND NON-UNIMODAL TREE at Gen {gen}, N={n_new}!")
                print(f"Poly: {A_new}")
                result = {
                    "found": True,
                    "n": n_new,
                    "poly": A_new,
                    "lc_ratio": lc
                }
                with open(args.out, "w") as f:
                    json.dump(result, f, indent=2)
                return
            
            if lc > best_lc_ratio:
                best_lc_ratio = lc
                best_lc_poly = A_new
                best_n = n_new
                
            new_items.append((A_new, B_new, n_new, lc))
            
        # Pruning
        # Sort by LC ratio
        new_items.sort(key=lambda x: x[3], reverse=True)
        
        # Keep top unique
        # Use hashing of poly A to dedupe
        seen = set()
        next_pool = []
        for item in new_items:
            sig = tuple(item[0])
            if sig in seen:
                continue
            seen.add(sig)
            # Store back in pool format (A, B, n)
            next_pool.append((item[0], item[1], item[2]))
            if len(next_pool) >= args.beam_width:
                break
        
        pool = next_pool
        print(f"Gen {gen}: Best LC={best_lc_ratio:.6f} (N={best_n}), Pool size={len(pool)}")
        
        # Inject fresh blood? 
        # Base case always available? 
        # Actually, since we combine from current pool, pool grows in N.
        # We should keep some small trees in pool to allow "unbalanced" growth (broom-like).
        # Add the base K1 back.
        pool.append(([1, 1], [1], 1))
        
        # Add a random broom-like structure to pool every gen
        # Broom(n, s): Path P_{n-s} attached to Star S_s center
        # We can simulate this by combining a Path and Star from basic components?
        # Simpler: just inject small paths/stars if we had them.
        # For now, base K1 injection helps.

    print("Finished search.")
    print(f"Best LC Ratio: {best_lc_ratio} at N={best_n}")
    
    with open(args.out, "w") as f:
        json.dump({
            "found": False,
            "max_n": best_n,
            "max_lc_ratio": best_lc_ratio,
            "max_lc_poly": best_lc_poly
        }, f, indent=2)

if __name__ == "__main__":
    main()
