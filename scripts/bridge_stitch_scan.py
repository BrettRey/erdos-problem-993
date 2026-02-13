#!/usr/bin/env python3
"""
Bridge Stitching Search:
Constructs synthetic trees by joining two rooted trees T1 and T2 via an edge between their roots.
The independence polynomial of the joined tree T is given by:
I(T) = I(T1-r1) * I(T2) + x * I(T1-N[r1]) * I(T2-r2)
     = A0 * (B0 + x*B1) + x * A1 * B0
     = A0*B0 + x*A0*B1 + x*A1*B0

where:
  A0 = I(T1-r1)
  A1 = I(T1-N[r1])
  B0 = I(T2-r2)
  B1 = I(T2-N[r2])

This effectively explores the space of trees formed by bridging two arbitrary trees.
"""

import argparse
import json
import os
import sys
import time
import random
from typing import Any, Tuple, List

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from trees import trees_geng
from indpoly import independence_poly, is_unimodal, _polymul, _polyadd, log_concavity_ratio

def get_rooted_signatures(n: int) -> List[Tuple[List[int], List[int]]]:
    """
    Generates signatures for all rooted trees of size n.
    Signature = (I(T-r), I(T-N[r]))
    """
    signatures = []
    # We iterate all trees. For each tree, we consider EACH node as a possible root.
    # Note: trees_geng yields unrooted graphs.
    for n_nodes, adj in trees_geng(n):
        for root in range(n_nodes):
            # Compute I(T-r)
            # T-r is simply the graph without vertex 'root'.
            # We can construct the adj of T-r by removing 'root'.
            # But simpler to just use a modified adjacency list or a helper.
            # However, independence_poly takes adj.
            
            # Construct T-r adj
            # Mapping: old node -> new node index
            # We skip 'root'.
            mapping = {}
            new_idx = 0
            for i in range(n_nodes):
                if i != root:
                    mapping[i] = new_idx
                    new_idx += 1
            
            adj_minus_r = [[] for _ in range(n_nodes - 1)]
            for u in range(n_nodes):
                if u == root: continue
                for v in adj[u]:
                    if v == root: continue
                    adj_minus_r[mapping[u]].append(mapping[v])
            
            poly_minus_r = independence_poly(n_nodes - 1, adj_minus_r)
            
            # Construct T-N[r]
            # Remove root AND its neighbors
            neighbors = set(adj[root])
            to_remove = neighbors | {root}
            
            mapping_nr = {}
            new_idx_nr = 0
            for i in range(n_nodes):
                if i not in to_remove:
                    mapping_nr[i] = new_idx_nr
                    new_idx_nr += 1
            
            adj_minus_nr = [[] for _ in range(n_nodes - len(to_remove))]
            for u in range(n_nodes):
                if u in to_remove: continue
                for v in adj[u]:
                    if v in to_remove: continue
                    adj_minus_nr[mapping_nr[u]].append(mapping_nr[v])
            
            poly_minus_nr = independence_poly(n_nodes - len(to_remove), adj_minus_nr)
            
            signatures.append((poly_minus_r, poly_minus_nr))
            
    return signatures

def solve_bridge_poly(sig1: Tuple[List[int], List[int]], sig2: Tuple[List[int], List[int]]) -> List[int]:
    """
    Computes I(T) where T is formed by bridging T1(root r1) and T2(root r2).
    sig1 = (A0, A1) = (I(T1-r1), I(T1-N[r1]))
    sig2 = (B0, B1) = (I(T2-r2), I(T2-N[r2]))
    
    I(T) = A0*B0 + x*(A0*B1 + A1*B0)
    """
    A0, A1 = sig1
    B0, B1 = sig2
    
    # Term 1: A0 * B0
    T1 = _polymul(A0, B0)
    
    # Term 2: x * (A0*B1 + A1*B0)
    # Inner: A0*B1 + A1*B0
    inner = _polyadd(_polymul(A0, B1), _polymul(A1, B0))
    # Multiply by x (shift right)
    T2 = [0] + inner
    
    return _polyadd(T1, T2)

def main():
    parser = argparse.ArgumentParser(description="Bridge Stitching Search")
    parser.add_argument("--max-n", type=int, default=12, help="Max size of component trees (total size ~2*max-n)")
    parser.add_argument("--random-samples", type=int, default=100000, help="Number of random pairs checked per large (n1,n2) batch")
    parser.add_argument("--seed", type=int, default=993, help="Random seed for sampling")
    parser.add_argument("--out", type=str, default="results/bridge_stitch_scan.json")
    args = parser.parse_args()
    random.seed(args.seed)
    
    print(f"Generating signatures up to n={args.max_n}...")
    signatures_by_n = {}
    total_signatures = 0
    
    for n in range(1, args.max_n + 1):
        sigs = get_rooted_signatures(n)
        # Deduplicate signatures
        # Tuples of lists aren't hashable, convert to tuple of tuples
        unique_sigs = list(set( (tuple(a), tuple(b)) for a, b in sigs ))
        # Convert back to list of lists for processing
        signatures_by_n[n] = [ (list(a), list(b)) for a, b in unique_sigs ]
        print(f"n={n}: {len(unique_sigs)} unique rooted signatures")
        total_signatures += len(unique_sigs)
        
    print(f"Total unique rooted signatures: {total_signatures}")
    
    # Search phase
    print("Starting search...")
    start_time = time.time()
    checked = 0
    failures = []
    
    max_lc_ratio = 0.0
    max_lc_poly = []
    
    # We want to check pairs (n1, n2) such that n1 + n2 <= 2 * args.max_n
    # Or just iterate all pairs.
    # To cover larger n, we probably want n1 near max_n and n2 near max_n.
    
    # Exhaustive for small sums, random for large?
    # Let's try exhaustive for all combinations up to limit if feasible.
    # With n=12, we might have ~5000 signatures. 5000^2 = 25M pairs. Feasible.
    
    # Iterate n1 from 1 to max_n
    # Iterate n2 from n1 to max_n (symmetry)
    
    for n1 in range(1, args.max_n + 1):
        for n2 in range(n1, args.max_n + 1):
            sigs1 = signatures_by_n[n1]
            sigs2 = signatures_by_n[n2]
            
            # Determine if we should sample
            num_s1 = len(sigs1)
            num_s2 = len(sigs2)
            total_pairs = num_s1 * num_s2
            
            use_sampling = False
            if total_pairs > 500000:
                use_sampling = True
                
            print(f"Checking pairs n1={n1} ({num_s1}) x n2={n2} ({num_s2}) - Total {total_pairs} {'(Sampling)' if use_sampling else ''}...")
            
            # If sampling, pick random indices
            # Otherwise iterate all
            
            iterations = args.random_samples if use_sampling else total_pairs
            
            # Create iterators or random samplers
            if use_sampling:
                # Reservoir sampling or just random choice? Random choice with replacement is fine for large space
                pairs_iter = ((random.choice(sigs1), random.choice(sigs2)) for _ in range(iterations))
            else:
                pairs_iter = ((s1, s2) for s1 in sigs1 for s2 in sigs2)

            best_lc_this_batch = 0.0
            
            for i, (s1, s2) in enumerate(pairs_iter):
                checked += 1
                poly = solve_bridge_poly(s1, s2)
                
                if not is_unimodal(poly):
                    print(f"FOUND NON-UNIMODAL! n={n1+n2}")
                    print(f"Poly: {poly}")
                    failures.append({
                        "n": n1+n2,
                        "n1": n1,
                        "n2": n2,
                        "sig1": s1,
                        "sig2": s2,
                        "poly": poly
                    })
                    with open(args.out, "w") as f:
                        json.dump(failures, f, indent=2)
                    return 
                
                lc, _ = log_concavity_ratio(poly)
                if lc > max_lc_ratio:
                    max_lc_ratio = lc
                    max_lc_poly = poly
                    # print(f"New Max LC Ratio: {lc} at n={n1+n2}")
                if lc > best_lc_this_batch:
                    best_lc_this_batch = lc

                if i % 100000 == 0 and i > 0:
                     print(f"  ...checked {i} pairs. Max LC so far: {max_lc_ratio:.6f}")
            
            print(f"  Batch finished. Best LC: {best_lc_this_batch:.6f}")
                    
    print(f"Finished. Checked {checked} pairs.")
    print(f"Failures found: {len(failures)}")
    print(f"Global Max LC Ratio found: {max_lc_ratio}")
    if max_lc_poly:
        print(f"Poly with max LC: {max_lc_poly}")
    
    if not failures:
        # Save empty result to indicate completion
        with open(args.out, "w") as f:
            json.dump({
                "status": "clean", 
                "checked": checked, 
                "max_n": args.max_n,
                "random_samples": args.random_samples,
                "seed": args.seed,
                "max_lc_ratio": max_lc_ratio,
                "max_lc_poly": max_lc_poly
            }, f, indent=2)

if __name__ == "__main__":
    main()
