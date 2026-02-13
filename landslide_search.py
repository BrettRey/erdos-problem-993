"""
Landslide Search: A physics-inspired heuristic for finding unimodality counterexamples.

Strategy:
1. Construct "Base" trees that have steep coefficient slopes (Brooms).
2. Identify "Weak Layers" where the coefficient descent is fragile.
3. "Graft" small "debris" trees onto the base to disrupt the slope and trigger a "landslide" (second peak).
"""

import argparse
import random
import sys
from typing import Iterator

# Import from existing modules (assuming they are in the python path)
try:
    from trees import trees
    from indpoly import independence_poly
except ImportError:
    # Fallback if running from a different directory depth without path setup
    sys.path.append(".")
    from trees import trees
    from indpoly import independence_poly


def get_poly(adj: list[list[int]]) -> list[int]:
    """Compute independence polynomial for a tree adjacency list."""
    n = len(adj)
    return independence_poly(n, adj)


def is_unimodal(seq: list[int]) -> bool:
    """Check if sequence is unimodal."""
    if len(seq) <= 2:
        return True
    
    # Find peak
    peak_idx = 0
    for i in range(1, len(seq)):
        if seq[i] >= seq[i-1]:
            peak_idx = i
        else:
            break
            
    # Check for any increase after peak
    for i in range(peak_idx + 1, len(seq)):
        if seq[i] > seq[i-1]:
            return False
            
    return True


def analyze_slope(poly: list[int]) -> list[float]:
    """Return ratios r_k = poly[k] / poly[k-1]."""
    ratios = []
    for k in range(1, len(poly)):
        if poly[k-1] == 0:
            ratios.append(float('inf'))
        else:
            ratios.append(poly[k] / poly[k-1])
    return ratios


def graft(base_adj: list[list[int]], base_node: int, 
          debris_adj: list[list[int]], debris_node: int) -> list[list[int]]:
    """
    Graft debris tree onto base tree by adding an edge between base_node and debris_node.
    Returns new adjacency list.
    """
    n_base = len(base_adj)
    n_debris = len(debris_adj)
    n_total = n_base + n_debris
    
    new_adj = [[] for _ in range(n_total)]
    
    # Copy base
    for u in range(n_base):
        new_adj[u] = list(base_adj[u])
        
    # Copy debris (shifted by n_base)
    for u in range(n_debris):
        u_mapped = u + n_base
        for v in debris_adj[u]:
            v_mapped = v + n_base
            new_adj[u_mapped].append(v_mapped)
            
    # Add bridge edge
    u_bridge = base_node
    v_bridge = debris_node + n_base
    
    new_adj[u_bridge].append(v_bridge)
    new_adj[v_bridge].append(u_bridge)
    
    # Sort for consistency
    for u in range(n_total):
        new_adj[u].sort()
        
    return new_adj


def generate_brooms(n_range: range) -> Iterator[tuple[int, list[list[int]]]]:
    """Yield Broom trees B(n, s).
    A broom B(n, s) consists of a path P_{n-s} and a star K_{1, s} 
    identified at a leaf of the path and the center of the star.
    Actually, usually B(n, s) is path of length n-s attached to center of star size s.
    Let's stick to the definition: Path of length L, one end attached to center of Star size S.
    Total vertices n = L + S + 1 (center).
    """
    for n in n_range:
        # Try various handle lengths.
        # Handle length 1 to n-2 (leaving at least 1 leaf in star)
        for handle_len in range(1, n - 1):
            star_leaves = n - 1 - handle_len
            if star_leaves < 1: continue
            
            adj = [[] for _ in range(n)]
            
            # Create handle (path 0 to handle_len)
            # 0 is the "root" of the handle, handle_len is the center of the star
            for i in range(handle_len):
                adj[i].append(i+1)
                adj[i+1].append(i)
                
            # Create star leaves attached to handle_len
            center = handle_len
            for i in range(star_leaves):
                leaf = center + 1 + i
                adj[center].append(leaf)
                adj[leaf].append(center)
                
            yield (n, adj)


def generate_random_trees(n_range: range, count: int) -> Iterator[tuple[int, list[list[int]]]]:
    """Yield random trees using networkx."""
    import networkx as nx
    for n in n_range:
        for _ in range(count):
            T = nx.random_labeled_tree(n)
            adj = [[] for _ in range(n)]
            for u, v in T.edges():
                adj[u].append(v)
                adj[v].append(u)
            yield (n, adj)


def main():
    parser = argparse.ArgumentParser(description="Landslide Search")
    parser.add_argument("--base-type", type=str, default="broom", choices=["broom", "random"], help="Type of base tree")
    parser.add_argument("--base-n", type=int, default=100, help="Size of base tree")
    parser.add_argument("--base-count", type=int, default=100, help="Number of random base trees to try")
    parser.add_argument("--debris-n", type=int, default=12, help="Max size of debris trees")
    parser.add_argument("--samples", type=int, default=1000, help="Number of debris samples")
    args = parser.parse_args()
    
    print(f"Starting Landslide Search: Base ({args.base_type}) N={args.base_n}, Debris N<={args.debris_n}")
    
    # 1. Generate Debris Catalog
    print("Building debris catalog...")
    debris_catalog = []
    # Mix of small random trees and exact enumeration for very small n
    for n in range(4, args.debris_n + 1):
        if n <= 10:
             # Exact
             count = 0
             # Use trees() wrapper to handle backend selection (geng/networkx)
             for _, adj in trees(n):
                 debris_catalog.append(adj)
                 count += 1
             # print(f"  n={n}: {count} trees")
        else:
             # Random sample (using networkx for simplicity if available, or just skip)
             # For now, let's just use what we have. If n > 10 is requested, 
             # we might want random trees.
             # Let's stick to n<=12 exact if feasible, or just limit sample.
             # geng is fast. n=12 is ~500k trees, maybe too many for RAM if we keep all.
             # Let's just keep a sample if n is large.
             count = 0
             limit = 1000
             for _, adj in trees(n):
                 debris_catalog.append(adj)
                 count += 1
                 if count >= limit: break
    
    print(f"Catalog size: {len(debris_catalog)} trees")
    
    # 2. Iterate Base Trees
    if args.base_type == "broom":
        base_gen = generate_brooms(range(args.base_n, args.base_n + 1))
    else:
        base_gen = generate_random_trees(range(args.base_n, args.base_n + 1), args.base_count)
    
    best_ratio = 0.0
    
    for n_base, base_adj in base_gen:
        base_poly = get_poly(base_adj)
        
        # Analyze Base
        if not is_unimodal(base_poly):
            print(f"Found counterexample in BASE! {base_poly}")
            return

        # Find "Weak Layer"
        # We look for the first descent, then look for a "terrace"
        # i.e., after peak, where ratio i_{k+1}/i_k is largest (closest to 1)
        peak = base_poly.index(max(base_poly))
        
        candidates = []
        for k in range(peak, len(base_poly) - 1):
            if base_poly[k] == 0: continue
            ratio = base_poly[k+1] / base_poly[k]
            candidates.append((ratio, k))
            
        # Sort by ratio descending (steepest/closest to flat first)
        candidates.sort(key=lambda x: x[0], reverse=True)
        
        # Try grafting onto the top 3 weak spots
        # We need to map k (index in poly) to a node? No, we don't know which node contributes to which k.
        # We just graft to *various nodes* and hope it hits the mechanics.
        # Heuristic: Grafting to the end of the handle vs center of star.
        
        # Nodes to try:
        # 0 (start of handle)
        # mid handle
        # center of star (max degree node)
        # a leaf of the star
        
        degrees = [len(base_adj[u]) for u in range(n_base)]
        max_deg_node = degrees.index(max(degrees))
        leaves = [u for u, d in enumerate(degrees) if d == 1]
        leaf_handle = 0 # by construction
        leaf_star = leaves[-1] if leaves else 0
        
        targets = [0, n_base // 2, max_deg_node, leaf_star]
        targets = list(set(targets)) # unique
        
        # Optimization: Don't check every debris on every base.
        # Just check a random subset of debris.
        
        sample_debris = random.sample(debris_catalog, min(len(debris_catalog), args.samples))
        
        for debris_adj in sample_debris:
            # Try grafting debris root (0) to target
            for t_node in targets:
                # We can graft debris node 0.
                # Ideally we should try all debris nodes, but 0 is usually a leaf or specific if generated by geng?
                # Geng canonical labeling might not put specific node at 0.
                # Let's just try node 0 of debris.
                
                grafted_adj = graft(base_adj, t_node, debris_adj, 0)
                g_poly = get_poly(grafted_adj)
                
                if not is_unimodal(g_poly):
                    print(f"\n!!! COUNTEREXAMPLE FOUND !!!")
                    print(f"Base N: {n_base}, Debris N: {len(debris_adj)}")
                    print(f"Graft Node on Base: {t_node}")
                    print(f"Poly: {g_poly}")
                    return
                
                # Check "near miss"
                # If we made a "dip", report it.
                # A dip is a local min.
                # Or just check max ratio in tail.
                
                # Manually check near miss ratio in tail
                curr_peak = g_poly.index(max(g_poly))
                local_worst = 0.0
                for k in range(curr_peak, len(g_poly) - 1):
                    if g_poly[k] == 0: continue
                    r = g_poly[k+1] / g_poly[k]
                    if r > local_worst:
                        local_worst = r
                        
                if local_worst > best_ratio:
                    best_ratio = local_worst
                    print(f"New best ratio: {best_ratio:.5f} (Base len {n_base})")
                    if best_ratio > 0.99:
                         print(f"  Poly: {g_poly}")
                         print(f"  Graft: Base Node {t_node} + Debris (size {len(debris_adj)})")

if __name__ == "__main__":
    main()
