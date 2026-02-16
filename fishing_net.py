#!/usr/bin/env python3
"""
Fishing-Net Search: Market-Inspired Approach to Unimodality

Inspired by the Kerala fishermen story:
- Markets = subtrees of a tree
- Prices = independent set counts i_k
- Cellphone connectivity = structural information flow

Hypothesis: Trees with "isolated markets" (multiple far-apart branch vertices
with different-sized subtrees) are more likely to have unimodality violations.

The "fishing net" captures the idea of creating connections that allow
information (independent set counts) to flow between subtrees.
"""

import json
import random
import time
from collections import defaultdict
from itertools import combinations
from typing import List, Tuple, Dict, Any, Optional
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from indpoly import independence_poly, is_unimodal, is_log_concave, near_miss_ratio, log_concavity_ratio
# from trees import TreeGenerator  # Not needed for this approach


def compute_market_isolation_metrics(G) -> Dict[str, float]:
    """
    Compute metrics inspired by market information flow.
    High values suggest "isolated markets" that may lead to unimodality violations.
    """
    import networkx as nx
    import numpy as np
    
    n = len(G.nodes())
    
    # Branch vertices (degree >= 3)
    branch_vertices = [v for v in G.nodes() if G.degree(v) >= 3]
    num_branches = len(branch_vertices)
    
    # Distance between branch vertices
    if len(branch_vertices) >= 2:
        branch_distances = []
        for i, b1 in enumerate(branch_vertices):
            for b2 in branch_vertices[i+1:]:
                try:
                    branch_distances.append(nx.shortest_path_length(G, b1, b2))
                except nx.NetworkXNoPath:
                    branch_distances.append(n)
        avg_branch_dist = np.mean(branch_distances)
        max_branch_dist = max(branch_distances)
    else:
        avg_branch_dist = n
        max_branch_dist = n
    
    # Subtree size variance at branch vertices
    branch_subtree_sizes = []
    for b in branch_vertices:
        children = list(G.neighbors(b))
        for c in children:
            subtree = nx.descendants(G, c)
            branch_subtree_sizes.append(len(subtree) + 1)
    
    if branch_subtree_sizes:
        subtree_variance = np.var(branch_subtree_sizes)
        subtree_max_min_ratio = max(branch_subtree_sizes) / max(min(branch_subtree_sizes), 1)
    else:
        subtree_variance = 0
        subtree_max_min_ratio = 1
    
    # Core size (degree >= 2)
    core_vertices = [v for v in G.nodes() if G.degree(v) >= 2]
    core_size = len(core_vertices)
    
    # Path bottleneck
    max_bottleneck = 0
    for v in G.nodes():
        if G.degree(v) == 2:
            # Count consecutive degree-2 vertices
            length = 1
            queue = [v]
            visited = {v}
            while queue:
                curr = queue.pop(0)
                for nbr in G.neighbors(curr):
                    if nbr not in visited and G.degree(nbr) == 2:
                        visited.add(nbr)
                        queue.append(nbr)
                        length += 1
            max_bottleneck = max(max_bottleneck, length)
    
    # Combined "market isolation score"
    # High score = more isolated markets = potentially more unimodality issues
    market_isolation_score = (
        num_branches * 2.0 + 
        (avg_branch_dist / n) * 1.0 +
        (subtree_variance / n) * 0.5 +
        (max_bottleneck / n) * 0.5
    )
    
    return {
        'num_branches': num_branches,
        'avg_branch_dist': avg_branch_dist,
        'max_branch_dist': max_branch_dist,
        'core_size': core_size,
        'subtree_variance': subtree_variance,
        'subtree_max_min_ratio': subtree_max_min_ratio,
        'max_bottleneck': max_bottleneck,
        'market_isolation_score': market_isolation_score,
        'n': n
    }


def generate_multi_spider(n: int, branch_sizes: List[int], 
                          path_length: int = 0) -> 'nx.Graph':
    """
    Generate a "multi-spider" - multiple branches attached to a central path.
    This creates "isolated markets" that may not communicate well.
    
    Args:
        n: total nodes
        branch_sizes: list of branch lengths (from center to leaf)
        path_length: length of central spine (0 = star-like)
    """
    import networkx as nx
    
    G = nx.Graph()
    
    if not branch_sizes:
        # Just create a path if no branch sizes specified
        for i in range(n - 1):
            G.add_edge(i, i + 1)
        return G
    
    if path_length == 0:
        # Star-like: all branches from single center
        center = 0
        node_idx = 1
        for i, s in enumerate(branch_sizes):
            # Create path of length s from center
            prev = center
            for j in range(s):
                G.add_edge(prev, node_idx)
                prev = node_idx
                node_idx += 1
    else:
        # Path-based: branches from different points on a path
        spine = list(range(path_length))
        node_idx = path_length
        for i, s in enumerate(branch_sizes):
            # Attach branch to spine[i % len(spine)]
            attach_point = spine[i % len(spine)]
            prev = attach_point
            for j in range(s):
                G.add_edge(prev, node_idx)
                prev = node_idx
                node_idx += 1
    
    # Add remaining nodes to make up n
    while node_idx < n:
        # Attach to a random existing node
        target = random.choice(list(G.nodes()))
        G.add_edge(target, node_idx)
        node_idx += 1
    
    return G


def generate_bottleneck_tree(n: int, num_bottlenecks: int = 2) -> 'nx.Graph':
    """
    Generate a tree with clear bottlenecks - multiple large components
    connected through narrow paths.
    """
    import networkx as nx
    
    G = nx.Graph()
    
    # Create main spine
    spine_len = num_bottlenecks + 1
    spine = list(range(spine_len))
    for i in range(spine_len - 1):
        G.add_edge(spine[i], spine[i+1])
    
    node_idx = spine_len
    
    # Attach large subtrees at each spine endpoint
    remaining = n - spine_len
    per_end = remaining // (num_bottlenecks + 1)
    
    for i, s in enumerate(spine):
        # Attach a star or path
        if i % 2 == 0:
            # Star
            for _ in range(per_end):
                G.add_edge(s, node_idx)
                node_idx += 1
        else:
            # Path
            prev = s
            for _ in range(per_end):
                G.add_edge(prev, node_idx)
                prev = node_idx
                node_idx += 1
    
    return G


def generate_asymmetric_broom(n: int, p: int, q: int) -> 'nx.Graph':
    """
    Asymmetric broom: a path of length p, with a star of q leaves at one end.
    Creates market asymmetry.
    """
    import networkx as nx
    
    G = nx.Graph()
    
    # Main path of length p
    for i in range(p):
        G.add_edge(i, i+1)
    
    # Star at endpoint
    for i in range(q):
        G.add_edge(p, p + 1 + i)
    
    # Fill remaining with path extension
    node = p + q
    while node < n:
        G.add_edge(node, node + 1)
        node += 1
    
    return G


def generate_market_isolation_tree(n: int, num_markets: int = 3) -> 'nx.Graph':
    """
    Generate a tree with multiple "markets" (large subtrees) that are
    intentionally isolated from each other.
    """
    import networkx as nx
    
    G = nx.Graph()
    
    # Create num_markets independent clusters, connected by single edges
    # This is like having isolated fishing villages
    
    # Distribute remaining nodes
    nodes_per_market = (n - num_markets + 1) // num_markets
    nodes = [0]
    
    for m in range(num_markets):
        # Create a market (small tree)
        market_size = nodes_per_market
        if m == num_markets - 1:
            market_size = n - sum([nodes_per_market] * (m)) - 1
        
        # Start from current node
        for i in range(1, market_size + 1):
            if i == 1:
                # Connect to previous market
                G.add_edge(nodes[-1], nodes[-1] + 1)
            else:
                G.add_edge(nodes[-1] + i - 1, nodes[-1] + i)
        
        nodes.append(nodes[-1] + market_size)
    
    return G


def search_fishing_net(n_start: int = 10, n_end: int = 40,
                       samples_per_n: int = 100,
                       output_file: str = "results/fishing_net_search.json") -> Dict:
    """
    Main search using the fishing-net (market isolation) approach.
    """
    import networkx as nx
    from indpoly import independence_poly, is_unimodal, is_log_concave, near_miss_ratio, log_concavity_ratio
    
    results = {
        'approach': 'fishing_net_market_isolation',
        'description': 'Search for trees with market isolation patterns',
        'n_start': n_start,
        'n_end': n_end,
        'samples_per_n': samples_per_n,
        'total_trees_tested': 0,
        'non_unimodal_count': 0,
        'lc_failure_count': 0,
        'non_unimodal': [],
        'lc_failures': [],
        'near_misses': [],
        'market_isolation_analysis': defaultdict(list),
        'best_nm': 1.0,
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    
    for n in range(n_start, n_end + 1):
        print(f"Testing n={n}...")
        
        for sample in range(samples_per_n):
            # Generate tree with market isolation features
            tree_type = random.choice([
                'multi_spider',
                'bottleneck', 
                'asymmetric_broom',
                'market_isolation',
                'random_branch_heavy'
            ])
            
            if tree_type == 'multi_spider':
                # Multiple branches from center/path
                num_branches = random.randint(3, 6)
                branch_sizes = [random.randint(2, max(3, n//8)) for _ in range(num_branches)]
                path_length = random.randint(0, 3)
                G = generate_multi_spider(n, branch_sizes, path_length)
            
            elif tree_type == 'bottleneck':
                num_bottlenecks = random.randint(2, 4)
                G = generate_bottleneck_tree(n, num_bottlenecks)
            
            elif tree_type == 'asymmetric_broom':
                p = random.randint(3, n//2)
                q = random.randint(3, n//3)
                G = generate_asymmetric_broom(n, p, q)
            
            elif tree_type == 'market_isolation':
                num_markets = random.randint(3, 5)
                G = generate_market_isolation_tree(n, num_markets)
            
            else:
                # Random tree using networkx
                G = nx.random_labeled_tree(n, seed=random.randint(0, 1000000))
            
            # Ensure we have exactly n nodes and it's connected
            if len(G.nodes()) != n or not nx.is_connected(G):
                continue
            
            # Compute metrics
            metrics = compute_market_isolation_metrics(G)
            
            # Compute independence polynomial
            poly = independence_poly(len(G.nodes()), [list(G.neighbors(v)) for v in G.nodes()])
            results['total_trees_tested'] += 1
            
            # Check unimodality
            is_uni = is_unimodal(poly)
            is_lc = is_log_concave(poly)
            lc_ratio, lc_pos = log_concavity_ratio(poly)
            
            # Store market isolation analysis
            results['market_isolation_analysis'][n].append({
                'isolation_score': metrics['market_isolation_score'],
                'num_branches': metrics['num_branches'],
                'nm_ratio': near_miss_ratio(poly)[0] if poly else 1.0,
                'is_unimodal': is_uni,
                'is_log_concave': is_lc
            })
            
            if not is_uni:
                results['non_unimodal_count'] += 1
                results['non_unimodal'].append({
                    'n': n,
                    'poly': poly,
                    'metrics': metrics,
                    'tree_type': tree_type
                })
                print(f"  FOUND NON-UNIMODAL at n={n}! {tree_type}")
            
            if not is_lc:
                results['lc_failure_count'] += 1
                results['lc_failures'].append({
                    'n': n,
                    'poly': poly,
                    'lc_pos': lc_pos,
                    'lc_ratio': lc_ratio,
                    'metrics': metrics,
                    'tree_type': tree_type
                })
            
            # Track near misses
            nm = near_miss_ratio(poly)[0]
            if nm < results['best_nm']:
                results['best_nm'] = nm
                results['near_misses'].append({
                    'n': n,
                    'nm_ratio': nm,
                    'metrics': metrics,
                    'tree_type': tree_type
                })
                results['near_misses'].sort(key=lambda x: x['nm_ratio'])
                results['near_misses'] = results['near_misses'][:20]
        
        # Progress update
        print(f"  n={n}: tested {samples_per_n}, total {results['total_trees_tested']}, "
              f"LC failures: {results['lc_failure_count']}, best nm: {results['best_nm']:.4f}")
    
    # Analyze correlation between market isolation and unimodality
    print("\n=== Market Isolation Analysis ===")
    for n, analyses in results['market_isolation_analysis'].items():
        avg_isolation = sum(a['isolation_score'] for a in analyses) / len(analyses)
        avg_nm = sum(a['nm_ratio'] for a in analyses) / len(analyses)
        print(f"n={n}: avg isolation={avg_isolation:.3f}, avg nm={avg_nm:.4f}")
    
    # Save results
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {output_file}")
    print(f"Total trees tested: {results['total_trees_tested']}")
    print(f"Non-unimodal found: {results['non_unimodal_count']}")
    print(f"LC failures found: {results['lc_failure_count']}")
    print(f"Best near-miss ratio: {results['best_nm']:.6f}")
    
    return results


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Fishing-net search for unimodality')
    parser.add_argument('--n-start', type=int, default=10, help='Start n')
    parser.add_argument('--n-end', type=int, default=30, help='End n')
    parser.add_argument('--samples', type=int, default=100, help='Samples per n')
    parser.add_argument('--output', type=str, default='results/fishing_net_search.json')
    
    args = parser.parse_args()
    
    search_fishing_net(args.n_start, args.n_end, args.samples, args.output)
