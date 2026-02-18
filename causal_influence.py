#!/usr/bin/env python3
"""
Causal Influence Analysis for Independence Polynomials

This module analyzes which vertices in a tree can influence which coefficients
in the independence polynomial. The intuition is that if we understand how
information flows during the DP computation, we might understand what causes
unimodality violations.

Key questions:
1. For coefficient i_k, which vertices can contribute to it?
2. What's the "influence distance" - how far must information travel?
3. Do LC failures have distinctive influence patterns?
"""

import json
import random
import sys
from collections import defaultdict
from typing import Dict, List, Set, Tuple

sys.path.insert(0, '.')

from indpoly import independence_poly, is_log_concave, log_concavity_ratio


def compute_influence_cones(n: int, adj: List[List[int]]) -> Dict[int, Set[int]]:
    """
    For each coefficient position k, compute which vertices can influence it.
    
    A vertex v can influence coefficient k if there exists an independent set
    of size k that includes v. Equivalently, in the DP, v's contribution
    flows to position k in the polynomial.
    
    We trace this through the DP computation:
    - dp1[v] includes v itself (adds 1 to the size)
    - dp0[v] excludes v (no size increase)
    - When combining children, we add polynomials
    """
    # Build rooted tree structure
    parent = [-1] * n
    children = [[] for _ in range(n)]
    
    # BFS to root the tree at 0
    queue = [0]
    parent[0] = 0
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        for u in adj[v]:
            if parent[u] == -1:
                parent[u] = v
                children[v].append(u)
                queue.append(u)
    
    # Compute dp polynomials with influence tracking
    # dp0[v][k] = coefficient k when v is excluded
    # dp1[v][k] = coefficient k when v is included (contributes +1 to set size)
    dp0 = [None] * n
    dp1 = [None] * n
    
    # Post-order
    order = queue[::-1]  # reversed BFS order = post-order
    
    for v in order:
        if not children[v]:
            # Leaf
            dp0[v] = {0: 1}  # empty set
            dp1[v] = {1: 1}  # set containing only v
        else:
            # Combine children
            # dp0[v] = product over children of (dp0[c] + dp1[c])
            # dp1[v] = x * product over children of dp0[c]
            
            prod = {0: 1}
            for c in children[v]:
                combined = dict(dp0[c])
                for k, v1 in dp1[c].items():
                    if k in combined:
                        combined[k] += v1
                    else:
                        combined[k] = v1
                # Convolution
                new_prod = defaultdict(int)
                for k1, v1 in prod.items():
                    for k2, v2 in combined.items():
                        new_prod[k1 + k2] += v1 * v2
                prod = dict(new_prod)
            dp0[v] = prod
            
            # dp1[v]: exclude all children, multiply by x
            prod = {0: 1}
            for c in children[v]:
                new_prod = defaultdict(int)
                for k1, v1 in prod.items():
                    for k2, v2 in dp0[c].items():
                        new_prod[k1 + k2] += v1 * v2
                prod = dict(new_prod)
            # Shift by 1 (include v)
            dp1[v] = {k+1: v for k, v in prod.items()}
    
    # Now compute influence cones: for each vertex, which coefficients can it affect?
    # A vertex v contributes to coefficient k if:
    # - It appears in dp1 of some ancestor (including itself)
    # - And the polynomial combination reaches position k
    
    # Simpler approach: trace contributions
    # For each vertex, track which coefficient positions it can reach
    influence_cones = {v: set() for v in range(n)}
    
    # Process in post-order: for each vertex, determine what it contributes to
    for v in order:
        # Vertex v itself can influence position 1 (when included)
        influence_cones[v].add(1)
        
        # When v is excluded, its children contribute as normal
        # When v is included, its children must be excluded (dp0 only)
        
        # Let's also track: which positions can be reached from v's subtree
        # when considering different inclusion/exclusion patterns
    
    # More detailed: trace actual contributions through the DP
    # This is expensive but more accurate
    
    def trace_influence(v: int, included: bool, offset: int) -> Set[int]:
        """Trace which coefficients are affected by decisions at vertex v."""
        positions = set()
        
        if included:
            # v is in the independent set
            positions.add(offset + 1)  # v itself contributes 1
            # All children must be excluded
            for c in children[v]:
                child_positions = trace_influence(c, False, offset + 1)
                positions.update(child_positions)
        else:
            # v is not in the independent set
            # Children can be either included or excluded
            for c in children[v]:
                # Either include c
                inc_positions = trace_influence(c, True, offset)
                positions.update(inc_positions)
                # Or exclude c
                exc_positions = trace_influence(c, False, offset)
                positions.update(exc_positions)
        
        return positions
    
    # For each vertex, compute full influence cone
    for v in range(n):
        # Including v
        influence_cones[v].update(trace_influence(v, True, 0))
        # Excluding v
        influence_cones[v].update(trace_influence(v, False, 0))
    
    return influence_cones


def compute_influence_depth(n: int, adj: List[List[int]]) -> Dict[int, List[int]]:
    """
    For each coefficient position k, what's the minimum distance
    information must travel to affect it?
    """
    # BFS from each vertex to compute distances
    from collections import deque
    
    all_distances = []
    for start in range(n):
        dist = [-1] * n
        dist[start] = 0
        q = deque([start])
        while q:
            v = q.popleft()
            for u in adj[v]:
                if dist[u] == -1:
                    dist[u] = dist[v] + 1
                    q.append(u)
        all_distances.append(dist)
    
    # For each coefficient position, compute max/avg/min influence distance
    influence_cones = compute_influence_cones(n, adj)
    
    # For each k, find how far information must travel
    depth_by_k = defaultdict(list)
    for v in range(n):
        for k in influence_cones[v]:
            # Minimum distance from any leaf to affect position k?
            # This is approximate
            min_dist = min(all_distances[v][u] for u in range(n) if adj[u])
            depth_by_k[k].append(min_dist)
    
    result = {}
    for k, depths in depth_by_k.items():
        result[k] = {
            'min': min(depths) if depths else 0,
            'max': max(depths) if depths else 0,
            'avg': sum(depths) / len(depths) if depths else 0
        }
    
    return result


def analyze_tree_influence(n: int, adj: List[List[int]], label: str = "") -> Dict:
    """Analyze influence patterns for a single tree."""
    poly = independence_poly(n, adj)
    is_lc = is_log_concave(poly)
    lc_viol, lc_viol_pos = log_concavity_ratio(poly)
    
    # Compute influence cones
    influence_cones = compute_influence_cones(n, adj)
    
    # Statistics
    cone_sizes = [len(cones) for cones in influence_cones.values()]
    
    # For each coefficient, how many vertices can influence it?
    k_influencers = defaultdict(int)
    for v, cones in influence_cones.items():
        for k in cones:
            k_influencers[k] += 1
    
    # Correlation: does being influenced by many vertices correlate with LC issues?
    # If k has many influencers and there's a violation at k, that's interesting
    
    result = {
        'label': label,
        'n': n,
        'poly': poly,
        'log_concave': is_lc,
        'lc_violation': lc_viol,
        'lc_viol_pos': lc_viol_pos,
        'max_cone_size': max(cone_sizes) if cone_sizes else 0,
        'avg_cone_size': sum(cone_sizes) / len(cone_sizes) if cone_sizes else 0,
        'k_influencers': dict(k_influencers),
    }
    
    return result


def analyze_lc_failure_patterns():
    """Analyze influence patterns for known LC failures."""
    # Load the n=26 LC failures
    with open('results/analysis_n26.json', 'r') as f:
        data = json.load(f)
    
    failures = data.get('lc_failures', [])
    print(f"Analyzing {len(failures)} known LC failures...")
    
    results = []
    for failure in failures:
        # The graph6 strings are in the results but we need to parse them
        # For now, just report what we know
        results.append({
            'n': failure['n'],
            'lc_ratio': failure['lc_ratio'],
            'lc_pos': failure['lc_pos'],
            'poly': failure['poly'][:10],  # First 10 coefficients
        })
    
    return results


def run_small_survey(min_n: int = 5, max_n: int = 15, samples: int = 100):
    """Run a small survey of random trees to find influence patterns."""
    import networkx as nx
    
    results = []
    rng = random.Random(42)
    
    for n in range(min_n, max_n + 1):
        for _ in range(samples):
            G = nx.random_labeled_tree(n, seed=rng.randint(0, 100000))
            adj = [list(G.neighbors(v)) for v in range(n)]
            
            result = analyze_tree_influence(n, adj, f"random_{n}")
            results.append(result)
    
    # Find patterns
    lc_failures = [r for r in results if not r['log_concave']]
    lc_pass = [r for r in results if r['log_concave']]
    
    print(f"\nSurvey Results (n={min_n}-{max_n}, {samples} each):")
    print(f"Total trees: {len(results)}")
    print(f"LC failures: {len(lc_failures)}")
    print(f"LC passes: {len(lc_pass)}")
    
    if lc_failures:
        print(f"\nLC Failure patterns:")
        for r in lc_failures[:5]:
            print(f"  n={r['n']}, lc_pos={r['lc_pos']}, max_cone={r['max_cone_size']}")
    
    return results


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Causal influence analysis')
    parser.add_argument('--analyze-lc', action='store_true', help='Analyze known LC failures')
    parser.add_argument('--survey', action='store_true', help='Run small survey')
    parser.add_argument('--n', type=int, default=10, help='Test specific n')
    
    args = parser.parse_args()
    
    if args.analyze_lc:
        analyze_lc_failure_patterns()
    elif args.survey:
        run_small_survey()
    else:
        # Quick test on a small tree
        import networkx as nx
        G = nx.balanced_tree(2, 3)  # Binary tree depth 3
        n = len(G.nodes())
        adj = [list(G.neighbors(v)) for v in range(n)]
        result = analyze_tree_influence(n, adj, "balanced_tree")
        print(json.dumps(result, indent=2))
