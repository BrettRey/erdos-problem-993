#!/usr/bin/env python3
"""
Evolutionary search for unimodality counterexamples specifically within the set of
Homeomorphically Irreducible (HI) trees.

Algorithm:
1. Population of HI trees.
2. Mutate (add/remove/move leaves, etc).
3. Repair: Apply series-reduction to remove any degree-2 vertices.
4. Objective: Maximize near-miss ratio.
"""

import sys
import random
import time
import networkx as nx
import numpy as np
from indpoly import independence_poly, near_miss_ratio, is_unimodal

def random_hi_tree(n):
    """Generate random HI tree on n vertices."""
    # Start with random tree
    # Series reduce until HI
    # This might reduce size.
    # To get size n, we might need to start larger?
    # Or: Constructive approach (random_labeled_tree on kernel + leaves)
    # Re-use the logic from search_hi_trees.py
    
    if n < 4: return None
    max_k = (n - 2) // 2
    if max_k < 1: return None
    k = random.randint(1, max_k)
    T_I = nx.random_labeled_tree(k)
    degrees = dict(T_I.degree())
    required_leaves = {}
    total_req = 0
    for v in range(k):
        req = max(0, 3 - degrees[v])
        required_leaves[v] = req
        total_req += req
    if n - k < total_req: return None
    
    leaves_count = {v: required_leaves[v] for v in range(k)}
    remaining = (n - k) - total_req
    for _ in range(remaining):
        v = random.randint(0, k-1)
        leaves_count[v] += 1
        
    G = nx.Graph()
    G.add_edges_from(T_I.edges())
    leaf_idx = k
    for v in range(k):
        for _ in range(leaves_count[v]):
            G.add_edge(v, leaf_idx)
            leaf_idx += 1
            
    # Relabel to 0..n-1
    return nx.convert_node_labels_to_integers(G)

def series_reduce(G):
    """Remove degree-2 vertices by merging neighbors."""
    # Repeat until no degree-2 vertices
    while True:
        deg2 = [n for n, d in G.degree() if d == 2]
        if not deg2:
            break
        
        # Pick one
        v = deg2[0]
        nbrs = list(G.neighbors(v))
        if len(nbrs) != 2:
            # Should not happen for degree 2
            break
            
        u, w = nbrs
        G.remove_node(v)
        G.add_edge(u, w)
        
    return G

def mutate(G):
    """Apply a mutation and return a new HI tree."""
    # Mutation strategies:
    # 1. Leaf move: move a leaf to another vertex.
    # 2. Leaf add/remove: change size? No, keep size n fixed?
    # The user wants to search HI trees.
    # If we change size, we drift. Let's try to keep size fixed.
    
    G = G.copy()
    n = G.number_of_nodes()
    
    # Strategy 1: Move a leaf
    leaves = [n for n, d in G.degree() if d == 1]
    if leaves:
        l = random.choice(leaves)
        # Neighbor
        nbr = list(G.neighbors(l))[0]
        # Remove edge
        G.remove_edge(l, nbr)
        # Add edge to random other vertex (not l)
        options = [v for v in G.nodes() if v != l]
        if options:
            new_nbr = random.choice(options)
            G.add_edge(l, new_nbr)
            
    # Strategy 2: SPR (Subtree Prune Regraft) - simple version
    # Pick edge, remove, reconnect components
    # ...
    
    # Repair
    series_reduce(G)
    
    # If size changed (reduced), we need to add vertices back?
    # Series reduction reduces size.
    # We want to explore HI trees of size N.
    # So if reduced, we add leaves back?
    while G.number_of_nodes() < n:
        # Add leaf to random node
        nodes = list(G.nodes())
        target = random.choice(nodes)
        new_node = max(G.nodes()) + 1 # simplistic ID
        # Find a free ID
        while new_node in G: new_node += 1
        G.add_edge(target, new_node)
        
    return nx.convert_node_labels_to_integers(G)

def fitness(G):
    n = G.number_of_nodes()
    adj = [[] for _ in range(n)]
    for u, v in G.edges():
        adj[u].append(v)
        adj[v].append(u)
        
    poly = independence_poly(n, adj)
    
    if not is_unimodal(poly):
        return 1000.0, poly # HUGE fitness
        
    nm, _ = near_miss_ratio(poly)
    return nm, poly

def main():
    print("Evolutionary Search for HI Trees...")
    target_size = 50
    pop_size = 50
    generations = 50
    
    # Init population
    population = []
    for _ in range(pop_size):
        G = None
        while G is None:
            G = random_hi_tree(target_size)
        fit, _ = fitness(G)
        population.append((fit, G))
        
    population.sort(key=lambda x: x[0], reverse=True)
    print(f"Gen 0: Best nm = {population[0][0]:.5f}")
    
    for gen in range(generations):
        # Elite
        new_pop = population[:5]
        
        # Breed
        while len(new_pop) < pop_size:
            parent = random.choice(population[:20])[1]
            child = mutate(parent)
            fit, poly = fitness(child)
            
            if fit > 1.0:
                print(f"FOUND COUNTEREXAMPLE! Gen {gen}")
                print(nx.to_dict_of_lists(child))
                print(poly)
                sys.exit(0)
                
            new_pop.append((fit, child))
            
        population = new_pop
        population.sort(key=lambda x: x[0], reverse=True)
        print(f"Gen {gen+1}: Best nm = {population[0][0]:.5f}")
        
    print("Done.")

if __name__ == "__main__":
    main()
