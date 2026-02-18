"""
Landslide Evolution: Evolve the 'debris' tree to trigger a landslide on a fixed base.
"""

import argparse
import random
import sys
import copy
import networkx as nx
from typing import Any

# Import from existing modules
try:
    from indpoly import independence_poly, is_unimodal
except ImportError:
    sys.path.append(".")
    from indpoly import independence_poly, is_unimodal

def get_poly(adj: list[list[int]]) -> list[int]:
    n = len(adj)
    return independence_poly(n, adj)

def adj_to_nx(adj: list[list[int]]) -> nx.Graph:
    n = len(adj)
    g = nx.Graph()
    g.add_nodes_from(range(n))
    for u, neigh in enumerate(adj):
        for v in neigh:
            if u < v:
                g.add_edge(u, v)
    return g

def nx_to_adj(g: nx.Graph) -> list[list[int]]:
    mapping = {node: i for i, node in enumerate(g.nodes())}
    h = nx.relabel_nodes(g, mapping)
    n = h.number_of_nodes()
    adj = [[] for _ in range(n)]
    for u, v in h.edges():
        adj[u].append(v)
        adj[v].append(u)
    for u in range(n):
        adj[u].sort()
    return adj

def graft_nx(base: nx.Graph, base_node: int, debris: nx.Graph, debris_node: int) -> nx.Graph:
    """Graft debris onto base using networkx disjoint union."""
    # Relabel debris to avoid collisions
    base_n = base.number_of_nodes()
    debris_n = debris.number_of_nodes()
    
    mapping = {n: n + base_n for n in debris.nodes()}
    debris_relabeled = nx.relabel_nodes(debris, mapping)
    debris_root = debris_node + base_n
    
    U = nx.compose(base, debris_relabeled)
    U.add_edge(base_node, debris_root)
    return U

class Individual:
    def __init__(self, graph: nx.Graph, target_node: int = 0):
        self.graph = graph
        self.target_node = target_node
        self.fitness = -1.0
        self.poly = []
        self.ratio = 0.0
        
    def clone(self):
        ind = Individual(self.graph.copy(), self.target_node)
        ind.fitness = self.fitness
        ind.poly = list(self.poly)
        ind.ratio = self.ratio
        return ind

def evaluate(ind: Individual, base: nx.Graph):
    # Graft ind.graph (debris) onto base at ind.target_node
    debris_root = 0
    if debris_root not in ind.graph:
        debris_root = list(ind.graph.nodes())[0]
        
    # Ensure target_node is valid for base
    if ind.target_node >= base.number_of_nodes():
        ind.target_node = 0
        
    full_g = graft_nx(base, ind.target_node, ind.graph, debris_root)
    
    adj = nx_to_adj(full_g)
    poly = get_poly(adj)
    ind.poly = poly
    
    # Check unimodality and ratio
    if not is_unimodal(poly):
        ind.fitness = 1000.0 # SUCCESS!
        ind.ratio = 999.0
        return

    # Calculate worst ratio in tail
    peak_idx = poly.index(max(poly))
    worst = 0.0
    for k in range(peak_idx, len(poly) - 1):
        if poly[k] == 0: continue
        r = poly[k+1] / poly[k]
        if r > worst:
            worst = r
    
    ind.ratio = worst
    ind.fitness = worst # Maximize ratio

def mutate(ind: Individual, rng: random.Random, max_nodes: int, base_nodes: int):
    # Mutate graph
    if rng.random() < 0.8:
        g = ind.graph.copy()
        op = rng.choice(["add", "remove", "rewire"])
        
        if op == "add" and g.number_of_nodes() < max_nodes:
            # Add leaf
            parent = rng.choice(list(g.nodes()))
            new_node = max(g.nodes()) + 1 if g.number_of_nodes() else 0
            g.add_edge(parent, new_node)
            
        elif op == "remove" and g.number_of_nodes() > 2:
            # Remove leaf
            leaves = [v for v in g.nodes() if g.degree(v) == 1]
            if leaves:
                leaf = rng.choice(leaves)
                g.remove_node(leaf)
                    
        elif op == "rewire" and g.number_of_nodes() > 3:
            leaves = [v for v in g.nodes() if g.degree(v) == 1]
            if leaves:
                leaf = rng.choice(leaves)
                neighbors = list(g.neighbors(leaf))
                if neighbors:
                    old_parent = neighbors[0]
                    others = [n for n in g.nodes() if n != leaf and n != old_parent]
                    if others:
                        new_parent = rng.choice(others)
                        g.remove_edge(leaf, old_parent)
                        g.add_edge(leaf, new_parent)
        
        ind.graph = g
    
    # Mutate target node
    if rng.random() < 0.3:
        ind.target_node = rng.randint(0, base_nodes - 1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-n", type=int, default=60)
    parser.add_argument("--debris-n", type=int, default=12)
    parser.add_argument("--pop-size", type=int, default=100)
    parser.add_argument("--gens", type=int, default=200)
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()
    
    rng = random.Random(args.seed)
    
    # 1. Setup Base (Random Tree)
    # Using networkx random_labeled_tree
    base = nx.random_labeled_tree(args.base_n, seed=args.seed)
    print(f"Base: Random Tree(N={args.base_n}), Seed={args.seed}")
    
    # 2. Init Population (Debris + Target)
    pop = []
    for _ in range(args.pop_size):
        # Random tree of size ~debris_n
        size = rng.randint(4, args.debris_n)
        T = nx.random_labeled_tree(size)
        target = rng.randint(0, args.base_n - 1)
        pop.append(Individual(T, target))
        
    best_ever = 0.0
    
    for gen in range(args.gens):
        # Evaluate
        for ind in pop:
            evaluate(ind, base)
            
        pop.sort(key=lambda x: x.fitness, reverse=True)
        
        if gen % 10 == 0:
            print(f"Gen {gen}: Current Best {pop[0].fitness:.5f}")

        if pop[0].fitness > best_ever:
            best_ever = pop[0].fitness
            print(f"*** New Record at Gen {gen}: {best_ever:.5f} (Debris N={pop[0].graph.number_of_nodes()}, Target={pop[0].target_node}) ***")
            if best_ever >= 999.0:
                print("VIOLATION FOUND!")
                print(f"Poly: {pop[0].poly}")
                return
                
        # Selection / Breeding
        next_pop = []
        elite = pop[:10]
        next_pop.extend([x.clone() for x in elite])
        
        while len(next_pop) < args.pop_size:
            parent = rng.choice(elite) # Elitism selection
            child = parent.clone()
            mutate(child, rng, args.debris_n, args.base_n)
            next_pop.append(child)
            
        pop = next_pop

if __name__ == "__main__":
    main()
