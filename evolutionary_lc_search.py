import random
import time
import copy
import math
import networkx as nx
import argparse
import json
import os
import sys

# Import local modules
from graph6 import parse_graph6
from indpoly import independence_poly, log_concavity_ratio, is_unimodal

def adj_to_nx(n, adj):
    """Convert adjacency list to NetworkX graph."""
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for u, neighbors in enumerate(adj):
        for v in neighbors:
            if u < v:
                G.add_edge(u, v)
    return G

def nx_to_adj(G):
    """Convert NetworkX graph to adjacency list (relabeling nodes to 0..n-1)."""
    mapping = {node: i for i, node in enumerate(G.nodes())}
    G_relabeled = nx.relabel_nodes(G, mapping)
    n = G_relabeled.number_of_nodes()
    adj = [[] for _ in range(n)]
    for u, v in G_relabeled.edges():
        adj[u].append(v)
        adj[v].append(u)
    return n, adj

class Individual:
    def __init__(self, graph):
        self.graph = graph
        self.poly = None
        self.lc_ratio = None
        self.lc_pos = None
        self.is_unimodal = None
        self.n = graph.number_of_nodes()
        self.fitness = -1.0

    def evaluate(self):
        if self.poly is None:
            n, adj = nx_to_adj(self.graph)
            # Basic sanity check for connectivity
            if not nx.is_connected(self.graph):
                 # Disconnected graphs are not trees, but their independence poly 
                 # is product of components. The conjecture is for trees.
                 # We can just penalize disconnected graphs or extract the largest component.
                 # Let's extract largest component to be safe.
                 largest_cc = max(nx.connected_components(self.graph), key=len)
                 self.graph = self.graph.subgraph(largest_cc).copy()
                 # Re-convert
                 n, adj = nx_to_adj(self.graph)
                 self.n = n
            
            self.poly = independence_poly(n, adj)
            self.lc_ratio, self.lc_pos = log_concavity_ratio(self.poly)
            self.is_unimodal = is_unimodal(self.poly)
            
            # Fitness: Primary is LC ratio. 
            # Bonus for being non-unimodal (which is the ultimate goal).
            # We want to maximize this ratio.
            self.fitness = self.lc_ratio
            if not self.is_unimodal:
                self.fitness += 100.0  # Massive bonus found it!

    def clone(self):
        new_ind = Individual(self.graph.copy())
        # Don't copy poly/fitness to force re-evaluation if mutated, 
        # but if we just clone without mutation we could copy. 
        # For safety in GA, assume we might mutate copy.
        return new_ind

# --- Mutation Operators ---

def mutate_add_leaf(ind):
    """Attach a new leaf to a random existing node."""
    G = ind.graph.copy()
    if G.number_of_nodes() >= 60: # Limit size to prevent explosion
        return ind
    
    nodes = list(G.nodes())
    parent = random.choice(nodes)
    new_node = max(nodes) + 1 if nodes else 0
    G.add_edge(parent, new_node)
    return Individual(G)

def mutate_remove_leaf(ind):
    """Remove a random leaf."""
    G = ind.graph.copy()
    if G.number_of_nodes() <= 10: # Don't shrink too much
        return ind
    
    leaves = [x for x in G.nodes() if G.degree(x) == 1]
    if not leaves:
        return ind
    
    leaf = random.choice(leaves)
    G.remove_node(leaf)
    if not nx.is_connected(G):
        # Should not happen if we remove a leaf from a tree, 
        # but just in case of weird state
        largest_cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_cc).copy()
        
    return Individual(G)

def mutate_graft(ind):
    """Cut a branch and reattach it to a different node."""
    G = ind.graph.copy()
    n = G.number_of_nodes()
    if n < 4:
        return ind
    
    # 1. Pick a random edge to cut
    edges = list(G.edges())
    u, v = random.choice(edges)
    G.remove_edge(u, v)
    
    # 2. Identify the two components
    comps = list(nx.connected_components(G))
    if len(comps) != 2:
        # Should be 2 for a tree
        return ind
    
    comp1 = list(comps[0])
    comp2 = list(comps[1])
    
    # 3. Pick a new attachment point
    # We want to attach a node from comp1 to a node from comp2
    # But NOT the original edge (u,v).
    new_u = random.choice(comp1)
    new_v = random.choice(comp2)
    
    G.add_edge(new_u, new_v)
    return Individual(G)

def mutate_rewire_leaf(ind):
    """Move a leaf from one parent to another."""
    G = ind.graph.copy()
    leaves = [x for x in G.nodes() if G.degree(x) == 1]
    if not leaves:
        return ind
    
    leaf = random.choice(leaves)
    neighbor = list(G.neighbors(leaf))[0]
    
    G.remove_edge(leaf, neighbor)
    
    # Pick new parent
    possibilities = [n for n in G.nodes() if n != leaf and n != neighbor]
    if not possibilities:
        G.add_edge(leaf, neighbor) # Put it back
        return ind
        
    new_parent = random.choice(possibilities)
    G.add_edge(leaf, new_parent)
    return Individual(G)

def mutate_inject_broom_features(ind):
    """Tendency to make nodes high degree (broom hubs)."""
    G = ind.graph.copy()
    # Pick a high degree node and add a leaf to it, 
    # OR move a leaf to it
    degrees = sorted(G.degree, key=lambda x: x[1], reverse=True)
    hub = degrees[0][0] # Highest degree node
    
    if random.random() < 0.5:
        # Add leaf to hub
        new_node = max(G.nodes()) + 1
        G.add_edge(hub, new_node)
    else:
        # Move random leaf to hub
        leaves = [x for x in G.nodes() if G.degree(x) == 1 and list(G.neighbors(x))[0] != hub]
        if leaves:
            leaf = random.choice(leaves)
            old_parent = list(G.neighbors(leaf))[0]
            G.remove_edge(leaf, old_parent)
            G.add_edge(leaf, hub)
            
    return Individual(G)


# --- Crossover ---
def crossover_subtree_exchange(parent1, parent2):
    """Swap subtrees between two trees."""
    # This is tricky because we need to maintain validity.
    # Simplified approach: 
    # 1. Cut edge in P1 -> (Main1, Sub1)
    # 2. Cut edge in P2 -> (Main2, Sub2)
    # 3. Form Child1 = Main1 + Sub2
    # 4. Form Child2 = Main2 + Sub1
    
    def cut_random(G):
        if G.number_of_nodes() < 2:
            return None, None, None
        edges = list(G.edges())
        if not edges: return None, None, None
        u, v = random.choice(edges)
        G_cut = G.copy()
        G_cut.remove_edge(u, v)
        comps = list(nx.connected_components(G_cut))
        if len(comps) != 2: return None, None, None
        # Return the larger component as Main (usually root-side conceptually)
        # But for trees it's arbitrary. Let's just return both components sets.
        return comps[0], comps[1], (u, v)

    c1a, c1b, edge1 = cut_random(parent1.graph)
    c2a, c2b, edge2 = cut_random(parent2.graph)
    
    if not c1a or not c2a:
        return parent1.clone(), parent2.clone()
    
    # Relabel nodes to avoid collision
    # Create disjoint unions
    
    def stitch(nodes_main, nodes_sub):
        # Create a graph with these nodes
        # We need to preserve internal structure of Main and Sub
        # And add ONE bridge edge.
        
        # This requires extracting the subgraphs.
        # It's computationally expensive to do full isomorphism checks or relabeling safely.
        # Let's try a simpler strategy: Just graft a copy of Sub2 onto P1 (removing Sub1).
        pass

    # Actually, simpler crossover:
    # Child = Parent1 with a random subtree replaced by a random subtree from Parent2.
    # To implement this easily with NetworkX:
    # 1. Cut P1 at edge (u, v). Keep the larger part P1_keep.
    # 2. Cut P2 at edge (x, y). Keep the smaller part P2_sub.
    # 3. Disjoint union P1_keep and P2_sub (relabeling P2_sub).
    # 4. Add edge between u (in P1_keep) and y (in P2_sub).
    
    def extract_cut(G):
        if G.number_of_nodes() < 2: return None
        edges = list(G.edges())
        u, v = random.choice(edges)
        G_cut = G.copy()
        G_cut.remove_edge(u, v)
        comps = list(nx.connected_components(G_cut))
        if len(comps) != 2: return None
        
        # Heuristic: keep the part containing node 'u' as base, 'v' as sub
        base_nodes = comps[0] if u in comps[0] else comps[1]
        sub_nodes = comps[1] if u in comps[0] else comps[0]
        
        return G.subgraph(base_nodes).copy(), G.subgraph(sub_nodes).copy(), u, v
    
    res1 = extract_cut(parent1.graph)
    res2 = extract_cut(parent2.graph)
    
    if not res1 or not res2:
        return parent1.clone(), parent2.clone()
        
    base1, sub1, u1, v1 = res1
    base2, sub2, u2, v2 = res2
    
    # Child 1: Base1 + Sub2
    # We need to relabel Sub2 to not conflict with Base1
    max_id1 = max(base1.nodes())
    mapping2 = {n: n + max_id1 + 10 for n in sub2.nodes()}
    sub2_relabeled = nx.relabel_nodes(sub2, mapping2)
    
    C1 = nx.compose(base1, sub2_relabeled)
    # Add bridge edge. u1 is in base1. v2 mapped is in sub2_relabeled.
    C1.add_edge(u1, mapping2[v2])
    
    return Individual(C1), Individual(C1) # Just return one distinct child type for now

def generate_seeds(lc_failures):
    seeds = []
    for fail in lc_failures:
        if "graph6" in fail:
            n, adj = parse_graph6(fail["graph6"].encode('utf-8'))
            G = adj_to_nx(n, adj)
            ind = Individual(G)
            seeds.append(ind)
    return seeds

def main():
    parser = argparse.ArgumentParser(description="Evolutionary search for non-unimodal trees.")
    parser.add_argument("--generations", type=int, default=1000)
    parser.add_argument("--pop-size", type=int, default=100)
    parser.add_argument("--seeds-file", type=str, default="results/analysis_n26.json")
    args = parser.parse_args()

    # Load seeds
    with open(args.seeds_file, 'r') as f:
        data = json.load(f)
        
    failures = data.get("lc_failures", [])
    near_misses = data.get("top_near_misses", [])
    
    population = generate_seeds(failures)
    
    # Add near misses as well
    population.extend(generate_seeds(near_misses))
    
    # Fill remaining population with mutated clones of seeds
    initial_seeds = copy.deepcopy(population)
    while len(population) < args.pop_size:
        parent = random.choice(initial_seeds)
        # Apply multi-mutation to initial population to spread them out
        child = parent
        for _ in range(3):
            child = mutate_add_leaf(child)
        population.append(child)
        
    print(f"Initialized population: {len(population)}")
    
    # Evaluation loop
    best_ever_ratio = 0.0
    best_ever_ind = None
    
    for gen in range(args.generations):
        # Evaluate
        for ind in population:
            if ind.fitness < 0:
                ind.evaluate()
        
        # Sort
        population.sort(key=lambda x: x.fitness, reverse=True)
        
        best = population[0]
        if best.fitness > best_ever_ratio:
            best_ever_ratio = best.fitness
            best_ever_ind = best
            
        if gen % 50 == 0:
            diversity = len(set([ind.fitness for ind in population]))
            print(f"Gen {gen}: Best Ratio = {best.fitness:.6f} (n={best.n}, pos={best.lc_pos}) Unimodal={best.is_unimodal}, PopDiv={diversity}")
            if not best.is_unimodal:
                print("FOUND NON-UNIMODAL TREE!")
                print(f"Poly: {best.poly}")
                nx_to_adj(best.graph) # just to verify
                break

        # Selection (Elitism + Tournament)
        elite_count = int(args.pop_size * 0.1)
        new_pop = population[:elite_count]
        
        # Tournament selection pool
        pool = population[:int(args.pop_size * 0.5)] # Top 50%
        
        while len(new_pop) < args.pop_size:
            if random.random() < 0.2: # Crossover
                p1 = random.choice(pool)
                p2 = random.choice(pool)
                c1, c2 = crossover_subtree_exchange(p1, p2)
                new_pop.append(c1)
                if len(new_pop) < args.pop_size:
                    new_pop.append(c2)
            else: # Mutation
                parent = random.choice(pool)
                
                # Multi-mutation: apply 1 to 3 mutations
                child = parent.clone()
                num_mutations = random.randint(1, 3)
                for _ in range(num_mutations):
                    r = random.random()
                    if r < 0.4: # Increased add probability
                        child = mutate_add_leaf(child)
                    # elif r < 0.5: # REMOVED remove_leaf to force growth
                    #     child = mutate_remove_leaf(child)
                    elif r < 0.7:
                        child = mutate_graft(child)
                    elif r < 0.8:
                        child = mutate_rewire_leaf(child)
                    else:
                        child = mutate_inject_broom_features(child)
                
                new_pop.append(child)
        
        population = new_pop
        
    print("Search finished.")
    print(f"Best Ratio: {best_ever_ratio}")
    if best_ever_ind:
        n, adj = nx_to_adj(best_ever_ind.graph)
        print(f"Best Tree (n={n}): {adj}")
        
        # Save result
        result = {
            "n": n,
            "adj": adj,
            "poly": best_ever_ind.poly,
            "lc_ratio": best_ever_ind.lc_ratio,
            "lc_pos": best_ever_ind.lc_pos,
            "unimodal": best_ever_ind.is_unimodal
        }
        os.makedirs("results", exist_ok=True)
        with open("results/best_evolutionary_tree.json", "w") as f:
            json.dump(result, f, indent=2)

if __name__ == "__main__":
    main()
