#!/usr/bin/env python3
"""Adversarial ROOTS optimizer for Erdős Problem #993.

Optimizes for "High Frequency, Low Damping" roots in the independence polynomial.
Fitness = max(Imag(r)) / |max(Real(r))|

Hypothesis: Counterexamples (unimodality failures) live where roots are
highly oscillatory (high Imag) and close to the imaginary axis (low Damping).
"""

import argparse
import json
import os
import random
import time
from copy import deepcopy
import numpy as np

from indpoly import independence_poly, is_unimodal, near_miss_ratio
from targeted import make_broom, make_caterpillar, make_spider

# ── Tree utilities ──────────────────────────────────────────────────

def _validate_tree(n: int, adj: list[list[int]]) -> bool:
    """Check that adj represents a valid tree on n vertices."""
    if len(adj) != n:
        return False
    edge_count = sum(len(adj[v]) for v in range(n)) // 2
    if edge_count != n - 1:
        return False
    visited = [False] * n
    queue = [0]
    visited[0] = True
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                queue.append(u)
    return all(visited)

def _leaves(adj: list[list[int]]) -> list[int]:
    return [v for v in range(len(adj)) if len(adj[v]) == 1]

def _degree_sequence(adj: list[list[int]]) -> list[int]:
    return sorted((len(adj[v]) for v in range(len(adj))), reverse=True)

def _tree_signature(adj: list[list[int]]) -> str:
    n = len(adj)
    degs = sorted(len(adj[v]) for v in range(n))
    leaves = sum(1 for d in degs if d == 1)
    max_deg = degs[-1] if degs else 0
    branch_verts = sum(1 for d in degs if d >= 3)
    return f"n{n}_l{leaves}_d{max_deg}_b{branch_verts}"

# ── Mutations (Same as nm_optimizer) ────────────────────────────────

def mutate_leaf_move(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    adj = deepcopy(adj)
    lvs = _leaves(adj)
    if not lvs: return n, adj
    leaf = rng.choice(lvs)
    neighbor = adj[leaf][0]
    adj[neighbor].remove(leaf)
    adj[leaf].remove(neighbor)
    candidates = [v for v in range(n) if v != leaf and v != neighbor]
    if not candidates:
        adj[neighbor].append(leaf)
        adj[leaf].append(neighbor)
        return n, adj
    target = rng.choice(candidates)
    adj[target].append(leaf)
    adj[leaf].append(target)
    return n, adj

def mutate_spr(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    if n <= 3: return mutate_leaf_move(n, adj, rng)
    adj = deepcopy(adj)
    parent = [-1] * n
    order = [0]
    visited = [False] * n
    visited[0] = True
    head = 0
    while head < len(order):
        v = order[head]
        head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                order.append(u)
    pruned_root = rng.choice(order[1:])
    pruned_parent = parent[pruned_root]
    subtree = set()
    stack = [pruned_root]
    while stack:
        v = stack.pop()
        subtree.add(v)
        for u in adj[v]:
            if u != pruned_parent and u not in subtree:
                stack.append(u)
    adj[pruned_parent].remove(pruned_root)
    adj[pruned_root].remove(pruned_parent)
    remaining = [v for v in range(n) if v not in subtree]
    if not remaining:
        adj[pruned_parent].append(pruned_root)
        adj[pruned_root].append(pruned_parent)
        return n, adj
    target = rng.choice(remaining)
    adj[target].append(pruned_root)
    adj[pruned_root].append(target)
    return n, adj

def mutate_pendant_concentrate(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    adj = deepcopy(adj)
    lvs = _leaves(adj)
    if len(lvs) < 2: return n, adj
    max_deg_v = max(range(n), key=lambda v: len(adj[v]))
    candidates = [lf for lf in lvs if adj[lf][0] != max_deg_v]
    if not candidates: return n, adj
    leaf = rng.choice(candidates)
    old_neighbor = adj[leaf][0]
    adj[old_neighbor].remove(leaf)
    adj[leaf].remove(old_neighbor)
    if not adj[old_neighbor]:
        adj[max_deg_v].append(old_neighbor)
        adj[old_neighbor].append(max_deg_v)
        adj[max_deg_v].append(leaf)
        adj[leaf].append(max_deg_v)
    else:
        adj[max_deg_v].append(leaf)
        adj[leaf].append(max_deg_v)
    return n, adj

def mutate_edge_subdivide_contract(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    if n <= 3: return mutate_leaf_move(n, adj, rng)
    adj = deepcopy(adj)
    deg2 = [v for v in range(n) if len(adj[v]) == 2]
    if not deg2: return mutate_leaf_move(n, adj, rng)
    contract_v = rng.choice(deg2)
    u, w = adj[contract_v]
    adj[u].remove(contract_v)
    adj[w].remove(contract_v)
    adj[contract_v] = []
    if w not in adj[u]:
        adj[u].append(w)
        adj[w].append(u)
    edges = []
    for v in range(n):
        for nb in adj[v]:
            if v < nb:
                edges.append((v, nb))
    if not edges: return n, adj
    sub_u, sub_w = rng.choice(edges)
    adj[sub_u].remove(sub_w)
    adj[sub_w].remove(sub_u)
    adj[sub_u].append(contract_v)
    adj[contract_v].append(sub_u)
    adj[sub_w].append(contract_v)
    adj[contract_v].append(sub_w)
    return n, adj

MUTATIONS = [
    (mutate_leaf_move, 0.35),
    (mutate_spr, 0.30),
    (mutate_pendant_concentrate, 0.20),
    (mutate_edge_subdivide_contract, 0.15),
]

def mutate(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    r = rng.random()
    cumulative = 0.0
    for func, weight in MUTATIONS:
        cumulative += weight
        if r < cumulative:
            return func(n, adj, rng)
    return MUTATIONS[-1][0](n, adj, rng)

# ── Seed generation (Same as nm_optimizer) ──────────────────────────

def generate_seeds(n: int, count: int, rng: random.Random) -> list[tuple[int, list[list[int]], str]]:
    seeds = []
    for p in range(2, min(n - 1, 50)):
        s = n - p
        if s < 1: continue
        n_t, adj = make_broom(p, s)
        seeds.append((n_t, adj, f"broom({p},{s})"))
        if len(seeds) >= count // 2: break
    for k in [1, 2, 3, 5, 8]:
        spine = n // (k + 1)
        if spine < 2: continue
        pendants = [k] * spine
        total = spine + sum(pendants)
        if total > n: continue
        pendants[0] += n - total
        n_t, adj = make_caterpillar(spine, pendants)
        seeds.append((n_t, adj, f"cat_conc({spine},k={k})"))
    for num_legs in [3, 5, 8, 12]:
        leg_len = max(1, (n - 1) // num_legs)
        actual_n = 1 + num_legs * leg_len
        if actual_n <= n:
            legs = [leg_len] * num_legs
            n_t, adj = make_spider(legs)
            seeds.append((n_t, adj, f"spider({num_legs}x{leg_len})"))
    while len(seeds) < count:
        prufer = [rng.randrange(n) for _ in range(n - 2)]
        adj_r = [[] for _ in range(n)]
        degree = [1] * n
        for v in prufer: degree[v] += 1
        for v in prufer:
            for u in range(n):
                if degree[u] == 1:
                    adj_r[u].append(v)
                    adj_r[v].append(u)
                    degree[u] -= 1
                    degree[v] -= 1
                    break
        last = [u for u in range(n) if degree[u] == 1]
        if len(last) == 2:
            adj_r[last[0]].append(last[1])
            adj_r[last[1]].append(last[0])
        seeds.append((n, adj_r, f"random_prufer_{len(seeds)}"))
    return seeds[:count]

# ── Evaluation (Modified for Roots) ─────────────────────────────────

def evaluate(n: int, adj: list[list[int]]) -> dict:
    """Compute root statistics and near-miss ratio."""
    poly = independence_poly(n, adj)
    uni = is_unimodal(poly)
    nm, nm_pos = near_miss_ratio(poly)
    
    # Root analysis
    # numpy roots expects highest power first
    try:
        roots = np.roots(poly[::-1])
        real_parts = roots.real
        imag_parts = roots.imag
        
        # Max Real part (usually negative, so closest to 0 from left)
        max_re = np.max(real_parts)
        # Max Imag part (oscillation freq)
        max_im = np.max(np.abs(imag_parts))
        
        # Fitness: Maximize Near Miss Ratio (for this run)
        fitness = nm
    except Exception:
        # Fallback
        max_re = -100.0
        max_im = 0.0
        fitness = 0.0

    return {
        "nm": nm,
        "nm_pos": nm_pos,
        "unimodal": uni,
        "alpha": len(poly) - 1,
        "poly": poly,
        "max_re": max_re,
        "max_im": max_im,
        "fitness": fitness
    }

# ── Main optimizer ──────────────────────────────────────────────────

def run_optimizer(n: int, pop_size: int = 60, generations: int = 500, elite_frac: float = 0.2, mutations_per_ind: int = 3, seed: int = 42, verbose: bool = True) -> dict:
    rng = random.Random(seed)
    elite_count = max(2, int(pop_size * elite_frac))

    seeds = generate_seeds(n, pop_size, rng)
    population = []
    
    for n_t, adj, label in seeds:
        if n_t != n: continue
        ev = evaluate(n_t, adj)
        population.append({
            "n": n_t, "adj": adj, "label": label,
            "nm": ev["nm"], "unimodal": ev["unimodal"],
            "max_re": ev["max_re"], "max_im": ev["max_im"],
            "fitness": ev["fitness"],
            "poly": ev["poly"],
            "generation": 0,
        })
        if not ev["unimodal"]:
            return _counterexample_result(n_t, adj, label, ev, 0)

    # Sort by NM (Near Miss Ratio)
    population.sort(key=lambda x: -x["nm"])

    best_ever = population[0]["nm"]
    history = []

    if verbose:
        print(f"Gen   0: Best NM={population[0]['nm']:.4f} (Im={population[0]['max_im']:.2f}, Re={population[0]['max_re']:.4f})")

    for gen in range(1, generations + 1):
        elite = population[:elite_count]
        offspring = []
        for parent in elite:
            for _ in range(mutations_per_ind):
                child_n, child_adj = parent["n"], deepcopy(parent["adj"])
                num_muts = rng.randint(1, 3)
                for _ in range(num_muts):
                    child_n, child_adj = mutate(child_n, child_adj, rng)
                if not _validate_tree(child_n, child_adj): continue
                
                ev = evaluate(child_n, child_adj)
                child_label = f"mut_g{gen}"
                offspring.append({
                    "n": child_n, "adj": child_adj, "label": child_label,
                    "nm": ev["nm"], "unimodal": ev["unimodal"],
                    "max_re": ev["max_re"], "max_im": ev["max_im"],
                    "fitness": ev["fitness"],
                    "poly": ev["poly"],
                    "generation": gen,
                })
                if not ev["unimodal"]:
                    return _counterexample_result(child_n, child_adj, child_label, ev, gen)
        
        # Inject random
        for _ in range(max(2, pop_size // 10)):
            # Same random injection logic...
            prufer = [rng.randrange(n) for _ in range(n - 2)]
            adj_r = [[] for _ in range(n)]
            degree = [1] * n
            for v in prufer: degree[v] += 1
            for v in prufer:
                for u in range(n):
                    if degree[u] == 1:
                        adj_r[u].append(v)
                        adj_r[v].append(u)
                        degree[u] -= 1
                        degree[v] -= 1
                        break
            last = [u for u in range(n) if degree[u] == 1]
            if len(last) == 2:
                adj_r[last[0]].append(last[1])
                adj_r[last[1]].append(last[0])
            ev = evaluate(n, adj_r)
            offspring.append({
                "n": n, "adj": adj_r, "label": f"inject",
                "nm": ev["nm"], "unimodal": ev["unimodal"],
                "max_re": ev["max_re"], "max_im": ev["max_im"],
                "fitness": ev["fitness"],
                "poly": ev["poly"],
                "generation": gen,
            })

        combined = elite + offspring
        combined.sort(key=lambda x: -x["nm"])
        population = combined[:pop_size]

        if population[0]["nm"] > best_ever:
            best_ever = population[0]["nm"]

        if gen % 25 == 0 or gen == 1:
            best = population[0]
            if verbose:
                print(f"Gen {gen:>3}: Best NM={best['nm']:.4f} (Im={best['max_im']:.2f}, Re={best['max_re']:.4f})")

    best = population[0]
    return {
        "counterexample_found": False,
        "best_fitness": best["fitness"],
        "best_structure": _tree_signature(best["adj"]),
        "best_max_im": best["max_im"],
        "best_max_re": best["max_re"],
        "best_nm": best["nm"],
        "best_adj": best["adj"], # Return adjacency
        "best_poly": best["poly"]
    }

def _counterexample_result(n, adj, label, ev, gen):
    return {
        "counterexample_found": True,
        "n": n, "label": label, "generation": gen,
        "poly": ev["poly"], "nm": ev["nm"],
        "max_im": ev["max_im"], "max_re": ev["max_re"],
        "best_adj": adj
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=26)
    parser.add_argument("--generations", type=int, default=100)
    args = parser.parse_args()

    print(f"Starting Roots Optimizer for n={args.n}...")
    result = run_optimizer(args.n, generations=args.generations)
    
    if result["counterexample_found"]:
        print("\n*** COUNTEREXAMPLE FOUND! ***")
        print(f"Max Im: {result['max_im']}, Max Re: {result['max_re']}")
    else:
        print("\nOptimization finished.")
        print(f"Best Fitness (NM): {result['best_fitness']:.4f}")
        print(f"Max Im: {result['best_max_im']:.4f}")
        print(f"Max Re: {result['best_max_re']:.4f}")
        print(f"Near Miss Ratio: {result['best_nm']:.4f}")
        
    # Save best tree
    out_path = f"best_roots_tree_n{args.n}.json"
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"Saved best tree to {out_path}")

if __name__ == "__main__":
    main()
