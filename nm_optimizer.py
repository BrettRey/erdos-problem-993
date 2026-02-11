#!/usr/bin/env python3
"""Adversarial near-miss ratio optimizer for Erdős Problem #993.

Uses evolutionary strategies to search for trees that maximize the
near-miss ratio (nm). If nm > 1, we have a counterexample to the
unimodality conjecture. If the optimizer consistently converges to
broom-like structures, that's strong evidence brooms are the right
family.

Mutations:
  1. Leaf relocation: move a leaf to a different vertex.
  2. Subtree prune-and-regraft (SPR): detach a subtree and reattach
     at a random edge of the remaining tree.
  3. Pendant transfer: move a leaf from a high-degree to a low-degree
     vertex (biased toward creating broom-like concentration).

Seeds: brooms at various (p, s), caterpillars, SSTs, random trees.
"""

import argparse
import json
import os
import random
import time
from copy import deepcopy

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
    # BFS connectivity check
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
    """Quick structural fingerprint for dedup."""
    n = len(adj)
    degs = sorted(len(adj[v]) for v in range(n))
    leaves = sum(1 for d in degs if d == 1)
    max_deg = degs[-1] if degs else 0
    branch_verts = sum(1 for d in degs if d >= 3)
    return f"n{n}_l{leaves}_d{max_deg}_b{branch_verts}"


# ── Mutations ───────────────────────────────────────────────────────


def mutate_leaf_move(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    """Move a random leaf to a random other vertex."""
    adj = deepcopy(adj)
    lvs = _leaves(adj)
    if not lvs:
        return n, adj

    leaf = rng.choice(lvs)
    neighbor = adj[leaf][0]

    # Remove edge leaf--neighbor
    adj[neighbor].remove(leaf)
    adj[leaf].remove(neighbor)

    # Pick a new attachment point (not the leaf itself, not the old neighbor)
    candidates = [v for v in range(n) if v != leaf and v != neighbor]
    if not candidates:
        # Reattach to old neighbor (no-op)
        adj[neighbor].append(leaf)
        adj[leaf].append(neighbor)
        return n, adj

    target = rng.choice(candidates)
    adj[target].append(leaf)
    adj[leaf].append(target)
    return n, adj


def mutate_spr(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    """Subtree prune-and-regraft: detach a subtree, reattach elsewhere."""
    if n <= 3:
        return mutate_leaf_move(n, adj, rng)

    adj = deepcopy(adj)

    # Pick a random non-root edge to cut
    # First pick a random vertex v != 0 to define the edge v--parent(v)
    # We'll root at 0 for this purpose
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

    # Pick a random non-root vertex to define the pruned subtree
    pruned_root = rng.choice(order[1:])
    pruned_parent = parent[pruned_root]

    # Find all vertices in the pruned subtree
    subtree = set()
    stack = [pruned_root]
    while stack:
        v = stack.pop()
        subtree.add(v)
        for u in adj[v]:
            if u != pruned_parent and u not in subtree:
                stack.append(u)

    # Remove the edge pruned_root -- pruned_parent
    adj[pruned_parent].remove(pruned_root)
    adj[pruned_root].remove(pruned_parent)

    # Pick a random vertex NOT in the subtree to reattach to
    remaining = [v for v in range(n) if v not in subtree]
    if not remaining:
        # Shouldn't happen, but be safe
        adj[pruned_parent].append(pruned_root)
        adj[pruned_root].append(pruned_parent)
        return n, adj

    target = rng.choice(remaining)
    adj[target].append(pruned_root)
    adj[pruned_root].append(target)

    return n, adj


def mutate_pendant_concentrate(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    """Move a leaf from a low-degree vertex to the highest-degree vertex.

    Biased toward creating star/broom-like concentration.
    """
    adj = deepcopy(adj)
    lvs = _leaves(adj)
    if len(lvs) < 2:
        return n, adj

    # Find the vertex with highest degree
    max_deg_v = max(range(n), key=lambda v: len(adj[v]))

    # Pick a leaf NOT attached to max_deg_v
    candidates = [lf for lf in lvs if adj[lf][0] != max_deg_v]
    if not candidates:
        return n, adj

    leaf = rng.choice(candidates)
    old_neighbor = adj[leaf][0]

    # Detach
    adj[old_neighbor].remove(leaf)
    adj[leaf].remove(old_neighbor)

    # Check: would removing this leaf disconnect the old_neighbor's component?
    # In a tree, removing a leaf never disconnects. Safe.

    # But we need to make sure old_neighbor is still connected to the rest.
    # If old_neighbor had degree 1 (i.e., it was only connected to leaf),
    # then it becomes isolated. We need to handle this.
    # Actually in a tree, if old_neighbor had degree 1, it was a leaf too,
    # and we'd be creating an isolated vertex. Let's check:
    if not adj[old_neighbor]:
        # old_neighbor is now isolated -- reattach it to max_deg_v too
        adj[max_deg_v].append(old_neighbor)
        adj[old_neighbor].append(max_deg_v)
        # And attach the leaf to max_deg_v
        adj[max_deg_v].append(leaf)
        adj[leaf].append(max_deg_v)
    else:
        # Reattach leaf to max_deg_v
        adj[max_deg_v].append(leaf)
        adj[leaf].append(max_deg_v)

    return n, adj


def mutate_edge_subdivide_contract(n: int, adj: list[list[int]], rng: random.Random) -> tuple[int, list[list[int]]]:
    """Subdivide one edge and contract another, keeping n fixed.

    Explores degree-2 vertex placement, which matters for LC failures.
    """
    if n <= 3:
        return mutate_leaf_move(n, adj, rng)

    adj = deepcopy(adj)

    # Find a degree-2 vertex to contract (if any)
    deg2 = [v for v in range(n) if len(adj[v]) == 2]
    if not deg2:
        return mutate_leaf_move(n, adj, rng)

    # Contract: pick a random degree-2 vertex, merge it with one neighbor
    contract_v = rng.choice(deg2)
    u, w = adj[contract_v]

    # Remove edges contract_v--u and contract_v--w
    adj[u].remove(contract_v)
    adj[w].remove(contract_v)
    adj[contract_v] = []

    # Add edge u--w (if not already)
    if w not in adj[u]:
        adj[u].append(w)
        adj[w].append(u)

    # Now contract_v is isolated. Pick a random edge to subdivide with it.
    edges = []
    for v in range(n):
        for nb in adj[v]:
            if v < nb:
                edges.append((v, nb))

    if not edges:
        return n, adj  # degenerate

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
    """Apply a random mutation."""
    r = rng.random()
    cumulative = 0.0
    for func, weight in MUTATIONS:
        cumulative += weight
        if r < cumulative:
            return func(n, adj, rng)
    return MUTATIONS[-1][0](n, adj, rng)


# ── Seed generation ─────────────────────────────────────────────────


def generate_seeds(n: int, count: int, rng: random.Random) -> list[tuple[int, list[list[int]], str]]:
    """Generate seed trees at vertex count n."""
    seeds = []

    # Brooms: sweep path lengths
    for p in range(2, min(n - 1, 50)):
        s = n - p
        if s < 1:
            continue
        n_t, adj = make_broom(p, s)
        seeds.append((n_t, adj, f"broom({p},{s})"))
        if len(seeds) >= count // 2:
            break

    # Caterpillars with concentrated pendants
    for k in [1, 2, 3, 5, 8]:
        spine = n // (k + 1)
        if spine < 2:
            continue
        pendants = [k] * spine
        total = spine + sum(pendants)
        if total > n:
            continue
        # Pad remaining vertices as extra pendants on first vertex
        pendants[0] += n - total
        n_t, adj = make_caterpillar(spine, pendants)
        seeds.append((n_t, adj, f"cat_conc({spine},k={k})"))

    # Spiders with varying legs
    for num_legs in [3, 5, 8, 12]:
        leg_len = max(1, (n - 1) // num_legs)
        actual_n = 1 + num_legs * leg_len
        if actual_n <= n:
            legs = [leg_len] * num_legs
            n_t, adj = make_spider(legs)
            seeds.append((n_t, adj, f"spider({num_legs}x{leg_len})"))

    # Random trees via Prüfer codes
    while len(seeds) < count:
        prufer = [rng.randrange(n) for _ in range(n - 2)]
        adj_r: list[list[int]] = [[] for _ in range(n)]
        degree = [1] * n
        for v in prufer:
            degree[v] += 1
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


# ── Evaluation ──────────────────────────────────────────────────────


def evaluate(n: int, adj: list[list[int]]) -> dict:
    """Compute nm and check unimodality."""
    poly = independence_poly(n, adj)
    uni = is_unimodal(poly)
    nm, nm_pos = near_miss_ratio(poly)
    return {
        "nm": nm,
        "nm_pos": nm_pos,
        "unimodal": uni,
        "alpha": len(poly) - 1,
        "poly_len": len(poly),
    }


# ── Main optimizer ──────────────────────────────────────────────────


def run_optimizer(
    n: int,
    pop_size: int = 60,
    generations: int = 500,
    elite_frac: float = 0.2,
    mutations_per_ind: int = 3,
    seed: int = 42,
    verbose: bool = True,
) -> dict:
    """Run the evolutionary nm optimizer at fixed vertex count n."""
    rng = random.Random(seed)
    elite_count = max(2, int(pop_size * elite_frac))

    # Generate seeds
    seeds = generate_seeds(n, pop_size, rng)

    # Evaluate initial population
    population = []
    for n_t, adj, label in seeds:
        if n_t != n:
            continue
        ev = evaluate(n_t, adj)
        population.append({
            "n": n_t,
            "adj": adj,
            "label": label,
            "nm": ev["nm"],
            "unimodal": ev["unimodal"],
            "alpha": ev["alpha"],
            "generation": 0,
        })
        if not ev["unimodal"]:
            return _counterexample_result(n_t, adj, label, ev, 0)

    population.sort(key=lambda x: -x["nm"])

    best_ever = population[0]["nm"]
    best_ever_label = population[0]["label"]
    best_ever_gen = 0
    best_ever_adj = deepcopy(population[0]["adj"])
    stagnation = 0

    history = [{
        "generation": 0,
        "best_nm": population[0]["nm"],
        "mean_nm": sum(p["nm"] for p in population) / len(population),
        "best_label": population[0]["label"],
        "signature": _tree_signature(population[0]["adj"]),
    }]

    if verbose:
        print(f"Gen   0: best_nm={population[0]['nm']:.10f}  "
              f"mean={history[0]['mean_nm']:.6f}  "
              f"best={population[0]['label']}  "
              f"sig={_tree_signature(population[0]['adj'])}")

    for gen in range(1, generations + 1):
        # Select elite
        elite = population[:elite_count]

        # Generate offspring via mutation
        offspring = []
        for parent in elite:
            for _ in range(mutations_per_ind):
                child_n, child_adj = parent["n"], deepcopy(parent["adj"])
                # Apply 1-3 mutations
                num_muts = rng.randint(1, 3)
                for _ in range(num_muts):
                    child_n, child_adj = mutate(child_n, child_adj, rng)

                if not _validate_tree(child_n, child_adj):
                    continue

                ev = evaluate(child_n, child_adj)
                child_label = f"mut_g{gen}_{len(offspring)}"
                offspring.append({
                    "n": child_n,
                    "adj": child_adj,
                    "label": child_label,
                    "nm": ev["nm"],
                    "unimodal": ev["unimodal"],
                    "alpha": ev["alpha"],
                    "generation": gen,
                })

                if not ev["unimodal"]:
                    return _counterexample_result(
                        child_n, child_adj, child_label, ev, gen
                    )

        # Also inject a few fresh random individuals to avoid premature convergence
        for _ in range(max(2, pop_size // 10)):
            prufer = [rng.randrange(n) for _ in range(n - 2)]
            adj_r: list[list[int]] = [[] for _ in range(n)]
            degree = [1] * n
            for v in prufer:
                degree[v] += 1
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
                "n": n,
                "adj": adj_r,
                "label": f"inject_g{gen}",
                "nm": ev["nm"],
                "unimodal": ev["unimodal"],
                "alpha": ev["alpha"],
                "generation": gen,
            })

        # Merge and select
        combined = elite + offspring
        combined.sort(key=lambda x: -x["nm"])
        population = combined[:pop_size]

        if population[0]["nm"] > best_ever:
            best_ever = population[0]["nm"]
            best_ever_label = population[0]["label"]
            best_ever_gen = gen
            best_ever_adj = deepcopy(population[0]["adj"])
            stagnation = 0
        else:
            stagnation += 1

        if gen % 25 == 0 or gen == 1:
            sig = _tree_signature(population[0]["adj"])
            mean_nm = sum(p["nm"] for p in population) / len(population)
            history.append({
                "generation": gen,
                "best_nm": population[0]["nm"],
                "mean_nm": mean_nm,
                "best_label": population[0]["label"],
                "signature": sig,
            })
            if verbose:
                print(f"Gen {gen:>3}: best_nm={population[0]['nm']:.10f}  "
                      f"mean={mean_nm:.6f}  sig={sig}  "
                      f"stag={stagnation}")

    # Analyze best-ever tree
    best_degs = _degree_sequence(best_ever_adj)
    best_leaves = sum(1 for d in best_degs if d == 1)
    best_max_deg = best_degs[0] if best_degs else 0
    best_branches = sum(1 for d in best_degs if d >= 3)

    return {
        "counterexample_found": False,
        "n": n,
        "generations": generations,
        "pop_size": pop_size,
        "seed": seed,
        "best_nm": best_ever,
        "best_label": best_ever_label,
        "best_generation": best_ever_gen,
        "best_structure": {
            "leaves": best_leaves,
            "max_degree": best_max_deg,
            "branch_vertices": best_branches,
            "degree_top5": best_degs[:5],
            "signature": _tree_signature(best_ever_adj),
        },
        "history": history,
    }


def _counterexample_result(n, adj, label, ev, gen):
    poly = independence_poly(n, adj)
    return {
        "counterexample_found": True,
        "n": n,
        "label": label,
        "generation": gen,
        "nm": ev["nm"],
        "alpha": ev["alpha"],
        "poly": poly,
        "adj": adj,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Adversarial nm optimizer for Erdős Problem #993"
    )
    parser.add_argument(
        "--min-n", type=int, default=50,
        help="Minimum vertex count (default: 50)",
    )
    parser.add_argument(
        "--max-n", type=int, default=200,
        help="Maximum vertex count (default: 200)",
    )
    parser.add_argument(
        "--step-n", type=int, default=25,
        help="Step between vertex counts (default: 25)",
    )
    parser.add_argument(
        "--pop-size", type=int, default=60,
        help="Population size (default: 60)",
    )
    parser.add_argument(
        "--generations", type=int, default=500,
        help="Generations per vertex count (default: 500)",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed (default: 42)",
    )
    args = parser.parse_args()

    print("Adversarial nm optimizer for Erdős Problem #993")
    print(f"Vertex range: {args.min_n} to {args.max_n} (step {args.step_n})")
    print(f"Population: {args.pop_size}, Generations: {args.generations}")
    print(f"{'='*70}\n")

    os.makedirs("results", exist_ok=True)
    all_results = []

    for n in range(args.min_n, args.max_n + 1, args.step_n):
        print(f"\n--- n = {n} ---")
        t0 = time.time()
        result = run_optimizer(
            n=n,
            pop_size=args.pop_size,
            generations=args.generations,
            seed=args.seed,
        )
        elapsed = time.time() - t0
        result["elapsed_s"] = round(elapsed, 2)

        if result["counterexample_found"]:
            print(f"\n*** COUNTEREXAMPLE FOUND at n={n}! ***")
            print(f"    nm = {result['nm']}")
            print(f"    poly = {result['poly']}")
            path = f"results/counterexample_optimizer_n{n}.json"
            with open(path, "w") as f:
                json.dump(result, f, indent=2)
            print(f"    Saved to {path}")
            return

        print(f"  Best nm: {result['best_nm']:.10f}")
        print(f"  Structure: {result['best_structure']}")
        print(f"  Time: {elapsed:.1f}s")
        all_results.append(result)

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"{'n':>5}  {'best_nm':>14}  {'leaves':>6}  {'max_deg':>7}  {'branches':>8}  {'sig':>25}  {'time':>6}")
    for r in all_results:
        s = r["best_structure"]
        print(f"{r['n']:>5}  {r['best_nm']:>14.10f}  {s['leaves']:>6}  "
              f"{s['max_degree']:>7}  {s['branch_vertices']:>8}  "
              f"{s['signature']:>25}  {r['elapsed_s']:>5.1f}s")

    counterexample = any(r["counterexample_found"] for r in all_results)
    print(f"\nCounterexample found: {counterexample}")

    # Check: do the best trees look like brooms?
    broom_like = sum(
        1 for r in all_results
        if r["best_structure"]["branch_vertices"] <= 2
    )
    print(f"Best trees with <= 2 branch vertices: {broom_like}/{len(all_results)}")

    path = f"results/nm_optimizer_{args.min_n}_{args.max_n}.json"
    with open(path, "w") as f:
        json.dump({
            "description": "Adversarial nm optimizer results",
            "parameters": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "step_n": args.step_n,
                "pop_size": args.pop_size,
                "generations": args.generations,
                "seed": args.seed,
            },
            "counterexample_found": counterexample,
            "results": all_results,
        }, f, indent=2)
    print(f"\nSaved to {path}")


if __name__ == "__main__":
    main()
