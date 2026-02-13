#!/usr/bin/env python3
"""Multi-island evolutionary LC-breaker search.

Goal: maximize log-concavity violation ratio
    max_k (i_{k-1} i_{k+1} / i_k^2)
with a large bonus for non-unimodality hits.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import time
from dataclasses import dataclass
from typing import Any

import networkx as nx

from graph6 import parse_graph6
from indpoly import independence_poly, is_unimodal, log_concavity_ratio


def adj_to_nx(n: int, adj: list[list[int]]) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from(range(n))
    for u, neigh in enumerate(adj):
        for v in neigh:
            if u < v:
                g.add_edge(u, v)
    return g


def nx_to_adj(g: nx.Graph) -> tuple[int, list[list[int]]]:
    mapping = {node: i for i, node in enumerate(g.nodes())}
    h = nx.relabel_nodes(g, mapping)
    n = h.number_of_nodes()
    adj = [[] for _ in range(n)]
    for u, v in h.edges():
        adj[u].append(v)
        adj[v].append(u)
    for u in range(n):
        adj[u].sort()
    return n, adj


@dataclass
class Individual:
    graph: nx.Graph
    poly: list[int] | None = None
    lc_ratio: float | None = None
    lc_pos: int | None = None
    unimodal: bool | None = None
    fitness: float | None = None

    def clone(self, keep_eval: bool = False) -> "Individual":
        if keep_eval:
            return Individual(
                graph=self.graph.copy(),
                poly=None if self.poly is None else self.poly[:],
                lc_ratio=self.lc_ratio,
                lc_pos=self.lc_pos,
                unimodal=self.unimodal,
                fitness=self.fitness,
            )
        return Individual(graph=self.graph.copy())

    def evaluate(self) -> None:
        if self.fitness is not None:
            return
        if not nx.is_tree(self.graph):
            # Penalize invalid graphs hard.
            self.poly = [1]
            self.lc_ratio = 0.0
            self.lc_pos = -1
            self.unimodal = False
            self.fitness = -1e9
            return
        n, adj = nx_to_adj(self.graph)
        self.poly = independence_poly(n, adj)
        self.lc_ratio, self.lc_pos = log_concavity_ratio(self.poly)
        self.unimodal = is_unimodal(self.poly)
        self.fitness = self.lc_ratio + (100.0 if not self.unimodal else 0.0)


def mutate_add_leaf(ind: Individual, rng: random.Random, max_nodes: int) -> Individual:
    g = ind.graph.copy()
    if g.number_of_nodes() >= max_nodes:
        return ind.clone(keep_eval=True)
    parent = rng.choice(list(g.nodes()))
    new_node = max(g.nodes()) + 1 if g.number_of_nodes() else 0
    g.add_edge(parent, new_node)
    return Individual(g)


def mutate_remove_leaf(
    ind: Individual, rng: random.Random, min_nodes: int = 10
) -> Individual:
    g = ind.graph.copy()
    if g.number_of_nodes() <= min_nodes:
        return ind.clone(keep_eval=True)
    leaves = [v for v in g.nodes() if g.degree(v) == 1]
    if not leaves:
        return ind.clone(keep_eval=True)
    g.remove_node(rng.choice(leaves))
    return Individual(g) if nx.is_tree(g) else ind.clone(keep_eval=True)


def mutate_graft(ind: Individual, rng: random.Random) -> Individual:
    g = ind.graph.copy()
    if g.number_of_nodes() < 4:
        return ind.clone(keep_eval=True)
    u, v = rng.choice(list(g.edges()))
    g.remove_edge(u, v)
    comps = list(nx.connected_components(g))
    if len(comps) != 2:
        return ind.clone(keep_eval=True)
    c1 = list(comps[0])
    c2 = list(comps[1])
    a = rng.choice(c1)
    b = rng.choice(c2)
    g.add_edge(a, b)
    return Individual(g) if nx.is_tree(g) else ind.clone(keep_eval=True)


def mutate_rewire_leaf(ind: Individual, rng: random.Random) -> Individual:
    g = ind.graph.copy()
    leaves = [v for v in g.nodes() if g.degree(v) == 1]
    if not leaves:
        return ind.clone(keep_eval=True)
    leaf = rng.choice(leaves)
    old_parent = next(iter(g.neighbors(leaf)))
    candidates = [v for v in g.nodes() if v != leaf and v != old_parent]
    if not candidates:
        return ind.clone(keep_eval=True)
    new_parent = rng.choice(candidates)
    g.remove_edge(leaf, old_parent)
    g.add_edge(leaf, new_parent)
    return Individual(g) if nx.is_tree(g) else ind.clone(keep_eval=True)


def mutate_inject_broom(ind: Individual, rng: random.Random, max_nodes: int) -> Individual:
    g = ind.graph.copy()
    hub = max(g.nodes(), key=lambda x: g.degree(x))
    if rng.random() < 0.55 and g.number_of_nodes() < max_nodes:
        new_node = max(g.nodes()) + 1
        g.add_edge(hub, new_node)
        return Individual(g)
    leaves = [v for v in g.nodes() if g.degree(v) == 1 and next(iter(g.neighbors(v))) != hub]
    if not leaves:
        return ind.clone(keep_eval=True)
    leaf = rng.choice(leaves)
    old_parent = next(iter(g.neighbors(leaf)))
    g.remove_edge(leaf, old_parent)
    g.add_edge(leaf, hub)
    return Individual(g) if nx.is_tree(g) else ind.clone(keep_eval=True)


def _extract_cut(
    g: nx.Graph, rng: random.Random
) -> tuple[nx.Graph, nx.Graph, int, int] | None:
    if g.number_of_edges() == 0:
        return None
    u, v = rng.choice(list(g.edges()))
    h = g.copy()
    h.remove_edge(u, v)
    comps = list(nx.connected_components(h))
    if len(comps) != 2:
        return None
    cu, cv = (comps[0], comps[1]) if u in comps[0] else (comps[1], comps[0])
    # Use larger side as base to keep child size stable-ish.
    if len(cu) >= len(cv):
        base_nodes, sub_nodes = cu, cv
        attach_base, attach_sub = u, v
    else:
        base_nodes, sub_nodes = cv, cu
        attach_base, attach_sub = v, u
    base = g.subgraph(base_nodes).copy()
    sub = g.subgraph(sub_nodes).copy()
    return base, sub, attach_base, attach_sub


def _build_child(
    base: nx.Graph, sub: nx.Graph, attach_base: int, attach_sub: int
) -> Individual | None:
    if base.number_of_nodes() == 0 or sub.number_of_nodes() == 0:
        return None
    base_map = {node: i for i, node in enumerate(base.nodes())}
    base_r = nx.relabel_nodes(base, base_map)
    attach_base_r = base_map[attach_base]

    off = base_r.number_of_nodes()
    sub_map = {node: off + i for i, node in enumerate(sub.nodes())}
    sub_r = nx.relabel_nodes(sub, sub_map)
    attach_sub_r = sub_map[attach_sub]

    child = nx.compose(base_r, sub_r)
    child.add_edge(attach_base_r, attach_sub_r)
    if not nx.is_tree(child):
        return None
    return Individual(child)


def crossover_subtree_exchange(
    p1: Individual, p2: Individual, rng: random.Random
) -> tuple[Individual, Individual]:
    c1 = _extract_cut(p1.graph, rng)
    c2 = _extract_cut(p2.graph, rng)
    if c1 is None or c2 is None:
        return p1.clone(keep_eval=True), p2.clone(keep_eval=True)
    b1, s1, a1, x1 = c1
    b2, s2, a2, x2 = c2
    ch1 = _build_child(b1, s2, a1, x2)
    ch2 = _build_child(b2, s1, a2, x1)
    if ch1 is None or ch2 is None:
        return p1.clone(keep_eval=True), p2.clone(keep_eval=True)
    return ch1, ch2


def mutate_individual(
    ind: Individual,
    rng: random.Random,
    max_nodes: int,
    allow_remove_leaf: bool,
) -> Individual:
    r = rng.random()
    if r < 0.38:
        return mutate_add_leaf(ind, rng, max_nodes)
    if allow_remove_leaf and r < 0.48:
        return mutate_remove_leaf(ind, rng)
    if r < 0.72:
        return mutate_graft(ind, rng)
    if r < 0.86:
        return mutate_rewire_leaf(ind, rng)
    return mutate_inject_broom(ind, rng, max_nodes)


def load_seeds(path: str) -> list[Individual]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    out: list[Individual] = []
    seen: set[str] = set()
    for key in ("lc_failures", "top_near_misses"):
        for item in data.get(key, []):
            g6 = item.get("graph6")
            if not g6 or g6 in seen:
                continue
            seen.add(g6)
            n, adj = parse_graph6(g6.encode("ascii"))
            g = adj_to_nx(n, adj)
            if nx.is_tree(g):
                out.append(Individual(g))
    if not out:
        raise RuntimeError(f"No seeds loaded from {path}")
    return out


def init_islands(
    seeds: list[Individual],
    pop_size: int,
    islands: int,
    rng: random.Random,
    max_nodes: int,
    allow_remove_leaf: bool,
) -> list[list[Individual]]:
    islands = max(1, min(islands, pop_size))
    sizes = [pop_size // islands] * islands
    for i in range(pop_size % islands):
        sizes[i] += 1

    pops: list[list[Individual]] = []
    for island_idx, size in enumerate(sizes):
        pop: list[Individual] = []
        # Keep at least one exact seed anchor per island.
        pop.append(seeds[island_idx % len(seeds)].clone(keep_eval=True))
        while len(pop) < size:
            if rng.random() < 0.2:
                pop.append(rng.choice(seeds).clone(keep_eval=True))
                continue
            parent = rng.choice(seeds).clone()
            child = parent
            for _ in range(rng.randint(1, 3)):
                child = mutate_individual(child, rng, max_nodes, allow_remove_leaf)
            pop.append(child)
        pops.append(pop)
    return pops


def migrate_ring(
    islands: list[list[Individual]], migrants: int
) -> None:
    if len(islands) < 2 or migrants <= 0:
        return
    outgoing: list[list[Individual]] = []
    for pop in islands:
        pop.sort(key=lambda x: x.fitness if x.fitness is not None else -1e12, reverse=True)
        m = min(migrants, max(1, len(pop) // 3))
        outgoing.append([ind.clone(keep_eval=True) for ind in pop[:m]])
    for i in range(len(islands)):
        target = islands[(i + 1) % len(islands)]
        incoming = outgoing[i]
        target.sort(key=lambda x: x.fitness if x.fitness is not None else -1e12)
        del target[: len(incoming)]
        target.extend(incoming)


def main() -> None:
    ap = argparse.ArgumentParser(description="Multi-island evolutionary LC-breaker search.")
    ap.add_argument("--generations", type=int, default=600)
    ap.add_argument("--pop-size", type=int, default=120)
    ap.add_argument("--seeds-file", default="results/analysis_n26.json")
    ap.add_argument("--out", default="results/best_evolutionary_tree.json")
    ap.add_argument("--seed", type=int, default=993)
    ap.add_argument("--max-nodes", type=int, default=60)
    ap.add_argument("--islands", type=int, default=6)
    ap.add_argument("--elite-frac", type=float, default=0.15)
    ap.add_argument("--crossover-prob", type=float, default=0.25)
    ap.add_argument("--allow-remove-leaf", action="store_true")
    ap.add_argument("--migration-interval", type=int, default=20)
    ap.add_argument("--migrants", type=int, default=2)
    ap.add_argument("--report-every", type=int, default=25)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    t0 = time.time()
    seeds = load_seeds(args.seeds_file)
    islands = init_islands(
        seeds,
        args.pop_size,
        args.islands,
        rng,
        args.max_nodes,
        args.allow_remove_leaf,
    )

    best: Individual | None = None
    history: list[dict[str, Any]] = []
    counterexample: dict[str, Any] | None = None

    for gen in range(1, args.generations + 1):
        for pop in islands:
            for ind in pop:
                ind.evaluate()
            pop.sort(key=lambda x: x.fitness if x.fitness is not None else -1e12, reverse=True)
            if best is None or (pop[0].fitness or -1e12) > (best.fitness or -1e12):
                best = pop[0].clone(keep_eval=True)
            if not pop[0].unimodal:
                counterexample = {
                    "generation": gen,
                    "lc_ratio": pop[0].lc_ratio,
                    "lc_pos": pop[0].lc_pos,
                    "poly": pop[0].poly,
                }
                best = pop[0].clone(keep_eval=True)
                break
        if counterexample is not None:
            break

        next_islands: list[list[Individual]] = []
        for pop in islands:
            size = len(pop)
            elite_count = max(2, int(size * args.elite_frac))
            pool_count = max(elite_count, int(size * 0.5))
            elite = [x.clone(keep_eval=True) for x in pop[:elite_count]]
            pool = pop[:pool_count]
            nxt = elite[:]
            while len(nxt) < size:
                if rng.random() < args.crossover_prob and len(pool) >= 2:
                    p1, p2 = rng.sample(pool, 2)
                    c1, c2 = crossover_subtree_exchange(p1, p2, rng)
                    nxt.append(c1)
                    if len(nxt) < size:
                        nxt.append(c2)
                else:
                    parent = rng.choice(pool).clone(keep_eval=True)
                    child = parent
                    for _ in range(rng.randint(1, 3)):
                        child = mutate_individual(
                            child, rng, args.max_nodes, args.allow_remove_leaf
                        )
                    nxt.append(child)
            next_islands.append(nxt)
        islands = next_islands

        if args.migration_interval > 0 and gen % args.migration_interval == 0:
            # Evaluate before migration so ring uses known fitness values.
            for pop in islands:
                for ind in pop:
                    ind.evaluate()
            migrate_ring(islands, args.migrants)

        if gen % args.report_every == 0 or gen == 1:
            island_best = []
            for pop in islands:
                for ind in pop:
                    ind.evaluate()
                pop.sort(key=lambda x: x.fitness if x.fitness is not None else -1e12, reverse=True)
                island_best.append(pop[0].fitness or -1e12)
            row = {
                "generation": gen,
                "best_lc_ratio": best.lc_ratio if best else None,
                "best_unimodal": best.unimodal if best else None,
                "island_best_mean": sum(island_best) / len(island_best),
            }
            history.append(row)
            print(
                f"gen={gen:4d} best_lc={row['best_lc_ratio']:.12f} "
                f"best_unimodal={row['best_unimodal']} "
                f"island_mean={row['island_best_mean']:.6f}"
            )

    if best is None:
        raise RuntimeError("No individuals evaluated.")
    best.evaluate()
    n_best, adj_best = nx_to_adj(best.graph)

    result = {
        "generations": args.generations,
        "pop_size": args.pop_size,
        "islands": args.islands,
        "seed": args.seed,
        "max_nodes": args.max_nodes,
        "allow_remove_leaf": args.allow_remove_leaf,
        "migration_interval": args.migration_interval,
        "migrants": args.migrants,
        "elapsed_s": time.time() - t0,
        "counterexample_found": counterexample is not None,
        "counterexample": counterexample,
        "best": {
            "n": n_best,
            "adj": adj_best,
            "poly": best.poly,
            "lc_ratio": best.lc_ratio,
            "lc_pos": best.lc_pos,
            "unimodal": best.unimodal,
            "fitness": best.fitness,
        },
        "history": history,
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
