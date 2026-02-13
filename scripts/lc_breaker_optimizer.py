#!/usr/bin/env python3
"""Evolutionary search for extreme log-concavity violations at fixed n.

Seeds from known n=26 LC-failure trees in results/analysis_n26.json and
optimizes max_k (i_{k-1} i_{k+1} / i_k^2).
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from copy import deepcopy
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from graph6 import parse_graph6
from indpoly import independence_poly, is_unimodal, log_concavity_ratio, near_miss_ratio
from nm_optimizer import mutate


def validate_tree(n: int, adj: list[list[int]]) -> bool:
    """Check that adjacency encodes a connected tree."""
    if len(adj) != n:
        return False
    edge_count = sum(len(adj[v]) for v in range(n)) // 2
    if edge_count != n - 1:
        return False
    seen = [False] * n
    q = [0]
    seen[0] = True
    head = 0
    while head < len(q):
        v = q[head]
        head += 1
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True
                q.append(u)
    return all(seen)


def evaluate(n: int, adj: list[list[int]]) -> dict[str, Any]:
    poly = independence_poly(n, adj)
    lc_ratio, lc_pos = log_concavity_ratio(poly)
    nm_ratio, nm_pos = near_miss_ratio(poly)
    return {
        "lc_ratio": lc_ratio,
        "lc_pos": lc_pos,
        "nm_ratio": nm_ratio,
        "nm_pos": nm_pos,
        "unimodal": is_unimodal(poly),
        "poly": poly,
    }


def random_tree(n: int, rng: random.Random) -> list[list[int]]:
    """Generate a random labeled tree via Prüfer code."""
    prufer = [rng.randrange(n) for _ in range(n - 2)]
    deg = [1] * n
    for v in prufer:
        deg[v] += 1
    adj: list[list[int]] = [[] for _ in range(n)]
    for v in prufer:
        u = min(i for i in range(n) if deg[i] == 1)
        adj[u].append(v)
        adj[v].append(u)
        deg[u] -= 1
        deg[v] -= 1
    leaves = [i for i in range(n) if deg[i] == 1]
    a, b = leaves
    adj[a].append(b)
    adj[b].append(a)
    return adj


def load_seed_graphs(analysis_path: str) -> list[str]:
    with open(analysis_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return [item["graph6"] for item in data.get("lc_failures", [])]


def run(
    n: int,
    analysis_path: str,
    out_path: str,
    pop_size: int,
    generations: int,
    elite_frac: float,
    seed: int,
    verbose_every: int,
) -> dict[str, Any]:
    rng = random.Random(seed)
    elite_count = max(2, int(pop_size * elite_frac))
    seed_g6 = load_seed_graphs(analysis_path)
    if not seed_g6:
        raise RuntimeError(f"No lc_failures found in {analysis_path}")

    population: list[dict[str, Any]] = []

    for g6 in seed_g6:
        n0, adj0 = parse_graph6(g6.encode("ascii"))
        if n0 != n:
            continue
        ev = evaluate(n, adj0)
        population.append(
            {
                "adj": adj0,
                "origin": "seed_lc_failure",
                **ev,
            }
        )

    if not population:
        raise RuntimeError(f"No n={n} seed graphs available in {analysis_path}")

    # Fill initial population with local mutations + random diversity.
    while len(population) < pop_size:
        if rng.random() < 0.35:
            child_adj = random_tree(n, rng)
            ev = evaluate(n, child_adj)
            population.append(
                {
                    "adj": child_adj,
                    "origin": "random_seed",
                    **ev,
                }
            )
            continue
        base = rng.choice(population)
        child_adj = deepcopy(base["adj"])
        steps = rng.randint(1, 4)
        for _ in range(steps):
            _, child_adj = mutate(n, child_adj, rng)
        if not validate_tree(n, child_adj):
            continue
        ev = evaluate(n, child_adj)
        population.append(
            {
                "adj": child_adj,
                "origin": "seed_mutation",
                **ev,
            }
        )

    population.sort(key=lambda x: x["lc_ratio"], reverse=True)
    best = deepcopy(population[0])
    history: list[dict[str, Any]] = [
        {
            "generation": 0,
            "best_lc_ratio": best["lc_ratio"],
            "best_nm_ratio": best["nm_ratio"],
            "best_unimodal": best["unimodal"],
        }
    ]

    t0 = time.time()
    for gen in range(1, generations + 1):
        elite = population[:elite_count]
        offspring: list[dict[str, Any]] = []

        while len(offspring) < pop_size:
            if rng.random() < 0.15:
                child_adj = random_tree(n, rng)
                ev = evaluate(n, child_adj)
                offspring.append(
                    {
                        "adj": child_adj,
                        "origin": f"inject_random_g{gen}",
                        **ev,
                    }
                )
                continue
            parent = rng.choice(elite)
            child_adj = deepcopy(parent["adj"])
            steps = rng.randint(1, 4)
            for _ in range(steps):
                _, child_adj = mutate(n, child_adj, rng)
            if not validate_tree(n, child_adj):
                continue
            ev = evaluate(n, child_adj)
            offspring.append(
                {
                    "adj": child_adj,
                    "origin": f"gen_{gen}",
                    **ev,
                }
            )

        merged = elite + offspring
        merged.sort(key=lambda x: x["lc_ratio"], reverse=True)
        population = merged[:pop_size]
        if population[0]["lc_ratio"] > best["lc_ratio"]:
            best = deepcopy(population[0])

        if gen % verbose_every == 0 or gen == 1:
            row = {
                "generation": gen,
                "best_lc_ratio": population[0]["lc_ratio"],
                "best_nm_ratio": population[0]["nm_ratio"],
                "best_unimodal": population[0]["unimodal"],
            }
            history.append(row)
            print(
                f"gen={gen:4d} lc={row['best_lc_ratio']:.12f} "
                f"nm={row['best_nm_ratio']:.12f} unimodal={row['best_unimodal']}"
            )

    result = {
        "n": n,
        "analysis_path": analysis_path,
        "seed": seed,
        "pop_size": pop_size,
        "generations": generations,
        "elite_frac": elite_frac,
        "elapsed_s": time.time() - t0,
        "best": {
            "lc_ratio": best["lc_ratio"],
            "lc_pos": best["lc_pos"],
            "nm_ratio": best["nm_ratio"],
            "nm_pos": best["nm_pos"],
            "unimodal": best["unimodal"],
            "poly": best["poly"],
        },
        "history": history,
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"wrote {out_path}")
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=26)
    ap.add_argument("--analysis", default="results/analysis_n26.json")
    ap.add_argument("--out", default="results/lc_breaker_evo_n26.json")
    ap.add_argument("--pop-size", type=int, default=96)
    ap.add_argument("--generations", type=int, default=600)
    ap.add_argument("--elite-frac", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=993)
    ap.add_argument("--verbose-every", type=int, default=25)
    args = ap.parse_args()

    run(
        n=args.n,
        analysis_path=args.analysis,
        out_path=args.out,
        pop_size=args.pop_size,
        generations=args.generations,
        elite_frac=args.elite_frac,
        seed=args.seed,
        verbose_every=args.verbose_every,
    )


if __name__ == "__main__":
    main()
