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


def _adj_fingerprint(adj: list[list[int]]) -> str:
    """Canonical edge-list fingerprint for exact dedup."""
    edges = []
    for u, neighbors in enumerate(adj):
        for v in neighbors:
            if u < v:
                edges.append(f"{u}-{v}")
    return ";".join(edges)


def _selection_score(lc_ratio: float, prospect_lc_ratio: float, prospect_bonus: float) -> float:
    """Score a tree by its current LC ratio plus observed mutation upside."""
    return lc_ratio + prospect_bonus * max(0.0, prospect_lc_ratio - lc_ratio)


def _population_sort_key(candidate: dict[str, Any]) -> tuple:
    return (
        -candidate["score"],
        -candidate["lc_ratio"],
        -candidate["nm_ratio"],
        candidate["generation"],
        candidate["origin"],
    )


def _archive_sort_key(candidate: dict[str, Any]) -> tuple:
    return (
        -candidate["lc_ratio"],
        -candidate["score"],
        -candidate["nm_ratio"],
        candidate["generation"],
        candidate["origin"],
    )


def _dedupe_candidates(
    candidates: list[dict[str, Any]],
    *,
    archive_mode: bool,
) -> list[dict[str, Any]]:
    """Keep the strongest representative for each exact tree."""
    best_by_fp: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        fingerprint = candidate["fingerprint"]
        incumbent = best_by_fp.get(fingerprint)
        if incumbent is None:
            best_by_fp[fingerprint] = candidate
            continue
        better = (
            _archive_sort_key(candidate) < _archive_sort_key(incumbent)
            if archive_mode
            else _population_sort_key(candidate) < _population_sort_key(incumbent)
        )
        if better:
            best_by_fp[fingerprint] = candidate
    return list(best_by_fp.values())


def _update_archive(
    archive: list[dict[str, Any]],
    candidates: list[dict[str, Any]],
    archive_size: int,
) -> list[dict[str, Any]]:
    """Retain the best distinct exact trees by raw LC ratio."""
    merged = _dedupe_candidates(archive + [deepcopy(c) for c in candidates], archive_mode=True)
    merged.sort(key=_archive_sort_key)
    return merged[:archive_size]


def _make_individual(
    adj: list[list[int]],
    origin: str,
    ev: dict[str, Any],
    generation: int,
    prospect_bonus: float,
    prospect_lc_ratio: float | None = None,
) -> dict[str, Any]:
    """Package a tree with objective and selection metadata."""
    if prospect_lc_ratio is None:
        prospect_lc_ratio = ev["lc_ratio"]
    return {
        "adj": adj,
        "origin": origin,
        "generation": generation,
        "score": _selection_score(ev["lc_ratio"], prospect_lc_ratio, prospect_bonus),
        "prospect_lc_ratio": prospect_lc_ratio,
        "prospect_gain": max(0.0, prospect_lc_ratio - ev["lc_ratio"]),
        "fingerprint": _adj_fingerprint(adj),
        **ev,
    }


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
    archive_size: int,
    prospect_bonus: float,
    seed: int,
    verbose_every: int,
) -> dict[str, Any]:
    rng = random.Random(seed)
    elite_count = max(2, int(pop_size * elite_frac))
    archive_size = max(2, archive_size)
    seed_g6 = load_seed_graphs(analysis_path)
    if not seed_g6:
        raise RuntimeError(f"No lc_failures found in {analysis_path}")

    population: list[dict[str, Any]] = []

    for g6 in seed_g6:
        n0, adj0 = parse_graph6(g6.encode("ascii"))
        if n0 != n:
            continue
        ev = evaluate(n, adj0)
        population.append(_make_individual(adj0, "seed_lc_failure", ev, 0, prospect_bonus))

    if not population:
        raise RuntimeError(f"No n={n} seed graphs available in {analysis_path}")

    # Fill initial population with local mutations + random diversity.
    while len(population) < pop_size:
        if rng.random() < 0.35:
            child_adj = random_tree(n, rng)
            ev = evaluate(n, child_adj)
            population.append(_make_individual(child_adj, "random_seed", ev, 0, prospect_bonus))
            continue
        base = rng.choice(population)
        child_adj = deepcopy(base["adj"])
        steps = rng.randint(1, 4)
        for _ in range(steps):
            _, child_adj = mutate(n, child_adj, rng)
        if not validate_tree(n, child_adj):
            continue
        ev = evaluate(n, child_adj)
        population.append(_make_individual(child_adj, "seed_mutation", ev, 0, prospect_bonus))

    population = _dedupe_candidates(population, archive_mode=False)
    population.sort(key=_population_sort_key)
    population = population[:pop_size]
    archive = _update_archive([], population, archive_size)
    best = deepcopy(population[0])
    history: list[dict[str, Any]] = [
        {
            "generation": 0,
            "best_lc_ratio": best["lc_ratio"],
            "best_score": best["score"],
            "best_nm_ratio": best["nm_ratio"],
            "best_unimodal": best["unimodal"],
            "archive_best_lc_ratio": archive[0]["lc_ratio"],
            "archive_size": len(archive),
        }
    ]

    t0 = time.time()
    for gen in range(1, generations + 1):
        elite = population[:elite_count]
        offspring: list[dict[str, Any]] = []
        parent_best_child: dict[str, float] = {
            parent["fingerprint"]: parent["lc_ratio"] for parent in elite
        }

        while len(offspring) < pop_size:
            if rng.random() < 0.15:
                child_adj = random_tree(n, rng)
                ev = evaluate(n, child_adj)
                offspring.append(
                    _make_individual(
                        child_adj,
                        f"inject_random_g{gen}",
                        ev,
                        gen,
                        prospect_bonus,
                    )
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
            parent_best_child[parent["fingerprint"]] = max(
                parent_best_child[parent["fingerprint"]],
                ev["lc_ratio"],
            )
            offspring.append(
                _make_individual(
                    child_adj,
                    f"gen_{gen}",
                    ev,
                    gen,
                    prospect_bonus,
                )
            )

        refreshed_elite = [
            _make_individual(
                deepcopy(parent["adj"]),
                parent["origin"],
                {
                    "lc_ratio": parent["lc_ratio"],
                    "lc_pos": parent["lc_pos"],
                    "nm_ratio": parent["nm_ratio"],
                    "nm_pos": parent["nm_pos"],
                    "unimodal": parent["unimodal"],
                    "poly": parent["poly"],
                },
                parent["generation"],
                prospect_bonus,
                prospect_lc_ratio=parent_best_child[parent["fingerprint"]],
            )
            for parent in elite
        ]

        archive = _update_archive(archive, refreshed_elite + offspring, archive_size)
        merged = _dedupe_candidates(refreshed_elite + offspring + archive, archive_mode=False)
        merged.sort(key=_population_sort_key)
        population = merged[:pop_size]
        best_candidate = max(archive + population, key=lambda x: x["lc_ratio"])
        if best_candidate["lc_ratio"] > best["lc_ratio"]:
            best = deepcopy(best_candidate)

        if gen % verbose_every == 0 or gen == 1:
            row = {
                "generation": gen,
                "best_lc_ratio": population[0]["lc_ratio"],
                "best_score": population[0]["score"],
                "best_nm_ratio": population[0]["nm_ratio"],
                "best_unimodal": population[0]["unimodal"],
                "archive_best_lc_ratio": archive[0]["lc_ratio"],
                "archive_size": len(archive),
            }
            history.append(row)
            print(
                f"gen={gen:4d} lc={row['best_lc_ratio']:.12f} "
                f"score={row['best_score']:.12f} "
                f"archive_lc={row['archive_best_lc_ratio']:.12f} "
                f"nm={row['best_nm_ratio']:.12f} unimodal={row['best_unimodal']}"
            )

    result = {
        "n": n,
        "analysis_path": analysis_path,
        "seed": seed,
        "pop_size": pop_size,
        "generations": generations,
        "elite_frac": elite_frac,
        "archive_size": archive_size,
        "prospect_bonus": prospect_bonus,
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
    ap.add_argument("--archive-size", type=int, default=16)
    ap.add_argument("--prospect-bonus", type=float, default=0.5)
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
        archive_size=args.archive_size,
        prospect_bonus=args.prospect_bonus,
        seed=args.seed,
        verbose_every=args.verbose_every,
    )


if __name__ == "__main__":
    main()
