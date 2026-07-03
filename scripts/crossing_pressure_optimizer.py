#!/usr/bin/env python3
"""Evolutionary optimizer for post-descent crossing pressure.

This is a direct search for the obstruction isolated in
``notes/literature/separation_invariant_lemma_schema_2026-07-03.md``:
after the first strict descent in the independence sequence, maximize
``a[j+1] / a[j]``.  A value above 1 is a unimodality counterexample.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from copy import deepcopy
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import independence_poly, is_log_concave, is_unimodal, log_concavity_ratio  # noqa: E402
from nm_optimizer import (  # noqa: E402
    _adj_fingerprint,
    _degree_sequence,
    _tree_signature,
    _validate_tree,
    generate_seeds,
    mutate,
)
from scripts.analyze_prufer_corpus import graph6_from_adj, mode_interval  # noqa: E402


def first_descent(poly: list[int]) -> int | None:
    for k in range(1, len(poly)):
        if poly[k] < poly[k - 1]:
            return k
    return None


def max_post_descent_ratio(poly: list[int]) -> tuple[float, int | None, bool]:
    descent = first_descent(poly)
    if descent is None:
        return 0.0, None, False

    best_ratio = 0.0
    best_j: int | None = None
    strict_crossing = False
    for j in range(descent, len(poly) - 1):
        if poly[j] == 0:
            continue
        if poly[j + 1] > poly[j]:
            strict_crossing = True
        ratio = poly[j + 1] / poly[j]
        if ratio > best_ratio:
            best_ratio = ratio
            best_j = j
    return best_ratio, best_j, strict_crossing


def lc_bump_right_pressure(poly: list[int], start: int | None) -> tuple[float, int | None]:
    if start is None:
        return 0.0, None
    best_ratio = 0.0
    best_k: int | None = None
    for k in range(max(1, start), len(poly) - 1):
        if poly[k - 1] * poly[k + 1] <= poly[k] * poly[k]:
            continue
        if poly[k] == 0:
            continue
        ratio = poly[k + 1] / poly[k]
        if ratio > best_ratio:
            best_ratio = ratio
            best_k = k
    return best_ratio, best_k


def evaluate(adj: list[list[int]]) -> dict[str, Any]:
    n = len(adj)
    poly = independence_poly(n, adj)
    mode_first, mode_last = mode_interval(poly)
    descent = first_descent(poly)
    pressure, pressure_pos, crossing = max_post_descent_ratio(poly)
    lc_bump_ratio, lc_bump_pos = lc_bump_right_pressure(poly, descent)
    lc_ratio, lc_pos = log_concavity_ratio(poly)
    return {
        "n": n,
        "alpha": len(poly) - 1,
        "mode_first": mode_first,
        "mode_last": mode_last,
        "first_descent": descent,
        "crossing_pressure": pressure,
        "crossing_pos": pressure_pos,
        "crossing_reserve": 1.0 - pressure,
        "strict_crossing": crossing,
        "unimodal": is_unimodal(poly),
        "log_concave": is_log_concave(poly),
        "lc_ratio": lc_ratio,
        "lc_pos": lc_pos,
        "lc_bump_right_ratio": lc_bump_ratio,
        "lc_bump_pos": lc_bump_pos,
        "graph6": graph6_from_adj(adj),
    }


def selection_score(ev: dict[str, Any], prospect_pressure: float, prospect_bonus: float) -> float:
    return ev["crossing_pressure"] + prospect_bonus * max(
        0.0, prospect_pressure - ev["crossing_pressure"]
    )


def make_individual(
    adj: list[list[int]],
    origin: str,
    ev: dict[str, Any],
    generation: int,
    prospect_bonus: float,
    prospect_pressure: float | None = None,
) -> dict[str, Any]:
    if prospect_pressure is None:
        prospect_pressure = ev["crossing_pressure"]
    return {
        "adj": adj,
        "origin": origin,
        "generation": generation,
        "score": selection_score(ev, prospect_pressure, prospect_bonus),
        "prospect_pressure": prospect_pressure,
        "prospect_gain": max(0.0, prospect_pressure - ev["crossing_pressure"]),
        "fingerprint": _adj_fingerprint(adj),
        **ev,
    }


def population_key(candidate: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -candidate["score"],
        -candidate["crossing_pressure"],
        -candidate["lc_bump_right_ratio"],
        candidate["generation"],
        candidate["origin"],
    )


def archive_key(candidate: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -candidate["crossing_pressure"],
        -candidate["score"],
        -candidate["lc_bump_right_ratio"],
        candidate["generation"],
        candidate["origin"],
    )


def dedupe(candidates: list[dict[str, Any]], *, archive_mode: bool) -> list[dict[str, Any]]:
    best_by_fp: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        fingerprint = candidate["fingerprint"]
        incumbent = best_by_fp.get(fingerprint)
        if incumbent is None:
            best_by_fp[fingerprint] = candidate
            continue
        key = archive_key if archive_mode else population_key
        if key(candidate) < key(incumbent):
            best_by_fp[fingerprint] = candidate
    return list(best_by_fp.values())


def update_archive(
    archive: list[dict[str, Any]],
    candidates: list[dict[str, Any]],
    archive_size: int,
) -> list[dict[str, Any]]:
    merged = dedupe(archive + [deepcopy(candidate) for candidate in candidates], archive_mode=True)
    merged.sort(key=archive_key)
    return merged[:archive_size]


def random_tree(n: int, rng: random.Random) -> list[list[int]]:
    prufer = [rng.randrange(n) for _ in range(n - 2)]
    degree = [1] * n
    for v in prufer:
        degree[v] += 1
    adj: list[list[int]] = [[] for _ in range(n)]
    for v in prufer:
        leaf = min(i for i in range(n) if degree[i] == 1)
        adj[leaf].append(v)
        adj[v].append(leaf)
        degree[leaf] -= 1
        degree[v] -= 1
    leaves = [i for i in range(n) if degree[i] == 1]
    if len(leaves) == 2:
        u, v = leaves
        adj[u].append(v)
        adj[v].append(u)
    return adj


def compact(candidate: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "n",
        "alpha",
        "origin",
        "generation",
        "score",
        "crossing_pressure",
        "crossing_reserve",
        "crossing_pos",
        "strict_crossing",
        "unimodal",
        "log_concave",
        "lc_ratio",
        "lc_pos",
        "lc_bump_right_ratio",
        "lc_bump_pos",
        "mode_first",
        "first_descent",
        "prospect_pressure",
        "prospect_gain",
        "graph6",
    ]
    return {key: candidate[key] for key in keys}


def structure(adj: list[list[int]]) -> dict[str, Any]:
    degrees = _degree_sequence(adj)
    return {
        "leaves": sum(1 for value in degrees if value == 1),
        "max_degree": degrees[0] if degrees else 0,
        "branch_vertices": sum(1 for value in degrees if value >= 3),
        "degree_top5": degrees[:5],
        "signature": _tree_signature(adj),
    }


def run_for_n(
    *,
    n: int,
    pop_size: int,
    generations: int,
    elite_frac: float,
    mutations_per_ind: int,
    archive_size: int,
    prospect_bonus: float,
    random_inject_frac: float,
    seed: int,
    verbose_every: int,
) -> dict[str, Any]:
    rng = random.Random(seed)
    elite_count = max(2, int(pop_size * elite_frac))
    archive_size = max(2, archive_size)

    population: list[dict[str, Any]] = []
    for n0, adj, label in generate_seeds(n, pop_size, rng):
        if n0 != n:
            continue
        ev = evaluate(adj)
        population.append(make_individual(adj, label, ev, 0, prospect_bonus))
        if ev["strict_crossing"]:
            return counterexample_result(population[-1])

    while len(population) < pop_size:
        adj = random_tree(n, rng)
        ev = evaluate(adj)
        population.append(make_individual(adj, f"random_fill_{len(population)}", ev, 0, prospect_bonus))
        if ev["strict_crossing"]:
            return counterexample_result(population[-1])

    population = dedupe(population, archive_mode=False)
    population.sort(key=population_key)
    population = population[:pop_size]
    archive = update_archive([], population, archive_size)
    best = deepcopy(archive[0])
    history = []
    start = time.time()

    def record(gen: int) -> None:
        history.append(
            {
                "generation": gen,
                "best_pressure": population[0]["crossing_pressure"],
                "archive_best_pressure": archive[0]["crossing_pressure"],
                "archive_best_reserve": archive[0]["crossing_reserve"],
                "best_lc_bump_right_ratio": population[0]["lc_bump_right_ratio"],
                "archive_size": len(archive),
                "signature": _tree_signature(population[0]["adj"]),
            }
        )

    record(0)
    if verbose_every:
        print(
            f"n={n} gen=0 pressure={population[0]['crossing_pressure']:.12f} "
            f"reserve={population[0]['crossing_reserve']:.6g} "
            f"sig={_tree_signature(population[0]['adj'])}",
            flush=True,
        )

    for gen in range(1, generations + 1):
        elite = population[:elite_count]
        offspring: list[dict[str, Any]] = []
        refreshed_elite: list[dict[str, Any]] = []

        for parent in elite:
            best_child_pressure = parent["crossing_pressure"]
            for _ in range(mutations_per_ind):
                child_adj = deepcopy(parent["adj"])
                for _ in range(rng.randint(1, 3)):
                    _, child_adj = mutate(n, child_adj, rng)
                if not _validate_tree(n, child_adj):
                    continue
                ev = evaluate(child_adj)
                child = make_individual(
                    child_adj,
                    f"mut_g{gen}_{len(offspring)}",
                    ev,
                    gen,
                    prospect_bonus,
                )
                best_child_pressure = max(best_child_pressure, ev["crossing_pressure"])
                offspring.append(child)
                if ev["strict_crossing"]:
                    return counterexample_result(child)

            refreshed_elite.append(
                make_individual(
                    deepcopy(parent["adj"]),
                    parent["origin"],
                    parent,
                    parent["generation"],
                    prospect_bonus,
                    prospect_pressure=best_child_pressure,
                )
            )

        for _ in range(max(1, int(pop_size * random_inject_frac))):
            adj = random_tree(n, rng)
            ev = evaluate(adj)
            child = make_individual(adj, f"inject_g{gen}_{len(offspring)}", ev, gen, prospect_bonus)
            offspring.append(child)
            if ev["strict_crossing"]:
                return counterexample_result(child)

        archive = update_archive(archive, refreshed_elite + offspring, archive_size)
        merged = dedupe(refreshed_elite + offspring + archive, archive_mode=False)
        merged.sort(key=population_key)
        population = merged[:pop_size]
        if archive_key(archive[0]) < archive_key(best):
            best = deepcopy(archive[0])

        if gen == 1 or (verbose_every and gen % verbose_every == 0):
            record(gen)
            if verbose_every:
                print(
                    f"n={n} gen={gen} pressure={population[0]['crossing_pressure']:.12f} "
                    f"archive={archive[0]['crossing_pressure']:.12f} "
                    f"reserve={archive[0]['crossing_reserve']:.6g} "
                    f"lc_bump_right={population[0]['lc_bump_right_ratio']:.6g} "
                    f"sig={_tree_signature(population[0]['adj'])}",
                    flush=True,
                )

    elapsed = time.time() - start
    return {
        "counterexample_found": False,
        "n": n,
        "seed": seed,
        "generations": generations,
        "pop_size": pop_size,
        "best": compact(best),
        "best_structure": structure(best["adj"]),
        "archive": [compact(candidate) for candidate in archive],
        "history": history,
        "elapsed_s": round(elapsed, 2),
    }


def counterexample_result(candidate: dict[str, Any]) -> dict[str, Any]:
    return {
        "counterexample_found": True,
        "n": candidate["n"],
        "candidate": compact(candidate),
        "poly": independence_poly(candidate["n"], candidate["adj"]),
        "adj": candidate["adj"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-n", type=int, default=60)
    parser.add_argument("--max-n", type=int, default=180)
    parser.add_argument("--step-n", type=int, default=40)
    parser.add_argument("--pop-size", type=int, default=70)
    parser.add_argument("--generations", type=int, default=120)
    parser.add_argument("--elite-frac", type=float, default=0.25)
    parser.add_argument("--mutations-per-ind", type=int, default=3)
    parser.add_argument("--archive-size", type=int, default=16)
    parser.add_argument("--prospect-bonus", type=float, default=0.5)
    parser.add_argument("--random-inject-frac", type=float, default=0.12)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--verbose-every", type=int, default=25)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/crossing_pressure_optimizer_2026-07-03.json"),
    )
    args = parser.parse_args()

    if args.min_n < 2 or args.max_n < args.min_n or args.step_n < 1:
        parser.error("invalid n range")
    if args.pop_size < 4:
        parser.error("--pop-size must be at least 4")
    if args.generations < 0:
        parser.error("--generations must be nonnegative")
    if not 0.0 < args.elite_frac <= 1.0:
        parser.error("--elite-frac must be in (0, 1]")
    if args.mutations_per_ind < 1:
        parser.error("--mutations-per-ind must be positive")
    if args.archive_size < 1:
        parser.error("--archive-size must be positive")

    results = []
    found = False
    for n in range(args.min_n, args.max_n + 1, args.step_n):
        result = run_for_n(
            n=n,
            pop_size=args.pop_size,
            generations=args.generations,
            elite_frac=args.elite_frac,
            mutations_per_ind=args.mutations_per_ind,
            archive_size=args.archive_size,
            prospect_bonus=args.prospect_bonus,
            random_inject_frac=args.random_inject_frac,
            seed=args.seed + n,
            verbose_every=args.verbose_every,
        )
        results.append(result)
        if result["counterexample_found"]:
            found = True
            break

    summary = {
        "description": "Direct optimizer for max post-first-descent coefficient ratio",
        "parameters": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "step_n": args.step_n,
            "pop_size": args.pop_size,
            "generations": args.generations,
            "elite_frac": args.elite_frac,
            "mutations_per_ind": args.mutations_per_ind,
            "archive_size": args.archive_size,
            "prospect_bonus": args.prospect_bonus,
            "random_inject_frac": args.random_inject_frac,
            "seed": args.seed,
        },
        "counterexample_found": found,
        "best_overall": max(
            (result["best"] for result in results if not result["counterexample_found"]),
            key=lambda row: row["crossing_pressure"],
            default=None,
        ),
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.out}; counterexample_found={found}", flush=True)
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
