#!/usr/bin/env python3
"""Evolve known tree LC failures toward a genuine unimodality failure.

The usual near-miss ratio can reward an extremely broad but log-concave peak.
Here the primary objective is stricter.  For adjacent ratios

    r_k = i_{k+1} / i_k,

start at the first descent and retain an edge only if ``r_k`` is larger than
an earlier tail ratio.  Such an edge is a genuine log-concavity rebound.  Its
absolute ratio must cross one to disprove unimodality.

All 21 exhaustive order-26/order-28 LC failures seed the search.  They are
grown to the requested order by exact tree-preserving leaf additions and edge
subdivisions, then evolved with leaf moves, SPR, pendant concentration, and
degree-2 relocation.  Integer arithmetic decides every rebound and every
counterexample; floats are used only to rank the evolutionary population.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from copy import deepcopy
from pathlib import Path

from graph6 import parse_graph6
from indpoly import independence_poly, is_unimodal, near_miss_ratio
from nm_optimizer import (
    _adj_fingerprint,
    _random_tree,
    _validate_tree,
    mutate,
)


def load_lc_bases() -> list[tuple[str, list[list[int]]]]:
    sources = [
        Path("results/analysis_n26.json"),
        Path("results/analysis_n28_modal_lc_nm.json"),
    ]
    out: list[tuple[str, list[list[int]]]] = []
    for source in sources:
        data = json.loads(source.read_text())
        records = data.get("lc_failures") or data.get("top_lc_failures") or []
        for index, record in enumerate(records):
            n, adj = parse_graph6(record["graph6"].encode("ascii"))
            replay = independence_poly(n, adj)
            if record.get("poly") is not None and replay != record["poly"]:
                raise AssertionError(f"base replay mismatch in {source}")
            out.append((f"{source.name}:{index}", adj))
    return out


def add_leaf(adj: list[list[int]], support: int) -> list[list[int]]:
    out = deepcopy(adj)
    leaf = len(out)
    out.append([support])
    out[support].append(leaf)
    return out


def subdivide_edge(
    adj: list[list[int]], left: int, right: int,
) -> list[list[int]]:
    out = deepcopy(adj)
    vertex = len(out)
    out.append([left, right])
    out[left].remove(right)
    out[right].remove(left)
    out[left].append(vertex)
    out[right].append(vertex)
    return out


def grow_to_order(
    adj: list[list[int]], target: int, rng: random.Random,
) -> list[list[int]]:
    out = deepcopy(adj)
    while len(out) < target:
        if len(out) > 1 and rng.random() < 0.55:
            edges = [
                (left, right)
                for left, neighbors in enumerate(out)
                for right in neighbors if left < right
            ]
            left, right = rng.choice(edges)
            out = subdivide_edge(out, left, right)
        else:
            # Bias some additions toward existing branch/support vertices but
            # retain uniform choices for structural diversity.
            if rng.random() < 0.6:
                weights = [max(1, len(neighbors)) for neighbors in out]
                support = rng.choices(range(len(out)), weights=weights, k=1)[0]
            else:
                support = rng.randrange(len(out))
            out = add_leaf(out, support)
    return out


def ratio_metrics(poly: list[int]) -> dict[str, object]:
    first = next(
        (edge for edge in range(len(poly) - 1) if poly[edge] > poly[edge + 1]),
        -1,
    )
    if first < 0:
        return {
            "first_descent": -1,
            "rebound_k": -1,
            "rebound_num": 0,
            "rebound_den": 1,
            "rebound_ratio": 0.0,
            "rebound_gain": 0.0,
            "minimum_slope_k": -1,
        }

    minimum_num = poly[first + 1]
    minimum_den = poly[first]
    minimum_k = first
    best_num, best_den = 0, 1
    best_gain_num, best_gain_den = 0, 1
    best_k = -1
    best_minimum_k = -1
    for edge in range(first + 1, len(poly) - 1):
        numerator = poly[edge + 1]
        denominator = poly[edge]
        if numerator * minimum_den > minimum_num * denominator:
            gain_num = numerator * minimum_den
            gain_den = denominator * minimum_num
            better_ratio = numerator * best_den > best_num * denominator
            tie_ratio = numerator * best_den == best_num * denominator
            better_gain = gain_num * best_gain_den > best_gain_num * gain_den
            if better_ratio or (tie_ratio and better_gain):
                best_num, best_den = numerator, denominator
                best_gain_num, best_gain_den = gain_num, gain_den
                best_k = edge
                best_minimum_k = minimum_k
        if numerator * minimum_den < minimum_num * denominator:
            minimum_num, minimum_den = numerator, denominator
            minimum_k = edge
    return {
        "first_descent": first,
        "rebound_k": best_k,
        "rebound_num": best_num,
        "rebound_den": best_den,
        "rebound_ratio": best_num / best_den,
        "rebound_gain": best_gain_num / best_gain_den,
        "minimum_slope_k": best_minimum_k,
    }


def evaluate(adj: list[list[int]], label: str, generation: int) -> dict[str, object]:
    poly = independence_poly(len(adj), adj)
    metrics = ratio_metrics(poly)
    nm, nm_k = near_miss_ratio(poly)
    return {
        "adj": adj,
        "label": label,
        "generation": generation,
        "poly": poly,
        "alpha": len(poly) - 1,
        "unimodal": is_unimodal(poly),
        "near_miss_ratio": nm,
        "near_miss_k": nm_k,
        "fingerprint": _adj_fingerprint(adj),
        **metrics,
    }


def dedupe(candidates: list[dict[str, object]]) -> list[dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    for candidate in candidates:
        fingerprint = str(candidate["fingerprint"])
        old = best.get(fingerprint)
        if old is None or primary_key(candidate) < primary_key(old):
            best[fingerprint] = candidate
    return list(best.values())


def primary_key(candidate: dict[str, object]) -> tuple[float, float, float, int]:
    return (
        -float(candidate["rebound_ratio"]),
        -float(candidate["rebound_gain"]),
        -float(candidate["near_miss_ratio"]),
        int(candidate["generation"]),
    )


def nm_key(candidate: dict[str, object]) -> tuple[float, float, int]:
    return (
        -float(candidate["near_miss_ratio"]),
        -float(candidate["rebound_ratio"]),
        int(candidate["generation"]),
    )


def gain_key(candidate: dict[str, object]) -> tuple[float, float, int]:
    return (
        -float(candidate["rebound_gain"]),
        -float(candidate["rebound_ratio"]),
        int(candidate["generation"]),
    )


def select(candidates: list[dict[str, object]], population_size: int) -> list[dict[str, object]]:
    candidates = dedupe(candidates)
    quota = max(1, population_size // 3)
    merged: dict[str, dict[str, object]] = {}
    for key in (primary_key, nm_key, gain_key):
        for candidate in sorted(candidates, key=key)[:quota]:
            merged[str(candidate["fingerprint"])] = candidate
    if len(merged) < population_size:
        for candidate in sorted(candidates, key=primary_key):
            merged.setdefault(str(candidate["fingerprint"]), candidate)
            if len(merged) >= population_size:
                break
    return sorted(merged.values(), key=primary_key)[:population_size]


def encode(candidate: dict[str, object], include_poly: bool) -> dict[str, object]:
    adj = candidate["adj"]
    assert isinstance(adj, list)
    return {
        key: value
        for key, value in candidate.items()
        if key not in {"adj", "fingerprint", "poly"}
    } | {
        "n": len(adj),
        "edges": [
            [left, right]
            for left, neighbors in enumerate(adj)
            for right in neighbors if left < right
        ],
        "polynomial": candidate["poly"] if include_poly else None,
    }


def run(
    target: int,
    population_size: int,
    generations: int,
    offspring_per_parent: int,
    seed: int,
) -> dict[str, object]:
    rng = random.Random(seed)
    bases = load_lc_bases()
    population: list[dict[str, object]] = []
    expansions_per_base = max(2, population_size // len(bases) + 1)
    for base_label, base_adj in bases:
        if len(base_adj) > target:
            continue
        for expansion in range(expansions_per_base):
            adj = grow_to_order(base_adj, target, rng)
            population.append(evaluate(
                adj, f"{base_label}:grow{expansion}", 0,
            ))
    while len(population) < population_size:
        adj = _random_tree(target, rng)
        population.append(evaluate(adj, f"random{len(population)}", 0))
    population = select(population, population_size)
    evaluated = len(population)
    history: list[dict[str, object]] = []
    best = min(population, key=primary_key)

    for generation in range(generations + 1):
        if generation == 0 or generation % 25 == 0:
            leader = min(population, key=primary_key)
            row = {
                "generation": generation,
                "evaluated": evaluated,
                "rebound_ratio": leader["rebound_ratio"],
                "rebound_gain": leader["rebound_gain"],
                "rebound_k": leader["rebound_k"],
                "first_descent": leader["first_descent"],
                "near_miss_ratio": leader["near_miss_ratio"],
                "unimodal": leader["unimodal"],
            }
            history.append(row)
            print(json.dumps({"event": "generation", "n": target, **row}), flush=True)
        counterexample = next(
            (candidate for candidate in population if not candidate["unimodal"]),
            None,
        )
        if counterexample is not None:
            return {
                "counterexample_found": True,
                "n": target,
                "evaluated": evaluated,
                "seed": seed,
                "history": history,
                "counterexample": encode(counterexample, True),
            }
        leader = min(population, key=primary_key)
        if primary_key(leader) < primary_key(best):
            best = deepcopy(leader)
        if generation == generations:
            break

        parents = select(population, max(12, population_size // 3))
        offspring: list[dict[str, object]] = []
        for parent_index, parent in enumerate(parents):
            parent_adj = parent["adj"]
            assert isinstance(parent_adj, list)
            for child_index in range(offspring_per_parent):
                adj = deepcopy(parent_adj)
                for _ in range(rng.randint(1, 4)):
                    _, adj = mutate(target, adj, rng)
                if not _validate_tree(target, adj):
                    continue
                child = evaluate(
                    adj,
                    f"g{generation + 1}_p{parent_index}_c{child_index}",
                    generation + 1,
                )
                evaluated += 1
                if not child["unimodal"]:
                    return {
                        "counterexample_found": True,
                        "n": target,
                        "evaluated": evaluated,
                        "seed": seed,
                        "history": history,
                        "counterexample": encode(child, True),
                    }
                offspring.append(child)
        # Freshly grown bases prevent the population from becoming trapped in
        # one labelled basin.
        for injection in range(max(2, population_size // 20)):
            base_label, base_adj = rng.choice(bases)
            adj = grow_to_order(base_adj, target, rng)
            offspring.append(evaluate(
                adj, f"inject_g{generation + 1}_{base_label}_{injection}",
                generation + 1,
            ))
            evaluated += 1
        population = select(population + offspring, population_size)

    return {
        "counterexample_found": False,
        "n": target,
        "evaluated": evaluated,
        "seed": seed,
        "history": history,
        "best": encode(best, False),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=30)
    parser.add_argument("--max-n", type=int, default=100)
    parser.add_argument("--step-n", type=int, default=10)
    parser.add_argument("--population", type=int, default=180)
    parser.add_argument("--generations", type=int, default=500)
    parser.add_argument("--offspring-per-parent", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260808)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/rebound_optimizer_20260808.json"),
    )
    args = parser.parse_args()
    started = time.time()
    results = []
    for target in range(args.min_n, args.max_n + 1, args.step_n):
        result = run(
            target,
            args.population,
            args.generations,
            args.offspring_per_parent,
            args.seed + target,
        )
        results.append(result)
        if result["counterexample_found"]:
            break
    payload = {
        "claim_scope": "exact evolutionary search; evidence only unless a counterexample is present",
        "objective": "largest later adjacent ratio that exceeds an earlier tail ratio",
        "parameters": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "step_n": args.step_n,
            "population": args.population,
            "generations": args.generations,
            "offspring_per_parent": args.offspring_per_parent,
            "seed": args.seed,
        },
        "counterexample_found": any(row["counterexample_found"] for row in results),
        "results": results,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "counterexample_found": payload["counterexample_found"],
        "out": str(args.out),
        "elapsed_seconds": payload["elapsed_seconds"],
    }), flush=True)
    raise SystemExit(1 if payload["counterexample_found"] else 0)


if __name__ == "__main__":
    main()
