#!/usr/bin/env python3
"""Multi-basin exact evolutionary search against EADS and Erdos #993.

The existential adjacent-mode split statement (EADS) asks whether every tree
has a vertex ``v`` for which the modal intervals of

    A_v = I(T-v),              B_v = x I(T-N[v])

are at distance at most one.  Earlier scalar searches found two different
adversarial basins: trees with only three EADS-good vertices, and trees with
no EADS-good leaves.  Optimizing only one statistic rapidly loses the other
basin.  This search retains simultaneous elites for total-good, leaf-good,
internal-good, and split-distance objectives.

All independence polynomials, modal intervals, and unimodality tests use exact
integers.  An EADS failure would refute only the proposed proof route.  A
nonunimodal independence polynomial would be a counterexample to Erdos #993.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from copy import deepcopy
from pathlib import Path
from typing import Callable

from nm_optimizer import _adj_fingerprint, _random_tree, mutate_leaf_move, mutate_spr
from scratch_existential_split_mode_20260807 import tree_vertex_splits


Evaluation = dict[str, object]
Individual = dict[str, object]


def is_unimodal(coefficients: list[int]) -> bool:
    """Return whether a finite sequence is weakly unimodal."""
    descending = False
    for left, right in zip(coefficients, coefficients[1:]):
        if right < left:
            descending = True
        elif descending and right > left:
            return False
    return True


def adjacency(n: int, edges: list[list[int]]) -> list[list[int]]:
    """Build an adjacency list from a JSON edge list."""
    adj: list[list[int]] = [[] for _ in range(n)]
    for left, right in edges:
        adj[left].append(right)
        adj[right].append(left)
    return adj


def load_seeds(paths: list[Path]) -> list[tuple[str, list[list[int]]]]:
    """Load both optimizer ``best`` records and exhaustive ``failures``."""
    seeds: list[tuple[str, list[list[int]]]] = []
    for path in paths:
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        records: list[dict[str, object]] = []
        best = data.get("best")
        if isinstance(best, dict):
            records.append(best)
        leaders = data.get("leaders")
        if isinstance(leaders, dict):
            records.extend(
                record for record in leaders.values()
                if isinstance(record, dict)
            )
        failures = data.get("failures")
        if isinstance(failures, list):
            records.extend(record for record in failures if isinstance(record, dict))
        for index, record in enumerate(records):
            edges = record.get("edges")
            n = record.get("n")
            if isinstance(n, int) and isinstance(edges, list):
                seeds.append((f"{path.stem}:{index}", adjacency(n, edges)))
    return seeds


def evaluate(adj: list[list[int]]) -> Evaluation:
    """Compute the exact polynomial and every EADS split statistic."""
    polynomial_tuple, splits = tree_vertex_splits(adj)
    polynomial = list(polynomial_tuple)
    distances = [split.distance for split in splits]
    leaves = [vertex for vertex, neighbors in enumerate(adj) if len(neighbors) == 1]
    leaf_set = set(leaves)
    good_vertices = [split.vertex for split in splits if split.distance <= 1]
    good_leaves = [vertex for vertex in good_vertices if vertex in leaf_set]
    good_internal = [vertex for vertex in good_vertices if vertex not in leaf_set]
    ordered = sorted(distances)
    potentials = [
        separation_potential(list(split.a_poly), list(split.b_poly))
        for split in splits
    ]
    return {
        "good_count": len(good_vertices),
        "good_leaf_count": len(good_leaves),
        "good_internal_count": len(good_internal),
        "good_vertices": good_vertices,
        "good_leaves": good_leaves,
        "good_internal": good_internal,
        "leaves": leaves,
        "distances": distances,
        "ordered_distances": ordered,
        "separation_potentials": potentials,
        "zero_count": sum(distance == 0 for distance in distances),
        "polynomial": polynomial,
        "unimodal": is_unimodal(polynomial),
    }


def log_ratio(numerator: int, denominator: int) -> float:
    """Stable log of a positive exact-integer ratio."""
    if numerator <= 0 or denominator <= 0:
        return float("-inf")
    return math.log(numerator) - math.log(denominator)


def oriented_separation_potential(left: list[int], right: list[int]) -> float:
    """Measure proximity to all modes of ``left`` lying two before ``right``.

    At a cut ``k``, compare the largest left mass at or before ``k`` with its
    later maximum, and the largest right mass at or after ``k+2`` with its
    earlier maximum.  Both log ratios are positive exactly when the two mode
    intervals lie strictly on the desired sides of the cut.  Their minimum is
    a continuous pressure toward that event; zero is the exact phase boundary.
    """
    size = max(len(left), len(right))
    a = left + [0] * (size - len(left))
    b = right + [0] * (size - len(right))
    prefix_a: list[int] = []
    running = 0
    for value in a:
        running = max(running, value)
        prefix_a.append(running)
    suffix_a = [0] * size
    running = 0
    for index in range(size - 1, -1, -1):
        running = max(running, a[index])
        suffix_a[index] = running
    prefix_b: list[int] = []
    running = 0
    for value in b:
        running = max(running, value)
        prefix_b.append(running)
    suffix_b = [0] * size
    running = 0
    for index in range(size - 1, -1, -1):
        running = max(running, b[index])
        suffix_b[index] = running
    best = float("-inf")
    for cut in range(size - 2):
        left_pressure = log_ratio(prefix_a[cut], suffix_a[cut + 1])
        right_pressure = log_ratio(suffix_b[cut + 2], prefix_b[cut + 1])
        best = max(best, min(left_pressure, right_pressure))
    return best


def separation_potential(a: list[int], b: list[int]) -> float:
    """Best log-pressure toward either orientation of a distance-two split."""
    return max(
        oriented_separation_potential(a, b),
        oriented_separation_potential(b, a),
    )


def total_key(individual: Individual) -> tuple[object, ...]:
    ev = individual["evaluation"]
    assert isinstance(ev, dict)
    ordered = ev["ordered_distances"]
    potentials = ev["separation_potentials"]
    good = ev["good_vertices"]
    assert isinstance(ordered, list) and isinstance(potentials, list)
    assert isinstance(good, list)
    good_potentials = sorted(float(potentials[int(vertex)]) for vertex in good)
    return (
        ev["good_count"], tuple(-value for value in good_potentials),
        ev["zero_count"],
        tuple(-int(value) for value in ordered), individual["n"],
    )


def leaf_key(individual: Individual) -> tuple[object, ...]:
    ev = individual["evaluation"]
    assert isinstance(ev, dict)
    distances = ev["distances"]
    leaves = ev["leaves"]
    potentials = ev["separation_potentials"]
    good_leaves = ev["good_leaves"]
    assert isinstance(distances, list) and isinstance(leaves, list)
    assert isinstance(potentials, list) and isinstance(good_leaves, list)
    profile = sorted(int(distances[int(leaf)]) for leaf in leaves)
    good_potentials = sorted(
        float(potentials[int(vertex)]) for vertex in good_leaves
    )
    return (
        ev["good_leaf_count"], tuple(-value for value in good_potentials),
        ev["good_count"], ev["zero_count"],
        tuple(-value for value in profile), individual["n"],
    )


def internal_key(individual: Individual) -> tuple[object, ...]:
    ev = individual["evaluation"]
    assert isinstance(ev, dict)
    distances = ev["distances"]
    leaves = set(ev["leaves"])
    potentials = ev["separation_potentials"]
    good_internal = ev["good_internal"]
    assert isinstance(distances, list) and isinstance(potentials, list)
    assert isinstance(good_internal, list)
    profile = sorted(
        int(value) for vertex, value in enumerate(distances) if vertex not in leaves
    )
    good_potentials = sorted(
        float(potentials[int(vertex)]) for vertex in good_internal
    )
    return (
        ev["good_internal_count"], tuple(-value for value in good_potentials),
        ev["good_count"], ev["zero_count"],
        tuple(-value for value in profile), individual["n"],
    )


def balanced_key(individual: Individual) -> tuple[object, ...]:
    ev = individual["evaluation"]
    assert isinstance(ev, dict)
    leaf = int(ev["good_leaf_count"])
    internal = int(ev["good_internal_count"])
    return (max(leaf, internal), leaf + internal, abs(leaf - internal), total_key(individual))


def scaled_total_key(individual: Individual) -> tuple[object, ...]:
    """Retain structural margin gains that are not just 1/n flattening."""
    ev = individual["evaluation"]
    assert isinstance(ev, dict)
    potentials = ev["separation_potentials"]
    good = ev["good_vertices"]
    assert isinstance(potentials, list) and isinstance(good, list)
    n = int(individual["n"])
    scaled = sorted(n * float(potentials[int(vertex)]) for vertex in good)
    return (ev["good_count"], tuple(-value for value in scaled), total_key(individual))


KEYS: tuple[Callable[[Individual], tuple[object, ...]], ...] = (
    total_key, scaled_total_key, leaf_key, internal_key, balanced_key,
)


def make_individual(adj: list[list[int]], label: str) -> Individual:
    return {
        "n": len(adj),
        "adj": adj,
        "label": label,
        "fingerprint": _adj_fingerprint(adj),
        "evaluation": evaluate(adj),
    }


def dedupe(individuals: list[Individual]) -> list[Individual]:
    """Deduplicate labeled adjacency lists, retaining the first copy."""
    seen: set[str] = set()
    unique: list[Individual] = []
    for individual in individuals:
        fingerprint = str(individual["fingerprint"])
        if fingerprint not in seen:
            seen.add(fingerprint)
            unique.append(individual)
    return unique


def select(individuals: list[Individual], size: int) -> list[Individual]:
    """Round-robin elites from complementary objective orderings."""
    unique = dedupe(individuals)
    rankings = [sorted(unique, key=key) for key in KEYS]
    selected: list[Individual] = []
    selected_ids: set[int] = set()
    rank = 0
    while len(selected) < min(size, len(unique)):
        added = False
        for ranking in rankings:
            while rank < len(ranking) and id(ranking[rank]) in selected_ids:
                # Do not advance a shared cursor here: each ranking can have a
                # selected item at this rank but a useful item immediately after.
                break
            for individual in ranking[rank:]:
                marker = id(individual)
                if marker not in selected_ids:
                    selected.append(individual)
                    selected_ids.add(marker)
                    added = True
                    break
                if len(selected) >= size:
                    break
            if len(selected) >= size:
                break
        if not added:
            break
        rank += 1
    return selected


def add_leaf(adj: list[list[int]], rng: random.Random, target: int | None = None) -> list[list[int]]:
    out = deepcopy(adj)
    if target is None:
        target = rng.randrange(len(out))
    leaf = len(out)
    out.append([target])
    out[target].append(leaf)
    return out


def delete_leaf(adj: list[list[int]], rng: random.Random) -> list[list[int]]:
    """Delete a leaf and compact vertex labels."""
    leaves = [vertex for vertex, neighbors in enumerate(adj) if len(neighbors) == 1]
    if len(adj) <= 3 or not leaves:
        return deepcopy(adj)
    removed = rng.choice(leaves)
    mapping = {}
    next_vertex = 0
    for old in range(len(adj)):
        if old != removed:
            mapping[old] = next_vertex
            next_vertex += 1
    edges = [
        [mapping[left], mapping[right]]
        for left, neighbors in enumerate(adj)
        for right in neighbors
        if left < right and left != removed and right != removed
    ]
    return adjacency(len(adj) - 1, edges)


def subdivide_edge(adj: list[list[int]], rng: random.Random) -> list[list[int]]:
    out = deepcopy(adj)
    edges = [
        (left, right)
        for left, neighbors in enumerate(out)
        for right in neighbors
        if left < right
    ]
    if not edges:
        return out
    left, right = rng.choice(edges)
    out[left].remove(right)
    out[right].remove(left)
    middle = len(out)
    out.append([left, right])
    out[left].append(middle)
    out[right].append(middle)
    return out


def guided_leaf_move(adj: list[list[int]], ev: Evaluation, rng: random.Random) -> list[list[int]]:
    """Move a leaf with preference for the current EADS-good region."""
    out = deepcopy(adj)
    leaves = list(ev["leaves"])
    good = list(ev["good_vertices"])
    good_leaves = list(ev["good_leaves"])
    if not leaves:
        return out
    leaf = rng.choice(good_leaves if good_leaves and rng.random() < 0.7 else leaves)
    old = out[leaf][0]
    candidates = [vertex for vertex in range(len(out)) if vertex not in (leaf, old)]
    if not candidates:
        return out
    if good and rng.random() < 0.5:
        preferred = [vertex for vertex in good if vertex not in (leaf, old)]
        target = rng.choice(preferred) if preferred else rng.choice(candidates)
    else:
        target = rng.choice(candidates)
    out[old].remove(leaf)
    out[leaf].remove(old)
    out[target].append(leaf)
    out[leaf].append(target)
    return out


def guided_spr(adj: list[list[int]], ev: Evaluation, rng: random.Random) -> list[list[int]]:
    """Prune a branch incident with a difficult EADS-good vertex."""
    out = deepcopy(adj)
    good = list(ev["good_vertices"])
    potentials = list(ev["separation_potentials"])
    if not good:
        _, return_value = mutate_spr(len(out), out, rng)
        return return_value
    # The smallest potential is furthest from crossing its modal boundary and
    # is therefore the bottleneck once only a few good vertices remain.
    weakest = min(float(potentials[int(vertex)]) for vertex in good)
    pool = [
        int(vertex) for vertex in good
        if float(potentials[int(vertex)]) <= weakest + 0.03
    ]
    pivot = rng.choice(pool)
    if not out[pivot]:
        return out
    branch_root = rng.choice(out[pivot])
    out[pivot].remove(branch_root)
    out[branch_root].remove(pivot)
    branch = {branch_root}
    stack = [branch_root]
    while stack:
        vertex = stack.pop()
        for neighbor in out[vertex]:
            if neighbor not in branch:
                branch.add(neighbor)
                stack.append(neighbor)
    targets = [vertex for vertex in range(len(out)) if vertex not in branch]
    if not targets:
        out[pivot].append(branch_root)
        out[branch_root].append(pivot)
        return out
    target = rng.choice(targets)
    out[target].append(branch_root)
    out[branch_root].append(target)
    return out


def mutate(individual: Individual, rng: random.Random, min_n: int, max_n: int) -> list[list[int]]:
    adj = individual["adj"]
    ev = individual["evaluation"]
    assert isinstance(adj, list) and isinstance(ev, dict)
    move = rng.random()
    if move < 0.22:
        _, out = mutate_leaf_move(len(adj), adj, rng)
    elif move < 0.42:
        _, out = mutate_spr(len(adj), adj, rng)
    elif move < 0.62:
        out = guided_spr(adj, ev, rng)
    elif move < 0.78:
        out = guided_leaf_move(adj, ev, rng)
    elif move < 0.87 and len(adj) < max_n:
        good = list(ev["good_vertices"])
        target = rng.choice(good) if good and rng.random() < 0.65 else None
        out = add_leaf(adj, rng, target)
    elif move < 0.94 and len(adj) > min_n:
        out = delete_leaf(adj, rng)
    elif len(adj) < max_n:
        out = subdivide_edge(adj, rng)
    else:
        _, out = mutate_spr(len(adj), adj, rng)
    return out


def edges(adj: list[list[int]]) -> list[list[int]]:
    return [
        [left, right]
        for left, neighbors in enumerate(adj)
        for right in neighbors
        if left < right
    ]


def encode(individual: Individual) -> dict[str, object]:
    adj = individual["adj"]
    ev = individual["evaluation"]
    assert isinstance(adj, list) and isinstance(ev, dict)
    return {
        "label": individual["label"],
        "n": len(adj),
        "edges": edges(adj),
        **ev,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=10)
    parser.add_argument("--max-n", type=int, default=100)
    parser.add_argument("--population", type=int, default=96)
    parser.add_argument("--generations", type=int, default=500)
    parser.add_argument("--offspring", type=int, default=4)
    parser.add_argument("--random-seeds", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260808)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_pareto_optimizer_20260808.json"),
    )
    args = parser.parse_args()
    if args.min_n < 2 or args.max_n < args.min_n:
        parser.error("require 2 <= min-n <= max-n")
    rng = random.Random(args.seed)
    seed_paths = [
        Path("results/eads_seeded_optimizer_20260807.json"),
        Path("results/eads_seeded_n60_optimizer_20260807.json"),
        Path("results/eads_leaf_n18_20260808.json"),
        Path("results/eads_leaf_n19_20260808.json"),
        Path("results/eads_leaf_n20_20260808.json"),
        Path("results/eads_center_ball_r3_optimizer_smoke_20260808.json"),
        Path("results/eads_pareto_optimizer_20260808.json"),
    ]
    loaded = [
        (label, adj) for label, adj in load_seeds(seed_paths)
        if args.min_n <= len(adj) <= args.max_n
    ]
    population = [make_individual(adj, label) for label, adj in loaded]
    for index in range(args.random_seeds):
        n = rng.randint(args.min_n, args.max_n)
        population.append(make_individual(_random_tree(n, rng), f"random:{index}"))
    evaluated = len(population)
    population = select(population, args.population)
    archive = select(population, max(args.population, 32))
    main_counterexample: Individual | None = next(
        (item for item in population if not item["evaluation"]["unimodal"]),  # type: ignore[index]
        None,
    )
    eads_counterexample: Individual | None = next(
        (item for item in population if item["evaluation"]["good_count"] == 0),  # type: ignore[index]
        None,
    )

    for generation in range(args.generations + 1):
        if generation % 10 == 0 or eads_counterexample or main_counterexample:
            leaders = [min(population, key=key) for key in KEYS]
            print(json.dumps({
                "event": "generation",
                "generation": generation,
                "evaluated": evaluated,
                "population": len(population),
                "best_total_good": leaders[0]["evaluation"]["good_count"],
                "best_scaled_total_good": leaders[1]["evaluation"]["good_count"],
                "best_leaf_good": leaders[2]["evaluation"]["good_leaf_count"],
                "best_internal_good": leaders[3]["evaluation"]["good_internal_count"],
                "best_balanced": [
                    leaders[4]["evaluation"]["good_leaf_count"],
                    leaders[4]["evaluation"]["good_internal_count"],
                ],
            }), flush=True)
        if eads_counterexample or main_counterexample or generation == args.generations:
            break

        children: list[Individual] = []
        for parent_index, parent in enumerate(population):
            for child_index in range(args.offspring):
                child_adj = mutate(parent, rng, args.min_n, args.max_n)
                # Occasional compound mutations help cross the phase boundaries
                # that trap one-step scalar hill climbers.
                if rng.random() < 0.12:
                    temporary = {
                        "adj": child_adj,
                        "evaluation": evaluate(child_adj),
                    }
                    evaluated += 1
                    child_adj = mutate(temporary, rng, args.min_n, args.max_n)
                child = make_individual(
                    child_adj, f"g{generation + 1}:p{parent_index}:c{child_index}",
                )
                evaluated += 1
                children.append(child)
                child_ev = child["evaluation"]
                assert isinstance(child_ev, dict)
                if not child_ev["unimodal"]:
                    main_counterexample = child
                    break
                if child_ev["good_count"] == 0:
                    eads_counterexample = child
                    break
            if eads_counterexample or main_counterexample:
                break
        combined = population + children
        population = select(combined, args.population)
        archive = select(archive + children, max(args.population, 32))

    total_best = min(archive, key=total_key)
    scaled_total_best = min(archive, key=scaled_total_key)
    leaf_best = min(archive, key=leaf_key)
    internal_best = min(archive, key=internal_key)
    balanced_best = min(archive, key=balanced_key)
    result = {
        "claim_scope": "exact evolutionary computation; evidence only",
        "seed": args.seed,
        "evaluated": evaluated,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "eads_counterexample_found": eads_counterexample is not None,
        "erdos_993_counterexample": (
            encode(main_counterexample) if main_counterexample else None
        ),
        "eads_counterexample": (
            encode(eads_counterexample) if eads_counterexample else None
        ),
        "leaders": {
            "total_good": encode(total_best),
            "scaled_total_good": encode(scaled_total_best),
            "leaf_good": encode(leaf_best),
            "internal_good": encode(internal_best),
            "balanced": encode(balanced_best),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "evaluated": evaluated,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "eads_counterexample_found": eads_counterexample is not None,
        "best_total_good": total_best["evaluation"]["good_count"],  # type: ignore[index]
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if main_counterexample or eads_counterexample else 0)


if __name__ == "__main__":
    main()
