#!/usr/bin/env python3
"""Exact EADS search in the decorated-spider basin.

Delete all pendant leaves from the current two-certificate record and its core
is a three-arm spider.  Generic tree mutations change pendant multiplicities
one leaf at a time, although modal phases often require coordinated changes
of several leaves.  This optimizer works directly with the pendant load at
each core vertex and can also insert, delete, split, and merge core arms.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from scratch_eads_pareto_optimizer_20260808 import (
    encode,
    make_individual,
    scaled_total_key,
    select,
    total_key,
)


Spec = tuple[int, tuple[tuple[int, ...], ...]]


def canonical(center_load: int, arms: list[tuple[int, ...]]) -> Spec:
    if center_load < 0 or any(not arm for arm in arms):
        raise ValueError("invalid decorated-spider specification")
    return center_load, tuple(sorted(tuple(arm) for arm in arms))


def order(spec: Spec) -> int:
    center, arms = spec
    return 1 + sum(map(len, arms)) + center + sum(sum(arm) for arm in arms)


def build(spec: Spec) -> list[list[int]]:
    center_load, arms = spec
    core_order = 1 + sum(map(len, arms))
    adj: list[list[int]] = [[] for _ in range(core_order + center_load + sum(sum(a) for a in arms))]

    def edge(left: int, right: int) -> None:
        adj[left].append(right)
        adj[right].append(left)

    loads = [center_load]
    vertex = 1
    for arm in arms:
        previous = 0
        for load in arm:
            edge(previous, vertex)
            previous = vertex
            loads.append(load)
            vertex += 1
    leaf = core_order
    for support, load in enumerate(loads):
        for _ in range(load):
            edge(support, leaf)
            leaf += 1
    if leaf != len(adj):
        raise AssertionError((leaf, len(adj)))
    return adj


def seed_specs() -> list[Spec]:
    return [
        # The order-23 two-certificate record.
        canonical(0, [(3, 0, 4), (4,), (4, 1)]),
        # The earlier three-certificate caterpillar.
        canonical(1, [(4, 0, 3), (5, 6)]),
        # The first all-leaves-bad tree (order 14).
        canonical(0, [(3,), (3,), (0, 3)]),
        # The all-leaves-bad double broom at order 18.
        canonical(0, [(0, 0, 5), (0, 0, 6)]),
        # Raw-margin leaders from the first long decorated-spider run.
        canonical(1, [(5, 6, 0, 6), (5, 7), (6,), (7,), (8,), (8,)]),
        canonical(0, [(5, 6, 1, 5), (7,), (7,), (7,), (7, 6), (8,), (14,)]),
    ]


def random_spec(rng: random.Random, max_load: int) -> Spec:
    arm_count = rng.randint(2, 5)
    arms = [
        tuple(rng.randint(0, max_load) for _ in range(rng.randint(1, 5)))
        for _ in range(arm_count)
    ]
    return canonical(rng.randint(0, max_load), arms)


def mutate(spec: Spec, rng: random.Random, max_load: int) -> Spec:
    center, frozen_arms = spec
    arms = [list(arm) for arm in frozen_arms]
    # Locations use (-1, 0) for the centre and (arm, position) otherwise.
    locations = [(-1, 0)] + [
        (arm, position)
        for arm, values in enumerate(arms)
        for position in range(len(values))
    ]
    move = rng.random()
    if move < 0.34:
        arm, position = rng.choice(locations)
        delta = rng.choice((-4, -3, -2, -1, 1, 2, 3, 4))
        if arm < 0:
            center = min(max_load, max(0, center + delta))
        else:
            arms[arm][position] = min(
                max_load, max(0, arms[arm][position] + delta),
            )
    elif move < 0.58:
        donor_arm, donor_position = rng.choice(locations)
        receiver_arm, receiver_position = rng.choice(locations)
        if (donor_arm, donor_position) != (receiver_arm, receiver_position):
            donor_value = center if donor_arm < 0 else arms[donor_arm][donor_position]
            receiver_value = center if receiver_arm < 0 else arms[receiver_arm][receiver_position]
            if donor_value and receiver_value < max_load:
                amount = min(
                    rng.choice((1, 1, 2, 3)), donor_value,
                    max_load - receiver_value,
                )
                if donor_arm < 0:
                    center -= amount
                else:
                    arms[donor_arm][donor_position] -= amount
                if receiver_arm < 0:
                    center += amount
                else:
                    arms[receiver_arm][receiver_position] += amount
    elif move < 0.74:
        arm = rng.randrange(len(arms))
        position = rng.randrange(len(arms[arm]) + 1)
        arms[arm].insert(position, rng.randint(0, min(6, max_load)))
    elif move < 0.86 and any(len(arm) > 1 for arm in arms):
        arm = rng.choice([index for index, values in enumerate(arms) if len(values) > 1])
        position = rng.randrange(len(arms[arm]))
        load = arms[arm].pop(position)
        if load:
            target = min(position, len(arms[arm]) - 1)
            arms[arm][target] = min(max_load, arms[arm][target] + load)
    elif move < 0.94 and len(arms) < 7:
        arms.append([rng.randint(0, min(8, max_load))])
    elif len(arms) > 2:
        arm = rng.randrange(len(arms))
        removed = arms.pop(arm)
        if removed and sum(removed):
            target = rng.randrange(len(arms))
            arms[target][-1] = min(max_load, arms[target][-1] + sum(removed))
    return canonical(center, [tuple(arm) for arm in arms])


def individual(spec: Spec, label: str) -> dict[str, object]:
    item = make_individual(build(spec), label)
    item["spec"] = spec
    return item


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=10)
    parser.add_argument("--max-n", type=int, default=180)
    parser.add_argument("--max-load", type=int, default=30)
    parser.add_argument("--population", type=int, default=128)
    parser.add_argument("--generations", type=int, default=1200)
    parser.add_argument("--offspring", type=int, default=6)
    parser.add_argument("--random-seeds", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260808)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_decorated_spider_optimizer_20260808.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)
    specs = [spec for spec in seed_specs() if args.min_n <= order(spec) <= args.max_n]
    while len(specs) < len(seed_specs()) + args.random_seeds:
        spec = random_spec(rng, args.max_load)
        if args.min_n <= order(spec) <= args.max_n:
            specs.append(spec)
    population = [individual(spec, f"seed:{index}") for index, spec in enumerate(specs)]
    evaluated = len(population)
    population = select(population, args.population)
    archive = select(population, args.population)
    eads_counterexample = None
    main_counterexample = None

    for generation in range(args.generations + 1):
        leader = min(population, key=total_key)
        scaled_leader = min(population, key=scaled_total_key)
        leader_ev = leader["evaluation"]
        scaled_ev = scaled_leader["evaluation"]
        assert isinstance(leader_ev, dict)
        assert isinstance(scaled_ev, dict)
        if generation % 10 == 0 or eads_counterexample or main_counterexample:
            print(json.dumps({
                "event": "generation",
                "generation": generation,
                "evaluated": evaluated,
                "best_good": leader_ev["good_count"],
                "best_potentials": [
                    leader_ev["separation_potentials"][vertex]  # type: ignore[index]
                    for vertex in leader_ev["good_vertices"]  # type: ignore[union-attr]
                ],
                "best_spec": leader["spec"],
                "best_n": leader["n"],
                "scaled_best_good": scaled_ev["good_count"],
                "scaled_best_potentials": [
                    scaled_ev["separation_potentials"][vertex]  # type: ignore[index]
                    for vertex in scaled_ev["good_vertices"]  # type: ignore[union-attr]
                ],
                "scaled_best_spec": scaled_leader["spec"],
                "scaled_best_n": scaled_leader["n"],
            }), flush=True)
        if generation == args.generations or eads_counterexample or main_counterexample:
            break
        children = []
        for parent_index, parent in enumerate(population):
            parent_spec = parent["spec"]
            assert isinstance(parent_spec, tuple)
            for child_index in range(args.offspring):
                child_spec = mutate(parent_spec, rng, args.max_load)
                if not args.min_n <= order(child_spec) <= args.max_n:
                    continue
                child = individual(
                    child_spec, f"g{generation + 1}:p{parent_index}:c{child_index}",
                )
                evaluated += 1
                children.append(child)
                ev = child["evaluation"]
                assert isinstance(ev, dict)
                if not ev["unimodal"]:
                    main_counterexample = child
                    break
                if ev["good_count"] == 0:
                    eads_counterexample = child
                    break
            if eads_counterexample or main_counterexample:
                break
        population = select(population + children, args.population)
        archive = select(archive + children, args.population)

    leader = min(archive, key=total_key)
    scaled_leader = min(archive, key=scaled_total_key)
    result = {
        "claim_scope": "exact decorated-spider evolutionary computation; evidence only",
        "seed": args.seed,
        "evaluated": evaluated,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "eads_counterexample_found": eads_counterexample is not None,
        "erdos_993_counterexample": encode(main_counterexample) if main_counterexample else None,
        "eads_counterexample": encode(eads_counterexample) if eads_counterexample else None,
        "best": {"spec": leader["spec"], **encode(leader)},
        "scaled_best": {
            "spec": scaled_leader["spec"], **encode(scaled_leader),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "evaluated": evaluated,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "eads_counterexample_found": eads_counterexample is not None,
        "best_good": leader["evaluation"]["good_count"],  # type: ignore[index]
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if eads_counterexample or main_counterexample else 0)


if __name__ == "__main__":
    main()
