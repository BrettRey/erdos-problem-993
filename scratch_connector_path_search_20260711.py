"""Exact heuristic search for two decorated hubs joined by a path.

This is a disposable July 11 search driver.  Every objective comparison uses
the exact integer independence polynomial and only coefficients strictly after
the first descent.  A reported hit is rebuilt as an adjacency list and checked
independently with ``independence_poly``.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import dataclass

from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly, is_unimodal
from trees import trees_geng_raw


Poly = tuple[int, ...]


def mul(a: Poly, b: Poly) -> Poly:
    return tuple(_polymul(list(a), list(b)))


def add(*polys: Poly) -> Poly:
    out: list[int] = []
    for poly in polys:
        out = _polyadd(out, list(poly))
    return tuple(out)


def rooted_state(adj: list[list[int]], root: int) -> tuple[Poly, Poly]:
    """Return (full polynomial, polynomial conditional on root excluded)."""
    parent = [-1] * len(adj)
    order = [root]
    for vertex in order:
        for child in adj[vertex]:
            if child != parent[vertex]:
                parent[child] = vertex
                order.append(child)
    excluded: list[Poly | None] = [None] * len(adj)
    included: list[Poly | None] = [None] * len(adj)
    for vertex in reversed(order):
        zero: Poly = (1,)
        one: Poly = (1,)
        for child in adj[vertex]:
            if parent[child] == vertex:
                assert excluded[child] is not None and included[child] is not None
                zero = mul(zero, add(excluded[child], included[child]))
                one = mul(one, excluded[child])
        excluded[vertex] = zero
        included[vertex] = (0,) + one
    assert excluded[root] is not None and included[root] is not None
    return add(excluded[root], included[root]), excluded[root]


@dataclass(frozen=True)
class State:
    full: Poly
    excluded: Poly
    graph6: str
    root: int
    order: int


def state_catalog(max_order: int = 12) -> list[State]:
    states: list[State] = []
    seen: set[tuple[Poly, Poly]] = set()
    for order in range(1, max_order + 1):
        trees = [(1, [[]], b"@")] if order == 1 else trees_geng_raw(order)
        for _, adj, graph6 in trees:
            for root in range(order):
                full, excluded = rooted_state(adj, root)
                key = (full, excluded)
                if key not in seen:
                    seen.add(key)
                    states.append(
                        State(full, excluded, graph6.decode("ascii"), root, order)
                    )
    return states


def product(states: list[State], indices: list[int], excluded: bool) -> Poly:
    out: Poly = (1,)
    for index in indices:
        state = states[index]
        out = mul(out, state.excluded if excluded else state.full)
    return out


def connector_poly(
    states: list[State], left: list[int], right: list[int], distance: int
) -> Poly:
    """Decorated endpoint hubs at graph distance ``distance`` (>=1)."""
    left_full = product(states, left, False)
    left_excluded = product(states, left, True)
    right_full = product(states, right, False)
    right_excluded = product(states, right, True)

    # State at left endpoint, according as that hub is excluded/included.
    zero = left_full
    one = (0,) + left_excluded
    # Add each undecorated internal path vertex.
    for _ in range(distance - 1):
        zero, one = add(zero, one), (0,) + zero
    # Add the decorated right endpoint.
    return add(mul(add(zero, one), right_full), (0,) + mul(zero, right_excluded))


def tail_score(poly: Poly) -> tuple[int, int, int, int]:
    """Return numerator, denominator, edge, first-descent endpoint."""
    first = next(
        (index for index in range(1, len(poly)) if poly[index] < poly[index - 1]),
        -1,
    )
    if first < 0 or first == len(poly) - 1:
        return 0, 1, -1, first
    edge = first
    for index in range(first + 1, len(poly) - 1):
        if poly[index + 1] * poly[edge] > poly[edge + 1] * poly[index]:
            edge = index
    return poly[edge + 1], poly[edge], edge, first


def ratio(score: tuple[int, int, int, int]) -> float:
    return score[0] / score[1]


def better(
    first: tuple[int, int, int, int], second: tuple[int, int, int, int]
) -> bool:
    return first[0] * second[1] > second[0] * first[1]


def build_tree(
    states: list[State], left: list[int], right: list[int], distance: int
) -> list[list[int]]:
    """Build the actual simple tree represented by a connector configuration."""
    # Connector vertices are 0,...,distance.
    adj: list[list[int]] = [[] for _ in range(distance + 1)]
    for vertex in range(distance):
        adj[vertex].append(vertex + 1)
        adj[vertex + 1].append(vertex)
    for hub, indices in ((0, left), (distance, right)):
        for index in indices:
            state = states[index]
            _, branch = parse_graph6(state.graph6.encode("ascii"))
            offset = len(adj)
            adj.extend([] for _ in branch)
            for vertex, neighbors in enumerate(branch):
                adj[offset + vertex].extend(offset + neighbor for neighbor in neighbors)
            adj[hub].append(offset + state.root)
            adj[offset + state.root].append(hub)
    for neighbors in adj:
        neighbors.sort()
    return adj


def exact_hit(
    states: list[State], left: list[int], right: list[int], distance: int, poly: Poly
) -> dict[str, object] | None:
    if is_unimodal(list(poly)):
        return None
    adj = build_tree(states, left, right, distance)
    replay = tuple(independence_poly(len(adj), adj))
    assert replay == poly
    assert sum(map(len, adj)) == 2 * (len(adj) - 1)
    # Connectedness plus n-1 edges certifies acyclicity.
    reached = {0}
    queue = [0]
    for vertex in queue:
        for neighbor in adj[vertex]:
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    assert len(reached) == len(adj)
    return {
        "order": len(adj),
        "distance": distance,
        "left": [states[index].__dict__ for index in left],
        "right": [states[index].__dict__ for index in right],
        "adjacency": adj,
        "coefficients": list(poly),
        "tail_score": tail_score(poly),
    }


def coordinate_ascent(
    states: list[State],
    left: list[int],
    right: list[int],
    distance: int,
    candidate_indices: list[int],
    rounds: int = 8,
) -> tuple[list[int], list[int], Poly, tuple[int, int, int, int]]:
    poly = connector_poly(states, left, right, distance)
    score = tail_score(poly)
    for _ in range(rounds):
        changed = False
        positions = [(left, index) for index in range(len(left))] + [
            (right, index) for index in range(len(right))
        ]
        random.shuffle(positions)
        for side, position in positions:
            old_index = side[position]
            best_index, best_poly, best_score = old_index, poly, score
            for candidate in candidate_indices:
                side[position] = candidate
                trial_poly = connector_poly(states, left, right, distance)
                trial_score = tail_score(trial_poly)
                if trial_score[0] > trial_score[1]:
                    hit = exact_hit(states, left, right, distance, trial_poly)
                    assert hit is not None
                    print("HIT", json.dumps(hit))
                    raise SystemExit(0)
                if better(trial_score, best_score):
                    best_index, best_poly, best_score = candidate, trial_poly, trial_score
            side[position] = best_index
            if best_index != old_index:
                changed = True
                poly, score = best_poly, best_score
                print(
                    "ASCENT",
                    distance,
                    len(left),
                    len(right),
                    ratio(score),
                    score,
                    len(build_tree(states, left, right, distance)),
                    states[best_index].graph6,
                    states[best_index].root,
                    flush=True,
                )
        if not changed:
            break
    return left, right, poly, score


def evolutionary_search(
    states: list[State],
    distance: int,
    min_branches: int,
    max_branches: int,
    population_size: int,
    generations: int,
    seed: int,
) -> list[tuple[list[int], list[int], Poly, tuple[int, int, int, int]]]:
    rng = random.Random(seed)
    usable = list(range(1, len(states)))  # Exclude the one-vertex state.

    def random_config() -> tuple[list[int], list[int]]:
        total = rng.randint(min_branches, max_branches)
        left_count = rng.randint(1, total - 1)
        return (
            [rng.choice(usable) for _ in range(left_count)],
            [rng.choice(usable) for _ in range(total - left_count)],
        )

    population: list[tuple[list[int], list[int], Poly, tuple[int, int, int, int]]] = []
    for _ in range(population_size):
        left, right = random_config()
        poly = connector_poly(states, left, right, distance)
        population.append((left, right, poly, tail_score(poly)))

    for generation in range(generations):
        children = list(population)
        for _ in range(population_size * 8):
            parent = rng.choice(population)
            left, right = list(parent[0]), list(parent[1])
            mutation = rng.randrange(7)
            if mutation <= 2:
                side = rng.choice((left, right))
                side[rng.randrange(len(side))] = rng.choice(usable)
            elif mutation == 3 and len(left) + len(right) < max_branches:
                rng.choice((left, right)).append(rng.choice(usable))
            elif mutation == 4 and len(left) + len(right) > min_branches:
                side = rng.choice([part for part in (left, right) if len(part) > 1])
                del side[rng.randrange(len(side))]
            elif mutation == 5:
                source, target = (left, right) if rng.random() < 0.5 else (right, left)
                if len(source) > 1:
                    target.append(source.pop(rng.randrange(len(source))))
            else:
                # Multiplicity mutation: duplicate a currently useful branch.
                if len(left) + len(right) < max_branches:
                    side = rng.choice((left, right))
                    side.append(rng.choice(side))
            poly = connector_poly(states, left, right, distance)
            score = tail_score(poly)
            if score[0] > score[1]:
                hit = exact_hit(states, left, right, distance, poly)
                assert hit is not None
                print("HIT", json.dumps(hit))
                raise SystemExit(0)
            children.append((left, right, poly, score))
        # Deduplicate and retain both high score and structural diversity.
        unique: dict[tuple[tuple[int, ...], tuple[int, ...]], tuple] = {}
        for item in children:
            key = (tuple(sorted(item[0])), tuple(sorted(item[1])))
            old = unique.get(key)
            if old is None or better(item[3], old[3]):
                unique[key] = item
        population = sorted(
            unique.values(), key=lambda item: ratio(item[3]), reverse=True
        )[:population_size]
        if generation % 10 == 0:
            best = population[0]
            print(
                "GEN",
                generation,
                ratio(best[3]),
                best[3],
                len(best[0]),
                len(best[1]),
                len(build_tree(states, best[0], best[1], distance)),
                flush=True,
            )
    return population


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--distance", type=int, default=2)
    parser.add_argument("--min-branches", type=int, default=3)
    parser.add_argument("--max-branches", type=int, default=12)
    parser.add_argument("--population", type=int, default=160)
    parser.add_argument("--generations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--top-coordinate", type=int, default=8)
    args = parser.parse_args()
    states = state_catalog()
    print("CATALOG", len(states), flush=True)
    population = evolutionary_search(
        states,
        args.distance,
        args.min_branches,
        args.max_branches,
        args.population,
        args.generations,
        args.seed,
    )
    candidates = list(range(1, len(states)))
    best = population[0]
    for item in population[: args.top_coordinate]:
        result = coordinate_ascent(
            states, list(item[0]), list(item[1]), args.distance, candidates
        )
        if better(result[3], best[3]):
            best = result
    left, right, poly, score = best
    print(
        "FINAL",
        ratio(score),
        score,
        args.distance,
        [(states[index].graph6, states[index].root, states[index].order) for index in left],
        [(states[index].graph6, states[index].root, states[index].order) for index in right],
        list(poly),
    )


if __name__ == "__main__":
    main()
