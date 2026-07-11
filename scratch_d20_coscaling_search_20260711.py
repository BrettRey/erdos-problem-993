#!/usr/bin/env python3
"""One-pivot D20 search: co-scale TG defects with nested-star plateaux.

The fixed-species grammar is eventually smoothed by a binomial factor.  This
driver lets the Bautista--Ramos parameter t grow with the nested-star leaf
count s, so the hard polynomial is no longer a fixed core.

For G=T_(3,t,1), TG=TG_(m,t), and a half containing a TG copies plus b
nested s-leaf stars,

    E_TG=(1+x) P_G^m,             S_TG=x E_G^m,
    E_half=(E_TG+S_TG)^a ((1+x)^s+x)^b,
    S_half=x E_TG^a (1+x)^(sb).

Two halves are joined by one edge.  Floats rank; every alleged full or
selected-state crossing is replayed exactly, with treehood for a true hit.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import random
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.special import gammaln
from sympy.polys.domains import ZZ
from sympy.polys.rings import ring

from scratch_d17_checkpointed_search_20260711 import (
    append_record,
    exact_valley,
    strict_post_descent_pressure,
    utc_now,
)
from scratch_d17_decorated_path_search_20260711 import ratio_bump_after_descent
from scratch_d17_spherical_search_20260711 import (
    Scaled,
    add,
    multiply,
    normalize,
    power,
    spherical_states,
)
from scratch_d20_state_grammar_search_20260711 import (
    Feature,
    feature,
    scaled_log_mass,
    shift_x,
)
from scripts.analyze_prufer_corpus import make_bautista_ramos_tree


@dataclass(frozen=True, order=True)
class CoSpec:
    m: int
    t: int
    tg_count: int
    nested_star_leaves: int
    nested_star_count: int


@dataclass(frozen=True)
class CoCandidate:
    score: float
    left: CoSpec
    right: CoSpec
    order: int
    first_descent_at: int
    first_descent_ratio: float
    best_later_at: int
    best_later_ratio: float
    left_selected: Feature
    right_selected: Feature
    left_log_selected_to_excluded_mass: float
    right_log_selected_to_excluded_mass: float


def tg_order(m: int, t: int) -> int:
    return 2 + m * (4 + 6 * t)


def half_order(spec: CoSpec) -> int:
    return (
        1
        + spec.tg_count * tg_order(spec.m, spec.t)
        + spec.nested_star_count * (spec.nested_star_leaves + 1)
    )


def canonical_pair(left: CoSpec, right: CoSpec) -> tuple[CoSpec, CoSpec]:
    return (left, right) if left <= right else (right, left)


class CoBuilder:
    def __init__(self) -> None:
        self.binomial_cache: dict[int, Scaled] = {}
        self.tg_cache: dict[tuple[int, int], tuple[Scaled, Scaled]] = {}
        self.half_cache: dict[CoSpec, tuple[Scaled, Scaled]] = {}

    def binomial(self, exponent: int) -> Scaled:
        if exponent not in self.binomial_cache:
            if exponent == 0:
                result = Scaled(np.array([1.0]), 0.0)
            else:
                k = np.arange(exponent + 1, dtype=float)
                logs = (
                    gammaln(exponent + 1.0)
                    - gammaln(k + 1.0)
                    - gammaln(exponent - k + 1.0)
                    - exponent * math.log(2.0)
                )
                values = np.exp(logs)
                values /= values.sum()
                result = normalize(Scaled(values, exponent * math.log(2.0)))
            self.binomial_cache[exponent] = result
        return self.binomial_cache[exponent]

    def tg_states(self, m: int, t: int) -> tuple[Scaled, Scaled]:
        key = (m, t)
        if key not in self.tg_cache:
            g_excluded, g_selected, _ = spherical_states((1, t, 3))
            g_total = add(g_excluded, g_selected)
            tg_excluded = multiply(self.binomial(1), power(g_total, m))
            tg_selected = shift_x(power(g_excluded, m))
            self.tg_cache[key] = (tg_excluded, tg_selected)
        return self.tg_cache[key]

    def half_states(self, spec: CoSpec) -> tuple[Scaled, Scaled]:
        if spec not in self.half_cache:
            tg_excluded, tg_selected = self.tg_states(spec.m, spec.t)
            tg_total = add(tg_excluded, tg_selected)
            star_excluded = self.binomial(spec.nested_star_leaves)
            star_total = add(star_excluded, Scaled(np.array([0.0, 1.0]), 0.0))
            excluded = multiply(
                power(tg_total, spec.tg_count),
                power(star_total, spec.nested_star_count),
            )
            selected_base = multiply(
                power(tg_excluded, spec.tg_count),
                self.binomial(
                    spec.nested_star_leaves * spec.nested_star_count
                ),
            )
            self.half_cache[spec] = (excluded, shift_x(selected_base))
        return self.half_cache[spec]

    def joined(self, left: CoSpec, right: CoSpec) -> np.ndarray:
        left_e, left_s = self.half_states(left)
        right_e, right_s = self.half_states(right)
        full = add(
            multiply(left_e, right_e),
            add(multiply(left_s, right_e), multiply(left_e, right_s)),
        )
        values = np.maximum(full.values, 0.0)
        values /= values.sum()
        return values


class CoEvaluator:
    def __init__(self, args, builder: CoBuilder) -> None:
        self.args = args
        self.builder = builder
        self.cache: dict[tuple[CoSpec, CoSpec], CoCandidate | None] = {}

    def evaluate(self, left: CoSpec, right: CoSpec) -> CoCandidate | None:
        left, right = canonical_pair(left, right)
        key = (left, right)
        if key in self.cache:
            return self.cache[key]
        values = self.builder.joined(left, right)
        pressure = strict_post_descent_pressure(
            values,
            relative_floor=self.args.relative_floor,
            descent_tolerance=self.args.descent_tolerance,
            min_separation=self.args.min_separation,
        )
        if pressure is None:
            self.cache[key] = None
            return None
        left_e, left_s = self.builder.half_states(left)
        right_e, right_s = self.builder.half_states(right)
        candidate = CoCandidate(
            score=float(pressure["best_later_ratio"]),
            left=left,
            right=right,
            order=half_order(left) + half_order(right),
            first_descent_at=int(pressure["first_descent_at"]),
            first_descent_ratio=float(pressure["first_descent_ratio"]),
            best_later_at=int(pressure["best_later_at"]),
            best_later_ratio=float(pressure["best_later_ratio"]),
            left_selected=feature(
                left_s,
                relative_floor=self.args.relative_floor,
                descent_tolerance=self.args.descent_tolerance,
                min_separation=self.args.min_separation,
            ),
            right_selected=feature(
                right_s,
                relative_floor=self.args.relative_floor,
                descent_tolerance=self.args.descent_tolerance,
                min_separation=self.args.min_separation,
            ),
            left_log_selected_to_excluded_mass=(
                scaled_log_mass(left_s) - scaled_log_mass(left_e)
            ),
            right_log_selected_to_excluded_mass=(
                scaled_log_mass(right_s) - scaled_log_mass(right_e)
            ),
        )
        self.cache[key] = candidate
        return candidate


def selected_objective(candidate: CoCandidate) -> tuple[float, float]:
    return (
        max(
            candidate.left_selected.best_later_ratio,
            candidate.right_selected.best_later_ratio,
        ),
        max(
            candidate.left_selected.ratio_bump,
            candidate.right_selected.ratio_bump,
        ),
    )


def choose_frontier(candidates: list[CoCandidate], population: int) -> list[CoCandidate]:
    selected: dict[tuple[CoSpec, CoSpec], CoCandidate] = {}
    objectives = (
        lambda item: item.score,
        lambda item: selected_objective(item)[0],
        lambda item: selected_objective(item)[1],
    )
    quotas = (population // 2, population // 4, population // 4)
    for objective, quota in zip(objectives, quotas):
        for item in sorted(candidates, key=objective, reverse=True)[: max(1, quota)]:
            selected[(item.left, item.right)] = item
    return list(selected.values())[:population]


def co_scaled_s(rng: random.Random, t: int, maximum: int) -> int:
    multiples = (1, 2, 4, 8, 16, 32, 64, 128)
    if rng.random() < 0.8:
        return max(1, min(maximum, t * rng.choice(multiples)))
    return rng.randint(1, maximum)


def random_spec(args, rng: random.Random, cap: int) -> CoSpec:
    for _ in range(10_000):
        t = rng.randint(args.min_t, args.max_t)
        spec = CoSpec(
            m=rng.randint(2, args.max_m),
            t=t,
            tg_count=rng.randint(1, args.max_tg_count),
            nested_star_leaves=co_scaled_s(
                rng, t, args.max_plateau_leaves
            ),
            nested_star_count=rng.randint(1, args.max_nested_stars),
        )
        if half_order(spec) <= cap:
            return spec
    raise RuntimeError("unable to sample co-scaling spec")


def mutate_spec(args, spec: CoSpec, rng: random.Random, cap: int) -> CoSpec:
    for _ in range(1_000):
        values = asdict(spec)
        field = rng.choice(tuple(values))
        if field == "m":
            values[field] = max(2, min(args.max_m, values[field] + rng.choice((-2, -1, 1, 2))))
        elif field == "t":
            values[field] = max(
                args.min_t,
                min(args.max_t, values[field] + rng.choice((-32, -8, -2, -1, 1, 2, 8, 32))),
            )
        elif field == "tg_count":
            values[field] = rng.randint(1, args.max_tg_count)
        elif field == "nested_star_count":
            values[field] = rng.randint(1, args.max_nested_stars)
        else:
            old = values[field]
            values[field] = max(
                1,
                min(
                    args.max_plateau_leaves,
                    old + rng.choice((-1024, -256, -64, -16, -1, 1, 16, 64, 256, 1024)),
                ),
            )
        candidate = CoSpec(**values)
        if half_order(candidate) <= cap:
            return candidate
    return spec


def exact_g_states(polynomial_ring, x, t: int):
    excluded = polynomial_ring.one
    selected = x
    for branch in (1, t, 3):
        old_excluded = excluded
        excluded = (excluded + selected) ** branch
        selected = x * old_excluded**branch
    return excluded, selected


def exact_half_states(spec: CoSpec, polynomial_ring, x):
    g_e, g_s = exact_g_states(polynomial_ring, x, spec.t)
    tg_e = (1 + x) * (g_e + g_s) ** spec.m
    tg_s = x * g_e**spec.m
    star_e = (1 + x) ** spec.nested_star_leaves
    half_e = (tg_e + tg_s) ** spec.tg_count * (
        star_e + x
    ) ** spec.nested_star_count
    half_s = x * tg_e**spec.tg_count * (1 + x) ** (
        spec.nested_star_leaves * spec.nested_star_count
    )
    return half_e, half_s


def exact_replay(candidate: CoCandidate) -> dict[str, Any]:
    polynomial_ring, x = ring("x", ZZ)
    left_e, left_s = exact_half_states(candidate.left, polynomial_ring, x)
    right_e, right_s = exact_half_states(candidate.right, polynomial_ring, x)
    full = left_e * right_e + left_s * right_e + left_e * right_s
    full_coefficients = [int(value) for value in reversed(full.to_dense())]
    left_coefficients = [int(value) for value in reversed(left_s.to_dense())]
    right_coefficients = [int(value) for value in reversed(right_s.to_dense())]
    return {
        "full_degree": len(full_coefficients) - 1,
        "full_valley": exact_valley(full_coefficients),
        "left_selected_valley": exact_valley(left_coefficients),
        "right_selected_valley": exact_valley(right_coefficients),
    }


def build_tree(candidate: CoCandidate) -> tuple[list[tuple[int, int]], int]:
    edges: list[tuple[int, int]] = []
    next_vertex = 0

    def new_vertex() -> int:
        nonlocal next_vertex
        vertex = next_vertex
        next_vertex += 1
        return vertex

    def attach_component(hub: int, adjacency: list[list[int]]) -> None:
        mapping = [new_vertex() for _ in adjacency]
        for u, neighbors in enumerate(adjacency):
            for v in neighbors:
                if u < v:
                    edges.append((mapping[u], mapping[v]))
        edges.append((hub, mapping[0]))

    def build_half(spec: CoSpec) -> int:
        hub = new_vertex()
        tg_adjacency = make_bautista_ramos_tree(spec.m, spec.t)
        for _ in range(spec.tg_count):
            attach_component(hub, tg_adjacency)
        for _ in range(spec.nested_star_count):
            center = new_vertex()
            edges.append((hub, center))
            for _ in range(spec.nested_star_leaves):
                edges.append((center, new_vertex()))
        return hub

    left_hub = build_half(candidate.left)
    right_hub = build_half(candidate.right)
    edges.append((left_hub, right_hub))
    return edges, next_vertex


def treehood_certificate(candidate: CoCandidate) -> dict[str, Any]:
    edges, order = build_tree(candidate)
    adjacency = [[] for _ in range(order)]
    normalized = set()
    for u, v in edges:
        normalized.add((min(u, v), max(u, v)))
        adjacency[u].append(v)
        adjacency[v].append(u)
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adjacency[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    edge_text = "\n".join(f"{u} {v}" for u, v in sorted(normalized))
    certificate = {
        "vertices": order,
        "edges": len(edges),
        "connected_vertices": len(seen),
        "simple": len(normalized) == len(edges),
        "acyclic_by_connected_edge_count": len(seen) == order and len(edges) == order - 1,
        "edge_list_sha256": hashlib.sha256(edge_text.encode()).hexdigest(),
    }
    if not (
        certificate["simple"]
        and certificate["connected_vertices"] == order
        and certificate["acyclic_by_connected_edge_count"]
    ):
        raise AssertionError(f"treehood failure: {certificate}")
    return certificate


def load_checkpoint(path: Path, run_id: str) -> dict[str, Any] | None:
    if not path.exists():
        return None
    last = None
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if record.get("run_id") == run_id and record.get("kind") == "coscale_checkpoint":
            last = record
    return last


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--iterations", type=int, default=5_000)
    parser.add_argument("--seed", type=int, default=20_260_720)
    parser.add_argument("--max-order", type=int, default=100_000)
    parser.add_argument("--population", type=int, default=60)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    parser.add_argument("--max-m", type=int, default=12)
    parser.add_argument("--min-t", type=int, default=3)
    parser.add_argument("--max-t", type=int, default=500)
    parser.add_argument("--max-tg-count", type=int, default=4)
    parser.add_argument("--max-nested-stars", type=int, default=6)
    parser.add_argument("--max-plateau-leaves", type=int, default=20_000)
    parser.add_argument("--relative-floor", type=float, default=1e-13)
    parser.add_argument("--descent-tolerance", type=float, default=1e-8)
    parser.add_argument("--ascent-tolerance", type=float, default=1e-7)
    parser.add_argument("--min-separation", type=int, default=5)
    args = parser.parse_args()

    builder = CoBuilder()
    evaluator = CoEvaluator(args, builder)
    checkpoint = load_checkpoint(args.output, args.run_id)
    if checkpoint is None:
        rng = random.Random(args.seed)
        population_pairs = []
        for _ in range(args.population * 2):
            left = random_spec(args, rng, args.max_order // 2)
            right = left if rng.random() < 0.35 else random_spec(
                args, rng, args.max_order - half_order(left)
            )
            population_pairs.append(canonical_pair(left, right))
        start_iteration = 0
        best_candidate = None
        append_record(
            args.output,
            {
                "kind": "coscale_start",
                "run_id": args.run_id,
                "at": utc_now(),
                "parameters": vars(args) | {"output": str(args.output)},
                "equations": {
                    "E_TG": "(1+x) P_G^m",
                    "S_TG": "x E_G^m",
                    "E_half": "P_TG^a ((1+x)^s+x)^b",
                    "S_half": "x E_TG^a (1+x)^(sb)",
                },
            },
        )
    else:
        rng = random.Random()
        rng.setstate(ast.literal_eval(checkpoint["rng_state_repr"]))
        population_pairs = [
            (CoSpec(**pair[0]), CoSpec(**pair[1]))
            for pair in checkpoint["population_pairs"]
        ]
        start_iteration = int(checkpoint["next_iteration"])
        raw_best = checkpoint.get("best")
        best_candidate = None
        if raw_best is not None:
            best_candidate = CoCandidate(
                **{
                    **raw_best,
                    "left": CoSpec(**raw_best["left"]),
                    "right": CoSpec(**raw_best["right"]),
                    "left_selected": Feature(**raw_best["left_selected"]),
                    "right_selected": Feature(**raw_best["right_selected"]),
                }
            )
        append_record(
            args.output,
            {
                "kind": "coscale_resume",
                "run_id": args.run_id,
                "at": utc_now(),
                "start_iteration": start_iteration,
            },
        )

    started = time.monotonic()
    for iteration in range(start_iteration, args.iterations):
        candidates = [
            candidate
            for left, right in set(population_pairs)
            if (candidate := evaluator.evaluate(left, right)) is not None
        ]
        frontier = choose_frontier(candidates, args.population)
        population_pairs = [(item.left, item.right) for item in frontier]
        champion = max(frontier, key=lambda item: item.score)
        selected_champion = max(frontier, key=selected_objective)
        if best_candidate is None or champion.score > best_candidate.score:
            best_candidate = champion
            append_record(
                args.output,
                {
                    "kind": "coscale_best",
                    "run_id": args.run_id,
                    "at": utc_now(),
                    "iteration": iteration,
                    "candidate": asdict(champion),
                },
            )
            print(json.dumps({"iteration": iteration, "best": asdict(champion)}, sort_keys=True), flush=True)

        full_crossing = champion.score > 1.0 + args.ascent_tolerance
        selected_crossing = selected_objective(selected_champion)[0] > 1.0 + args.ascent_tolerance
        if full_crossing or selected_crossing:
            replay_candidate = champion if full_crossing else selected_champion
            exact_started = time.monotonic()
            exact = exact_replay(replay_candidate)
            record: dict[str, Any] = {
                "kind": "coscale_exact_replay",
                "run_id": args.run_id,
                "at": utc_now(),
                "iteration": iteration,
                "candidate": asdict(replay_candidate),
                "exact_seconds": time.monotonic() - exact_started,
                "exact": exact,
            }
            if any(exact[key] is not None for key in (
                "full_valley", "left_selected_valley", "right_selected_valley"
            )):
                record["treehood"] = treehood_certificate(replay_candidate)
            append_record(args.output, record)
            if exact["full_valley"] is not None:
                print(json.dumps(record, sort_keys=True), flush=True)
                return

        parent = rng.choice(frontier)
        if rng.random() < 0.5:
            left = mutate_spec(
                args,
                parent.left,
                rng,
                args.max_order - half_order(parent.right),
            )
            child = canonical_pair(left, parent.right)
        else:
            right = mutate_spec(
                args,
                parent.right,
                rng,
                args.max_order - half_order(parent.left),
            )
            child = canonical_pair(parent.left, right)
        if half_order(child[0]) + half_order(child[1]) <= args.max_order:
            population_pairs.append(child)

        if iteration % 25 == 0:
            for _ in range(12):
                left = random_spec(args, rng, args.max_order // 2)
                right = left if rng.random() < 0.25 else random_spec(
                    args, rng, args.max_order - half_order(left)
                )
                population_pairs.append(canonical_pair(left, right))

        if (iteration + 1) % args.checkpoint_every == 0 or iteration + 1 == args.iterations:
            record = {
                "kind": "coscale_checkpoint",
                "run_id": args.run_id,
                "at": utc_now(),
                "next_iteration": iteration + 1,
                "elapsed_seconds_this_invocation": time.monotonic() - started,
                "population_pairs": [
                    [asdict(left), asdict(right)] for left, right in population_pairs
                ],
                "rng_state_repr": repr(rng.getstate()),
                "best": asdict(best_candidate) if best_candidate is not None else None,
                "iteration_champion": asdict(champion),
                "selected_state_champion": asdict(selected_champion),
                "evaluated_pairs": len(evaluator.cache),
                "cached_halves": len(builder.half_cache),
                "cached_tg_states": len(builder.tg_cache),
            }
            append_record(args.output, record)
            print(
                json.dumps(
                    {
                        "checkpoint": iteration + 1,
                        "best_ratio": best_candidate.score if best_candidate else None,
                        "selected_objective": selected_objective(selected_champion),
                        "evaluated_pairs": len(evaluator.cache),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )

    append_record(
        args.output,
        {
            "kind": "coscale_complete",
            "run_id": args.run_id,
            "at": utc_now(),
            "iterations": args.iterations,
            "best": asdict(best_candidate) if best_candidate is not None else None,
        },
    )


if __name__ == "__main__":
    main()
