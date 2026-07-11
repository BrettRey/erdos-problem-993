#!/usr/bin/env python3
"""Checkpointed D20 search over heterogeneous two-scale rooted halves.

Each rooted half has a hub with four child classes:

* copies of one Bautista--Ramos TG species (hard selected-state bumps),
* optional copies of one Galvin species (hard total-state control),
* nested stars with s leaves (near-plateau without root-weight collapse), and
* direct leaves (near-plateau with an exact 2^{-d} selected-weight penalty).

For child root states (E_i,S_i), P_i=E_i+S_i, the half states are

    E = P_TG^a P_G^g ((1+x)^s+x)^b (1+x)^d,
    S = x E_TG^a E_G^g (1+x)^(sb).

Two halves are joined by one edge.  Floats rank only; any full-polynomial or
selected-state crossing is replayed exactly in ZZ[x].  A full crossing is
reported only with an explicit treehood certificate.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import os
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
)
from scripts.analyze_prufer_corpus import (
    make_bautista_ramos_tree,
    make_galvin_tree,
)


def poly_add(a: list[int], b: list[int]) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return out


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        if not left:
            continue
        for j, right in enumerate(b):
            if right:
                out[i + j] += left * right
    return out


def rooted_states(adjacency: list[list[int]], root: int = 0) -> tuple[list[int], list[int]]:
    parent = [-2] * len(adjacency)
    parent[root] = -1
    order = [root]
    for vertex in order:
        for neighbor in adjacency[vertex]:
            if parent[neighbor] == -2:
                parent[neighbor] = vertex
                order.append(neighbor)
    excluded: dict[int, list[int]] = {}
    selected: dict[int, list[int]] = {}
    for vertex in reversed(order):
        e_state = [1]
        s_state = [0, 1]
        for child in adjacency[vertex]:
            if parent[child] == vertex:
                e_state = poly_mul(e_state, poly_add(excluded[child], selected[child]))
                s_state = poly_mul(s_state, excluded[child])
        excluded[vertex] = e_state
        selected[vertex] = s_state
    return excluded[root], selected[root]


def exact_to_scaled(coefficients: list[int]) -> Scaled:
    total = sum(coefficients)
    values = np.array([value / total for value in coefficients], dtype=float)
    return normalize(Scaled(values, math.log(total)))


def shift_x(polynomial: Scaled) -> Scaled:
    return Scaled(np.pad(polynomial.values, (1, 0)), polynomial.log_scale)


@dataclass(frozen=True)
class Species:
    label: str
    adjacency: list[list[int]]
    excluded_exact: list[int]
    selected_exact: list[int]
    excluded_float: Scaled
    selected_float: Scaled
    total_float: Scaled

    @property
    def order(self) -> int:
        return len(self.adjacency)


def make_species(label: str, adjacency: list[list[int]]) -> Species:
    excluded, selected = rooted_states(adjacency)
    excluded_float = exact_to_scaled(excluded)
    selected_float = exact_to_scaled(selected)
    return Species(
        label=label,
        adjacency=adjacency,
        excluded_exact=excluded,
        selected_exact=selected,
        excluded_float=excluded_float,
        selected_float=selected_float,
        total_float=add(excluded_float, selected_float),
    )


def species_libraries() -> tuple[list[Species], list[Species]]:
    tg_parameters = ((2, 5), (4, 6), (6, 7), (8, 8))
    galvin_parameters = ((3, 4), (6, 6), (14, 8), (21, 11))
    tg = [
        make_species(f"TG_{m}_{t}", make_bautista_ramos_tree(m, t))
        for m, t in tg_parameters
    ]
    galvin = [
        make_species(f"G_{m}_{t}", make_galvin_tree(m, t))
        for m, t in galvin_parameters
    ]
    return tg, galvin


@dataclass(frozen=True, order=True)
class HalfSpec:
    tg_index: int
    tg_count: int
    galvin_index: int
    galvin_count: int
    nested_star_leaves: int
    nested_star_count: int
    direct_leaves: int


@dataclass(frozen=True)
class Feature:
    first_descent_at: int | None
    first_descent_ratio: float | None
    best_later_at: int | None
    best_later_ratio: float
    ratio_bump: float
    ratio_bump_at: int | None


@dataclass(frozen=True)
class Candidate:
    score: float
    left: HalfSpec
    right: HalfSpec
    order: int
    first_descent_at: int
    first_descent_ratio: float
    best_later_at: int
    best_later_ratio: float
    left_selected: Feature
    right_selected: Feature
    left_log_selected_to_excluded_mass: float
    right_log_selected_to_excluded_mass: float


class FloatBuilder:
    def __init__(self, tg: list[Species], galvin: list[Species]) -> None:
        self.tg = tg
        self.galvin = galvin
        self.binomial_cache: dict[int, Scaled] = {}
        self.half_cache: dict[HalfSpec, tuple[Scaled, Scaled]] = {}

    def binomial(self, exponent: int) -> Scaled:
        if exponent not in self.binomial_cache:
            if exponent == 0:
                result = Scaled(np.array([1.0]), 0.0)
            else:
                k = np.arange(exponent + 1, dtype=float)
                log_values = (
                    gammaln(exponent + 1.0)
                    - gammaln(k + 1.0)
                    - gammaln(exponent - k + 1.0)
                    - exponent * math.log(2.0)
                )
                probabilities = np.exp(log_values)
                probabilities /= probabilities.sum()
                result = normalize(
                    Scaled(probabilities, exponent * math.log(2.0))
                )
            self.binomial_cache[exponent] = result
        return self.binomial_cache[exponent]

    def half_states(self, spec: HalfSpec) -> tuple[Scaled, Scaled]:
        if spec in self.half_cache:
            return self.half_cache[spec]
        tg = self.tg[spec.tg_index]
        excluded = power(tg.total_float, spec.tg_count)
        selected_base = power(tg.excluded_float, spec.tg_count)
        if spec.galvin_count:
            galvin = self.galvin[spec.galvin_index]
            excluded = multiply(excluded, power(galvin.total_float, spec.galvin_count))
            selected_base = multiply(
                selected_base,
                power(galvin.excluded_float, spec.galvin_count),
            )
        if spec.nested_star_count:
            star_excluded = self.binomial(spec.nested_star_leaves)
            star_selected = Scaled(np.array([0.0, 1.0]), 0.0)
            star_total = add(star_excluded, star_selected)
            excluded = multiply(
                excluded,
                power(star_total, spec.nested_star_count),
            )
            selected_base = multiply(
                selected_base,
                self.binomial(spec.nested_star_leaves * spec.nested_star_count),
            )
        if spec.direct_leaves:
            excluded = multiply(excluded, self.binomial(spec.direct_leaves))
        selected = shift_x(selected_base)
        self.half_cache[spec] = (excluded, selected)
        return excluded, selected

    def joined(self, left: HalfSpec, right: HalfSpec) -> np.ndarray:
        left_excluded, left_selected = self.half_states(left)
        right_excluded, right_selected = self.half_states(right)
        full = add(
            multiply(left_excluded, right_excluded),
            add(
                multiply(left_selected, right_excluded),
                multiply(left_excluded, right_selected),
            ),
        )
        values = np.maximum(full.values, 0.0)
        values /= values.sum()
        return values


def half_order(spec: HalfSpec, tg: list[Species], galvin: list[Species]) -> int:
    return (
        1
        + spec.tg_count * tg[spec.tg_index].order
        + spec.galvin_count * galvin[spec.galvin_index].order
        + spec.nested_star_count * (spec.nested_star_leaves + 1)
        + spec.direct_leaves
    )


def scaled_log_mass(polynomial: Scaled) -> float:
    return polynomial.log_scale + math.log(float(np.sum(polynomial.values)))


def feature(
    polynomial: Scaled,
    *,
    relative_floor: float,
    descent_tolerance: float,
    min_separation: int,
) -> Feature:
    values = np.maximum(polynomial.values, 0.0)
    values /= values.sum()
    pressure = strict_post_descent_pressure(
        values,
        relative_floor=relative_floor,
        descent_tolerance=descent_tolerance,
        min_separation=min_separation,
    )
    if pressure is None:
        return Feature(None, None, None, 0.0, 0.0, None)
    bump, bump_at = ratio_bump_after_descent(
        values,
        int(pressure["first_descent_at"]),
        relative_floor,
    )
    return Feature(
        first_descent_at=int(pressure["first_descent_at"]),
        first_descent_ratio=float(pressure["first_descent_ratio"]),
        best_later_at=int(pressure["best_later_at"]),
        best_later_ratio=float(pressure["best_later_ratio"]),
        ratio_bump=bump,
        ratio_bump_at=bump_at,
    )


def canonical_pair(left: HalfSpec, right: HalfSpec) -> tuple[HalfSpec, HalfSpec]:
    return (left, right) if left <= right else (right, left)


class Evaluator:
    def __init__(self, args: argparse.Namespace, builder: FloatBuilder, tg, galvin) -> None:
        self.args = args
        self.builder = builder
        self.tg = tg
        self.galvin = galvin
        self.cache: dict[tuple[HalfSpec, HalfSpec], Candidate | None] = {}

    def evaluate(self, left: HalfSpec, right: HalfSpec) -> Candidate | None:
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
        candidate = Candidate(
            score=float(pressure["best_later_ratio"]),
            left=left,
            right=right,
            order=half_order(left, self.tg, self.galvin)
            + half_order(right, self.tg, self.galvin),
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


def selected_objective(candidate: Candidate) -> tuple[float, float]:
    features = (candidate.left_selected, candidate.right_selected)
    return (
        max(item.best_later_ratio for item in features),
        max(item.ratio_bump for item in features),
    )


def choose_frontier(candidates: list[Candidate], population: int) -> list[Candidate]:
    selected: dict[tuple[HalfSpec, HalfSpec], Candidate] = {}
    by_full = sorted(candidates, key=lambda item: item.score, reverse=True)
    by_selected = sorted(candidates, key=selected_objective, reverse=True)
    by_bump = sorted(
        candidates,
        key=lambda item: max(
            item.left_selected.ratio_bump,
            item.right_selected.ratio_bump,
        ),
        reverse=True,
    )
    main_count = max(1, population // 2)
    secondary_count = max(1, population // 4)
    for item in by_full[:main_count] + by_selected[:secondary_count] + by_bump[:secondary_count]:
        selected[(item.left, item.right)] = item
    return list(selected.values())[:population]


def random_plateau_size(rng: random.Random, maximum: int) -> int:
    anchors = (0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1000, 2000, 5000)
    if rng.random() < 0.65:
        return min(maximum, rng.choice(anchors))
    return rng.randint(0, maximum)


def random_spec(args, rng, tg, galvin, cap: int) -> HalfSpec:
    for _ in range(10_000):
        nested_count = rng.randint(0, args.max_nested_stars)
        nested_leaves = (
            random_plateau_size(rng, args.max_plateau_leaves)
            if nested_count
            else 0
        )
        direct = random_plateau_size(rng, args.max_direct_leaves)
        if not nested_count and not direct:
            if args.max_direct_leaves:
                direct = 1
            else:
                nested_count = 1
                nested_leaves = max(
                    1, random_plateau_size(rng, args.max_plateau_leaves)
                )
        spec = HalfSpec(
            tg_index=rng.randrange(len(tg)),
            tg_count=rng.randint(1, args.max_hard_count),
            galvin_index=rng.randrange(len(galvin)),
            galvin_count=rng.randint(0, args.max_galvin_count),
            nested_star_leaves=nested_leaves,
            nested_star_count=nested_count,
            direct_leaves=direct,
        )
        if half_order(spec, tg, galvin) <= cap:
            return spec
    raise RuntimeError("unable to sample half under order cap")


def mutate_spec(args, spec: HalfSpec, rng, tg, galvin, cap: int) -> HalfSpec:
    for _ in range(1_000):
        values = asdict(spec)
        field = rng.choice(tuple(values))
        if field == "tg_index":
            values[field] = rng.randrange(len(tg))
        elif field == "tg_count":
            values[field] = rng.randint(1, args.max_hard_count)
        elif field == "galvin_index":
            values[field] = rng.randrange(len(galvin))
        elif field == "galvin_count":
            values[field] = rng.randint(0, args.max_galvin_count)
        elif field == "nested_star_count":
            values[field] = rng.randint(0, args.max_nested_stars)
            if values[field] == 0:
                values["nested_star_leaves"] = 0
        elif field == "nested_star_leaves":
            if values["nested_star_count"] == 0:
                values["nested_star_count"] = 1
            old = int(values[field])
            values[field] = max(
                0,
                min(
                    args.max_plateau_leaves,
                    old + rng.choice((-256, -64, -16, -4, -1, 1, 4, 16, 64, 256)),
                ),
            )
        else:
            old = int(values[field])
            values[field] = max(
                0,
                min(
                    args.max_direct_leaves,
                    old + rng.choice((-256, -64, -16, -4, -1, 1, 4, 16, 64, 256)),
                ),
            )
        if not values["nested_star_count"] and not values["direct_leaves"]:
            if args.max_direct_leaves:
                values["direct_leaves"] = 1
            else:
                values["nested_star_count"] = 1
                values["nested_star_leaves"] = max(
                    1, int(values["nested_star_leaves"])
                )
        candidate = HalfSpec(**values)
        if half_order(candidate, tg, galvin) <= cap:
            return candidate
    return spec


def exact_element(polynomial_ring, coefficients: list[int]):
    return polynomial_ring.from_dict(
        {(index,): value for index, value in enumerate(coefficients) if value}
    )


def exact_half_states(spec: HalfSpec, tg, galvin, polynomial_ring, x):
    tg_species = tg[spec.tg_index]
    tg_e = exact_element(polynomial_ring, tg_species.excluded_exact)
    tg_s = exact_element(polynomial_ring, tg_species.selected_exact)
    excluded = (tg_e + tg_s) ** spec.tg_count
    selected_base = tg_e**spec.tg_count
    if spec.galvin_count:
        g_species = galvin[spec.galvin_index]
        g_e = exact_element(polynomial_ring, g_species.excluded_exact)
        g_s = exact_element(polynomial_ring, g_species.selected_exact)
        excluded *= (g_e + g_s) ** spec.galvin_count
        selected_base *= g_e**spec.galvin_count
    if spec.nested_star_count:
        binomial = (1 + x) ** spec.nested_star_leaves
        excluded *= (binomial + x) ** spec.nested_star_count
        selected_base *= (1 + x) ** (
            spec.nested_star_leaves * spec.nested_star_count
        )
    if spec.direct_leaves:
        excluded *= (1 + x) ** spec.direct_leaves
    return excluded, x * selected_base


def exact_replay(candidate: Candidate, tg, galvin) -> dict[str, Any]:
    polynomial_ring, x = ring("x", ZZ)
    left_e, left_s = exact_half_states(candidate.left, tg, galvin, polynomial_ring, x)
    right_e, right_s = exact_half_states(candidate.right, tg, galvin, polynomial_ring, x)
    full = left_e * right_e + left_s * right_e + left_e * right_s
    coefficients = [int(value) for value in reversed(full.to_dense())]
    left_selected = [int(value) for value in reversed(left_s.to_dense())]
    right_selected = [int(value) for value in reversed(right_s.to_dense())]
    return {
        "full_degree": len(coefficients) - 1,
        "full_valley": exact_valley(coefficients),
        "left_selected_valley": exact_valley(left_selected),
        "right_selected_valley": exact_valley(right_selected),
    }


def build_tree(candidate: Candidate, tg, galvin) -> tuple[list[tuple[int, int]], int]:
    edges: list[tuple[int, int]] = []
    next_vertex = 0

    def new_vertex() -> int:
        nonlocal next_vertex
        vertex = next_vertex
        next_vertex += 1
        return vertex

    def attach_component(hub: int, species: Species) -> None:
        mapping = [new_vertex() for _ in species.adjacency]
        for u, neighbors in enumerate(species.adjacency):
            for v in neighbors:
                if u < v:
                    edges.append((mapping[u], mapping[v]))
        edges.append((hub, mapping[0]))

    def build_half(spec: HalfSpec) -> int:
        hub = new_vertex()
        for _ in range(spec.tg_count):
            attach_component(hub, tg[spec.tg_index])
        for _ in range(spec.galvin_count):
            attach_component(hub, galvin[spec.galvin_index])
        for _ in range(spec.nested_star_count):
            center = new_vertex()
            edges.append((hub, center))
            for _ in range(spec.nested_star_leaves):
                leaf = new_vertex()
                edges.append((center, leaf))
        for _ in range(spec.direct_leaves):
            leaf = new_vertex()
            edges.append((hub, leaf))
        return hub

    left_hub = build_half(candidate.left)
    right_hub = build_half(candidate.right)
    edges.append((left_hub, right_hub))
    return edges, next_vertex


def treehood_certificate(candidate: Candidate, tg, galvin) -> dict[str, Any]:
    edges, order = build_tree(candidate, tg, galvin)
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
        if record.get("run_id") == run_id and record.get("kind") == "d20_checkpoint":
            last = record
    return last


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--iterations", type=int, default=5_000)
    parser.add_argument("--seed", type=int, default=20_260_711)
    parser.add_argument("--max-order", type=int, default=100_000)
    parser.add_argument("--population", type=int, default=60)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    parser.add_argument("--max-hard-count", type=int, default=5)
    parser.add_argument("--max-galvin-count", type=int, default=3)
    parser.add_argument("--max-nested-stars", type=int, default=3)
    parser.add_argument("--max-plateau-leaves", type=int, default=5_000)
    parser.add_argument("--max-direct-leaves", type=int, default=5_000)
    parser.add_argument("--relative-floor", type=float, default=1e-13)
    parser.add_argument("--descent-tolerance", type=float, default=1e-8)
    parser.add_argument("--ascent-tolerance", type=float, default=1e-7)
    parser.add_argument("--min-separation", type=int, default=5)
    args = parser.parse_args()

    tg, galvin = species_libraries()
    builder = FloatBuilder(tg, galvin)
    evaluator = Evaluator(args, builder, tg, galvin)
    checkpoint = load_checkpoint(args.output, args.run_id)
    cap = args.max_order - 1
    if checkpoint is None:
        rng = random.Random(args.seed)
        population_pairs = []
        for _ in range(args.population * 2):
            left = random_spec(args, rng, tg, galvin, cap // 2)
            right = left if rng.random() < 0.35 else random_spec(
                args, rng, tg, galvin, args.max_order - half_order(left, tg, galvin)
            )
            if half_order(left, tg, galvin) + half_order(right, tg, galvin) <= args.max_order:
                population_pairs.append(canonical_pair(left, right))
        start_iteration = 0
        best = None
        append_record(
            args.output,
            {
                "kind": "d20_start",
                "run_id": args.run_id,
                "at": utc_now(),
                "parameters": vars(args) | {"output": str(args.output)},
                "tg_species": [item.label for item in tg],
                "galvin_species": [item.label for item in galvin],
                "half_equations": {
                    "E": "P_TG^a P_G^g ((1+x)^s+x)^b (1+x)^d",
                    "S": "x E_TG^a E_G^g (1+x)^(sb)",
                },
            },
        )
    else:
        rng = random.Random()
        rng.setstate(ast.literal_eval(checkpoint["rng_state_repr"]))
        population_pairs = [
            (HalfSpec(**pair[0]), HalfSpec(**pair[1]))
            for pair in checkpoint["population_pairs"]
        ]
        start_iteration = int(checkpoint["next_iteration"])
        best = checkpoint.get("best")
        append_record(
            args.output,
            {
                "kind": "d20_resume",
                "run_id": args.run_id,
                "at": utc_now(),
                "start_iteration": start_iteration,
            },
        )

    started = time.monotonic()
    best_candidate: Candidate | None = None
    if best is not None:
        best_candidate = Candidate(
            **{
                **best,
                "left": HalfSpec(**best["left"]),
                "right": HalfSpec(**best["right"]),
                "left_selected": Feature(**best["left_selected"]),
                "right_selected": Feature(**best["right_selected"]),
            }
        )

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
                    "kind": "d20_best",
                    "run_id": args.run_id,
                    "at": utc_now(),
                    "iteration": iteration,
                    "candidate": asdict(champion),
                },
            )
            print(json.dumps({"iteration": iteration, "best": asdict(champion)}, sort_keys=True), flush=True)

        full_crossing = champion.score > 1.0 + args.ascent_tolerance
        selected_crossing = (
            max(selected_objective(selected_champion)[0], 0.0)
            > 1.0 + args.ascent_tolerance
        )
        if full_crossing or selected_crossing:
            replay_candidate = champion if full_crossing else selected_champion
            exact_started = time.monotonic()
            exact = exact_replay(replay_candidate, tg, galvin)
            record: dict[str, Any] = {
                "kind": "d20_exact_replay",
                "run_id": args.run_id,
                "at": utc_now(),
                "iteration": iteration,
                "candidate": asdict(replay_candidate),
                "exact_seconds": time.monotonic() - exact_started,
                "exact": exact,
            }
            if (
                exact["full_valley"] is not None
                or exact["left_selected_valley"] is not None
                or exact["right_selected_valley"] is not None
            ):
                record["treehood"] = treehood_certificate(replay_candidate, tg, galvin)
            append_record(args.output, record)
            if exact["full_valley"] is not None:
                print(json.dumps(record, sort_keys=True), flush=True)
                return

        parent = rng.choice(frontier)
        mutate_left = rng.random() < 0.5
        if mutate_left:
            cap_left = args.max_order - half_order(parent.right, tg, galvin)
            left = mutate_spec(args, parent.left, rng, tg, galvin, cap_left)
            child = canonical_pair(left, parent.right)
        else:
            cap_right = args.max_order - half_order(parent.left, tg, galvin)
            right = mutate_spec(args, parent.right, rng, tg, galvin, cap_right)
            child = canonical_pair(parent.left, right)
        if half_order(child[0], tg, galvin) + half_order(child[1], tg, galvin) <= args.max_order:
            population_pairs.append(child)

        if iteration % 25 == 0:
            for _ in range(12):
                left = random_spec(args, rng, tg, galvin, cap // 2)
                right = left if rng.random() < 0.25 else random_spec(
                    args,
                    rng,
                    tg,
                    galvin,
                    args.max_order - half_order(left, tg, galvin),
                )
                population_pairs.append(canonical_pair(left, right))

        if (iteration + 1) % args.checkpoint_every == 0 or iteration + 1 == args.iterations:
            checkpoint_record = {
                "kind": "d20_checkpoint",
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
            }
            append_record(args.output, checkpoint_record)
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
            "kind": "d20_complete",
            "run_id": args.run_id,
            "at": utc_now(),
            "iterations": args.iterations,
            "best": asdict(best_candidate) if best_candidate is not None else None,
        },
    )


if __name__ == "__main__":
    main()
