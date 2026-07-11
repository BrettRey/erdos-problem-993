#!/usr/bin/env python3
"""Crash-safe D17 search over paths of rooted spherical blocks.

For disjoint rooted blocks with root-state polynomials (E_i,S_i), join their
roots in a path.  If O_i and N_i count configurations of the first i blocks
whose last root is respectively excluded and selected, then

    O_i = E_i (O_{i-1}+N_{i-1}),   N_i = S_i O_{i-1}.

This is the materially new m>=3 connector: unlike a single joined edge it has
several independently tunable forbidden root-status corners.  Floating point
only ranks candidates; every alleged crossing is replayed in ZZ[x] and its
constructed graph is certified to be a tree before being emitted as a hit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from sympy.polys.domains import ZZ
from sympy.polys.rings import ring

from scratch_d17_checkpointed_search_20260711 import (
    append_record,
    exact_valley,
    random_branching,
    rooted_order,
    strict_post_descent_pressure,
    utc_now,
)
from scratch_d17_spherical_search_20260711 import (
    Scaled,
    add,
    multiply,
    spherical_states,
)


Branching = tuple[int, ...]
DecoratedPath = tuple[Branching, ...]


@dataclass(frozen=True)
class PathCandidate:
    score: float
    blocks: DecoratedPath
    order: int
    first_descent_at: int
    first_descent_ratio: float
    best_later_at: int
    best_later_ratio: float
    rebound: float
    ratio_bump: float
    ratio_bump_at: int | None


class FloatPathEvaluator:
    def __init__(self) -> None:
        self.state_cache: dict[Branching, tuple[Scaled, Scaled]] = {}

    def states(self, branching: Branching) -> tuple[Scaled, Scaled]:
        if branching not in self.state_cache:
            excluded, selected, _ = spherical_states(branching)
            self.state_cache[branching] = (excluded, selected)
        return self.state_cache[branching]

    def polynomial(self, blocks: DecoratedPath) -> np.ndarray:
        previous_out: Scaled | None = None
        previous_total: Scaled | None = None
        for index, branching in enumerate(blocks):
            excluded, selected = self.states(branching)
            if index == 0:
                current_out = excluded
                current_on = selected
            else:
                assert previous_out is not None and previous_total is not None
                current_out = multiply(excluded, previous_total)
                current_on = multiply(selected, previous_out)
            previous_out = current_out
            previous_total = add(current_out, current_on)
        assert previous_total is not None
        values = np.maximum(previous_total.values, 0.0)
        values /= values.sum()
        return values


def ratio_bump_after_descent(
    poly: np.ndarray, first: int, floor: float
) -> tuple[float, int | None]:
    """Largest increase of an adjacent ratio after the first descent."""
    maximum = float(np.max(poly))
    prior: float | None = None
    best = 0.0
    best_at = None
    for k in range(first, len(poly) - 1):
        if poly[k] <= maximum * floor or poly[k + 1] <= maximum * floor:
            continue
        ratio = float(poly[k + 1] / poly[k])
        if prior is not None and ratio - prior > best:
            best = ratio - prior
            best_at = k
        prior = ratio
    return best, best_at


def path_order(blocks: DecoratedPath) -> int:
    return sum(rooted_order(block) for block in blocks)


def canonical_path(blocks: DecoratedPath) -> DecoratedPath:
    reverse = tuple(reversed(blocks))
    return min(blocks, reverse)


def retain_best(candidates: list[PathCandidate], limit: int) -> list[PathCandidate]:
    unique: dict[DecoratedPath, PathCandidate] = {}
    for candidate in candidates:
        key = canonical_path(candidate.blocks)
        prior = unique.get(key)
        if prior is None or (candidate.score, candidate.ratio_bump) > (
            prior.score,
            prior.ratio_bump,
        ):
            unique[key] = candidate
    return sorted(
        unique.values(),
        key=lambda item: (item.score, item.ratio_bump, item.rebound),
        reverse=True,
    )[:limit]


def retain_bump_best(candidates: list[PathCandidate], limit: int) -> list[PathCandidate]:
    unique: dict[DecoratedPath, PathCandidate] = {}
    for candidate in candidates:
        key = canonical_path(candidate.blocks)
        prior = unique.get(key)
        if prior is None or (candidate.ratio_bump, candidate.score) > (
            prior.ratio_bump,
            prior.score,
        ):
            unique[key] = candidate
    return sorted(
        unique.values(),
        key=lambda item: (item.ratio_bump, item.score),
        reverse=True,
    )[:limit]


def make_pool(args: argparse.Namespace) -> list[Branching]:
    rng = random.Random(args.seed ^ 0xD17DEC0)
    cap = args.max_order // args.blocks
    seed_rows: set[Branching] = {
        (1, 1, 2, 1, 2, 3, 2),  # D14 root-outward witness, reversed.
        (4, 4, 2, 1, 2, 2, 1, 2, 4, 4),
        (6, 6, 6, 3, 3, 3, 7),
        (8, 5, 5, 1, 2, 1, 5, 5, 8),
    }
    rows = {row for row in seed_rows if rooted_order(row) <= cap}
    for depth in range(args.min_depth, args.max_depth + 1):
        for branch in range(1, args.max_branch + 1):
            row = (branch,) * depth
            if rooted_order(row) <= cap:
                rows.add(row)
    while len(rows) < args.pool_size:
        rows.add(
            random_branching(
                rng,
                min_depth=args.min_depth,
                max_depth=args.max_depth,
                max_branch=args.max_branch,
                max_root_order=cap,
            )
        )
    return sorted(rows, key=lambda row: (rooted_order(row), row))


def sample_path(args: argparse.Namespace, pool: list[Branching], rng: random.Random) -> DecoratedPath:
    if args.symmetry == "palindromic":
        half = tuple(rng.choice(pool) for _ in range((args.blocks + 1) // 2))
        if args.blocks % 2:
            blocks = half + tuple(reversed(half[:-1]))
        else:
            blocks = half + tuple(reversed(half))
    else:
        blocks = tuple(rng.choice(pool) for _ in range(args.blocks))
    if path_order(blocks) > args.max_order:
        raise AssertionError("pool cap failed to enforce total order")
    return canonical_path(blocks)


def exact_path_polynomial(blocks: DecoratedPath) -> list[int]:
    polynomial_ring, x = ring("x", ZZ)

    def states(branching: Branching):
        excluded = polynomial_ring.one
        selected = x
        for branch in branching:
            old_excluded = excluded
            excluded = (excluded + selected) ** branch
            selected = x * old_excluded**branch
        return excluded, selected

    previous_out = None
    previous_total = None
    for index, branching in enumerate(blocks):
        excluded, selected = states(branching)
        if index == 0:
            current_out = excluded
            current_on = selected
        else:
            current_out = excluded * previous_total
            current_on = selected * previous_out
        previous_out = current_out
        previous_total = current_out + current_on
    return [int(value) for value in reversed(previous_total.to_dense())]


def build_decorated_path(blocks: DecoratedPath) -> tuple[list[tuple[int, int]], int]:
    edges: list[tuple[int, int]] = []
    next_vertex = 0

    def build_spherical(branching: Branching) -> int:
        nonlocal next_vertex
        root = next_vertex
        next_vertex += 1
        if branching:
            for _ in range(branching[-1]):
                child = build_spherical(branching[:-1])
                edges.append((root, child))
        return root

    roots = [build_spherical(block) for block in blocks]
    edges.extend(zip(roots, roots[1:]))
    return edges, next_vertex


def treehood_certificate(blocks: DecoratedPath) -> dict[str, Any]:
    edges, order = build_decorated_path(blocks)
    adjacency = [[] for _ in range(order)]
    normalized_edges = set()
    for u, v in edges:
        normalized_edges.add((min(u, v), max(u, v)))
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
    edge_text = "\n".join(f"{u} {v}" for u, v in sorted(normalized_edges))
    certificate = {
        "vertices": order,
        "edges": len(edges),
        "connected_vertices": len(seen),
        "simple": len(normalized_edges) == len(edges),
        "acyclic_by_connected_edge_count": len(seen) == order and len(edges) == order - 1,
        "edge_list_sha256": hashlib.sha256(edge_text.encode()).hexdigest(),
    }
    if not (
        certificate["simple"]
        and certificate["connected_vertices"] == order
        and certificate["acyclic_by_connected_edge_count"]
    ):
        raise AssertionError(f"treehood replay failed: {certificate}")
    return certificate


def completed_batches(path: Path, run_id: str) -> set[int]:
    if not path.exists():
        return set()
    completed = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if record.get("run_id") == run_id and record.get("kind") == "path_batch":
            completed.add(int(record["batch"]))
    return completed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--blocks", type=int, choices=(3, 4, 5, 6), required=True)
    parser.add_argument("--symmetry", choices=("palindromic", "asymmetric"), required=True)
    parser.add_argument("--batches", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=1_000)
    parser.add_argument("--pool-size", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20_260_711)
    parser.add_argument("--min-depth", type=int, default=3)
    parser.add_argument("--max-depth", type=int, default=13)
    parser.add_argument("--max-branch", type=int, default=8)
    parser.add_argument("--max-order", type=int, default=100_000)
    parser.add_argument("--top", type=int, default=15)
    parser.add_argument("--relative-floor", type=float, default=1e-13)
    parser.add_argument("--descent-tolerance", type=float, default=1e-8)
    parser.add_argument("--ascent-tolerance", type=float, default=1e-7)
    parser.add_argument("--min-separation", type=int, default=5)
    parser.add_argument("--exact-replay-best-bump", action="store_true")
    args = parser.parse_args()

    pool = make_pool(args)
    evaluator = FloatPathEvaluator()
    done = completed_batches(args.output, args.run_id)
    append_record(
        args.output,
        {
            "kind": "path_resume" if done else "path_start",
            "run_id": args.run_id,
            "at": utc_now(),
            "completed_batches": sorted(done),
            "parameters": vars(args) | {"output": str(args.output)},
            "pool": [
                {"branching_leaves_to_root": row, "order": rooted_order(row)} for row in pool
            ],
            "exact_recurrence": "O_i=E_i(O_(i-1)+N_(i-1)); N_i=S_i O_(i-1)",
        },
    )

    global_best: list[PathCandidate] = []
    global_bump_best: list[PathCandidate] = []
    tested_total = 0
    for batch in range(args.batches):
        if batch in done:
            continue
        started = time.monotonic()
        rng = random.Random(args.seed + 1_000_003 * batch)
        batch_best: list[PathCandidate] = []
        batch_bump_best: list[PathCandidate] = []
        seen: set[DecoratedPath] = set()
        attempts = 0
        while len(seen) < args.batch_size and attempts < 10 * args.batch_size:
            attempts += 1
            blocks = sample_path(args, pool, rng)
            if blocks in seen:
                continue
            seen.add(blocks)
            poly = evaluator.polynomial(blocks)
            pressure = strict_post_descent_pressure(
                poly,
                relative_floor=args.relative_floor,
                descent_tolerance=args.descent_tolerance,
                min_separation=args.min_separation,
            )
            if pressure is None:
                continue
            ratio_bump, ratio_bump_at = ratio_bump_after_descent(
                poly,
                int(pressure["first_descent_at"]),
                args.relative_floor,
            )
            candidate = PathCandidate(
                score=float(pressure["best_later_ratio"]),
                blocks=blocks,
                order=path_order(blocks),
                first_descent_at=int(pressure["first_descent_at"]),
                first_descent_ratio=float(pressure["first_descent_ratio"]),
                best_later_at=int(pressure["best_later_at"]),
                best_later_ratio=float(pressure["best_later_ratio"]),
                rebound=float(pressure["rebound"]),
                ratio_bump=ratio_bump,
                ratio_bump_at=ratio_bump_at,
            )
            batch_best = retain_best([*batch_best, candidate], args.top)
            batch_bump_best = retain_bump_best([*batch_bump_best, candidate], args.top)
            if candidate.score > 1.0 + args.ascent_tolerance:
                exact_started = time.monotonic()
                coefficients = exact_path_polynomial(blocks)
                valley = exact_valley(coefficients)
                exact_record: dict[str, Any] = {
                    "kind": "path_exact_replay",
                    "run_id": args.run_id,
                    "batch": batch,
                    "at": utc_now(),
                    "candidate": asdict(candidate),
                    "exact_seconds": time.monotonic() - exact_started,
                    "exact_degree": len(coefficients) - 1,
                    "exact_valley": valley,
                }
                if valley is not None:
                    exact_record["treehood"] = treehood_certificate(blocks)
                append_record(args.output, exact_record)
                if valley is not None:
                    print(json.dumps(exact_record, sort_keys=True), flush=True)
                    return
        tested_total += len(seen)
        global_best = retain_best([*global_best, *batch_best], args.top)
        global_bump_best = retain_bump_best(
            [*global_bump_best, *batch_bump_best], args.top
        )
        checkpoint = {
            "kind": "path_batch",
            "run_id": args.run_id,
            "batch": batch,
            "at": utc_now(),
            "tested": len(seen),
            "attempts": attempts,
            "tested_total_this_invocation": tested_total,
            "elapsed_seconds": time.monotonic() - started,
            "cached_states": len(evaluator.state_cache),
            "batch_best": [asdict(item) for item in batch_best],
            "batch_bump_best": [asdict(item) for item in batch_bump_best],
            "global_best_this_invocation": [asdict(item) for item in global_best],
            "global_bump_best_this_invocation": [
                asdict(item) for item in global_bump_best
            ],
        }
        append_record(args.output, checkpoint)
        print(json.dumps(checkpoint, sort_keys=True), flush=True)

    if args.exact_replay_best_bump and global_bump_best:
        candidate = global_bump_best[0]
        exact_started = time.monotonic()
        coefficients = exact_path_polynomial(candidate.blocks)
        bump_at = candidate.ratio_bump_at
        exact_bump = None
        if bump_at is not None and 1 <= bump_at < len(coefficients) - 1:
            exact_bump = {
                "at": bump_at,
                "left_ratio": [coefficients[bump_at], coefficients[bump_at - 1]],
                "right_ratio": [coefficients[bump_at + 1], coefficients[bump_at]],
                "lc_bump_numerator": (
                    coefficients[bump_at - 1] * coefficients[bump_at + 1]
                    - coefficients[bump_at] ** 2
                ),
                "right_endpoint_below_one": (
                    coefficients[bump_at + 1] < coefficients[bump_at]
                ),
            }
        append_record(
            args.output,
            {
                "kind": "path_exact_bump_replay",
                "run_id": args.run_id,
                "at": utc_now(),
                "candidate": asdict(candidate),
                "exact_seconds": time.monotonic() - exact_started,
                "exact_degree": len(coefficients) - 1,
                "exact_valley": exact_valley(coefficients),
                "exact_bump": exact_bump,
                "treehood": treehood_certificate(candidate.blocks),
            },
        )

    append_record(
        args.output,
        {
            "kind": "path_complete",
            "run_id": args.run_id,
            "at": utc_now(),
            "tested_total_this_invocation": tested_total,
            "best_this_invocation": [asdict(item) for item in global_best],
            "bump_best_this_invocation": [
                asdict(item) for item in global_bump_best
            ],
        },
    )


if __name__ == "__main__":
    main()
