#!/usr/bin/env python3
"""Crash-safe D17 search over edge-joined rooted spherical trees.

The floating-point calculation is only a candidate ranker.  A hit is emitted
as a crossing only after exact integer replay and an explicit treehood check.
Each completed batch is appended, flushed, and fsynced to a JSONL ledger, so a
crash can lose at most the currently running batch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import time
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from sympy.polys.domains import ZZ
from sympy.polys.rings import ring

from scratch_d17_spherical_search_20260711 import joined_polynomial


Branching = tuple[int, ...]


@dataclass(frozen=True)
class RankedCandidate:
    score: float
    left: Branching
    right: Branching
    order: int
    first_descent_at: int
    first_descent_ratio: float
    best_later_at: int
    best_later_ratio: float
    rebound: float


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


def rooted_order(branching: Branching) -> int:
    order = 1
    for branch in branching:
        order = 1 + branch * order
    return order


def append_record(path: Path, record: dict[str, Any]) -> None:
    """Append one durable checkpoint; a torn final line is safe to ignore."""
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(encoded)
        handle.flush()
        os.fsync(handle.fileno())


def strict_post_descent_pressure(
    poly: np.ndarray,
    *,
    relative_floor: float,
    descent_tolerance: float,
    min_separation: int,
) -> dict[str, float | int] | None:
    """Rank only a later edge after a prior strict coefficient descent.

    The score is the largest later adjacent ratio, never a peak tie or a
    log-concavity ratio.  Thus score > 1 is precisely the requested float
    pattern: a strict descent followed by a later strict ascent.
    """
    maximum = float(np.max(poly))
    support = np.flatnonzero(poly > maximum * relative_floor)
    if len(support) < 3:
        return None
    start = int(support[0])
    stop = int(support[-1])
    first: int | None = None
    first_ratio = math.inf
    best_at: int | None = None
    best_ratio = -math.inf
    for k in range(start, stop):
        if poly[k] <= 0.0:
            continue
        ratio = float(poly[k + 1] / poly[k])
        if first is None:
            if ratio < 1.0 - descent_tolerance:
                first = k
                first_ratio = ratio
            continue
        if k >= first + min_separation and ratio > best_ratio:
            best_at = k
            best_ratio = ratio
    if first is None or best_at is None:
        return None
    return {
        "first_descent_at": first,
        "first_descent_ratio": first_ratio,
        "best_later_at": best_at,
        "best_later_ratio": best_ratio,
        "rebound": best_ratio - first_ratio,
    }


def random_branching(
    rng: random.Random,
    *,
    min_depth: int,
    max_depth: int,
    max_branch: int,
    max_root_order: int,
) -> Branching:
    """Sample both smooth and deliberately inhomogeneous spherical profiles."""
    for _ in range(10_000):
        depth = rng.randint(min_depth, max_depth)
        style = rng.randrange(6)
        if style == 0:
            row = tuple(rng.randint(1, max_branch) for _ in range(depth))
        elif style == 1:
            cut = rng.randint(1, depth - 1)
            a, b = rng.randint(1, max_branch), rng.randint(1, max_branch)
            row = (a,) * cut + (b,) * (depth - cut)
        elif style == 2:
            a, b = rng.randint(1, max_branch), rng.randint(1, max_branch)
            row = tuple(a if j % 2 == 0 else b for j in range(depth))
        elif style == 3:
            # Long unary corridors separated by one or two high-branch layers.
            values = [1] * depth
            for j in rng.sample(range(depth), k=rng.randint(1, min(2, depth))):
                values[j] = rng.randint(2, max_branch)
            row = tuple(values)
        elif style == 4:
            half = [rng.randint(1, max_branch) for _ in range((depth + 1) // 2)]
            row = tuple((half + half[-2::-1])[:depth])
        elif depth >= 3:
            # Three constant shells can align distinct conditional shoulders.
            c1 = rng.randint(1, depth - 2)
            c2 = rng.randint(c1 + 1, depth - 1)
            a, b, c = (rng.randint(1, max_branch) for _ in range(3))
            row = (a,) * c1 + (b,) * (c2 - c1) + (c,) * (depth - c2)
        else:
            row = tuple(rng.randint(1, max_branch) for _ in range(depth))
        if rooted_order(row) <= max_root_order:
            return row
    raise RuntimeError("could not sample a branching row under the order cap")


def mutate_branching(
    row: Branching,
    rng: random.Random,
    *,
    max_branch: int,
    max_root_order: int,
) -> Branching:
    """Make a nearby asymmetric phase profile while preserving the size cap."""
    for _ in range(1_000):
        values = list(row)
        operation = rng.randrange(4)
        if operation == 0 and len(values) > 2:
            del values[rng.randrange(len(values))]
        elif operation == 1:
            values.insert(rng.randrange(len(values) + 1), rng.randint(1, max_branch))
        else:
            j = rng.randrange(len(values))
            values[j] = rng.randint(1, max_branch)
        candidate = tuple(values)
        if rooted_order(candidate) <= max_root_order:
            return candidate
    return row


def sample_pair(args: argparse.Namespace, rng: random.Random) -> tuple[Branching, Branching]:
    left_cap = args.max_order // 2 if args.mode == "symmetric" else args.max_order - 1
    left = random_branching(
        rng,
        min_depth=args.min_depth,
        max_depth=args.max_depth,
        max_branch=args.max_branch,
        max_root_order=left_cap,
    )
    if args.mode == "symmetric":
        return left, left
    for _ in range(10_000):
        if rng.random() < args.correlated_probability:
            right = mutate_branching(
                left,
                rng,
                max_branch=args.max_branch,
                max_root_order=args.max_order - rooted_order(left),
            )
        else:
            right = random_branching(
                rng,
                min_depth=args.min_depth,
                max_depth=args.max_depth,
                max_branch=args.max_branch,
                max_root_order=args.max_order - rooted_order(left),
            )
        if rooted_order(left) + rooted_order(right) <= args.max_order:
            return left, right
    raise RuntimeError("could not sample an asymmetric pair under the order cap")


def exact_joined_polynomial(left: Branching, right: Branching) -> list[int]:
    """Replay the spherical recurrence in ZZ[x] using SymPy's exact ring."""
    polynomial_ring, x = ring("x", ZZ)

    def states(branching: Branching):
        excluded = polynomial_ring.one
        selected = x
        for branch in branching:
            old_excluded = excluded
            excluded = (excluded + selected) ** branch
            selected = x * old_excluded**branch
        return excluded, selected

    left_excluded, left_selected = states(left)
    right_excluded, right_selected = states(right)
    full = (
        left_excluded * right_excluded
        + left_selected * right_excluded
        + left_excluded * right_selected
    )
    return [int(value) for value in reversed(full.to_dense())]


def exact_valley(coefficients: list[int]) -> dict[str, Any] | None:
    """Return the first exact descent with any later exact ascent."""
    first: int | None = None
    for k in range(len(coefficients) - 1):
        if first is None and coefficients[k] > coefficients[k + 1]:
            first = k
            continue
        if first is not None and coefficients[k] < coefficients[k + 1]:
            return {
                "first_descent_at": first,
                "later_ascent_at": k,
                "descent_pair": [coefficients[first], coefficients[first + 1]],
                "ascent_pair": [coefficients[k], coefficients[k + 1]],
                "coefficient_window": coefficients[max(0, first - 2) : min(len(coefficients), k + 4)],
                "window_start": max(0, first - 2),
            }
    return None


def build_joined_tree(left: Branching, right: Branching) -> tuple[list[tuple[int, int]], int]:
    edges: list[tuple[int, int]] = []
    next_vertex = 0

    def build(branching: Branching) -> int:
        nonlocal next_vertex
        root = next_vertex
        next_vertex += 1
        if not branching:
            return root
        for _ in range(branching[-1]):
            child = build(branching[:-1])
            edges.append((root, child))
        return root

    left_root = build(left)
    right_root = build(right)
    edges.append((left_root, right_root))
    return edges, next_vertex


def treehood_certificate(left: Branching, right: Branching) -> dict[str, Any]:
    edges, order = build_joined_tree(left, right)
    adjacency = [[] for _ in range(order)]
    for u, v in edges:
        if u == v:
            raise AssertionError("self-loop in constructed tree")
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
    normalized = "\n".join(f"{min(u, v)} {max(u, v)}" for u, v in sorted(edges))
    certificate = {
        "vertices": order,
        "edges": len(edges),
        "connected_vertices": len(seen),
        "simple": len({tuple(sorted(edge)) for edge in edges}) == len(edges),
        "acyclic_by_connected_edge_count": len(seen) == order and len(edges) == order - 1,
        "edge_list_sha256": hashlib.sha256(normalized.encode()).hexdigest(),
    }
    if not (
        certificate["simple"]
        and certificate["connected_vertices"] == order
        and certificate["acyclic_by_connected_edge_count"]
    ):
        raise AssertionError(f"treehood replay failed: {certificate}")
    return certificate


def retain_best(candidates: Iterable[RankedCandidate], limit: int) -> list[RankedCandidate]:
    unique: dict[tuple[Branching, Branching], RankedCandidate] = {}
    for candidate in candidates:
        key = (candidate.left, candidate.right)
        prior = unique.get(key)
        if prior is None or candidate.score > prior.score:
            unique[key] = candidate
    return sorted(unique.values(), key=lambda item: (item.score, item.rebound), reverse=True)[:limit]


def parse_completed_batches(path: Path, run_id: str) -> set[int]:
    if not path.exists():
        return set()
    completed: set[int] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if record.get("run_id") == run_id and record.get("kind") == "batch":
            completed.add(int(record["batch"]))
    return completed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("symmetric", "asymmetric"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=1_000)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--min-depth", type=int, default=3)
    parser.add_argument("--max-depth", type=int, default=14)
    parser.add_argument("--max-branch", type=int, default=12)
    parser.add_argument("--max-order", type=int, default=100_000)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--relative-floor", type=float, default=1e-13)
    parser.add_argument("--descent-tolerance", type=float, default=1e-8)
    parser.add_argument("--ascent-tolerance", type=float, default=1e-7)
    parser.add_argument("--min-separation", type=int, default=1)
    parser.add_argument("--correlated-probability", type=float, default=0.6)
    parser.add_argument("--exact-max-order", type=int, default=20_000)
    args = parser.parse_args()

    if args.min_depth < 2 or args.max_depth < args.min_depth:
        parser.error("require 2 <= min-depth <= max-depth")
    if args.max_order > args.exact_max_order:
        # A float crossing above this cap is still replayed; this flag is only
        # a warning threshold recorded in the ledger, never permission to skip.
        exact_warning = "crossings above exact-max-order may be slow but are still replayed"
    else:
        exact_warning = None

    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    completed = parse_completed_batches(args.output, args.run_id)
    append_record(
        args.output,
        {
            "kind": "resume" if completed else "start",
            "run_id": args.run_id,
            "at": utc_now(),
            "completed_batches": sorted(completed),
            "parameters": {
                key: str(value) if isinstance(value, Path) else value
                for key, value in vars(args).items()
            },
            "script_sha256": script_hash,
            "exact_warning": exact_warning,
        },
    )

    global_best: list[RankedCandidate] = []
    tested_total = 0
    for batch in range(args.batches):
        if batch in completed:
            continue
        started = time.monotonic()
        rng = random.Random(args.seed + 1_000_003 * batch)
        batch_best: list[RankedCandidate] = []
        tested = 0
        rejected = 0
        for _ in range(args.batch_size):
            left, right = sample_pair(args, rng)
            try:
                poly, order = joined_polynomial(left, right)
            except (MemoryError, ValueError):
                rejected += 1
                continue
            ranked = strict_post_descent_pressure(
                poly,
                relative_floor=args.relative_floor,
                descent_tolerance=args.descent_tolerance,
                min_separation=args.min_separation,
            )
            tested += 1
            if ranked is None:
                continue
            candidate = RankedCandidate(
                score=float(ranked["best_later_ratio"]),
                left=left,
                right=right,
                order=order,
                first_descent_at=int(ranked["first_descent_at"]),
                first_descent_ratio=float(ranked["first_descent_ratio"]),
                best_later_at=int(ranked["best_later_at"]),
                best_later_ratio=float(ranked["best_later_ratio"]),
                rebound=float(ranked["rebound"]),
            )
            batch_best = retain_best((*batch_best, candidate), args.top)
            if candidate.score > 1.0 + args.ascent_tolerance:
                exact_started = time.monotonic()
                coefficients = exact_joined_polynomial(left, right)
                valley = exact_valley(coefficients)
                exact_record = {
                    "kind": "exact_replay",
                    "run_id": args.run_id,
                    "batch": batch,
                    "at": utc_now(),
                    "candidate": asdict(candidate),
                    "exact_seconds": time.monotonic() - exact_started,
                    "exact_degree": len(coefficients) - 1,
                    "exact_valley": valley,
                }
                if valley is not None:
                    exact_record["treehood"] = treehood_certificate(left, right)
                    append_record(args.output, exact_record)
                    print(json.dumps(exact_record, sort_keys=True), flush=True)
                    return
                append_record(args.output, exact_record)
        tested_total += tested
        global_best = retain_best((*global_best, *batch_best), args.top)
        checkpoint = {
            "kind": "batch",
            "run_id": args.run_id,
            "batch": batch,
            "at": utc_now(),
            "tested": tested,
            "rejected": rejected,
            "tested_total_this_invocation": tested_total,
            "elapsed_seconds": time.monotonic() - started,
            "batch_best": [asdict(item) for item in batch_best],
            "global_best_this_invocation": [asdict(item) for item in global_best],
        }
        append_record(args.output, checkpoint)
        print(json.dumps(checkpoint, sort_keys=True), flush=True)

    append_record(
        args.output,
        {
            "kind": "complete",
            "run_id": args.run_id,
            "at": utc_now(),
            "tested_total_this_invocation": tested_total,
            "best_this_invocation": [asdict(item) for item in global_best],
        },
    )


if __name__ == "__main__":
    main()
