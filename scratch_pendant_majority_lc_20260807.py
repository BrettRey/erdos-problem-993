#!/usr/bin/env python3
"""Stress-test a pendant-leaf-majority log-concavity conjecture.

Let ``Delta`` be a downset on ``b`` labelled core vertices and attach ``q_i``
private pendant leaves to vertex ``i``.  The resulting independence enumerator
is

    F(x) = sum_{S in Delta} x^|S| (1+x)^(ell-q(S)),

where ``ell=sum(q_i)`` and ``q(S)=sum_{i in S} q_i``.  For a graph core,
``Delta`` is its independence complex.  This script deliberately tests the
stronger statement for arbitrary downsets containing every singleton.

The motivating conjecture is that ``ell >= b+2`` forces ``F`` to be
log-concave.  A counterexample is emitted with the complete face set, load
vector, polynomial, and failed minor, so it can be replayed independently.
The calculations use exact Python integers.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path


Poly = tuple[int, ...]


@lru_cache(maxsize=None)
def labelled_downsets(b: int) -> tuple[int, ...]:
    """Return all downsets on ``b`` elements as bitmasks of subset masks."""
    if b == 0:
        # The empty family and {empty set}; both are needed by the recurrence.
        return (0, 1)
    previous = labelled_downsets(b - 1)
    offset = 1 << (b - 1)
    out: list[int] = []
    for lower in previous:
        for upper in previous:
            if upper & ~lower == 0:
                out.append(lower | (upper << offset))
    return tuple(out)


def has_all_vertices(downset: int, b: int) -> bool:
    """Whether the downset contains the empty face and every singleton."""
    required = 1
    for vertex in range(b):
        required |= 1 << (1 << vertex)
    return downset & required == required


def weak_compositions(total: int, parts: int):
    """Yield all ordered weak compositions of ``total`` into ``parts``."""
    if parts == 1:
        yield (total,)
        return
    for cuts in combinations(range(total + parts - 1), parts - 1):
        boundaries = (-1,) + cuts + (total + parts - 1,)
        yield tuple(
            boundaries[j + 1] - boundaries[j] - 1
            for j in range(parts)
        )


def faces_from_mask(downset: int) -> list[int]:
    """Decode a downset bitmask into its subset masks."""
    faces: list[int] = []
    while downset:
        low = downset & -downset
        faces.append(low.bit_length() - 1)
        downset ^= low
    return faces


def pendant_polynomial(faces: list[int], q: tuple[int, ...]) -> Poly:
    """Compute the exact weighted-downset polynomial from the displayed sum."""
    ell = sum(q)
    degree = ell + len(q)
    out = [0] * (degree + 1)
    subset_load = [0] * (1 << len(q))
    for subset in range(1, 1 << len(q)):
        low = subset & -subset
        vertex = low.bit_length() - 1
        subset_load[subset] = subset_load[subset ^ low] + q[vertex]
    for face in faces:
        rank = face.bit_count()
        exponent = ell - subset_load[face]
        for j in range(exponent + 1):
            out[rank + j] += comb(exponent, j)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def first_failed_minor(poly: Poly) -> tuple[int, int, int] | None:
    """Return ``(k, a_k^2, a_{k-1}a_{k+1})`` for the first LC failure."""
    for k in range(1, len(poly) - 1):
        lhs = poly[k] * poly[k]
        rhs = poly[k - 1] * poly[k + 1]
        if lhs < rhs:
            return k, lhs, rhs
    return None


def worst_minor(poly: Poly) -> tuple[int, int, int]:
    """Return ``(k, lhs, rhs)`` maximizing the exact ratio ``rhs/lhs``."""
    if len(poly) < 3:
        return 0, 1, 0
    worst = (1, poly[1] * poly[1], poly[0] * poly[2])
    for k in range(2, len(poly) - 1):
        lhs = poly[k] * poly[k]
        rhs = poly[k - 1] * poly[k + 1]
        if rhs * worst[1] > worst[2] * lhs:
            worst = (k, lhs, rhs)
    return worst


def encode_faces(faces: list[int], b: int) -> list[list[int]]:
    return [
        [vertex for vertex in range(b) if face & (1 << vertex)]
        for face in faces
    ]


def random_downset(b: int, rng: random.Random) -> int:
    """Generate a vertex-complete downset from random maximal-face seeds."""
    mask = 1
    for vertex in range(b):
        mask |= 1 << (1 << vertex)
    seed_count = rng.randint(1, max(2, 2 * b))
    for _ in range(seed_count):
        probability = rng.uniform(0.05, 0.95)
        facet = 0
        for vertex in range(b):
            if rng.random() < probability:
                facet |= 1 << vertex
        subset = facet
        while True:
            mask |= 1 << subset
            if subset == 0:
                break
            subset = (subset - 1) & facet
    return mask


def random_composition(total: int, parts: int, rng: random.Random) -> tuple[int, ...]:
    cuts = sorted(rng.sample(range(total + parts - 1), parts - 1))
    boundaries = (-1,) + tuple(cuts) + (total + parts - 1,)
    return tuple(
        boundaries[j + 1] - boundaries[j] - 1
        for j in range(parts)
    )


def make_failure(
    b: int,
    ell: int,
    downset: int,
    q: tuple[int, ...],
    poly: Poly,
    failed: tuple[int, int, int],
) -> dict[str, object]:
    faces = faces_from_mask(downset)
    k, lhs, rhs = failed
    return {
        "b": b,
        "ell": ell,
        "q": list(q),
        "faces": encode_faces(faces, b),
        "polynomial": list(poly),
        "failed_minor": {"k": k, "lhs": lhs, "rhs": rhs},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-exhaustive-b", type=int, default=1)
    parser.add_argument("--exhaustive-b", type=int, default=5)
    parser.add_argument("--random-trials", type=int, default=200_000)
    parser.add_argument("--random-max-b", type=int, default=16)
    parser.add_argument("--extra-leaves", type=int, default=2)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/pendant_majority_downsets_20260807.json"),
    )
    args = parser.parse_args()

    started = time.time()
    tested = 0
    failure: dict[str, object] | None = None
    worst: dict[str, object] | None = None
    worst_lhs = 1
    worst_rhs = 0
    exhaustive_orders: list[dict[str, int]] = []

    for b in range(args.min_exhaustive_b, args.exhaustive_b + 1):
        ell = b + args.extra_leaves
        if ell < 0:
            raise ValueError("b + extra-leaves must be nonnegative")
        downsets = [
            family for family in labelled_downsets(b)
            if has_all_vertices(family, b)
        ]
        compositions = list(weak_compositions(ell, b))
        order_tests = 0
        for downset in downsets:
            faces = faces_from_mask(downset)
            for q in compositions:
                poly = pendant_polynomial(faces, q)
                failed = first_failed_minor(poly)
                k, lhs, rhs = worst_minor(poly)
                if rhs * worst_lhs > worst_rhs * lhs:
                    worst_lhs, worst_rhs = lhs, rhs
                    worst = make_failure(b, ell, downset, q, poly, (k, lhs, rhs))
                tested += 1
                order_tests += 1
                if failed is not None:
                    failure = make_failure(b, ell, downset, q, poly, failed)
                    break
            if failure is not None:
                break
        exhaustive_orders.append({
            "b": b,
            "downsets": len(downsets),
            "compositions": len(compositions),
            "tests": order_tests,
        })
        print(exhaustive_orders[-1], flush=True)
        if failure is not None:
            break

    random_tests = 0
    if failure is None:
        rng = random.Random(args.seed)
        for _ in range(args.random_trials):
            b = rng.randint(max(1, args.exhaustive_b + 1), args.random_max_b)
            ell = b + args.extra_leaves + rng.randint(0, 3 * b)
            downset = random_downset(b, rng)
            q = random_composition(ell, b, rng)
            faces = faces_from_mask(downset)
            poly = pendant_polynomial(faces, q)
            failed = first_failed_minor(poly)
            k, lhs, rhs = worst_minor(poly)
            if rhs * worst_lhs > worst_rhs * lhs:
                worst_lhs, worst_rhs = lhs, rhs
                worst = make_failure(b, ell, downset, q, poly, (k, lhs, rhs))
            tested += 1
            random_tests += 1
            if failed is not None:
                failure = make_failure(b, ell, downset, q, poly, failed)
                break

    summary: dict[str, object] = {
        "claim_scope": "exact computation; evidence only",
        "statement_tested": (
            "for every tested vertex-complete downset Delta on b vertices "
            "and nonnegative q with sum(q)=b+extra_leaves, the polynomial "
            "sum_{S in Delta} x^|S|(1+x)^(sum(q)-q(S)) is log-concave"
        ),
        "exhaustive_extra_leaves": args.extra_leaves,
        "exhaustive_orders": exhaustive_orders,
        "random_seed": args.seed,
        "random_tests": random_tests,
        "total_tests": tested,
        "failure": failure,
        "worst_minor": worst,
        "worst_ratio": worst_rhs / worst_lhs,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    raise SystemExit(1 if failure is not None else 0)


if __name__ == "__main__":
    main()
