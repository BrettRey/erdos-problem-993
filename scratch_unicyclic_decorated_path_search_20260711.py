#!/usr/bin/env python3
"""Search one-chord closures of decorated spherical paths.

The block roots form a path.  A chord between two nonadjacent roots makes one
cycle.  Root-state assignments are summed directly (at most 2^6 terms), so
the float filter and exact replay use the same finite state grammar.
"""

from __future__ import annotations

import argparse
import random
from itertools import combinations
from types import SimpleNamespace

import numpy as np
from sympy.polys.domains import ZZ
from sympy.polys.rings import ring

from scratch_d17_checkpointed_search_20260711 import (
    exact_valley,
    rooted_order,
    strict_post_descent_pressure,
)
from scratch_d17_decorated_path_search_20260711 import make_pool
from scratch_d17_spherical_search_20260711 import Scaled, add, multiply, spherical_states


Closure = tuple[tuple[int, int], ...]


def allowed(mask: int, count: int, closure: Closure) -> bool:
    if any(((mask >> i) & 3) == 3 for i in range(count - 1)):
        return False
    return not any((mask >> u & 1) and (mask >> v & 1) for u, v in closure)


class Evaluator:
    def __init__(self) -> None:
        self.cache: dict[tuple[int, ...], tuple[Scaled, Scaled]] = {}

    def states(self, block: tuple[int, ...]) -> tuple[Scaled, Scaled]:
        if block not in self.cache:
            excluded, selected, _ = spherical_states(block)
            self.cache[block] = excluded, selected
        return self.cache[block]

    def polynomials(
        self,
        blocks: tuple[tuple[int, ...], ...],
        closures: list[Closure],
    ) -> dict[Closure, np.ndarray]:
        states = [self.states(block) for block in blocks]
        terms = {0: Scaled(np.array([1.0]), 0.0)}
        for i, pair in enumerate(states):
            terms = {
                mask | (selected << i): multiply(term, pair[selected])
                for mask, term in terms.items()
                for selected in (0, 1)
            }
        out = {}
        for closure in closures:
            total: Scaled | None = None
            for mask, term in terms.items():
                if allowed(mask, len(blocks), closure):
                    total = term if total is None else add(total, term)
            assert total is not None
            values = np.maximum(total.values, 0.0)
            out[closure] = values / values.sum()
        return out


def exact_polynomial(
    blocks: tuple[tuple[int, ...], ...], closure: Closure
) -> list[int]:
    polynomial_ring, x = ring("x", ZZ)

    def states(branching: tuple[int, ...]):
        excluded = polynomial_ring.one
        selected = x
        for branch in branching:
            old_excluded = excluded
            excluded = (excluded + selected) ** branch
            selected = x * old_excluded**branch
        return excluded, selected

    rows = [states(block) for block in blocks]
    total = polynomial_ring.zero
    for mask in range(1 << len(blocks)):
        if not allowed(mask, len(blocks), closure):
            continue
        term = polynomial_ring.one
        for i, pair in enumerate(rows):
            term *= pair[bool(mask >> i & 1)]
        total += term
    return [int(value) for value in reversed(total.to_dense())]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=5000)
    parser.add_argument("--blocks", type=int, choices=(4, 5, 6), default=6)
    parser.add_argument("--pool-size", type=int, default=400)
    parser.add_argument("--min-depth", type=int, default=3)
    parser.add_argument("--max-depth", type=int, default=13)
    parser.add_argument("--max-branch", type=int, default=8)
    parser.add_argument("--max-order", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--extra-edges", type=int, choices=(1, 2), default=1)
    args = parser.parse_args()
    pool_args = SimpleNamespace(**vars(args))
    pool = make_pool(pool_args)
    rng = random.Random(args.seed ^ 0xC1C1E)
    evaluator = Evaluator()
    chords = [pair for pair in combinations(range(args.blocks), 2) if pair[1] > pair[0] + 1]
    closures = [tuple(row) for row in combinations(chords, args.extra_edges)]
    best: tuple[float, object] | None = None
    for sample in range(1, args.samples + 1):
        half = tuple(rng.choice(pool) for _ in range((args.blocks + 1) // 2))
        blocks = (
            half + tuple(reversed(half[:-1]))
            if args.blocks % 2
            else half + tuple(reversed(half))
        )
        order = sum(rooted_order(block) for block in blocks)
        for closure, poly in evaluator.polynomials(blocks, closures).items():
            pressure = strict_post_descent_pressure(
                poly,
                relative_floor=1e-13,
                descent_tolerance=1e-8,
                min_separation=1,
            )
            if pressure is None:
                continue
            score = float(pressure["best_later_ratio"])
            row = {"blocks": blocks, "closure": closure, "order": order, **pressure}
            if best is None or score > best[0]:
                best = score, row
            if score > 1.000001:
                coefficients = exact_polynomial(blocks, closure)
                valley = exact_valley(coefficients)
                print(
                    {
                        "status": "exact_hit" if valley else "float_artifact",
                        **row,
                        "exact_valley": valley,
                        "coefficients": coefficients,
                    },
                    flush=True,
                )
                if valley:
                    return
        if sample % 100 == 0:
            print({"status": "progress", "sample": sample, "best": best}, flush=True)
    print({"status": "passed", "samples": args.samples, "best": best}, flush=True)


if __name__ == "__main__":
    main()
