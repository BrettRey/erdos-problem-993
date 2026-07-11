#!/usr/bin/env python3
"""Exact probes for maximum-matching compression in D13.

A perfect-matching tree can be obtained from a quotient tree Q by replacing
each quotient vertex with an edge (a two-vertex block), then attaching every
quotient edge to one endpoint in each incident block.  Concentrating all
external incidences at one endpoint of every block gives the corona Q o K1.

This script tests two candidate extremal steps for the prefix-GSB ratio:

  (orientation)  R_r(any lift of Q) <= R_r(Q o K1),
  (quotient)     R_r(Q o K1) <= R_r(K_1,a-1 o K1).

All comparisons use exact integers.  The candidate is expected only for
r>=2; rank one has the ordinary star as an exact obstruction.
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
from fractions import Fraction
from math import ceil

from indpoly import independence_poly
from trees import trees_geng_raw

from scratch_d13_block_extremal_20260711 import (
    corona_star_poly,
    gsb_ratio,
    prefix_last,
)


Adj = list[list[int]]


def from_edges(n: int, edges: list[tuple[int, int]]) -> Adj:
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def edge_list(adj: Adj) -> list[tuple[int, int]]:
    return [(u, v) for u, ns in enumerate(adj) for v in ns if u < v]


def lift_quotient(adj: Adj, choices: tuple[int, ...]) -> Adj:
    """Replace every quotient vertex by an edge, using 2 bits/Q-edge."""
    qedges = edge_list(adj)
    assert len(choices) == 2 * len(qedges)
    edges = [(2 * u, 2 * u + 1) for u in range(len(adj))]
    for j, (u, v) in enumerate(qedges):
        edges.append((2 * u + choices[2 * j], 2 * v + choices[2 * j + 1]))
    return from_edges(2 * len(adj), edges)


def canonical_lift_choices(adj: Adj):
    """Yield one incidence-bit representative per block-flip isomorphism.

    Flipping the two vertices inside one matching block complements all
    incidence bits at that quotient vertex and does not change the lifted
    tree up to relabelling.  Fixing the first incidence at every block to zero
    leaves sum_v(deg(v)-1)=|V|-2 free bits.
    """
    qedges = edge_list(adj)
    positions: list[list[int]] = [[] for _ in adj]
    for j, (u, v) in enumerate(qedges):
        positions[u].append(2 * j)
        positions[v].append(2 * j + 1)
    free = [pos for ps in positions for pos in ps[1:]]
    for bits in itertools.product((0, 1), repeat=len(free)):
        choices = [0] * (2 * len(qedges))
        for pos, bit in zip(free, bits):
            choices[pos] = bit
        yield tuple(choices)


def corona_poly(base_adj: Adj) -> list[int]:
    """I(base o K1)=(1+x)^a I_base(x/(1+x)), coefficientwise."""
    from math import comb

    a = len(base_adj)
    p = independence_poly(a, base_adj)
    out = [0] * (a + 1)
    for s, count in enumerate(p):
        for j in range(s, a + 1):
            out[j] += count * comb(a - s, j - s)
    return out


def star_adj(a: int) -> Adj:
    return from_edges(a, [(0, v) for v in range(1, a)])


def scan_quotients(max_alpha: int) -> dict:
    failures = []
    closest = []
    for alpha in range(3, max_alpha + 1):
        target = corona_star_poly(alpha)
        last = prefix_last(alpha)
        count = 0
        for _, adj, raw in trees_geng_raw(alpha):
            count += 1
            p = corona_poly(adj)
            for rank in range(2, last + 1):
                got = gsb_ratio(p, rank)
                bound = gsb_ratio(target, rank)
                if got > bound:
                    failures.append({
                        "alpha": alpha, "rank": rank,
                        "ratio": [got.numerator, got.denominator],
                        "bound": [bound.numerator, bound.denominator],
                        "quotient_graph6": raw.decode(),
                    })
                    return {"failures": failures, "closest": closest}
                closest.append((got / bound, alpha, rank, raw.decode(), got, bound))
        closest.sort(reverse=True, key=lambda row: row[0])
        del closest[20:]
        print(json.dumps({"alpha": alpha, "quotients": count}), flush=True)
    return {
        "failures": failures,
        "closest": [
            {
                "relative": [row[0].numerator, row[0].denominator],
                "alpha": row[1], "rank": row[2], "graph6": row[3],
                "ratio": [row[4].numerator, row[4].denominator],
                "bound": [row[5].numerator, row[5].denominator],
            }
            for row in closest
        ],
    }


def scan_lifts(
    exhaustive_alpha: int, random_samples: int, random_max_alpha: int, seed: int
) -> dict:
    orientation_failures = []
    global_failures = []
    closest = []

    def check(base: Adj, base_g6: str, choices: tuple[int, ...]) -> None:
        alpha = len(base)
        lift = lift_quotient(base, choices)
        p = independence_poly(2 * alpha, lift)
        corona = corona_poly(base)
        global_bound_poly = corona_star_poly(alpha)
        last = prefix_last(alpha)
        for rank in range(2, last + 1):
            got = gsb_ratio(p, rank)
            bound = gsb_ratio(corona, rank)
            global_bound = gsb_ratio(global_bound_poly, rank)
            if got > bound:
                if len(orientation_failures) < 20:
                    orientation_failures.append({
                        "alpha": alpha, "rank": rank,
                        "ratio": [got.numerator, got.denominator],
                        "corona_ratio": [bound.numerator, bound.denominator],
                        "quotient_graph6": base_g6,
                        "choices": choices,
                        "lift_edges": edge_list(lift),
                    })
            if got > global_bound:
                global_failures.append({
                    "alpha": alpha, "rank": rank,
                    "ratio": [got.numerator, got.denominator],
                    "global_bound": [global_bound.numerator, global_bound.denominator],
                    "quotient_graph6": base_g6,
                    "choices": choices,
                    "lift_edges": edge_list(lift),
                })
            closest.append((
                got / global_bound, alpha, rank, base_g6, choices, got,
                global_bound,
            ))
            closest.sort(reverse=True, key=lambda row: row[0])
            del closest[20:]

    for alpha in range(3, exhaustive_alpha + 1):
        lift_count = 0
        for _, base, raw in trees_geng_raw(alpha):
            for choices in itertools.product((0, 1), repeat=2 * (alpha - 1)):
                lift_count += 1
                check(base, raw.decode(), choices)
        print(json.dumps({"alpha": alpha, "lifts": lift_count}), flush=True)

    rng = random.Random(seed)
    bases: dict[int, list[tuple[Adj, str]]] = {}
    for alpha in range(exhaustive_alpha + 1, random_max_alpha + 1):
        # Random labelled quotient trees via a Prüfer decoder local to this file.
        bases[alpha] = []
    for sample in range(random_samples):
        alpha = rng.randint(exhaustive_alpha + 1, random_max_alpha)
        seq = [rng.randrange(alpha) for _ in range(alpha - 2)]
        degree = [1] * alpha
        for v in seq:
            degree[v] += 1
        leaves = [v for v in range(alpha) if degree[v] == 1]
        edges = []
        for v in seq:
            u = min(leaves)
            leaves.remove(u)
            edges.append((u, v))
            degree[u] -= 1
            degree[v] -= 1
            if degree[v] == 1:
                leaves.append(v)
        edges.append(tuple(leaves))
        base = from_edges(alpha, edges)
        choices = tuple(rng.randrange(2) for _ in range(2 * (alpha - 1)))
        check(base, "random-labelled", choices)
        if (sample + 1) % 1000 == 0:
            print(json.dumps({"random_lifts": sample + 1}), flush=True)

    return {
        "orientation_failures": orientation_failures,
        "global_failures": global_failures,
        "closest": [
            {
                "relative": [row[0].numerator, row[0].denominator],
                "alpha": row[1], "rank": row[2], "quotient_graph6": row[3],
                "choices": row[4],
                "ratio": [row[5].numerator, row[5].denominator],
                "corona_ratio": [row[6].numerator, row[6].denominator],
            }
            for row in closest
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    q = sub.add_parser("quotients")
    q.add_argument("--max-alpha", type=int, default=18)
    l = sub.add_parser("lifts")
    l.add_argument("--exhaustive-alpha", type=int, default=7)
    l.add_argument("--random-samples", type=int, default=20000)
    l.add_argument("--random-max-alpha", type=int, default=40)
    l.add_argument("--seed", type=int, default=993)
    args = parser.parse_args()
    if args.mode == "quotients":
        result = scan_quotients(args.max_alpha)
    else:
        result = scan_lifts(
            args.exhaustive_alpha, args.random_samples,
            args.random_max_alpha, args.seed,
        )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
