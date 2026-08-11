#!/usr/bin/env python3
"""Exact rooted-pair falsification search for the two local EADS rules.

Every edge of a tree is obtained by joining the roots of two rooted trees.
For a rooted side U write

    E_U = I(U-root),       R_U = x I(U-N[root]).

After joining rooted sides U and V, the deletion split at the two endpoints is

    A_U = E_U (E_V+R_V),   B_U = R_U E_V,
    A_V = E_V (E_U+R_U),   B_V = R_V E_U.

This script enumerates exact rooted states from all trees up to a requested
order, deduplicates equal (E,R,root-degree) triples, and samples pairs directly.
It is a counterexample search for the empirical assertions that adjacent bad
vertices cannot both be right-oriented and that a left-oriented bad vertex
cannot have degree at most three.  A hit is replayed on the fully joined tree.
No negative result is a proof.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from dataclasses import dataclass
from pathlib import Path

from indpoly import _polyadd, _polymul
from scratch_eads_local_constraint_breaker_20260808 import signed_orientation
from scratch_existential_split_mode_20260807 import (
    mode_interval,
    tree_vertex_splits,
)
from trees import trees_geng_raw


Poly = tuple[int, ...]


def add(left: Poly, right: Poly) -> Poly:
    return tuple(_polyadd(list(left), list(right)))


def mul(left: Poly, right: Poly) -> Poly:
    return tuple(_polymul(list(left), list(right)))


def separation(left: Poly, right: Poly) -> tuple[str, int]:
    """Orient B=right relative to A=left, returning modal separation."""
    a_modes = mode_interval(left)
    b_modes = mode_interval(right)
    if b_modes[0] > a_modes[1]:
        return "right", b_modes[0] - a_modes[1]
    if a_modes[0] > b_modes[1]:
        return "left", a_modes[0] - b_modes[1]
    return "overlap", 0


@dataclass(frozen=True)
class RootedState:
    order: int
    graph6: str
    edges: tuple[tuple[int, int], ...]
    root: int
    root_degree: int
    excluded: Poly
    included: Poly


def enumerate_states(max_order: int) -> list[RootedState]:
    states: dict[tuple[Poly, Poly, int], RootedState] = {}
    for order in range(1, max_order + 1):
        rows = [(1, [[]], b"@")] if order == 1 else trees_geng_raw(order)
        tree_count = 0
        for _, adjacency, graph6 in rows:
            tree_count += 1
            _, splits = tree_vertex_splits(adjacency)
            edges = tuple(
                (u, v)
                for u, neighbors in enumerate(adjacency)
                for v in neighbors
                if u < v
            )
            for split in splits:
                key = (
                    split.a_poly,
                    split.b_poly,
                    len(adjacency[split.vertex]),
                )
                states.setdefault(key, RootedState(
                    order=order,
                    graph6=graph6.decode("ascii"),
                    edges=edges,
                    root=split.vertex,
                    root_degree=len(adjacency[split.vertex]),
                    excluded=split.a_poly,
                    included=split.b_poly,
                ))
        print(json.dumps({
            "event": "library_order",
            "order": order,
            "trees": tree_count,
            "distinct_states": len(states),
        }), flush=True)
    return list(states.values())


def endpoint_data(
    left: RootedState, right: RootedState,
) -> tuple[tuple[str, int], tuple[str, int]]:
    left_a = mul(left.excluded, add(right.excluded, right.included))
    left_b = mul(left.included, right.excluded)
    right_a = mul(right.excluded, add(left.excluded, left.included))
    right_b = mul(right.included, left.excluded)
    return separation(left_a, left_b), separation(right_a, right_b)


def joined_tree(left: RootedState, right: RootedState) -> list[list[int]]:
    n = left.order + right.order
    adjacency: list[list[int]] = [[] for _ in range(n)]
    for u, v in left.edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    offset = left.order
    for u, v in right.edges:
        adjacency[offset + u].append(offset + v)
        adjacency[offset + v].append(offset + u)
    u = left.root
    v = offset + right.root
    adjacency[u].append(v)
    adjacency[v].append(u)
    return adjacency


def encode_state(state: RootedState) -> dict[str, object]:
    return {
        "order": state.order,
        "graph6": state.graph6,
        "edges": [list(edge) for edge in state.edges],
        "root": state.root,
        "root_degree_before_join": state.root_degree,
        "excluded": list(state.excluded),
        "included": list(state.included),
    }


def encode_hit(
    kind: str, left: RootedState, right: RootedState,
    orientations: tuple[tuple[str, int], tuple[str, int]], sample: int,
) -> dict[str, object]:
    adjacency = joined_tree(left, right)
    full, splits = tree_vertex_splits(adjacency)
    endpoints = (left.root, left.order + right.root)
    replay = [signed_orientation(splits[vertex]) for vertex in endpoints]
    if replay != list(orientations):
        raise AssertionError((orientations, replay))
    return {
        "kind": kind,
        "sample": sample,
        "order": len(adjacency),
        "endpoints": list(endpoints),
        "endpoint_degrees": [len(adjacency[v]) for v in endpoints],
        "orientations": [list(value) for value in orientations],
        "left": encode_state(left),
        "right": encode_state(right),
        "joined_edges": [
            [u, v]
            for u, neighbors in enumerate(adjacency)
            for v in neighbors
            if u < v
        ],
        "independence_polynomial": list(full),
        "endpoint_splits": [{
            "vertex": vertex,
            "a_modes": list(splits[vertex].a_modes),
            "b_modes": list(splits[vertex].b_modes),
            "a_poly": list(splits[vertex].a_poly),
            "b_poly": list(splits[vertex].b_poly),
        } for vertex in endpoints],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rooted-order", type=int, default=15)
    parser.add_argument("--samples", type=int, default=5_000_000)
    parser.add_argument("--seed", type=int, default=202608085)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_rooted_pair_breaker_20260808.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)
    started = time.time()
    states = enumerate_states(args.max_rooted_order)
    # Low-degree candidates are oversampled because they are rare in the full
    # state bank and are the only possible left endpoints of the second rule.
    low_states = [state for state in states if state.root_degree <= 2]
    best_right = None
    best_right_score = -1
    best_low_left = None
    best_low_left_score = -1
    hit = None

    for sample in range(1, args.samples + 1):
        if sample % 2:
            left = rng.choice(states)
            right = rng.choice(states)
        else:
            left = rng.choice(low_states)
            right = rng.choice(states)
        orientations = endpoint_data(left, right)
        left_o, right_o = orientations
        right_score = (
            min(left_o[1], right_o[1])
            if left_o[0] == right_o[0] == "right" else -1
        )
        if right_score > best_right_score:
            best_right_score = right_score
            best_right = (left, right, orientations, sample)
        low_left_score = max(
            (
                value[1]
                for value, state in zip(orientations, (left, right))
                if value[0] == "left" and state.root_degree + 1 <= 3
            ),
            default=-1,
        )
        if low_left_score > best_low_left_score:
            best_low_left_score = low_left_score
            best_low_left = (left, right, orientations, sample)
        if right_score >= 2:
            hit = encode_hit(
                "adjacent_right_bad", left, right, orientations, sample,
            )
            break
        if low_left_score >= 2:
            hit = encode_hit(
                "low_degree_left_bad", left, right, orientations, sample,
            )
            break
        if sample % 250_000 == 0:
            print(json.dumps({
                "event": "samples",
                "sample": sample,
                "best_right_score": best_right_score,
                "best_low_left_score": best_low_left_score,
                "elapsed_seconds": time.time() - started,
            }), flush=True)

    def best_record(row):
        if row is None:
            return None
        left, right, orientations, sample = row
        return {
            "sample": sample,
            "orientations": [list(value) for value in orientations],
            "left": encode_state(left),
            "right": encode_state(right),
        }

    result = {
        "claim_scope": "exact rooted-state pair sampling; evidence only",
        "seed": args.seed,
        "max_rooted_order": args.max_rooted_order,
        "distinct_states": len(states),
        "low_degree_states": len(low_states),
        "samples_completed": sample,
        "counterexample": hit,
        "best_right_score": best_right_score,
        "best_low_degree_left_score": best_low_left_score,
        "best_right": best_record(best_right),
        "best_low_degree_left": best_record(best_low_left),
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "samples": sample,
        "counterexample": hit is not None,
        "best_right_score": best_right_score,
        "best_low_left_score": best_low_left_score,
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if hit is not None else 0)


if __name__ == "__main__":
    main()
