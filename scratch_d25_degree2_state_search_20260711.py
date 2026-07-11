#!/usr/bin/env python3
"""Exact rooted-state falsifier for the oriented D25 fork lemma at degree two."""

from __future__ import annotations

import argparse
from dataclasses import dataclass

from indpoly import _polyadd, _polymul
from trees import trees


@dataclass(frozen=True)
class State:
    total: tuple[int, ...]
    excluded: tuple[int, ...]
    reserve: tuple[int, ...]
    label: str
    order: int


def rooted_state(adj: list[list[int]], root: int, label: str) -> State:
    parent = [-2] * len(adj)
    parent[root] = -1
    order = [root]
    for vertex in order:
        for neighbor in adj[vertex]:
            if parent[neighbor] == -2:
                parent[neighbor] = vertex
                order.append(neighbor)
    total: dict[int, list[int]] = {}
    excluded: dict[int, list[int]] = {}
    reserve: dict[int, list[int]] = {}
    for vertex in reversed(order):
        children = [u for u in adj[vertex] if parent[u] == vertex]
        e = [1]
        r = [1]
        for child in children:
            e = _polymul(e, total[child])
            r = _polymul(r, excluded[child])
        p = _polyadd(e, [0] + r)
        total[vertex], excluded[vertex], reserve[vertex] = p, e, r
    return State(
        tuple(total[root]),
        tuple(excluded[root]),
        tuple(reserve[root]),
        label,
        len(adj),
    )


def coeff(poly: list[int], rank: int) -> int:
    return poly[rank] if rank < len(poly) else 0


def gap(left: State, right: State):
    p0, e0, r0 = map(list, (left.total, left.excluded, left.reserve))
    p1, e1, r1 = map(list, (right.total, right.excluded, right.reserve))
    whole = _polyadd(_polymul(p0, p1), [0] + _polymul(e0, e1))
    a_c_poly = _polymul(e0, e1)
    a_0_poly = _polymul(r0, p1)
    a_1_poly = _polymul(r1, p0)
    joint = _polymul(r0, r1)
    for rank, n_r in enumerate(whole):
        ac = coeff(a_c_poly, rank)
        a0 = coeff(a_0_poly, rank)
        a1 = coeff(a_1_poly, rank)
        j = coeff(joint, rank)
        value = a1 * a1 + 2 * a1 * (a0 + ac) - 2 * n_r * j
        if value < 0:
            return rank, value, n_r, ac, a0, a1, j, whole
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-branch-n", type=int, default=10)
    args = parser.parse_args()
    unique: dict[tuple, State] = {}
    for n in range(1, args.max_branch_n + 1):
        for graph6, adj in trees(n):
            for root in range(n):
                state = rooted_state(adj, root, f"{graph6}@{root}")
                key = (state.total, state.excluded, state.reserve)
                unique.setdefault(key, state)
    states = list(unique.values())
    pairs = 0
    for left in states:
        for right in states:
            pairs += 1
            witness = gap(left, right)
            if witness:
                print(
                    {
                        "left": left,
                        "right": right,
                        "center_order": left.order + right.order + 1,
                        "witness": witness,
                    },
                    flush=True,
                )
                raise SystemExit(1)
    print({"states": len(states), "ordered_pairs": pairs, "failure": None})


if __name__ == "__main__":
    main()
