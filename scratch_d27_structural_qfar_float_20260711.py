#!/usr/bin/env python3
"""Tilted structural-orbit qfar evaluation for a large symmetric tree.

This is a floating falsification filter, not a certificate.  It interns every
directed rooted-component shape, so the 9k-vertex D25 hard composition only
requires polynomial messages for a small number of structural types.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter
from functools import lru_cache

import numpy as np

from scratch_d17_spherical_search_20260711 import (
    Scaled,
    add,
    multiply,
    normalize,
    power,
)
from scripts.analyze_prufer_corpus import make_bautista_ramos_tree, make_galvin_tree


def shifted(item: Scaled, log_activity: float) -> Scaled:
    return Scaled(
        np.concatenate((np.asarray([0.0]), item.values)),
        item.log_scale + log_activity,
    )


def compose_branches(
    parent: list[list[int]], child: list[list[int]], multiplicity: int
) -> list[list[int]]:
    adjacency: list[list[int]] = [[]]

    def append_branch(branch: list[list[int]]) -> None:
        offset = len(adjacency)
        adjacency.extend([[offset + w for w in neighbors] for neighbors in branch])
        adjacency[0].append(offset)
        adjacency[offset].append(0)

    append_branch(parent)
    for _ in range(multiplicity):
        append_branch(child)
    return adjacency


def compose_tree(multiplicity: int) -> list[list[int]]:
    return compose_branches(
        make_galvin_tree(40, 20),
        make_bautista_ramos_tree(4, 6),
        multiplicity,
    )


class StructuralMessages:
    def __init__(self, adjacency: list[list[int]], log_activity: float):
        self.adjacency = adjacency
        self.log_activity = log_activity
        self.signatures: list[tuple[int, ...]] = []
        self.signature_ids: dict[tuple[int, ...], int] = {}

    def intern(self, signature: tuple[int, ...]) -> int:
        if signature not in self.signature_ids:
            self.signature_ids[signature] = len(self.signatures)
            self.signatures.append(signature)
        return self.signature_ids[signature]

    @lru_cache(maxsize=None)
    def directed_type(self, vertex: int, parent: int) -> int:
        children = tuple(
            sorted(
                self.directed_type(neighbor, vertex)
                for neighbor in self.adjacency[vertex]
                if neighbor != parent
            )
        )
        return self.intern(children)

    def full_type(self, vertex: int) -> int:
        return self.intern(
            tuple(
                sorted(
                    self.directed_type(neighbor, vertex)
                    for neighbor in self.adjacency[vertex]
                )
            )
        )

    @lru_cache(maxsize=None)
    def state(self, type_id: int) -> tuple[Scaled, Scaled, Scaled]:
        """Return total P, root-excluded E, and selected-tail R."""
        counts = Counter(self.signatures[type_id])
        excluded = Scaled(np.asarray([1.0]), 0.0)
        reserve = Scaled(np.asarray([1.0]), 0.0)
        for child_type, exponent in counts.items():
            child_total, child_excluded, _ = self.state(child_type)
            excluded = multiply(excluded, power(child_total, exponent))
            reserve = multiply(reserve, power(child_excluded, exponent))
        total = add(excluded, shifted(reserve, self.log_activity))
        return total, excluded, reserve

    @lru_cache(maxsize=None)
    def joint_two(self, type_id: int) -> Scaled:
        """Sum distance-two joint-deletion polynomials at this center."""
        zero = Scaled(np.asarray([1.0]), 0.0)
        once = Scaled(np.asarray([0.0]), 0.0)
        twice = Scaled(np.asarray([0.0]), 0.0)
        for child_type in self.signatures[type_id]:
            unrestricted, _, endpoint = self.state(child_type)
            twice = add(
                multiply(twice, unrestricted), multiply(once, endpoint)
            )
            once = add(
                multiply(once, unrestricted), multiply(zero, endpoint)
            )
            zero = multiply(zero, unrestricted)
        return twice


def coefficient_log(item: Scaled, rank: int) -> float:
    if rank < 0 or rank >= len(item.values) or item.values[rank] <= 0:
        return -math.inf
    return item.log_scale + math.log(float(item.values[rank]))


def logsum(rows: list[tuple[float, float]]) -> float:
    """log(sum weight*exp(logvalue)); rows have positive weights."""
    rows = [(logvalue, weight) for logvalue, weight in rows if weight and math.isfinite(logvalue)]
    if not rows:
        return -math.inf
    top = max(logvalue for logvalue, _ in rows)
    return top + math.log(sum(weight * math.exp(logvalue - top) for logvalue, weight in rows))


def evaluate_adjacency(
    adjacency: list[list[int]], rank: int, activity: float
) -> dict:
    messages = StructuralMessages(adjacency, math.log(activity))
    vertex_types = [messages.full_type(vertex) for vertex in range(len(adjacency))]
    type_counts = Counter(vertex_types)
    whole, _, _ = messages.state(vertex_types[0])
    log_n = coefficient_log(whole, rank)
    log_i2 = coefficient_log(whole, rank + 2)

    a_logs = {
        type_id: coefficient_log(messages.state(type_id)[2], rank)
        for type_id in type_counts
    }
    a_top = max(a_logs.values())
    a_scaled = {
        type_id: math.exp(value - a_top) if math.isfinite(value) else 0.0
        for type_id, value in a_logs.items()
    }
    d1 = sum(count * a_scaled[type_id] for type_id, count in type_counts.items())
    squares = sum(
        count * a_scaled[type_id] ** 2 for type_id, count in type_counts.items()
    )
    edges = sum(
        a_scaled[vertex_types[u]] * a_scaled[vertex_types[v]]
        for u, neighbors in enumerate(adjacency)
        for v in neighbors
        if u < v
    )
    distance_two = 0.0
    for center, neighbors in enumerate(adjacency):
        values = [a_scaled[vertex_types[v]] for v in neighbors]
        distance_two += (sum(values) ** 2 - sum(value * value for value in values)) / 2

    joint_logs = [
        (coefficient_log(messages.joint_two(type_id), rank), count)
        for type_id, count in type_counts.items()
    ]
    log_joint_two = logsum(joint_logs)

    # Scale every term in qfar by exp(2*a_top).  The all-joint term uses
    # i_(r+2), whose tilted coefficient carries two surplus powers of z.
    scale = 2 * a_top
    terms = {
        'N_joint_all': math.exp(
            log_n
            + log_i2
            + math.log((rank + 2) * (rank + 1) / 2)
            - 2 * math.log(activity)
            - scale
        ),
        'N_joint_two': math.exp(log_n + log_joint_two - scale),
        'all_pair_products': (d1 * d1 - squares) / 2,
        'edge_products': edges,
        'distance_two_products': distance_two,
    }
    qfar_scaled = (
        terms['N_joint_all']
        - terms['N_joint_two']
        - terms['all_pair_products']
        + terms['edge_products']
        + terms['distance_two_products']
    )
    return {
        'n': len(adjacency),
        'rank': rank,
        'activity': activity,
        'structural_types': len(messages.signatures),
        'vertex_orbits': len(type_counts),
        'qfar_scaled': qfar_scaled,
        'terms_scaled': terms,
    }


def evaluate(multiplicity: int, rank: int, activity: float) -> dict:
    row = evaluate_adjacency(compose_tree(multiplicity), rank, activity)
    row['m'] = multiplicity
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--multiplicity', type=int, default=48)
    parser.add_argument('--rank', type=int, default=3279)
    parser.add_argument('--activity', type=float, default=1.172169686408871)
    args = parser.parse_args()
    print(evaluate(args.multiplicity, args.rank, args.activity), flush=True)


if __name__ == '__main__':
    main()
