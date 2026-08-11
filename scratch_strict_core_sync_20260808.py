#!/usr/bin/env python3
"""Exact checks for rooted invariants in the strict-core tree class.

A strict-core tree has a unique maximum independent set N, and every vertex
outside N has degree at least four.  Rooting at an N-leaf gives two rooted
state polynomials at every vertex v:

    E_v = product(F_c),
    J_v = product(E_c),
    F_v = E_v + x J_v.

Here E_v excludes v and F_v is the full independence polynomial of the
rooted subtree.  This script tests corrected partial ratio-domination,
ratio-domination, and synchronization relations using exact integer
arithmetic.  Its output is computational evidence, not a proof.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter
from pathlib import Path

from indpoly import _polyadd, _polymul, is_log_concave
from trees import trees


Poly = list[int]


def trim(poly: Poly) -> Poly:
    """Remove trailing zero coefficients, retaining one coefficient."""
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def shift(poly: Poly, amount: int = 1) -> Poly:
    return [0] * amount + list(poly)


def partial_ratio_leq(left: Poly, right: Poly) -> bool:
    """Return left <= right in the corrected 2023 partial ratio order."""
    left = trim(left)
    right = trim(right)
    degree = max(len(left), len(right)) - 1
    a = left + [0] * (degree + 1 - len(left))
    b = right + [0] * (degree + 1 - len(right))
    delta_left = next(k for k, value in enumerate(a) if value)
    degree_right = max(k for k, value in enumerate(b) if value)
    if delta_left > degree_right + 1:
        return False
    return all(
        a[k] * b[k - 1] <= b[k] * a[k - 1]
        for k in range(1, degree + 1)
    )


def ratio_dominated(left: Poly, right: Poly) -> bool:
    """Return whether right is ratio-dominant over left."""
    return partial_ratio_leq(left, right) and partial_ratio_leq(
        right, shift(left),
    )


def synchronized(left: Poly, right: Poly) -> bool:
    """Return the Gross--Mansour--Tucker--Wang synchronization relation."""
    return (
        is_log_concave(left)
        and is_log_concave(right)
        and partial_ratio_leq(left, shift(right))
        and partial_ratio_leq(right, shift(left))
    )


def partially_synchronized(left: Poly, right: Poly) -> bool:
    """Return the Hu--Wang--Zhao--Zhao partial synchronization relation."""
    degree = max(len(left), len(right)) - 1
    a = list(left) + [0] * (degree + 3 - len(left))
    b = list(right) + [0] * (degree + 3 - len(right))

    def coefficient(poly: Poly, index: int) -> int:
        return 0 if index < 0 else poly[index]

    for lower in range(degree + 2):
        for upper in range(lower, degree + 2):
            lhs = (
                coefficient(a, upper + 1) * coefficient(b, lower - 1)
                + coefficient(a, lower - 1) * coefficient(b, upper + 1)
            )
            rhs = (
                coefficient(a, upper) * coefficient(b, lower)
                + coefficient(a, lower) * coefficient(b, upper)
            )
            if lhs > rhs:
                return False
    return True


def unique_maximum_independent_set(adj: list[list[int]]) -> set[int] | None:
    """Return the unique maximum independent set, or None if non-unique."""
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in adj]
    order = [0]
    for vertex in order:
        for neighbor in adj[vertex]:
            if neighbor == parent[vertex]:
                continue
            parent[neighbor] = vertex
            children[vertex].append(neighbor)
            order.append(neighbor)

    exclude_size = [0] * n
    exclude_count = [0] * n
    include_size = [0] * n
    include_count = [0] * n
    for vertex in reversed(order):
        include_size[vertex] = 1 + sum(
            exclude_size[child] for child in children[vertex]
        )
        include_count[vertex] = 1
        for child in children[vertex]:
            include_count[vertex] *= exclude_count[child]

        exclude_count[vertex] = 1
        for child in children[vertex]:
            if include_size[child] > exclude_size[child]:
                exclude_size[vertex] += include_size[child]
                exclude_count[vertex] *= include_count[child]
            elif exclude_size[child] > include_size[child]:
                exclude_size[vertex] += exclude_size[child]
                exclude_count[vertex] *= exclude_count[child]
            else:
                exclude_size[vertex] += exclude_size[child]
                exclude_count[vertex] *= (
                    include_count[child] + exclude_count[child]
                )

    if include_size[0] > exclude_size[0]:
        root_state = 1
        optimum_count = include_count[0]
    elif exclude_size[0] > include_size[0]:
        root_state = 0
        optimum_count = exclude_count[0]
    else:
        return None
    if optimum_count != 1:
        return None

    independent_set: set[int] = set()
    stack = [(0, root_state)]
    while stack:
        vertex, state = stack.pop()
        if state:
            independent_set.add(vertex)
            stack.extend((child, 0) for child in children[vertex])
            continue
        for child in children[vertex]:
            if include_size[child] > exclude_size[child]:
                stack.append((child, 1))
            elif exclude_size[child] > include_size[child]:
                stack.append((child, 0))
            else:
                raise AssertionError("unique optimum has a tied child state")
    return independent_set


def rooted_states(
    adj: list[list[int]], root: int,
) -> tuple[list[list[int]], list[Poly], list[Poly], list[Poly]]:
    """Return children and the exact E, J, F state polynomials."""
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in adj]
    order = [root]
    for vertex in order:
        for neighbor in adj[vertex]:
            if neighbor == parent[vertex]:
                continue
            parent[neighbor] = vertex
            children[vertex].append(neighbor)
            order.append(neighbor)

    excluded: list[Poly] = [[] for _ in adj]
    included_tail: list[Poly] = [[] for _ in adj]
    full: list[Poly] = [[] for _ in adj]
    for vertex in reversed(order):
        e_poly = [1]
        j_poly = [1]
        for child in children[vertex]:
            e_poly = _polymul(e_poly, full[child])
            j_poly = _polymul(j_poly, excluded[child])
        excluded[vertex] = e_poly
        included_tail[vertex] = j_poly
        full[vertex] = _polyadd(e_poly, shift(j_poly))
    return children, excluded, included_tail, full


def random_labelled_core(
    p_count: int, n_count: int, rng: random.Random,
) -> tuple[list[list[int]], list[str]]:
    """Build a random P/N-labelled tree with N independent."""
    adj: list[list[int]] = [[]]
    labels = ["P"]
    remaining = ["P"] * (p_count - 1) + ["N"] * n_count
    rng.shuffle(remaining)
    for label in remaining:
        allowed = [
            vertex
            for vertex, old_label in enumerate(labels)
            if label == "P" or old_label == "P"
        ]
        parent = rng.choice(allowed)
        vertex = len(adj)
        adj.append([parent])
        adj[parent].append(vertex)
        labels.append(label)
    return adj, labels


def complete_strict_core(
    adj: list[list[int]], labels: list[str], rng: random.Random,
) -> None:
    """Attach N-leaves so each P has degree >=4 and >=2 N-neighbors."""
    for support, label in list(enumerate(labels)):
        if label != "P":
            continue
        n_neighbors = sum(labels[neighbor] == "N" for neighbor in adj[support])
        add = max(0, 2 - n_neighbors, 4 - len(adj[support]))
        if rng.random() < 0.35:
            add += rng.randint(0, 3)
        for _ in range(add):
            leaf = len(adj)
            adj.append([support])
            adj[support].append(leaf)
            labels.append("N")


def encode_failure(
    source: str,
    adj: list[list[int]],
    labels: list[str],
    root: int,
    vertex: int,
    relation: str,
    excluded: Poly,
    included_tail: Poly,
    full: Poly,
) -> dict[str, object]:
    return {
        "source": source,
        "order": len(adj),
        "root": root,
        "vertex": vertex,
        "type": labels[vertex],
        "relation": relation,
        "edges": [
            [left, right]
            for left, neighbors in enumerate(adj)
            for right in neighbors
            if left < right
        ],
        "labels": labels,
        "E": excluded,
        "J": included_tail,
        "F": full,
    }


def analyze_tree(
    source: str,
    adj: list[list[int]],
    labels: list[str],
    root: int,
    stats: Counter[str],
    failures: dict[str, dict[str, object]],
) -> None:
    children, excluded, included_tail, full = rooted_states(adj, root)
    leaf_factor = [1, 1]
    for vertex in range(len(adj)):
        vertex_type = labels[vertex]
        e_poly = excluded[vertex]
        j_poly = included_tail[vertex]
        f_poly = full[vertex]
        n_children = sum(labels[child] == "N" for child in children[vertex])
        relations = {
            "F_sync_E": synchronized(f_poly, e_poly),
            "E_lc": is_log_concave(e_poly),
            "F_lc": is_log_concave(f_poly),
            "E_partial_sync_xJ": partially_synchronized(
                e_poly, shift(j_poly),
            ),
            "E_partial_sync_J": partially_synchronized(e_poly, j_poly),
        }
        if vertex_type == "N":
            relations["N_E_ratio_below_F"] = ratio_dominated(e_poly, f_poly)
        else:
            relations["P_F_ratio_below_leaf_times_E"] = ratio_dominated(
                f_poly, _polymul(leaf_factor, e_poly),
            )
            relations["P_J_partial_below_E"] = partial_ratio_leq(
                j_poly, e_poly,
            )
            relations["P_E_partial_below_shifted_J"] = partial_ratio_leq(
                e_poly, shift(j_poly, n_children),
            )

        stats[f"vertices_{vertex_type}"] += 1
        stats[f"{vertex_type}_n_children_{n_children}"] += 1
        for relation, holds in relations.items():
            stats[f"{relation}_{'pass' if holds else 'fail'}"] += 1
            if not holds and relation not in failures:
                failures[relation] = encode_failure(
                    source,
                    adj,
                    labels,
                    root,
                    vertex,
                    relation,
                    e_poly,
                    j_poly,
                    f_poly,
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=16)
    parser.add_argument("--backend", choices=("auto", "geng", "networkx"), default="auto")
    parser.add_argument("--random-trials", type=int, default=1_000)
    parser.add_argument("--max-p", type=int, default=45)
    parser.add_argument("--seed", type=int, default=20260808)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/strict_core_sync_20260808.json"),
    )
    args = parser.parse_args()

    stats: Counter[str] = Counter()
    failures: dict[str, dict[str, object]] = {}
    qualified_trees = 0
    for order in range(1, args.max_n + 1):
        order_qualified = 0
        for _, adj in trees(order, backend=args.backend):
            maximum_set = unique_maximum_independent_set(adj)
            if maximum_set is None:
                continue
            if any(
                len(adj[vertex]) < 4
                for vertex in range(order)
                if vertex not in maximum_set
            ):
                continue
            roots = [
                vertex
                for vertex in maximum_set
                if len(adj[vertex]) <= 1
            ]
            if not roots:
                continue
            labels = [
                "N" if vertex in maximum_set else "P"
                for vertex in range(order)
            ]
            analyze_tree(
                f"exhaustive_n{order}",
                adj,
                labels,
                min(roots),
                stats,
                failures,
            )
            qualified_trees += 1
            order_qualified += 1
        print(json.dumps({
            "event": "exhaustive_order",
            "order": order,
            "qualified": order_qualified,
        }), flush=True)

    rng = random.Random(args.seed)
    random_max_order = 0
    for trial in range(args.random_trials):
        p_count = rng.randint(1, args.max_p)
        initial_n_count = rng.randint(0, p_count - 1)
        adj, labels = random_labelled_core(p_count, initial_n_count, rng)
        complete_strict_core(adj, labels, rng)
        maximum_set = unique_maximum_independent_set(adj)
        labelled_n = {
            vertex for vertex, label in enumerate(labels) if label == "N"
        }
        if maximum_set != labelled_n:
            raise AssertionError("constructed N is not the unique maximum set")
        roots = [
            vertex
            for vertex in maximum_set
            if len(adj[vertex]) == 1
        ]
        analyze_tree(
            f"random_trial_{trial}",
            adj,
            labels,
            min(roots),
            stats,
            failures,
        )
        random_max_order = max(random_max_order, len(adj))

    result = {
        "claim_scope": "exact computation; evidence only",
        "corrected_partial_ratio_order": (
            "Bautista-Ramos--Guillen-Galvan--Gomez-Salgado (2019), "
            "with the added support condition from their 2023 corrigendum"
        ),
        "max_exhaustive_order": args.max_n,
        "qualified_exhaustive_trees": qualified_trees,
        "random_seed": args.seed,
        "random_trials": args.random_trials,
        "random_max_order": random_max_order,
        "stats": dict(sorted(stats.items())),
        "first_failures": failures,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "qualified_exhaustive_trees": qualified_trees,
        "random_trials": args.random_trials,
        "random_max_order": random_max_order,
        "failure_relations": sorted(failures),
        "out": str(args.out),
    }), flush=True)


if __name__ == "__main__":
    main()
