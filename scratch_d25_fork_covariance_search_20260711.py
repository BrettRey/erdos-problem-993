#!/usr/bin/env python3
"""Exact falsification search for the D25 fork-covariance inequality.

For a uniform rank-r independent set, X_v indicates that v is addable and
q2 is the unnormalised sum of Cov(X_u,X_v) over distance-two pairs.  The
candidate is

    2 q2 <= sum_v a_v^2 + 2 sum_(uv in E) a_u a_v,

where a_v is the number of rank-r sets leaving v addable.  Equivalently, the
joint distance-two mass is paid by product-marginal mass at distance at most
two.  This driver tests exact Python-integer arithmetic.
"""

from __future__ import annotations

import argparse
import random
from fractions import Fraction
from math import ceil

from indpoly import _polyadd, independence_poly
from scratch_d18_covariance_search_20260711 import (
    cavity_messages,
    coeff,
    distance_two_pairs,
    product,
    prufer_tree,
)
from targeted import make_T_m_t_1, make_T_m_t_d
from trees import trees


def rows(adj: list[list[int]], prefix_only: bool) -> tuple[int, int, list[dict]]:
    poly = independence_poly(len(adj), adj)
    alpha = len(poly) - 1
    limit = ceil((2 * alpha - 1) / 3)
    total, excluded = cavity_messages(adj)
    one = [product([excluded[(w, v)] for w in adj[v]]) for v in range(len(adj))]
    two_pairs = list(distance_two_pairs(adj))
    joint_two = [0]
    joint_by_center = [[0] for _ in adj]
    for u, center, v in two_pairs:
        factors = [excluded[(w, u)] for w in adj[u] if w != center]
        factors += [excluded[(w, v)] for w in adj[v] if w != center]
        factors += [
            total[(w, center)]
            for w in adj[center]
            if w != u and w != v
        ]
        joint = product(factors)
        joint_two = _polyadd(joint_two, joint)
        joint_by_center[center] = _polyadd(joint_by_center[center], joint)

    out = []
    parent = [-2] * len(adj)
    if adj:
        parent[0] = -1
        order = [0]
        for vertex in order:
            for neighbor in adj[vertex]:
                if parent[neighbor] == -2:
                    parent[neighbor] = vertex
                    order.append(neighbor)
    for rank, n_r in enumerate(poly):
        if prefix_only and not (rank >= 4 and rank <= limit - 2):
            continue
        a = [coeff(p, rank) for p in one]
        p_two = sum(a[u] * a[v] for u, _, v in two_pairs)
        q_two = n_r * coeff(joint_two, rank) - p_two
        sum_squares = sum(x * x for x in a)
        edge_products = sum(
            a[u] * a[v]
            for u, neighbors in enumerate(adj)
            for v in neighbors
            if u < v
        )
        budget = sum_squares + 2 * edge_products
        gap = budget - 2 * q_two
        center_gaps = []
        oriented_center_gaps = []
        all_parent_gaps = []
        inverse_degree_center_gaps = []
        pair_bridge_gaps = []
        neighbor_subpoisson_gaps = []
        for center, neighbors in enumerate(adj):
            closed_sum = a[center] + sum(a[u] for u in neighbors)
            incident_edge_squares = sum(
                (a[center] + a[u]) ** 2 for u in neighbors
            )
            # Twice the proposed centerwise allocation.  Summing and then
            # dividing by two gives exactly the global fork gap.
            center_gap2 = (
                2 * closed_sum * closed_sum
                - incident_edge_squares
                - 4 * n_r * coeff(joint_by_center[center], rank)
            )
            center_gaps.append(center_gap2)
            oriented_budget = closed_sum * closed_sum
            if parent[center] >= 0:
                oriented_budget -= (a[center] + a[parent[center]]) ** 2
            oriented_center_gaps.append(
                oriented_budget
                - 2 * n_r * coeff(joint_by_center[center], rank)
            )
            for proposed_parent in neighbors:
                all_parent_gaps.append(
                    (
                        closed_sum * closed_sum
                        - (a[center] + a[proposed_parent]) ** 2
                        - 2 * n_r * coeff(joint_by_center[center], rank)
                    )
                )
            product_pairs = sum(
                a[u] * a[v]
                for index, u in enumerate(neighbors)
                for v in neighbors[index + 1 :]
            )
            q_center = (
                n_r * coeff(joint_by_center[center], rank) - product_pairs
            )
            neighbor_sum = sum(a[u] for u in neighbors)
            neighbor_subpoisson_gaps.append(
                neighbor_sum * neighbor_sum
                - 2 * n_r * coeff(joint_by_center[center], rank)
            )
            inverse_budget = sum(
                Fraction(a[u] * a[u], len(adj[u])) for u in neighbors
            ) + a[center] * sum(a[u] for u in neighbors)
            inverse_degree_center_gaps.append(inverse_budget - 2 * q_center)
            for endpoint in neighbors:
                for other in neighbors:
                    if endpoint == other:
                        continue
                    pair_joint = product(
                        [
                            excluded[(w, endpoint)]
                            for w in adj[endpoint]
                            if w != center
                        ]
                        + [
                            excluded[(w, other)]
                            for w in adj[other]
                            if w != center
                        ]
                        + [
                            total[(w, center)]
                            for w in neighbors
                            if w not in (endpoint, other)
                        ]
                    )
                    pair_bridge_gaps.append(
                        a[endpoint] * (a[other] + a[center])
                        - n_r * coeff(pair_joint, rank)
                    )
        out.append(
            {
                "rank": rank,
                "N": n_r,
                "q2": q_two,
                "sum_squares": sum_squares,
                "edge_products": edge_products,
                "gap": gap,
                "ratio": Fraction(2 * q_two, budget) if budget else Fraction(0),
                "min_center_gap2": min(center_gaps, default=0),
                "min_center": min(range(len(adj)), key=center_gaps.__getitem__)
                if center_gaps
                else None,
                "min_oriented_center_gap": min(oriented_center_gaps, default=0),
                "min_oriented_center": min(
                    range(len(adj)), key=oriented_center_gaps.__getitem__
                )
                if oriented_center_gaps
                else None,
                "min_all_parent_gap": min(all_parent_gaps, default=0),
                "min_inverse_degree_center_gap": min(
                    inverse_degree_center_gaps, default=Fraction(0)
                ),
                "min_pair_bridge_gap": min(pair_bridge_gaps, default=0),
                "min_neighbor_subpoisson_gap": min(
                    neighbor_subpoisson_gaps, default=0
                ),
            }
        )
    return alpha, limit, out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exhaustive-max-n", type=int, default=0)
    parser.add_argument("--random", type=int, default=0)
    parser.add_argument("--min-random-n", type=int, default=20)
    parser.add_argument("--max-random-n", type=int, default=100)
    parser.add_argument("--galvin-max-n", type=int, default=0)
    parser.add_argument("--generalized-max-n", type=int, default=0)
    parser.add_argument("--all-ranks", action="store_true")
    parser.add_argument("--seed", type=int, default=251993)
    args = parser.parse_args()
    prefix_only = not args.all_ranks
    best = None

    def inspect(label: str, adj: list[list[int]]) -> None:
        nonlocal best
        alpha, limit, entries = rows(adj, prefix_only)
        for row in entries:
            item = (row["ratio"], label, len(adj), alpha, limit, row)
            if best is None or item[0] > best[0]:
                best = item
            if row["gap"] < 0:
                print({"failure": item, "adj": adj}, flush=True)
                raise SystemExit(1)
            if row["min_oriented_center_gap"] < 0:
                print({"oriented_center_failure": item, "adj": adj}, flush=True)
                raise SystemExit(2)
            if row["min_all_parent_gap"] < 0:
                print({"all_parent_failure": item, "adj": adj}, flush=True)
                raise SystemExit(3)
            if row["min_inverse_degree_center_gap"] < 0:
                print({"inverse_degree_failure": item, "adj": adj}, flush=True)
                raise SystemExit(4)
            if row["min_neighbor_subpoisson_gap"] < 0:
                print({"neighbor_subpoisson_failure": item, "adj": adj}, flush=True)
                raise SystemExit(6)

    for n in range(1, args.exhaustive_max_n + 1):
        count = 0
        for graph6, adj in trees(n):
            inspect(f"tree:{graph6}", adj)
            count += 1
        print({"scan_n": n, "trees": count, "best": best}, flush=True)

    rng = random.Random(args.seed)
    for trial in range(args.random):
        n = rng.randint(args.min_random_n, args.max_random_n)
        inspect(f"prufer:{trial}", prufer_tree(n, rng))
        if (trial + 1) % 100 == 0:
            print({"random": trial + 1, "best": best}, flush=True)

    if args.galvin_max_n:
        count = 0
        for t in range(2, 30):
            for m in range(1, 80):
                n, adj = make_T_m_t_1(m, t)
                if n > args.galvin_max_n:
                    break
                inspect(f"T({m},{t},1)", adj)
                count += 1
        print({"galvin": count, "best": best}, flush=True)

    if args.generalized_max_n:
        count = 0
        for d in range(2, 6):
            for t in range(2, 16):
                for m in range(1, 30):
                    n, adj = make_T_m_t_d(m, t, d)
                    if n > args.generalized_max_n:
                        break
                    inspect(f"T({m},{t},{d})", adj)
                    count += 1
        print({"generalized": count, "best": best}, flush=True)

    print({"certificate": "no obstruction", "best": best}, flush=True)


if __name__ == "__main__":
    main()
