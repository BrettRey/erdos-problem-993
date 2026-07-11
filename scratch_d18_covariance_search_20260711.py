#!/usr/bin/env python3
"""Exact D18 distance-class covariance profiles for uniform independent sets.

For a uniform independent r-set S, X_v indicates that v can be adjoined to S.
This script computes exact numerators for

    C2   = sum_{dist(u,v)=2} Cov(X_u, X_v),
    Cfar = sum_{dist(u,v)>=3} Cov(X_u, X_v).

The joint counts at distance two are obtained from directed cavity messages;
no subset enumeration is used.  This is a falsification/identity tool for the
candidate prefix bounds Cfar <= 0 and 2 C2 <= E[sum_v X_v].
"""

from __future__ import annotations

import argparse
import random
from functools import lru_cache
from math import ceil, comb

from indpoly import _polyadd, _polymul, independence_poly
from trees import trees


Poly = list[int]


def coeff(poly: Poly, rank: int) -> int:
    return poly[rank] if 0 <= rank < len(poly) else 0


def product(polys: list[Poly]) -> Poly:
    out = [1]
    for poly in polys:
        out = _polymul(out, poly)
    return out


def poly_subtract(left: Poly, right: Poly) -> Poly:
    out = [0] * max(len(left), len(right))
    for i, value in enumerate(left):
        out[i] += value
    for i, value in enumerate(right):
        out[i] -= value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def induced_poly(adj: list[list[int]], removed: set[int]) -> Poly:
    """Direct deletion implementation, used only to audit cavity formulas."""
    kept = [v for v in range(len(adj)) if v not in removed]
    index = {v: i for i, v in enumerate(kept)}
    induced = [
        [index[w] for w in adj[v] if w in index]
        for v in kept
    ]
    return independence_poly(len(kept), induced)


def cavity_messages(
    adj: list[list[int]],
) -> tuple[dict[tuple[int, int], Poly], dict[tuple[int, int], Poly]]:
    """Return P_{u->p} and E_{u->p} for every directed edge.

    P is the independence polynomial of the u-side component of T-up;
    E is the same polynomial conditioned on u being excluded.
    """

    @lru_cache(maxsize=None)
    def message(u: int, parent: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
        child_messages = [message(w, u) for w in adj[u] if w != parent]
        excluded = product([list(p) for p, _ in child_messages])
        selected_tail = product([list(e) for _, e in child_messages])
        total = _polyadd(excluded, [0] + selected_tail)
        return tuple(total), tuple(excluded)

    total: dict[tuple[int, int], Poly] = {}
    excluded: dict[tuple[int, int], Poly] = {}
    for u in range(len(adj)):
        for parent in adj[u]:
            p, e = message(u, parent)
            total[(u, parent)] = list(p)
            excluded[(u, parent)] = list(e)
    return total, excluded


def closed_neighborhood(adj: list[list[int]], vertex: int) -> set[int]:
    return {vertex, *adj[vertex]}


def distance_two_pairs(adj: list[list[int]]):
    """Yield each distance-two pair once, together with its unique center."""
    for center, neighbors in enumerate(adj):
        for i, u in enumerate(neighbors):
            for v in neighbors[i + 1 :]:
                yield u, center, v


def unique_path(adj: list[list[int]], source: int, target: int) -> list[int]:
    parent = [-1] * len(adj)
    parent[source] = source
    queue = [source]
    for vertex in queue:
        if vertex == target:
            break
        for neighbor in adj[vertex]:
            if parent[neighbor] == -1:
                parent[neighbor] = vertex
                queue.append(neighbor)
    path = [target]
    while path[-1] != source:
        path.append(parent[path[-1]])
    path.reverse()
    return path


def path_rayleigh_factor(
    adj: list[list[int]],
    path: list[int],
    total: dict[tuple[int, int], Poly],
    excluded: dict[tuple[int, int], Poly],
) -> Poly:
    """Return the positive factor W_P in the tree path identity.

    If the path has d edges and endpoints u,v, then

      A_u A_v - I(T) A_uv = (-1)^(d+1) x^(d-1) W_P,

    where A_u=I(T-N[u]) and A_uv=I(T-(N[u] union N[v])).
    """
    factors: list[Poly] = []
    for i, vertex in enumerate(path):
        previous = path[i - 1] if i else -1
        following = path[i + 1] if i + 1 < len(path) else -1
        side_neighbors = [
            neighbor
            for neighbor in adj[vertex]
            if neighbor != previous and neighbor != following
        ]
        a_i = product([total[(neighbor, vertex)] for neighbor in side_neighbors])
        b_i = product(
            [excluded[(neighbor, vertex)] for neighbor in side_neighbors]
        )
        factors.extend((a_i, b_i))
    return product(factors)


def audit_cavity_factorization(adj: list[list[int]]) -> None:
    total, excluded = cavity_messages(adj)
    for v in range(len(adj)):
        cavity = product([excluded[(w, v)] for w in adj[v]])
        direct = induced_poly(adj, closed_neighborhood(adj, v))
        if cavity != direct:
            raise AssertionError(("one", v, cavity, direct))

    for u, center, v in distance_two_pairs(adj):
        factors = [
            excluded[(w, u)] for w in adj[u] if w != center
        ]
        factors += [
            excluded[(w, v)] for w in adj[v] if w != center
        ]
        factors += [
            total[(w, center)]
            for w in adj[center]
            if w != u and w != v
        ]
        cavity = product(factors)
        removed = closed_neighborhood(adj, u) | closed_neighborhood(adj, v)
        direct = induced_poly(adj, removed)
        if cavity != direct:
            raise AssertionError(("two", u, center, v, cavity, direct))

    # Christoffel--Darboux/Cassini factorization along the unique tree path.
    tree_poly = independence_poly(len(adj), adj)
    one_polys = [
        product([excluded[(w, v)] for w in adj[v]])
        for v in range(len(adj))
    ]
    distances = distance_classes(adj)
    for u in range(len(adj)):
        for v in range(u + 1, len(adj)):
            distance = distances[u][v]
            if distance < 2:
                continue
            removed = closed_neighborhood(adj, u) | closed_neighborhood(adj, v)
            joint = induced_poly(adj, removed)
            left = poly_subtract(
                _polymul(one_polys[u], one_polys[v]),
                _polymul(tree_poly, joint),
            )
            path = unique_path(adj, u, v)
            factor = path_rayleigh_factor(adj, path, total, excluded)
            right = [0] * (distance - 1) + factor
            if (distance + 1) % 2:
                right = [-value for value in right]
            while len(right) > 1 and right[-1] == 0:
                right.pop()
            if left != right:
                raise AssertionError(
                    ("path", u, v, distance, left, right, path)
                )


def distance_classes(adj: list[list[int]]) -> list[list[int]]:
    n = len(adj)
    distances = [[-1] * n for _ in range(n)]
    for source in range(n):
        distances[source][source] = 0
        queue = [source]
        for vertex in queue:
            for neighbor in adj[vertex]:
                if distances[source][neighbor] == -1:
                    distances[source][neighbor] = distances[source][vertex] + 1
                    queue.append(neighbor)
    return distances


def covariance_profile(adj: list[list[int]]) -> dict:
    n = len(adj)
    poly = independence_poly(n, adj)
    alpha = len(poly) - 1
    limit = ceil((2 * alpha - 1) / 3)
    total, excluded = cavity_messages(adj)

    one_polys = [
        product([excluded[(w, v)] for w in adj[v]])
        for v in range(n)
    ]

    joint_two_poly = [0]
    for u, center, v in distance_two_pairs(adj):
        factors = [
            excluded[(w, u)] for w in adj[u] if w != center
        ]
        factors += [
            excluded[(w, v)] for w in adj[v] if w != center
        ]
        factors += [
            total[(w, center)]
            for w in adj[center]
            if w != u and w != v
        ]
        joint_two_poly = _polyadd(joint_two_poly, product(factors))

    distances = distance_classes(adj)
    rows = []
    for rank, count in enumerate(poly):
        a = [coeff(p, rank) for p in one_polys]
        d1 = sum(a)
        p_edge = 0
        p_two = 0
        p_far = 0
        for u in range(n):
            for v in range(u + 1, n):
                value = a[u] * a[v]
                distance = distances[u][v]
                if distance == 1:
                    p_edge += value
                elif distance == 2:
                    p_two += value
                else:
                    p_far += value

        j_two = coeff(joint_two_poly, rank)
        j_all = comb(rank + 2, 2) * coeff(poly, rank + 2)
        j_far = j_all - j_two
        q_two = count * j_two - p_two
        q_far = count * j_far - p_far
        rows.append(
            {
                "rank": rank,
                "N": count,
                "d1": d1,
                "q2": q_two,
                "qfar": q_far,
                "two_bound_gap": count * d1 - 2 * q_two,
                "prefix": rank <= limit - 2,
            }
        )
    return {"n": n, "alpha": alpha, "limit": limit, "rows": rows}


def prufer_tree(n: int, rng: random.Random) -> list[list[int]]:
    import networkx as nx

    graph = nx.from_prufer_sequence([rng.randrange(n) for _ in range(n - 2)])
    return [sorted(graph.neighbors(v)) for v in range(n)]


def known_obstruction_certificate() -> None:
    """Replay the exact tail and local-polarization obstructions."""
    from targeted import make_T_m_t_1

    # The smaller Galvin tree shows why the far statement must be windowed.
    _, small = make_T_m_t_1(3, 4)
    small_poly = independence_poly(len(small), small)
    small_profile = covariance_profile(small)
    tail = small_profile["rows"][13]
    assert small_profile["alpha"] == 15
    assert small_profile["limit"] == 10
    assert not tail["prefix"]
    assert tail["N"] == 5410
    assert tail["d1"] == 840
    assert tail["q2"] == 47262
    assert tail["qfar"] == 196002

    total, excluded = cavity_messages(small)
    one = [
        product([excluded[(w, v)] for w in small[v]])
        for v in range(len(small))
    ]
    distances = distance_classes(small)
    by_distance: dict[int, int] = {}
    for u in range(len(small)):
        for v in range(u + 1, len(small)):
            distance = distances[u][v]
            if distance < 3:
                continue
            joint = induced_poly(
                small,
                closed_neighborhood(small, u) | closed_neighborhood(small, v),
            )
            numerator = (
                small_poly[13] * coeff(joint, 13)
                - coeff(one[u], 13) * coeff(one[v], 13)
            )
            by_distance[distance] = by_distance.get(distance, 0) + numerator
    assert by_distance == {
        3: -23136,
        4: 123474,
        5: -39168,
        6: 134832,
    }

    # A rooted branch-ratio sign needed by naive bivariate shortening reverses
    # well inside the global prefix on the larger Galvin tree.
    _, large = make_T_m_t_1(14, 8)
    large_poly = independence_poly(len(large), large)
    alpha = len(large_poly) - 1
    limit = ceil((2 * alpha - 1) / 3)
    assert alpha == 126 and limit == 84 and 35 <= limit - 2
    total, excluded = cavity_messages(large)
    root = 0
    root_excluded = product([total[(w, root)] for w in large[root]])
    root_link = product([excluded[(w, root)] for w in large[root]])
    determinant = (
        root_link[35] * root_excluded[34]
        - root_excluded[35] * root_link[34]
    )
    assert determinant == int(
        "6335498169539288470952778004241096231168434068322913107924643509073457709056"
    )
    assert determinant > 0

    print(
        {
            "tail_far_obstruction": by_distance,
            "tail_qfar": tail["qfar"],
            "prefix_branch_polarization_obstruction": determinant,
            "certificate": "passed",
        },
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-max-n", type=int, default=0)
    parser.add_argument("--exhaustive-max-n", type=int, default=0)
    parser.add_argument("--random", type=int, default=0)
    parser.add_argument("--min-random-n", type=int, default=20)
    parser.add_argument("--max-random-n", type=int, default=80)
    parser.add_argument("--seed", type=int, default=180993)
    parser.add_argument("--known-obstructions", action="store_true")
    args = parser.parse_args()

    if args.known_obstructions:
        known_obstruction_certificate()

    if args.audit_max_n:
        for n in range(1, args.audit_max_n + 1):
            count = 0
            for _, adj in trees(n):
                audit_cavity_factorization(adj)
                count += 1
            print({"audit_n": n, "trees": count}, flush=True)

    if args.exhaustive_max_n:
        for n in range(1, args.exhaustive_max_n + 1):
            count = 0
            for _, adj in trees(n):
                profile = covariance_profile(adj)
                for row in profile["rows"]:
                    if row["prefix"] and row["rank"] >= 4:
                        if row["qfar"] > 0 or row["two_bound_gap"] < 0:
                            print({"failure": profile, "adj": adj})
                            raise SystemExit(1)
                count += 1
            print({"scan_n": n, "trees": count, "failures": 0}, flush=True)

    if args.random:
        rng = random.Random(args.seed)
        best_far = None
        best_two = None
        for trial in range(args.random):
            n = rng.randint(args.min_random_n, args.max_random_n)
            adj = prufer_tree(n, rng)
            profile = covariance_profile(adj)
            for row in profile["rows"]:
                if not row["prefix"] or row["rank"] < 4 or row["d1"] == 0:
                    continue
                far_key = (row["qfar"] / (row["N"] * row["d1"]),)
                two_key = (2 * row["q2"] / (row["N"] * row["d1"]),)
                payload = (profile["n"], profile["alpha"], row["rank"], adj)
                if best_far is None or far_key > best_far[0]:
                    best_far = (far_key, payload, row)
                if best_two is None or two_key > best_two[0]:
                    best_two = (two_key, payload, row)
                if row["qfar"] > 0 or row["two_bound_gap"] < 0:
                    print({"failure": payload, "row": row}, flush=True)
                    raise SystemExit(1)
            if (trial + 1) % 100 == 0:
                print(
                    {"trials": trial + 1, "best_far": best_far, "best_two": best_two},
                    flush=True,
                )
        print({"trials": args.random, "best_far": best_far, "best_two": best_two})


if __name__ == "__main__":
    main()
