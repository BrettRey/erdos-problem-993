#!/usr/bin/env python3
"""Verify the maximum-matching bag representation of a tree polynomial."""

import networkx as nx

from indpoly import independence_poly


def add(a: list[int], b: list[int]) -> list[int]:
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(n)]


def mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def bag_polynomial(tree: nx.Graph) -> tuple[list[int], int, int]:
    n = tree.number_of_nodes()
    matching = nx.max_weight_matching(tree, maxcardinality=True)
    mate: dict[int, int] = {}
    for u, v in matching:
        mate[u] = v
        mate[v] = u

    bags: list[tuple[int, ...]] = []
    bag_of: dict[int, int] = {}
    seen: set[int] = set()
    for v in tree.nodes():
        if v in seen:
            continue
        bag = (v, mate[v]) if v in mate else (v,)
        bag = tuple(sorted(bag))
        idx = len(bags)
        bags.append(bag)
        for w in bag:
            bag_of[w] = idx
            seen.add(w)

    skeleton = nx.Graph()
    skeleton.add_nodes_from(range(len(bags)))
    port: dict[tuple[int, int], tuple[int, int]] = {}
    for u, v in tree.edges():
        a, b = bag_of[u], bag_of[v]
        if a == b:
            continue
        skeleton.add_edge(a, b)
        port[(a, b)] = (u, v)
        port[(b, a)] = (v, u)
    assert nx.is_tree(skeleton)

    root = 0
    parent = {root: -1}
    order = [root]
    for v in order:
        for w in skeleton.neighbors(v):
            if w != parent[v]:
                parent[w] = v
                order.append(w)

    states = [bag + (None,) for bag in bags]
    dp: list[dict[int | None, list[int]]] = [{} for _ in bags]
    for v in reversed(order):
        for state in states[v]:
            poly = [1] if state is None else [0, 1]
            for w in skeleton.neighbors(v):
                if parent.get(w) != v:
                    continue
                pv, pw = port[(v, w)]
                child = [0]
                for child_state, child_poly in dp[w].items():
                    if state == pv and child_state == pw:
                        continue
                    child = add(child, child_poly)
                poly = mul(poly, child)
            dp[v][state] = poly

    total = [0]
    for poly in dp[root].values():
        total = add(total, poly)
    singles = sum(len(bag) == 1 for bag in bags)
    return total, len(bags), singles


def main() -> None:
    count = 0
    for n in range(2, 13):
        for tree in nx.generators.nonisomorphic_trees(n):
            adj = [list(tree.neighbors(v)) for v in range(n)]
            direct = independence_poly(n, adj)
            via_bags, alpha, deficiency = bag_polynomial(tree)
            assert via_bags == direct
            assert alpha == len(direct) - 1
            assert deficiency == 2 * alpha - n
            count += 1
    print(f"PASS: matching-bag polynomial and deficiency on {count} trees (n<=12)")


if __name__ == "__main__":
    main()
