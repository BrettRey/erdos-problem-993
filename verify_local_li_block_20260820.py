#!/usr/bin/env python3
"""Finite injectivity audit of the local Li block, including p=2."""

from itertools import product

import networkx as nx


def spider(lengths: tuple[int, ...]) -> tuple[nx.Graph, list[list[int]]]:
    graph = nx.Graph()
    graph.add_node(0)
    arms: list[list[int]] = []
    nxt = 1
    for length in lengths:
        arm: list[int] = []
        prev = 0
        for _ in range(length):
            graph.add_edge(prev, nxt)
            arm.append(nxt)
            prev = nxt
            nxt += 1
        arms.append(arm)
    return graph, arms


def clan_graph(graph: nx.Graph, alpha: tuple[int, ...]) -> nx.Graph:
    clan = nx.Graph()
    copies: dict[int, list[tuple[int, int]]] = {}
    for vertex, multiplicity in enumerate(alpha):
        copies[vertex] = [(vertex, i) for i in range(multiplicity)]
        clan.add_nodes_from(copies[vertex])
        for i in range(multiplicity):
            for j in range(i):
                clan.add_edge(copies[vertex][i], copies[vertex][j])
    for u, v in graph.edges():
        for x in copies[u]:
            for y in copies[v]:
                clan.add_edge(x, y)
    return clan


def local_transform(alpha: tuple[int, ...], arms: list[list[int]]) -> tuple[int, ...] | None:
    odd: list[tuple[int, int]] = []
    for arm_index, arm in enumerate(arms):
        prefix = 0
        for vertex in arm:
            if alpha[vertex] != 1:
                break
            prefix += 1
        if prefix % 2:
            odd.append((prefix, arm_index))
    if len(odd) < 2:
        return None

    odd.sort()
    length, first = odd[0]
    out = list(alpha)
    out[0] = 0
    if length == 1:
        out[arms[first][0]] = 2
    else:
        second = odd[1][1]
        for pos in range(length):
            out[arms[first][pos]] = 2 if pos % 2 == 0 else 0
        for pos in range(length - 1):
            out[arms[second][pos]] = 2 if pos % 2 == 0 else 0
    return tuple(out)


def main() -> None:
    checked = 0
    for lengths in ((1, 1), (2, 3, 3), (3, 3, 3), (1, 3, 3, 3)):
        graph, arms = spider(lengths)
        images: dict[tuple[int, ...], tuple[int, ...]] = {}
        for tail in product(range(3), repeat=graph.number_of_nodes() - 1):
            alpha = (1,) + tail
            if not nx.is_bipartite(clan_graph(graph, alpha)):
                continue
            image = local_transform(alpha, arms)
            if image is None:
                continue
            assert sum(alpha) == sum(image)
            assert nx.is_bipartite(clan_graph(graph, image))
            assert image not in images or images[image] == alpha
            images[image] = alpha
            checked += 1
    print(f"PASS: local transform injective and weight-preserving on {checked} clan maps")


if __name__ == "__main__":
    main()
