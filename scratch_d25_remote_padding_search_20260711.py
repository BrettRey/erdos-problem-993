#!/usr/bin/env python3
"""Exact remote-binomial padding search for a prefix fork obstruction.

Start from the 105-vertex deep-tail oriented-fork witness built from four
copies of the first n=26 non-log-concave tree.  Modify only the parent branch,
at a vertex at distance at least two from the central attachment root, by
adding either s leaves or s pendant P2 arms.  This is a tree-realizable way to
convolve the local cavity states with a large near-binomial factor and move the
known r=52 defect toward the prefix window.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import ceil

from graph6 import parse_graph6
from indpoly import _polyadd, _polymul
from scratch_d25_degree2_state_search_20260711 import rooted_state


def trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def power(poly: list[int], exponent: int) -> list[int]:
    out = [1]
    while exponent:
        if exponent & 1:
            out = _polymul(out, poly)
        exponent //= 2
        if exponent:
            poly = _polymul(poly, poly)
    return trim(out)


def add_padding(
    base: list[list[int]], site: int, count: int, kind: str
) -> list[list[int]]:
    adj = [list(neighbors) for neighbors in base]
    for _ in range(count):
        first = len(adj)
        adj.append([site])
        adj[site].append(first)
        if kind == "p2":
            leaf = len(adj)
            adj.append([first])
            adj[first].append(leaf)
    return adj


def distances(adj: list[list[int]], root: int) -> list[int]:
    out = [-1] * len(adj)
    out[root] = 0
    queue = [root]
    for vertex in queue:
        for neighbor in adj[vertex]:
            if out[neighbor] < 0:
                out[neighbor] = out[vertex] + 1
                queue.append(neighbor)
    return out


def state_lists(adj: list[list[int]], root: int):
    state = rooted_state(adj, root, "padding")
    return map(list, (state.total, state.excluded, state.reserve))


def compose_rows(
    parent_adj: list[list[int]],
    parent_root: int,
    child_state: tuple[list[int], list[int], list[int]],
    child_root_degree: int,
):
    p0, e0, r0 = state_lists(parent_adj, parent_root)
    p1, e1, r1 = child_state
    m = 3
    p1_m = power(p1, m)
    e1_m = power(e1, m)
    whole = trim(_polyadd(_polymul(p0, p1_m), [0] + _polymul(e0, e1_m)))
    ac_poly = _polymul(e0, e1_m)
    ap_poly = _polymul(r0, p1_m)
    au_poly = _polymul(_polymul(r1, p0), power(p1, m - 1))
    jpc_poly = _polymul(_polymul(r0, r1), power(p1, m - 1))
    jcc_poly = _polymul(_polymul(_polymul(r1, r1), p0), p1)
    alpha = len(whole) - 1
    limit = ceil((2 * alpha - 1) / 3)
    best_oriented = None
    best_inverse = None
    first_oriented_failure = None
    first_inverse_failure = None
    parent_degree = len(parent_adj[parent_root]) + 1
    for rank in range(4, max(4, limit - 1)):
        def c(poly):
            return poly[rank] if rank < len(poly) else 0

        n_r, ac, ap, au, jpc, jcc = map(
            c, (whole, ac_poly, ap_poly, au_poly, jpc_poly, jcc_poly)
        )
        j_ordered = 2 * m * jpc + m * (m - 1) * jcc
        child_sum = m * au
        oriented_rhs = child_sum * (child_sum + 2 * (ap + ac))
        oriented_lhs = n_r * j_ordered
        if oriented_rhs:
            ratio = Fraction(oriented_lhs, oriented_rhs)
            row = (ratio, rank, oriented_rhs - oriented_lhs)
            if best_oriented is None or row[0] > best_oriented[0]:
                best_oriented = row
            if ratio > 1 and first_oriented_failure is None:
                first_oriented_failure = row

        sum_a = ap + m * au
        sum_squares = ap * ap + m * au * au
        product_ordered = sum_a * sum_a - sum_squares
        q_ordered = n_r * j_ordered - product_ordered
        inverse_budget = (
            Fraction(ap * ap, parent_degree)
            + m * Fraction(au * au, child_root_degree)
            + ac * sum_a
        )
        if inverse_budget:
            ratio = Fraction(q_ordered, inverse_budget)
            row = (ratio, rank, inverse_budget - q_ordered)
            if best_inverse is None or row[0] > best_inverse[0]:
                best_inverse = row
            if ratio > 1 and first_inverse_failure is None:
                first_inverse_failure = row
    return {
        "alpha": alpha,
        "limit": limit,
        "best_oriented": best_oriented,
        "best_inverse": best_inverse,
        "oriented_failure": first_oriented_failure,
        "inverse_failure": first_inverse_failure,
    }


def main() -> None:
    data = json.load(open("results/analysis_n26.json"))
    _, base = parse_graph6(data["lc_failures"][0]["graph6"].encode())
    child_state = tuple(state_lists(base, 23))
    child_root_degree = len(base[23]) + 1
    dist = distances(base, 0)
    sites = [vertex for vertex, d in enumerate(dist) if d >= 2]
    counts = (60, 80, 100, 120, 150, 180, 220, 280, 360)
    best_oriented = None
    best_inverse = None
    tested = 0
    for kind in ("leaves", "p2"):
        for site in sites:
            for count in counts:
                adj = add_padding(base, site, count, kind)
                result = compose_rows(
                    adj, 0, child_state, child_root_degree
                )
                tested += 1
                payload = {
                    "kind": kind,
                    "site": site,
                    "site_distance": dist[site],
                    "count": count,
                    "branch_order": len(adj),
                    **result,
                }
                for key, current in (
                    ("best_oriented", best_oriented),
                    ("best_inverse", best_inverse),
                ):
                    value = result[key]
                    if value is not None and (
                        current is None or value[0] > current[0]
                    ):
                        if key == "best_oriented":
                            best_oriented = (value[0], payload)
                        else:
                            best_inverse = (value[0], payload)
                if result["oriented_failure"] or result["inverse_failure"]:
                    print({"prefix_failure": payload}, flush=True)
                    raise SystemExit(1)
        print(
            {
                "completed_kind": kind,
                "tested": tested,
                "best_oriented": best_oriented,
                "best_inverse": best_inverse,
            },
            flush=True,
        )
    print(
        {
            "tested": tested,
            "prefix_failure": None,
            "best_oriented": best_oriented,
            "best_inverse": best_inverse,
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
