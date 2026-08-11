#!/usr/bin/env python3
"""Exact kill test for the existential adjacent deletion-split conjecture.

For a forest F and vertex v, put

    A_v(x) = I(F-v; x),
    B_v(x) = x I(F-N[v]; x).

The conjectural certificate is that some v has modes of A_v and B_v at
distance at most one (where either endpoint of a modal plateau may be used).
If this holds for every forest, a smallest non-unimodal forest cannot exist:
A_v and B_v are smaller-forest polynomials, and a sum of two nonnegative
unimodal sequences with adjacent modes is unimodal.

This script is evidence only.  It computes every directed rooted-tree message
exactly, checks the deletion identity A_v+B_v=I(F), and emits a complete JSON
summary.  A failure record contains graph6 component certificates and the full
integer polynomials needed for independent replay.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly
from trees import trees_geng_raw


Poly = tuple[int, ...]


def polyadd(a: Poly, b: Poly) -> Poly:
    return tuple(_polyadd(list(a), list(b)))


def polymul(a: Poly, b: Poly) -> Poly:
    return tuple(_polymul(list(a), list(b)))


def mode_interval(poly: Poly) -> tuple[int, int]:
    peak = max(poly)
    positions = [k for k, value in enumerate(poly) if value == peak]
    return positions[0], positions[-1]


def interval_distance(left: tuple[int, int], right: tuple[int, int]) -> int:
    """Minimum absolute distance between two integer intervals."""
    if left[1] < right[0]:
        return right[0] - left[1]
    if right[1] < left[0]:
        return left[0] - right[1]
    return 0


@dataclass(frozen=True)
class VertexSplit:
    vertex: int
    a_poly: Poly
    b_poly: Poly
    a_modes: tuple[int, int]
    b_modes: tuple[int, int]
    distance: int


def tree_vertex_splits(adj: list[list[int]]) -> tuple[Poly, list[VertexSplit]]:
    """Return I(T) and all exact deletion-split polynomials for one tree."""
    n = len(adj)
    if n == 1:
        full = (1, 1)
        return full, [VertexSplit(0, (1,), (0, 1), (0, 0), (1, 1), 1)]

    @lru_cache(maxsize=None)
    def message(root: int, parent: int) -> tuple[Poly, Poly]:
        # I is the full component polynomial.  E counts configurations with
        # root excluded, equivalently I(component-root).
        excluded: Poly = (1,)
        root_included_tail: Poly = (1,)
        for child in adj[root]:
            if child == parent:
                continue
            child_i, child_e = message(child, root)
            excluded = polymul(excluded, child_i)
            root_included_tail = polymul(root_included_tail, child_e)
        included = (0,) + root_included_tail
        return polyadd(excluded, included), excluded

    # Populate both orientations.  Each directed component is acyclic, so the
    # memoized recursion is exact and terminates without choosing a global root.
    for u in range(n):
        for v in adj[u]:
            message(u, v)

    splits: list[VertexSplit] = []
    full: Poly | None = None
    for v in range(n):
        a_poly: Poly = (1,)
        closed_tail: Poly = (1,)
        for u in adj[v]:
            branch_i, branch_e = message(u, v)
            a_poly = polymul(a_poly, branch_i)
            closed_tail = polymul(closed_tail, branch_e)
        b_poly = (0,) + closed_tail
        total = polyadd(a_poly, b_poly)
        if full is None:
            full = total
        elif total != full:
            raise AssertionError("vertex deletion identities disagree")
        a_modes = mode_interval(a_poly)
        b_modes = mode_interval(b_poly)
        splits.append(VertexSplit(
            vertex=v,
            a_poly=a_poly,
            b_poly=b_poly,
            a_modes=a_modes,
            b_modes=b_modes,
            distance=interval_distance(a_modes, b_modes),
        ))
    assert full is not None
    return full, splits


def convolve_component_splits(
    components: list[tuple[Poly, list[VertexSplit]]],
) -> tuple[Poly, list[VertexSplit]]:
    """Lift component splits to their disjoint-union forest."""
    prefix: list[Poly] = [(1,)]
    for poly, _ in components:
        prefix.append(polymul(prefix[-1], poly))
    suffix: list[Poly] = [(1,)] * (len(components) + 1)
    for j in range(len(components) - 1, -1, -1):
        suffix[j] = polymul(components[j][0], suffix[j + 1])

    out: list[VertexSplit] = []
    offset = 0
    for j, (poly, splits) in enumerate(components):
        other = polymul(prefix[j], suffix[j + 1])
        for split in splits:
            a_poly = polymul(other, split.a_poly)
            b_poly = polymul(other, split.b_poly)
            a_modes = mode_interval(a_poly)
            b_modes = mode_interval(b_poly)
            out.append(VertexSplit(
                vertex=offset + split.vertex,
                a_poly=a_poly,
                b_poly=b_poly,
                a_modes=a_modes,
                b_modes=b_modes,
                distance=interval_distance(a_modes, b_modes),
            ))
        offset += max(1, poly[1] if len(poly) > 1 else 0)
    return prefix[-1], out


def encode_split(split: VertexSplit, include_polys: bool) -> dict[str, object]:
    return {
        "vertex": split.vertex,
        "a_modes": list(split.a_modes),
        "b_modes": list(split.b_modes),
        "distance": split.distance,
        "a_poly": list(split.a_poly) if include_polys else None,
        "b_poly": list(split.b_poly) if include_polys else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=1)
    parser.add_argument("--max-n", type=int, default=17)
    parser.add_argument("--out", type=Path,
                        default=Path("results/existential_split_mode_n17_20260807.json"))
    parser.add_argument("--stop-on-first", action="store_true")
    args = parser.parse_args()

    started = time.time()
    summary: dict[str, object] = {
        "claim_scope": "exhaustive evidence, not a proof",
        "statement": (
            "some vertex v has adjacent modal intervals for "
            "I(F-v) and x I(F-N[v])"
        ),
        "min_n": args.min_n,
        "max_n": args.max_n,
        "orders": [],
        "failures": [],
    }
    total_trees = 0
    total_vertices = 0
    identity_failures = 0
    distance_hist: dict[int, int] = {}

    for n in range(args.min_n, args.max_n + 1):
        order_trees = 0
        order_vertices = 0
        order_failures = 0
        order_best_hist: dict[int, int] = {}
        iterator = [
            (1, [[]], b"@")
        ] if n == 1 else trees_geng_raw(n)
        for tree_index, (_, adj, graph6) in enumerate(iterator):
            order_trees += 1
            order_vertices += n
            full, splits = tree_vertex_splits(adj)
            replay = tuple(independence_poly(n, adj))
            if full != replay:
                identity_failures += 1
                raise AssertionError((n, tree_index, full, replay))
            best = min(split.distance for split in splits)
            order_best_hist[best] = order_best_hist.get(best, 0) + 1
            for split in splits:
                distance_hist[split.distance] = distance_hist.get(split.distance, 0) + 1
            if best > 1:
                order_failures += 1
                failure = {
                    "n": n,
                    "tree_index": tree_index,
                    "graph6": graph6.decode("ascii"),
                    "independence_polynomial": list(full),
                    "best_distance": best,
                    "splits": [encode_split(split, True) for split in splits],
                }
                summary["failures"].append(failure)  # type: ignore[index]
                if args.stop_on_first:
                    break
        total_trees += order_trees
        total_vertices += order_vertices
        row = {
            "n": n,
            "trees": order_trees,
            "vertex_splits": order_vertices,
            "failures": order_failures,
            "best_distance_histogram": {
                str(k): v for k, v in sorted(order_best_hist.items())
            },
            "elapsed_seconds": time.time() - started,
        }
        summary["orders"].append(row)  # type: ignore[index]
        print(json.dumps({"event": "order", **row}), flush=True)
        if args.stop_on_first and order_failures:
            break

    summary.update({
        "total_trees": total_trees,
        "total_vertex_splits": total_vertices,
        "identity_failures": identity_failures,
        "distance_histogram": {str(k): v for k, v in sorted(distance_hist.items())},
        "failure_count": len(summary["failures"]),  # type: ignore[arg-type]
        "elapsed_seconds": time.time() - started,
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "trees": total_trees,
        "failures": len(summary["failures"]),  # type: ignore[arg-type]
        "out": str(args.out),
        "elapsed": time.time() - started,
    }), flush=True)


if __name__ == "__main__":
    main()
