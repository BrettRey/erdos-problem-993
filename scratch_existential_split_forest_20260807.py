#!/usr/bin/env python3
"""Exhaustive forest extension of the adjacent deletion-split kill test.

Every unlabeled forest is a multiset of unlabeled trees.  This driver builds
that multiset library with nauty ``geng``, lifts every exact component split
through the product of the other components, and checks whether some vertex
has adjacent modal intervals for ``I(F-v)`` and ``x I(F-N[v])``.

The output is computational evidence, not a proof.  A failure is a replayable
list of graph6 component certificates plus exact coefficient vectors.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from scratch_existential_split_mode_20260807 import (
    Poly,
    VertexSplit,
    convolve_component_splits,
    encode_split,
    tree_vertex_splits,
)
from trees import trees_geng_raw


@dataclass(frozen=True)
class Component:
    order: int
    graph6: str
    poly: Poly
    splits: tuple[VertexSplit, ...]


def component_library(max_n: int) -> list[Component]:
    out = [Component(1, "@", (1, 1), tuple(tree_vertex_splits([[]])[1]))]
    for n in range(2, max_n + 1):
        for _, adj, graph6 in trees_geng_raw(n):
            poly, splits = tree_vertex_splits(adj)
            out.append(Component(n, graph6.decode("ascii"), poly, tuple(splits)))
    return out


def forest_multisets(
    library: list[Component], total_order: int,
) -> Iterator[tuple[int, ...]]:
    """Yield nondecreasing component-index tuples of the requested order."""
    chosen: list[int] = []

    def visit(remaining: int, start: int) -> Iterator[tuple[int, ...]]:
        for index in range(start, len(library)):
            order = library[index].order
            if order > remaining:
                break
            chosen.append(index)
            if order == remaining:
                yield tuple(chosen)
            else:
                yield from visit(remaining - order, index)
            chosen.pop()

    yield from visit(total_order, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=15)
    parser.add_argument("--out", type=Path,
                        default=Path("results/existential_split_forest_n15_20260807.json"))
    parser.add_argument("--stop-on-first", action="store_true")
    args = parser.parse_args()

    started = time.time()
    library = component_library(args.max_n)
    print(json.dumps({
        "event": "library",
        "components": len(library),
        "elapsed": time.time() - started,
    }), flush=True)

    result: dict[str, object] = {
        "claim_scope": "exhaustive evidence, not a proof",
        "statement": (
            "some vertex v has adjacent modal intervals for "
            "I(F-v) and x I(F-N[v])"
        ),
        "max_n": args.max_n,
        "component_library_size": len(library),
        "orders": [],
        "failures": [],
    }
    total_forests = 0
    total_splits = 0
    best_distance_hist: dict[int, int] = {}

    for n in range(1, args.max_n + 1):
        order_forests = 0
        order_splits = 0
        order_failures = 0
        order_hist: dict[int, int] = {}
        for forest_index, indices in enumerate(forest_multisets(library, n)):
            order_forests += 1
            components = [
                (library[index].poly, list(library[index].splits))
                for index in indices
            ]
            full, splits = convolve_component_splits(components)
            order_splits += len(splits)
            best = min(split.distance for split in splits)
            order_hist[best] = order_hist.get(best, 0) + 1
            best_distance_hist[best] = best_distance_hist.get(best, 0) + 1
            if best > 1:
                order_failures += 1
                failure = {
                    "n": n,
                    "forest_index": forest_index,
                    "component_graph6": [library[index].graph6 for index in indices],
                    "independence_polynomial": list(full),
                    "best_distance": best,
                    "splits": [encode_split(split, True) for split in splits],
                }
                result["failures"].append(failure)  # type: ignore[index]
                if args.stop_on_first:
                    break
        total_forests += order_forests
        total_splits += order_splits
        row = {
            "n": n,
            "forests": order_forests,
            "vertex_splits": order_splits,
            "failures": order_failures,
            "best_distance_histogram": {
                str(k): v for k, v in sorted(order_hist.items())
            },
            "elapsed_seconds": time.time() - started,
        }
        result["orders"].append(row)  # type: ignore[index]
        print(json.dumps({"event": "order", **row}), flush=True)
        if args.stop_on_first and order_failures:
            break

    result.update({
        "total_forests": total_forests,
        "total_vertex_splits": total_splits,
        "best_distance_histogram": {
            str(k): v for k, v in sorted(best_distance_hist.items())
        },
        "failure_count": len(result["failures"]),  # type: ignore[arg-type]
        "elapsed_seconds": time.time() - started,
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "forests": total_forests,
        "failures": len(result["failures"]),  # type: ignore[arg-type]
        "out": str(args.out),
        "elapsed": time.time() - started,
    }), flush=True)


if __name__ == "__main__":
    main()
