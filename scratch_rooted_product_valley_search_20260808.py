#!/usr/bin/env python3
"""Search rooted tree products for a nonunimodal independence polynomial.

If ``G`` is a tree on ``g`` vertices and ``(H, r)`` is a rooted tree, the
rooted product identifies the root of one copy of ``H`` with every vertex of
``G``.  Put

    A = I(H-r; x),
    B = x I(H-N[r]; x).

The resulting tree has independence polynomial

    A(x)^g I(G; B(x) / A(x))
      = sum_k i_k(G) B(x)^k A(x)^(g-k).

This script tests that identity exactly against the known order-26 and
order-28 log-concavity failures.  A rooted gadget can therefore amplify a
known ratio reversal into a genuine coefficient valley without ever leaving
the class of trees.  Every reported counterexample includes enough graph
data for an independent reconstruction and direct DP replay.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly, is_unimodal
from trees import trees_geng_raw


def rooted_states(adj: list[list[int]], root: int) -> tuple[list[int], list[int]]:
    """Return the exact excluded and included-root polynomials ``(A, B)``."""
    n = len(adj)
    keep_a = [vertex != root for vertex in range(n)]
    keep_b = [vertex != root and vertex not in adj[root] for vertex in range(n)]

    def induced_poly(keep: list[bool]) -> list[int]:
        vertices = [vertex for vertex, retained in enumerate(keep) if retained]
        if not vertices:
            return [1]
        relabel = {old: new for new, old in enumerate(vertices)}
        sub_adj = [[] for _ in vertices]
        for old in vertices:
            sub_adj[relabel[old]] = [
                relabel[neighbor]
                for neighbor in adj[old]
                if keep[neighbor]
            ]
        return independence_poly(len(vertices), sub_adj)

    a_poly = induced_poly(keep_a)
    b_poly = [0] + induced_poly(keep_b)
    return a_poly, b_poly


def homogeneous_substitution(
    base_poly: list[int], base_order: int, a_poly: list[int], b_poly: list[int],
) -> list[int]:
    """Compute ``A^n I(G; B/A)`` by homogeneous Horner evaluation."""
    if len(base_poly) > base_order + 1:
        raise ValueError("base polynomial degree exceeds the base order")
    coefficients = base_poly + [0] * (base_order + 1 - len(base_poly))
    value = [coefficients[base_order]]
    a_power = [1]
    for index in range(base_order - 1, -1, -1):
        a_power = _polymul(a_power, a_poly)
        value = _polyadd(
            _polymul(value, b_poly),
            [coefficients[index] * coefficient for coefficient in a_power],
        )
    while len(value) > 1 and value[-1] == 0:
        value.pop()
    return value


def valley_witness(poly: list[int]) -> dict[str, int] | None:
    """Return the first adjacent down-then-up witness, if present."""
    first_descent = None
    for index in range(len(poly) - 1):
        if poly[index] > poly[index + 1]:
            first_descent = index
            break
    if first_descent is None:
        return None
    for index in range(first_descent + 1, len(poly) - 1):
        if poly[index] < poly[index + 1]:
            return {
                "first_descent_k": first_descent,
                "later_ascent_k": index,
                "descent_left": poly[first_descent],
                "descent_right": poly[first_descent + 1],
                "ascent_left": poly[index],
                "ascent_right": poly[index + 1],
            }
    return None


def near_miss(poly: list[int]) -> tuple[int, int, int]:
    """Return the largest genuine slope rebound after an initial descent.

    Write ``r_k = p[k+1]/p[k]``.  After the first edge with ``r_k < 1``, a
    candidate edge is counted only when its ratio is larger than an earlier
    tail ratio.  This excludes the slowly flattening first edge of a broad
    log-concave peak.  A nonunimodal polynomial has a counted edge with ratio
    greater than one; a useful near miss has a genuine ratio reversal whose
    later ratio is still just below one.
    """
    first_descent = None
    for index in range(len(poly) - 1):
        if poly[index] > poly[index + 1]:
            first_descent = index
            break
    first_later_edge = first_descent + 1 if first_descent is not None else None
    if first_later_edge is None or first_later_edge + 1 >= len(poly):
        return -1, 0, 1
    minimum_num = poly[first_descent + 1]
    minimum_den = poly[first_descent]
    best_k, best_num, best_den = -1, 0, 1
    for index in range(first_later_edge, len(poly) - 1):
        numerator = poly[index + 1]
        denominator = poly[index]
        if numerator * minimum_den > minimum_num * denominator:
            if numerator * best_den > best_num * denominator:
                best_k, best_num, best_den = index, numerator, denominator
        if numerator * minimum_den < minimum_num * denominator:
            minimum_num, minimum_den = numerator, denominator
    return best_k, best_num, best_den


def build_rooted_product(
    base_adj: list[list[int]], gadget_adj: list[list[int]], root: int,
) -> list[list[int]]:
    """Construct the rooted-product adjacency list for direct verification."""
    g = len(base_adj)
    h = len(gadget_adj)
    nonroots = [vertex for vertex in range(h) if vertex != root]
    out: list[list[int]] = [[] for _ in range(g * h)]

    # Put base vertices at 0, h, 2h, ... and give every copy its own block.
    def mapped(copy: int, gadget_vertex: int) -> int:
        if gadget_vertex == root:
            return copy * h
        position = nonroots.index(gadget_vertex)
        return copy * h + 1 + position

    for u, neighbors in enumerate(base_adj):
        for v in neighbors:
            if u < v:
                left, right = u * h, v * h
                out[left].append(right)
                out[right].append(left)
    for copy in range(g):
        for u, neighbors in enumerate(gadget_adj):
            for v in neighbors:
                if u < v:
                    left, right = mapped(copy, u), mapped(copy, v)
                    out[left].append(right)
                    out[right].append(left)
    return out


def load_bases() -> list[dict[str, object]]:
    """Load every known exhaustive log-concavity-failing tree certificate."""
    sources = [
        Path("results/analysis_n26.json"),
        Path("results/analysis_n28_modal_lc_nm.json"),
    ]
    bases: list[dict[str, object]] = []
    for source in sources:
        data = json.loads(source.read_text())
        records = data.get("lc_failures") or data.get("top_lc_failures") or []
        for record_index, record in enumerate(records):
            graph6 = record["graph6"]
            n, adj = parse_graph6(graph6.encode("ascii"))
            replay = independence_poly(n, adj)
            stored = record.get("poly")
            if stored is not None and replay != stored:
                raise AssertionError(f"base replay mismatch in {source}")
            bases.append({
                "source": str(source),
                "record_index": record_index,
                "graph6": graph6,
                "n": n,
                "adj": adj,
                "poly": replay,
            })
    return bases


def encode_candidate(
    base: dict[str, object], gadget_graph6: str, gadget_adj: list[list[int]],
    root: int, poly: list[int], include_poly: bool,
) -> dict[str, object]:
    k, numerator, denominator = near_miss(poly)
    return {
        "base_source": base["source"],
        "base_record_index": base["record_index"],
        "base_graph6": base["graph6"],
        "base_order": base["n"],
        "gadget_graph6": gadget_graph6,
        "gadget_order": len(gadget_adj),
        "gadget_root": root,
        "gadget_edges": [
            [u, v]
            for u, neighbors in enumerate(gadget_adj)
            for v in neighbors if u < v
        ],
        "product_order": int(base["n"]) * len(gadget_adj),
        "degree": len(poly) - 1,
        "near_miss": {
            "k": k,
            "numerator": numerator,
            "denominator": denominator,
            "ratio": numerator / denominator,
        },
        "valley": valley_witness(poly),
        "polynomial": poly if include_poly else None,
    }


def better(candidate: dict[str, object], incumbent: dict[str, object] | None) -> bool:
    if incumbent is None:
        return True
    c = candidate["near_miss"]
    i = incumbent["near_miss"]
    assert isinstance(c, dict) and isinstance(i, dict)
    return int(c["numerator"]) * int(i["denominator"]) > int(i["numerator"]) * int(c["denominator"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-gadget-order", type=int, default=10)
    parser.add_argument("--out", type=Path,
                        default=Path("results/rooted_product_valley_20260808.json"))
    parser.add_argument("--progress-every", type=int, default=100)
    args = parser.parse_args()

    started = time.time()
    bases = load_bases()
    evaluated = 0
    best: dict[str, object] | None = None
    counterexample: dict[str, object] | None = None

    for gadget_order in range(1, args.max_gadget_order + 1):
        iterator = [(1, [[]], b"@")]
        if gadget_order > 1:
            iterator = trees_geng_raw(gadget_order)
        rooted_gadgets = 0
        for _, gadget_adj, graph6_bytes in iterator:
            gadget_graph6 = graph6_bytes.decode("ascii")
            for root in range(gadget_order):
                rooted_gadgets += 1
                a_poly, b_poly = rooted_states(gadget_adj, root)
                for base in bases:
                    poly = homogeneous_substitution(
                        base["poly"], int(base["n"]), a_poly, b_poly,
                    )
                    evaluated += 1
                    candidate = encode_candidate(
                        base, gadget_graph6, gadget_adj, root, poly, False,
                    )
                    if better(candidate, best):
                        best = candidate
                    if not is_unimodal(poly):
                        product_adj = build_rooted_product(
                            base["adj"], gadget_adj, root,
                        )
                        replay = independence_poly(len(product_adj), product_adj)
                        if replay != poly:
                            raise AssertionError("rooted-product identity replay failed")
                        counterexample = encode_candidate(
                            base, gadget_graph6, gadget_adj, root, poly, True,
                        )
                        counterexample["product_edges"] = [
                            [u, v]
                            for u, neighbors in enumerate(product_adj)
                            for v in neighbors if u < v
                        ]
                        break
                if counterexample is not None:
                    break
            if counterexample is not None:
                break
        print(json.dumps({
            "event": "order",
            "gadget_order": gadget_order,
            "rooted_gadgets": rooted_gadgets,
            "evaluated": evaluated,
            "best_ratio": (
                best["near_miss"]["ratio"] if best is not None else None
            ),
            "counterexample_found": counterexample is not None,
            "elapsed_seconds": time.time() - started,
        }), flush=True)
        if counterexample is not None:
            break

    result = {
        "claim_scope": "exact rooted-product search",
        "identity": "I(G[H];x)=A(x)^|V(G)| I(G;B(x)/A(x))",
        "base_count": len(bases),
        "max_gadget_order_requested": args.max_gadget_order,
        "evaluated_products": evaluated,
        "counterexample_found": counterexample is not None,
        "counterexample": counterexample,
        "best_near_miss": best,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "evaluated": evaluated,
        "counterexample_found": counterexample is not None,
        "out": str(args.out),
        "elapsed_seconds": time.time() - started,
    }), flush=True)
    raise SystemExit(1 if counterexample is not None else 0)


if __name__ == "__main__":
    main()
