#!/usr/bin/env python3
"""Probe a Vatter-style residual-tail condition for tree IS polynomials.

For a tree T, write I(T;x) = sum_k i_k(T)x^k and let

    d(T) = min {k >= 1 : i_k(T) < i_{k-1}(T)}.

Fix a leaf-peeling order v_1, ..., v_n.  Just before v_j is removed, let
R_j be the remaining tree and put H_j = R_j - N_{R_j}[v_j].  The condition
tested here is

    i_k(H_j) <= i_{k-1}(H_j)  for every j and every k >= d(T).       (RTD)

The least-vertex decomposition gives

    i_{k+1}(T) = sum_j i_k(H_j),    i_k(T) = sum_j i_{k-1}(H_j).

Thus RTD implies i_{k+1}(T) <= i_k(T) for every k >= d(T), ruling out
every rebound after the first descent.  ``--order-scope all`` checks every
leaf-peeling order exactly by enumerating all reachable remaining subtrees.
The fixed-order modes support larger adversarial trees.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import shutil
import sys
import time
from typing import Callable, Iterable

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from graph6 import parse_graph6  # noqa: E402
from indpoly import independence_poly  # noqa: E402
from targeted import (  # noqa: E402
    make_T_m_t_1,
    make_broom,
    make_double_star,
    make_spider,
)
from trees import trees, trees_geng_raw  # noqa: E402


def first_strict_descent(poly: list[int] | tuple[int, ...]) -> int | None:
    """Return the first k with i_k < i_{k-1}, or None if absent."""
    return next(
        (k for k in range(1, len(poly)) if poly[k] < poly[k - 1]),
        None,
    )


def tail_ascent(
    poly: list[int] | tuple[int, ...],
    start: int,
) -> tuple[int, int, int] | None:
    """Return the first (k, i_{k-1}, i_k) ascent with k >= start."""
    for k in range(start, len(poly)):
        if poly[k] > poly[k - 1]:
            return k, poly[k - 1], poly[k]
    return None


def largest_tail_ratio(
    poly: list[int] | tuple[int, ...],
    start: int,
) -> dict | None:
    """Return the largest defined i_k/i_{k-1} for k >= ``start``."""
    best = None
    for k in range(start, len(poly)):
        denominator = poly[k - 1]
        if denominator == 0:
            continue
        candidate = {
            "k": k,
            "numerator": poly[k],
            "denominator": denominator,
            "ratio": poly[k] / denominator,
        }
        if _larger_ratio(candidate, best):
            best = candidate
    return best


def _larger_ratio(candidate: dict | None, incumbent: dict | None) -> bool:
    if candidate is None:
        return False
    if incumbent is None:
        return True
    return (
        candidate["numerator"] * incumbent["denominator"]
        > incumbent["numerator"] * candidate["denominator"]
    )


def induced_polynomial(adj: list[list[int]], vertices: Iterable[int]) -> list[int]:
    """Compute the IS polynomial of the subgraph induced by ``vertices``."""
    keep = list(vertices)
    position = {v: i for i, v in enumerate(keep)}
    sub_adj = [
        [position[u] for u in adj[v] if u in position]
        for v in keep
    ]
    return independence_poly(len(keep), sub_adj)


def leaf_peeling_order(
    adj: list[list[int]],
    choose: Callable[[list[int]], int],
) -> list[int]:
    """Return a leaf-peeling order using ``choose`` to break leaf ties."""
    alive = set(range(len(adj)))
    degree = [len(neighbors) for neighbors in adj]
    order: list[int] = []
    while alive:
        leaves = [v for v in alive if degree[v] <= 1]
        if not leaves:
            raise ValueError("adjacency list is not a forest")
        v = choose(leaves)
        if v not in leaves:
            raise ValueError("leaf chooser returned a non-leaf")
        order.append(v)
        alive.remove(v)
        for u in adj[v]:
            if u in alive:
                degree[u] -= 1
    return order


def residual_polynomials(
    adj: list[list[int]],
    order: list[int],
) -> list[tuple[int, list[int], list[int]]]:
    """Return ``(v, residual vertices, I(residual))`` along an order."""
    alive = set(range(len(adj)))
    out: list[tuple[int, list[int], list[int]]] = []
    for v in order:
        if v not in alive:
            raise ValueError("order repeats a vertex")
        residual = sorted(alive - {v} - set(adj[v]))
        out.append((v, residual, induced_polynomial(adj, residual)))
        alive.remove(v)
    if alive:
        raise ValueError("order omits a vertex")
    return out


def check_fixed_order(
    adj: list[list[int]],
    order: list[int],
) -> dict | None:
    """Return an RTD counterexample for one leaf order, if one exists."""
    parent = independence_poly(len(adj), adj)
    descent = first_strict_descent(parent)
    if descent is None:
        return None

    for step, (v, residual, poly) in enumerate(residual_polynomials(adj, order)):
        ascent = tail_ascent(poly, descent)
        if ascent is not None:
            k, previous, current = ascent
            return {
                "descent": descent,
                "k": k,
                "step": step,
                "vertex": v,
                "remaining_residual": residual,
                "parent_poly": parent,
                "residual_poly": poly,
                "previous": previous,
                "current": current,
                "order": order,
            }
    return None


class InducedPolynomialOracle:
    """Exact memoized IS polynomials for induced subgraphs of one tree."""

    def __init__(self, adj: list[list[int]]) -> None:
        self.adj_masks = [sum(1 << u for u in neighbors) for neighbors in adj]
        self.cache: dict[int, tuple[int, ...]] = {0: (1,)}

    def polynomial(self, mask: int) -> tuple[int, ...]:
        cached = self.cache.get(mask)
        if cached is not None:
            return cached
        vertex_bit = mask & -mask
        v = vertex_bit.bit_length() - 1
        without_v = self.polynomial(mask ^ vertex_bit)
        without_closed_neighborhood = self.polynomial(
            mask & ~vertex_bit & ~self.adj_masks[v]
        )
        out = list(without_v)
        needed = len(without_closed_neighborhood) + 1
        if len(out) < needed:
            out.extend([0] * (needed - len(out)))
        for k, value in enumerate(without_closed_neighborhood):
            out[k + 1] += value
        result = tuple(out)
        self.cache[mask] = result
        return result


def check_all_leaf_orders(adj: list[list[int]]) -> tuple[dict | None, dict]:
    """Check RTD for every leaf-peeling order exactly."""
    n = len(adj)
    oracle = InducedPolynomialOracle(adj)
    full_mask = (1 << n) - 1
    parent = oracle.polynomial(full_mask)
    descent = first_strict_descent(parent)
    if descent is None:
        return None, {
            "remaining_states": 1,
            "residual_checks": 0,
            "polynomials": len(oracle.cache),
            "worst_ratio": None,
        }

    seen = {full_mask}
    stack = [full_mask]
    residual_checks = 0
    worst_ratio = None
    while stack:
        mask = stack.pop()
        remaining_bits = mask
        while remaining_bits:
            vertex_bit = remaining_bits & -remaining_bits
            remaining_bits ^= vertex_bit
            v = vertex_bit.bit_length() - 1
            if (oracle.adj_masks[v] & mask).bit_count() > 1:
                continue

            residual_mask = mask & ~vertex_bit & ~oracle.adj_masks[v]
            residual_poly = oracle.polynomial(residual_mask)
            residual_checks += 1
            candidate_ratio = largest_tail_ratio(residual_poly, descent)
            if _larger_ratio(candidate_ratio, worst_ratio):
                worst_ratio = {
                    **candidate_ratio,
                    "remaining": _mask_vertices(mask),
                    "vertex": v,
                    "remaining_residual": _mask_vertices(residual_mask),
                    "residual_poly": list(residual_poly),
                }
            ascent = tail_ascent(residual_poly, descent)
            if ascent is not None:
                k, previous, current = ascent
                failure = {
                    "descent": descent,
                    "k": k,
                    "remaining": _mask_vertices(mask),
                    "vertex": v,
                    "remaining_residual": _mask_vertices(residual_mask),
                    "parent_poly": list(parent),
                    "residual_poly": list(residual_poly),
                    "previous": previous,
                    "current": current,
                }
                return failure, {
                    "remaining_states": len(seen),
                    "residual_checks": residual_checks,
                    "polynomials": len(oracle.cache),
                    "worst_ratio": worst_ratio,
                }

            next_mask = mask ^ vertex_bit
            if next_mask and next_mask not in seen:
                seen.add(next_mask)
                stack.append(next_mask)

    return None, {
        "remaining_states": len(seen),
        "residual_checks": residual_checks,
        "polynomials": len(oracle.cache),
        "worst_ratio": worst_ratio,
    }


def _mask_vertices(mask: int) -> list[int]:
    vertices: list[int] = []
    while mask:
        bit = mask & -mask
        mask ^= bit
        vertices.append(bit.bit_length() - 1)
    return vertices


def _tree_iterator(
    n: int,
    backend: str,
) -> Iterable[tuple[list[list[int]], str | None]]:
    resolved = backend
    if resolved == "auto":
        resolved = "geng" if shutil.which("geng") else "networkx"
    if resolved == "geng":
        for _, adj, raw in trees_geng_raw(n):
            yield adj, raw.decode("ascii")
        return
    for _, adj in trees(n, backend=resolved):
        yield adj, None


def run_exhaustive(args: argparse.Namespace) -> dict:
    totals = {
        "trees": 0,
        "remaining_states": 0,
        "residual_checks": 0,
        "polynomials": 0,
        "worst_ratio": None,
    }
    per_n = []
    failure = None
    for n in range(args.min_n, args.max_n + 1):
        started = time.monotonic()
        row = {
            "n": n,
            "trees": 0,
            "remaining_states": 0,
            "residual_checks": 0,
            "polynomials": 0,
            "worst_ratio": None,
        }
        for tree_index, (adj, graph6) in enumerate(
            _tree_iterator(n, args.backend),
            start=1,
        ):
            if args.order_scope == "all":
                bad, stats = check_all_leaf_orders(adj)
            else:
                bad, stats = check_selected_orders(
                    adj,
                    args.order_scope,
                    args.random_orders,
                    args.seed + totals["trees"],
                )
            row["trees"] += 1
            totals["trees"] += 1
            for key in ("remaining_states", "residual_checks", "polynomials"):
                row[key] += stats[key]
                totals[key] += stats[key]
            if _larger_ratio(stats.get("worst_ratio"), totals["worst_ratio"]):
                totals["worst_ratio"] = {
                    **stats["worst_ratio"],
                    "n": n,
                    "tree_index": tree_index,
                    "graph6": graph6,
                    "adj": adj if graph6 is None else None,
                }
            if _larger_ratio(stats.get("worst_ratio"), row["worst_ratio"]):
                row["worst_ratio"] = {
                    **stats["worst_ratio"],
                    "tree_index": tree_index,
                    "graph6": graph6,
                    "adj": adj if graph6 is None else None,
                }
            if bad is not None:
                failure = {
                    "n": n,
                    "tree_index": tree_index,
                    "graph6": graph6,
                    "adj": adj if graph6 is None else None,
                    "violation": bad,
                }
                break
        row["elapsed_seconds"] = time.monotonic() - started
        per_n.append(row)
        print(
            f"n={n}: trees={row['trees']:,}, "
            f"states={row['remaining_states']:,}, "
            f"checks={row['residual_checks']:,}, "
            f"elapsed={row['elapsed_seconds']:.2f}s",
            flush=True,
        )
        if failure is not None:
            break
    return {"totals": totals, "per_n": per_n, "failure": failure}


def check_selected_orders(
    adj: list[list[int]],
    scope: str,
    random_orders: int,
    seed: int,
) -> tuple[dict | None, dict]:
    """Check min/max and/or randomized leaf orders on one tree."""
    rng = random.Random(seed)
    orders: list[tuple[str, list[int]]] = []
    if scope in {"minmax", "selected"}:
        orders.append(("min", leaf_peeling_order(adj, min)))
        orders.append(("max", leaf_peeling_order(adj, max)))
    if scope in {"random", "selected"}:
        for index in range(random_orders):
            orders.append(
                (
                    f"random-{index}",
                    leaf_peeling_order(adj, rng.choice),
                )
            )
    for label, order in orders:
        failure = check_fixed_order(adj, order)
        if failure is not None:
            failure["order_label"] = label
            return failure, {
                "remaining_states": 0,
                "residual_checks": len(order),
                "polynomials": 0,
                "worst_ratio": None,
            }
    return None, {
        "remaining_states": 0,
        "residual_checks": len(adj) * len(orders),
        "polynomials": 0,
        "worst_ratio": None,
    }


def run_witnesses(args: argparse.Namespace) -> dict:
    source = _ROOT / "results" / "analysis_n26.json"
    with source.open() as handle:
        witnesses = json.load(handle)["lc_failures"]
    rows = []
    failure = None
    for index, witness in enumerate(witnesses, start=1):
        _, adj = parse_graph6(witness["graph6"].encode("ascii"))
        bad, stats = check_selected_orders(
            adj,
            "selected",
            args.witness_orders,
            args.seed + index,
        )
        row = {
            "index": index,
            "graph6": witness["graph6"],
            "orders": args.witness_orders + 2,
            "stats": stats,
            "failure": bad,
        }
        rows.append(row)
        print(
            f"n=26 LC witness {index}: orders={row['orders']:,}, "
            f"status={'FAIL' if bad else 'pass'}",
            flush=True,
        )
        if bad is not None:
            failure = row
            break
    return {"source": str(source), "rows": rows, "failure": failure}


def targeted_cases(max_n: int) -> Iterable[tuple[str, list[list[int]]]]:
    """Yield the bounded adversarial family suite used by the gate."""
    for path_len in range(2, 61, 3):
        for star_size in range(1, 61, 3):
            n, adj = make_broom(path_len, star_size)
            if n <= max_n:
                yield f"broom({path_len},{star_size})", adj
    for left in range(1, 61, 3):
        for right in range(left, 61, 5):
            n, adj = make_double_star(left, right)
            if n <= max_n:
                yield f"double-star({left},{right})", adj
    for legs in range(2, 61, 2):
        n, adj = make_spider([2] * legs)
        if n <= max_n:
            yield f"spider(2^{legs})", adj
        n, adj = make_spider([1] * (legs // 2) + [2] * legs)
        if n <= max_n:
            yield f"spider(1^{legs // 2},2^{legs})", adj
    for branches in range(2, 15):
        for width in range(1, 12):
            n, adj = make_T_m_t_1(branches, width)
            if n <= max_n:
                yield f"T({branches},{width},1)", adj


def run_targeted(args: argparse.Namespace) -> dict:
    rows = []
    failure = None
    started = time.monotonic()
    for index, (label, adj) in enumerate(targeted_cases(args.targeted_max_n), start=1):
        bad, stats = check_selected_orders(
            adj,
            "selected",
            args.targeted_random_orders,
            args.seed + index,
        )
        if bad is not None:
            failure = {
                "index": index,
                "label": label,
                "n": len(adj),
                "adj": adj,
                "violation": bad,
            }
            break
        rows.append({"label": label, "n": len(adj), "stats": stats})
        if index % 200 == 0:
            print(f"targeted: checked {index:,} cases", flush=True)
    return {
        "cases": len(rows) + (1 if failure else 0),
        "orders_per_case": args.targeted_random_orders + 2,
        "max_n": args.targeted_max_n,
        "elapsed_seconds": time.monotonic() - started,
        "failure": failure,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite",
        choices=("exhaustive", "witnesses", "targeted", "gate"),
        default="gate",
    )
    parser.add_argument("--min-n", type=int, default=1)
    parser.add_argument("--max-n", type=int, default=14)
    parser.add_argument(
        "--backend",
        choices=("auto", "geng", "networkx"),
        default="auto",
    )
    parser.add_argument(
        "--order-scope",
        choices=("all", "minmax", "random", "selected"),
        default="all",
    )
    parser.add_argument("--random-orders", type=int, default=10)
    parser.add_argument("--witness-orders", type=int, default=1000)
    parser.add_argument("--targeted-random-orders", type=int, default=1)
    parser.add_argument("--targeted-max-n", type=int, default=180)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--out", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.monotonic()
    result = {
        "candidate": (
            "For every tree T, every leaf-peeling order, every residual "
            "H_j=R_j-N_Rj[v_j], and every k>=d(T), "
            "i_k(H_j)<=i_{k-1}(H_j)."
        ),
        "argv": sys.argv,
        "seed": args.seed,
        "suite": args.suite,
    }
    if args.suite in {"exhaustive", "gate"}:
        result["exhaustive"] = run_exhaustive(args)
    if not _has_failure(result) and args.suite in {"witnesses", "gate"}:
        result["witnesses"] = run_witnesses(args)
    if not _has_failure(result) and args.suite in {"targeted", "gate"}:
        result["targeted"] = run_targeted(args)
    result["elapsed_seconds"] = time.monotonic() - started
    result["status"] = "fail" if _has_failure(result) else "pass"

    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        print(f"wrote {args.out}", flush=True)
    print(
        f"RTD gate: {result['status'].upper()} "
        f"({result['elapsed_seconds']:.2f}s)",
        flush=True,
    )
    return 1 if result["status"] == "fail" else 0


def _has_failure(result: dict) -> bool:
    return any(
        isinstance(value, dict) and value.get("failure") is not None
        for value in result.values()
    )


if __name__ == "__main__":
    raise SystemExit(main())
