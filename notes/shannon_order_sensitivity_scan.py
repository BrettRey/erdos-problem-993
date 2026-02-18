#!/usr/bin/env python3
"""Order-sensitivity scan for Shannon/Lamport boundary behavior.

For each rooted parent state, test multiple child-attachment orders and check
whether any order creates an interior tail violation (k > d).
"""

from __future__ import annotations

import argparse
import datetime as dt
import itertools
import json
import os
import random
import shutil
import sys
import time
from dataclasses import dataclass
from typing import Any

# Allow imports from repository root when running as `python notes/...`.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from trees import trees, trees_geng_raw


def _trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_add(a: list[int], b: list[int]) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += v
    return _trim(out)


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            out[i + j] += ai * bj
    return _trim(out)


def deltas(seq: list[int]) -> list[int]:
    return [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]


def first_descent(seq: list[int]) -> int | None:
    for k in range(len(seq) - 1):
        if seq[k + 1] < seq[k]:
            return k
    return None


def rooted_dp_all(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    n = len(adj)
    parent = [-1] * n
    parent[root] = root
    order = [root]
    for v in order:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v
                order.append(w)

    children: list[list[int]] = [[] for _ in range(n)]
    for v in order:
        if v != root:
            children[parent[v]].append(v)

    f: list[list[int]] = [[] for _ in range(n)]
    g: list[list[int]] = [[] for _ in range(n)]
    for v in reversed(order):
        if not children[v]:
            f[v] = [1]
            g[v] = [0, 1]
        else:
            prod_h = [1]
            prod_f = [1]
            for c in children[v]:
                prod_h = poly_mul(prod_h, poly_add(f[c], g[c]))
                prod_f = poly_mul(prod_f, f[c])
            f[v] = prod_h
            g[v] = [0] + prod_f
    return children, f, g


def violates_interior_tail(order: tuple[int, ...], f: list[list[int]], g: list[list[int]]) -> dict[str, Any] | None:
    P = [1]
    Q = [0, 1]

    for pos, c in enumerate(order):
        U = f[c]
        V = g[c]
        I = poly_add(P, Q)
        IU = poly_mul(I, U)
        PV = poly_mul(P, V)
        Iprime = poly_add(IU, PV)

        d = first_descent(IU)
        if d is not None:
            dIU = deltas(IU)
            dPV = deltas(PV)
            dIp = deltas(Iprime)
            L = max(len(dIU), len(dPV), len(dIp))
            dIU += [0] * (L - len(dIU))
            dPV += [0] * (L - len(dPV))
            dIp += [0] * (L - len(dIp))

            bad_dd = [k for k in range(d + 1, L) if dPV[k] > -dIU[k]]
            if bad_dd:
                return {
                    "type": "interior_dd_failure",
                    "step_pos": pos,
                    "d": d,
                    "k": bad_dd,
                    "delta_IU": dIU,
                    "delta_PV": dPV,
                    "delta_Iprime": dIp,
                }

            bad_tail = [k for k in range(d + 1, L) if dIp[k] > 0]
            if bad_tail:
                return {
                    "type": "positive_tail_after_boundary",
                    "step_pos": pos,
                    "d": d,
                    "k": bad_tail,
                    "delta_IU": dIU,
                    "delta_PV": dPV,
                    "delta_Iprime": dIp,
                }

        P = poly_mul(P, poly_add(U, V))
        Q = poly_mul(Q, U)

    return None


@dataclass
class Stats:
    trees_scanned: int = 0
    rooted_states_tested: int = 0
    orders_tested: int = 0
    exact_permutations_tested: int = 0
    sampled_orders_tested: int = 0
    interior_violations: int = 0
    first_violation: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "trees_scanned": self.trees_scanned,
            "rooted_states_tested": self.rooted_states_tested,
            "orders_tested": self.orders_tested,
            "exact_permutations_tested": self.exact_permutations_tested,
            "sampled_orders_tested": self.sampled_orders_tested,
            "interior_violations": self.interior_violations,
            "first_violation": self.first_violation,
        }


def _tree_iter(n: int, backend: str):
    if n == 1:
        yield 1, [[]], None
        return

    resolved = backend
    if resolved == "auto":
        resolved = "geng" if shutil.which("geng") else "networkx"

    if resolved == "geng":
        for n_out, adj, raw in trees_geng_raw(n):
            yield n_out, adj, raw.decode("ascii")
    else:
        for n_out, adj in trees(n, backend=resolved):
            yield n_out, adj, None


def scan(max_n: int, backend: str, exact_degree: int, sample_orders: int, seed: int) -> dict[str, Any]:
    rng = random.Random(seed)
    stats = Stats()
    by_n: list[dict[str, Any]] = []
    t0 = time.time()

    for n in range(1, max_n + 1):
        n0 = time.time()
        s0 = Stats()

        for tree_idx, (n_out, adj, g6) in enumerate(_tree_iter(n, backend), start=1):
            s0.trees_scanned += 1

            for root in range(n_out):
                children, f, g = rooted_dp_all(adj, root)

                for v in range(n_out):
                    ch = children[v]
                    d = len(ch)
                    if d <= 1:
                        continue

                    s0.rooted_states_tested += 1

                    orders: list[tuple[int, ...]] = []
                    if d <= exact_degree:
                        orders = list(itertools.permutations(ch))
                        s0.exact_permutations_tested += len(orders)
                    else:
                        seen: set[tuple[int, ...]] = set()
                        base = tuple(ch)
                        seen.add(base)
                        orders.append(base)
                        while len(orders) < sample_orders:
                            tmp = ch[:]
                            rng.shuffle(tmp)
                            tup = tuple(tmp)
                            if tup in seen:
                                continue
                            seen.add(tup)
                            orders.append(tup)
                        s0.sampled_orders_tested += len(orders)

                    for ord_ch in orders:
                        s0.orders_tested += 1
                        witness = violates_interior_tail(ord_ch, f, g)
                        if witness is None:
                            continue

                        s0.interior_violations += 1
                        if s0.first_violation is None:
                            s0.first_violation = {
                                "n": n_out,
                                "tree_index": tree_idx,
                                "graph6": g6,
                                "root": root,
                                "parent": v,
                                "order": list(ord_ch),
                                **witness,
                            }

        by_n.append(
            {
                "n": n,
                "elapsed_s": round(time.time() - n0, 3),
                **s0.to_dict(),
            }
        )

        stats.trees_scanned += s0.trees_scanned
        stats.rooted_states_tested += s0.rooted_states_tested
        stats.orders_tested += s0.orders_tested
        stats.exact_permutations_tested += s0.exact_permutations_tested
        stats.sampled_orders_tested += s0.sampled_orders_tested
        stats.interior_violations += s0.interior_violations
        if stats.first_violation is None and s0.first_violation is not None:
            stats.first_violation = s0.first_violation

        print(
            f"n={n} trees={s0.trees_scanned} states={s0.rooted_states_tested} "
            f"orders={s0.orders_tested} interior={s0.interior_violations} "
            f"elapsed={time.time()-n0:.2f}s",
            flush=True,
        )

    return {
        "range": {"min_n": 1, "max_n": max_n},
        "backend": backend,
        "exact_degree": exact_degree,
        "sample_orders": sample_orders,
        "seed": seed,
        "elapsed_total_s": round(time.time() - t0, 3),
        "by_n": by_n,
        "totals": stats.to_dict(),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Order-sensitivity scan for Lamport boundary behavior")
    ap.add_argument("--max-n", type=int, default=12)
    ap.add_argument("--backend", choices=["auto", "geng", "networkx"], default="auto")
    ap.add_argument("--exact-degree", type=int, default=7)
    ap.add_argument("--sample-orders", type=int, default=120)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    report = {
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "description": (
            "Order-sensitivity stress test for Lamport transition. "
            "Checks whether any child-ordering creates an interior (k>d) tail violation."
        ),
        "definitions": {
            "interior_dd_failure": "exists k>d with Delta(PV)_k > -Delta(IU)_k",
            "positive_tail_after_boundary": "exists k>d with Delta(I')_k > 0",
        },
        "scan": scan(
            max_n=args.max_n,
            backend=args.backend,
            exact_degree=args.exact_degree,
            sample_orders=args.sample_orders,
            seed=args.seed,
        ),
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
