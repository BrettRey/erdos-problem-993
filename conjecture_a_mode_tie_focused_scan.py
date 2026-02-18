#!/usr/bin/env python3
"""Focused mode-tie scan for independence polynomials of trees.

For each eligible tree T with independence polynomial

    I_T(x) = sum_k i_k x^k,

let m be the leftmost mode index at lambda=1, define

    lambda_m = i_{m-1} / i_m,

and check the focused tie inequality

    mu(lambda_m) >= m - 1,

where mu(lambda) is the hard-core mean IS size at fugacity lambda.

This script can run on all trees or restricted to d_leaf <= 1 trees.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import time
from collections import Counter
from dataclasses import asdict, dataclass
from heapq import heappush, heapreplace
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import independence_poly


@dataclass
class Witness:
    n: int
    g6: str
    mode: int
    lambda_mode: float
    mu_lambda_mode: float
    margin: float
    degree_signature: dict[str, int]
    is_balanced_len2_spider: bool
    spider_k: int


@dataclass
class PerN:
    seen: int
    considered: int
    checked: int
    fail: int
    min_margin: float | None
    min_witness: dict[str, Any] | None
    wall_s: float


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def mean_at_lambda(poly: list[int], lam: float) -> float:
    z = 0.0
    m = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        m += k * w
        p *= lam
    return m / z if z > 0.0 else float("nan")


def balanced_len2_spider_k(adj: list[list[int]]) -> int:
    n = len(adj)
    deg = [len(nb) for nb in adj]
    k = max(deg)
    if n != 2 * k + 1:
        return 0
    if deg.count(k) != 1 or deg.count(1) != k or deg.count(2) != k:
        return 0

    center = deg.index(k)
    if any(deg[u] != 2 for u in adj[center]):
        return 0

    for u in adj[center]:
        other = [v for v in adj[u] if v != center]
        if len(other) != 1:
            return 0
        leaf = other[0]
        if deg[leaf] != 1:
            return 0

    return k


def make_witness(n: int, g6: str, adj: list[list[int]], m: int, lam: float, mu: float) -> Witness:
    k = balanced_len2_spider_k(adj)
    degree_counts = Counter(len(nb) for nb in adj)
    return Witness(
        n=n,
        g6=g6,
        mode=m,
        lambda_mode=lam,
        mu_lambda_mode=mu,
        margin=mu - (m - 1),
        degree_signature={str(d): c for d, c in sorted(degree_counts.items())},
        is_balanced_len2_spider=(k > 0),
        spider_k=k,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Focused mode-tie scan: mu(lambda_mode) >= mode-1"
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument(
        "--all-trees",
        action="store_true",
        help="Check all trees (default checks only d_leaf<=1).",
    )
    parser.add_argument("--stop-on-first", action="store_true")
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0
    total_checked = 0
    total_fail = 0

    min_margin = math.inf
    min_witness: Witness | None = None

    # Keep top-k smallest margins as a max-heap by (-margin, idx, Witness).
    heap: list[tuple[float, int, Witness]] = []
    serial = 0

    per_n: dict[str, dict[str, Any]] = {}
    t_all = time.time()

    scope = "all trees" if args.all_trees else "d_leaf<=1 trees"
    print(f"Focused mode-tie scan on {scope}", flush=True)
    print("Property: mu(lambda_mode) >= mode - 1", flush=True)
    print("-" * 88, flush=True)

    should_stop = False
    first_fail: Witness | None = None

    for n in range(args.min_n, args.max_n + 1):
        if should_stop:
            break

        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_checked = 0
        n_fail = 0
        n_min_margin = math.inf
        n_min_witness: Witness | None = None

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1

            if not args.all_trees and not is_dleaf_le_1(n0, adj):
                continue

            n_considered += 1
            total_considered += 1

            poly = independence_poly(n0, adj)
            m = mode_index_leftmost(poly)
            if m == 0 or poly[m - 1] <= 0 or poly[m] <= 0:
                continue

            lam = poly[m - 1] / poly[m]
            mu = mean_at_lambda(poly, lam)
            margin = mu - (m - 1)

            n_checked += 1
            total_checked += 1

            g6 = line.decode("ascii").strip()
            wit = make_witness(n0, g6, adj, m, lam, mu)

            if margin < min_margin:
                min_margin = margin
                min_witness = wit

            if margin < n_min_margin:
                n_min_margin = margin
                n_min_witness = wit

            serial += 1
            key = -margin
            item = (key, serial, wit)
            if len(heap) < args.top_k:
                heappush(heap, item)
            elif key > heap[0][0]:
                heapreplace(heap, item)

            if margin < -args.tol:
                n_fail += 1
                total_fail += 1
                if first_fail is None:
                    first_fail = wit
                if args.stop_on_first:
                    should_stop = True
                    proc.kill()
                    break

        proc.wait()

        per_n[str(n)] = asdict(
            PerN(
                seen=n_seen,
                considered=n_considered,
                checked=n_checked,
                fail=n_fail,
                min_margin=None if n_min_witness is None else n_min_margin,
                min_witness=None if n_min_witness is None else asdict(n_min_witness),
                wall_s=time.time() - t0,
            )
        )

        if n_min_witness is None:
            msg = (
                f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
                f"checked={n_checked:8d} fail={n_fail:4d} (no eligible) "
                f"({time.time()-t0:.1f}s)"
            )
        else:
            msg = (
                f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
                f"checked={n_checked:8d} fail={n_fail:4d} "
                f"min_margin={n_min_margin:.12f} "
                f"spider={n_min_witness.is_balanced_len2_spider} "
                f"({time.time()-t0:.1f}s)"
            )
        print(msg, flush=True)

    top_worst = sorted((item[2] for item in heap), key=lambda w: w.margin)

    print("-" * 88, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"checked={total_checked:,} fail={total_fail:,} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(f"Global minimum margin: {min_margin}", flush=True)
    print(
        f"Global min witness: {None if min_witness is None else asdict(min_witness)}",
        flush=True,
    )
    if first_fail is not None:
        print(f"First failure witness: {asdict(first_fail)}", flush=True)

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "tol": args.tol,
            "all_trees": args.all_trees,
            "top_k": args.top_k,
            "stop_on_first": args.stop_on_first,
        },
        "summary": {
            "seen": total_seen,
            "considered": total_considered,
            "checked": total_checked,
            "fail": total_fail,
            "minimum_margin": min_margin,
            "minimum_witness": None if min_witness is None else asdict(min_witness),
            "first_failure": None if first_fail is None else asdict(first_fail),
            "top_worst": [asdict(w) for w in top_worst],
            "wall_s": time.time() - t_all,
        },
        "per_n": per_n,
    }

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
