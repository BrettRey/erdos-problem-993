#!/usr/bin/env python3
"""Scan for the "squeeze from both ends" condition on independence sequences.

For each tree T, compute the independence sequence i_k. Let alpha be the
independence number (degree of the polynomial), and let

  t = ceil((2*alpha - 1) / 3) = floor((2*alpha + 1) / 3).

Levit–Mandrescu (2006) show the tail is strictly decreasing for k >= t.
If i_k is nondecreasing up to k = t-1 (equivalently first descent >= t),
then unimodality follows immediately. This script measures how often that
condition holds in the exhaustive data.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from typing import Any, Iterable, Iterator

from indpoly import independence_poly
from trees import trees, trees_geng_raw


def first_descent(seq: list[int]) -> int | None:
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            return i
    return None


def peak_index(seq: list[int]) -> int:
    maxv = max(seq)
    return seq.index(maxv)  # earliest peak


def _min_positive_peak_ratio(seq: list[int]) -> float:
    alpha = len(seq) - 1
    if alpha <= 0:
        return 0.0
    maxv = max(seq)
    first = seq.index(maxv)
    last = len(seq) - 1 - seq[::-1].index(maxv)
    # avoid a ratio of 0 when the peak is at k=0 (rare but possible)
    if first == 0 and last > 0:
        return last / alpha
    return first / alpha


def tail_start(alpha: int) -> int:
    # ceil((2*alpha - 1)/3) == floor((2*alpha + 1)/3)
    return (2 * alpha + 1) // 3


def iter_trees(
    n: int, backend: str, res: int | None, mod: int | None
) -> Iterator[tuple[int, list[list[int]], str | None]]:
    if backend == "geng":
        if (res is None) ^ (mod is None):
            raise ValueError("res and mod must be provided together")
        for n_out, adj, g6 in trees_geng_raw(n, res=res, mod=mod):
            yield n_out, adj, g6.decode("ascii")
        return
    # networkx or auto
    for n_out, adj in trees(n, backend=backend):
        yield n_out, adj, None


def scan_range(
    min_n: int,
    max_n: int,
    backend: str,
    res: int | None,
    mod: int | None,
    stop_on_fail: bool,
    max_examples: int,
) -> dict[str, Any]:
    by_n: list[dict[str, Any]] = []
    examples: list[dict[str, Any]] = []
    total = ok = fail = 0
    worst_margin = None
    worst_peak_ratio = None
    t0 = time.time()

    for n in range(min_n, max_n + 1):
        n_total = n_ok = n_fail = 0
        n_worst_margin = None
        n_worst_peak_ratio = None
        n_min_first_descent = None
        n_min_peak_index = None
        start = time.time()

        for _, adj, g6 in iter_trees(n, backend, res, mod):
            poly = independence_poly(n, adj)
            alpha = len(poly) - 1
            t = tail_start(alpha)
            fd = first_descent(poly)
            pk = peak_index(poly)
            peak_ratio = _min_positive_peak_ratio(poly)

            if fd is None:
                margin = None
            else:
                margin = fd - t

            squeeze_ok = (fd is None) or (fd >= t)
            n_total += 1
            if squeeze_ok:
                n_ok += 1
            else:
                n_fail += 1
                if max_examples > 0 and len(examples) < max_examples:
                    examples.append(
                        {
                            "n": n,
                            "alpha": alpha,
                            "tail_start": t,
                            "first_descent": fd,
                            "peak_index": pk,
                            "peak_ratio": peak_ratio,
                            "graph6": g6,
                        }
                    )
                if stop_on_fail:
                    break

            if margin is not None and (n_worst_margin is None or margin < n_worst_margin):
                n_worst_margin = margin
                n_min_first_descent = fd
            if n_worst_peak_ratio is None or peak_ratio < n_worst_peak_ratio:
                n_worst_peak_ratio = peak_ratio
                n_min_peak_index = pk

        elapsed = time.time() - start
        by_n.append(
            {
                "n": n,
                "total_trees": n_total,
                "squeeze_ok": n_ok,
                "squeeze_fail": n_fail,
                "min_first_descent": n_min_first_descent,
                "min_peak_index": n_min_peak_index,
                "min_peak_ratio": n_worst_peak_ratio,
                "min_margin": n_worst_margin,
                "elapsed_s": round(elapsed, 3),
            }
        )

        total += n_total
        ok += n_ok
        fail += n_fail
        if n_worst_margin is not None and (worst_margin is None or n_worst_margin < worst_margin):
            worst_margin = n_worst_margin
        if worst_peak_ratio is None or (
            n_worst_peak_ratio is not None and n_worst_peak_ratio < worst_peak_ratio
        ):
            worst_peak_ratio = n_worst_peak_ratio

        if stop_on_fail and n_fail > 0:
            break

    return {
        "description": "Squeeze-from-both-ends scan: first descent vs Levit–Mandrescu tail start",
        "range": {"min_n": min_n, "max_n": max_n},
        "backend": backend,
        "partition": {"res": res, "mod": mod},
        "by_n": by_n,
        "totals": {
            "total_trees": total,
            "squeeze_ok": ok,
            "squeeze_fail": fail,
            "worst_margin": worst_margin,
            "worst_peak_ratio": worst_peak_ratio,
            "elapsed_total_s": round(time.time() - t0, 3),
        },
        "examples": examples,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-n", type=int, default=1)
    parser.add_argument("--max-n", type=int, required=True)
    parser.add_argument("--backend", type=str, default="auto")
    parser.add_argument("--res", type=int, default=None)
    parser.add_argument("--mod", type=int, default=None)
    parser.add_argument("--out", type=str, default=None)
    parser.add_argument("--stop-on-fail", action="store_true")
    parser.add_argument("--max-examples", type=int, default=3)
    args = parser.parse_args()

    data = scan_range(
        min_n=args.min_n,
        max_n=args.max_n,
        backend=args.backend,
        res=args.res,
        mod=args.mod,
        stop_on_fail=args.stop_on_fail,
        max_examples=args.max_examples,
    )

    if args.out:
        with open(args.out, "w") as f:
            json.dump(data, f, indent=2)
        print(f"Wrote {args.out}")
    else:
        print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
