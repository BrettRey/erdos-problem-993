#!/usr/bin/env python3
"""Exact tail scan of the Bautista--Ramos pattern-tree family.

For

    U_m = (T_{1,l}^v : S_{2,n}^w)_k^(m),

Theorem 4.11 of arXiv:2603.14204 gives

    I(U_m) = I(S_{2,l}) I(T_{k,n})^m
             + x I(P_2)^l I(S_{2,n})^(k m),

where

    I(P_2) = 1 + 2x,
    I(S_{2,t}) = (1 + 2x)^t + x(1 + x)^t,
    I(T_{k,n}) = I(S_{2,n})^k + x I(P_2)^(kn).

The two summands have the same degree, so reversal commutes with their sum.
This script computes a fixed number of *top* coefficients exactly, even for
very large m, and looks for a valley (which certifies non-unimodality).
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def add(a: list[int], b: list[int], limit: int | None = None) -> list[int]:
    size = max(len(a), len(b))
    if limit is not None:
        size = min(size, limit)
    out = [0] * size
    for i, value in enumerate(a[:size]):
        out[i] += value
    for i, value in enumerate(b[:size]):
        out[i] += value
    return out


def mul(a: list[int], b: list[int], limit: int | None = None) -> list[int]:
    if not a or not b:
        return []
    size = len(a) + len(b) - 1
    if limit is not None:
        size = min(size, limit)
    out = [0] * size
    for i, ai in enumerate(a):
        if ai == 0 or i >= size:
            continue
        stop = min(len(b), size - i)
        for j in range(stop):
            out[i + j] += ai * b[j]
    return out


def power(base: list[int], exponent: int, limit: int | None = None) -> list[int]:
    result = [1]
    factor = base[:limit]
    while exponent:
        if exponent & 1:
            result = mul(result, factor, limit)
        exponent >>= 1
        if exponent:
            factor = mul(factor, factor, limit)
    return result


def shift(poly: list[int]) -> list[int]:
    return [0] + poly


def s2(t: int) -> list[int]:
    return add(power([1, 2], t), shift(power([1, 1], t)))


def tree_t(k: int, n: int) -> list[int]:
    return add(power(s2(n), k), shift(power([1, 2], k * n)))


def reverse(poly: list[int]) -> list[int]:
    return poly[::-1]


def reversed_s2(t: int, limit: int) -> list[int]:
    """Top coefficients of S_{2,t}, in reverse order.

    Reversing at the common degree t+1 gives

        rev(S_{2,t}) = (1+x)^t + x(2+x)^t.

    Computing this expression directly avoids constructing all coefficients
    when t is much larger than the requested tail length.
    """
    return add(
        power([1, 1], t, limit),
        shift(power([2, 1], t, max(0, limit - 1))),
        limit,
    )


def reversed_tree_t(k: int, n: int, limit: int) -> list[int]:
    """Top coefficients of T_{k,n}, in reverse order.

    At degree k(n+1), the exact reversal identity is

        rev(T_{k,n}) = rev(S_{2,n})^k + x^(k-1)(2+x)^(kn).
    """
    first = power(reversed_s2(n, limit), k, limit)
    if k - 1 >= limit:
        return first
    correction = [0] * (k - 1) + power(
        [2, 1], k * n, limit - (k - 1)
    )
    return add(first, correction, limit)


def reversed_u_tail(k: int, ell: int, n: int, m: int, limit: int) -> list[int]:
    # Reversal is multiplicative and the summands have common degree.
    first = mul(
        reversed_s2(ell, limit),
        power(reversed_tree_t(k, n, limit), m, limit),
        limit,
    )
    second = mul(
        power([2, 1], ell, limit),
        power(reversed_s2(n, limit), k * m, limit),
        limit,
    )
    return add(first, second, limit)


def valleys(seq: list[int]) -> list[int]:
    """Return strict-valley indices; a plateau valley is checked separately."""
    return [
        i
        for i in range(1, len(seq) - 1)
        if seq[i - 1] > seq[i] and seq[i + 1] > seq[i]
    ]


def is_unimodal(seq: list[int]) -> bool:
    decreasing = False
    for left, right in zip(seq, seq[1:]):
        if right < left:
            decreasing = True
        elif right > left and decreasing:
            return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k-max", type=int, default=8)
    parser.add_argument("--ell-max", type=int, default=30)
    parser.add_argument("--n-max", type=int, default=30)
    parser.add_argument("--m-max", type=int, default=500)
    parser.add_argument("--degree-max", type=int, default=20_000)
    parser.add_argument("--tail", type=int, default=80)
    parser.add_argument("--progress-every", type=int, default=10_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/pattern_family_tail_scan_20260808.json"),
    )
    args = parser.parse_args()

    checked = 0
    closest: dict[str, object] | None = None
    certificate: dict[str, object] | None = None

    for k in range(1, args.k_max + 1):
        for n in range(1, args.n_max + 1):
            stride = k * (n + 1)
            for ell in range(args.ell_max + 1):
                max_m = min(args.m_max, (args.degree_max - ell - 1) // stride)
                for m in range(1, max_m + 1):
                    tail = reversed_u_tail(k, ell, n, m, args.tail)
                    checked += 1
                    if not is_unimodal(tail):
                        certificate = {
                            "k": k,
                            "ell": ell,
                            "n": n,
                            "m": m,
                            "degree": stride * m + ell + 1,
                            "tail": tail,
                            "strict_valleys": valleys(tail),
                        }
                        print(json.dumps({"event": "counterexample", **certificate}))
                        break

                    # Track the closest post-decrease rebound ratio in this tail.
                    decreasing = False
                    for j, (left, right) in enumerate(zip(tail, tail[1:])):
                        if right < left:
                            decreasing = True
                        elif decreasing and left:
                            candidate = (right / left, j)
                            if closest is None or candidate[0] > closest["ratio"]:
                                closest = {
                                    "ratio": candidate[0],
                                    "offset": j,
                                    "k": k,
                                    "ell": ell,
                                    "n": n,
                                    "m": m,
                                    "degree": stride * m + ell + 1,
                                }
                    if checked % args.progress_every == 0:
                        print(
                            json.dumps(
                                {
                                    "event": "progress",
                                    "checked": checked,
                                    "parameters": [k, ell, n, m],
                                    "closest": closest,
                                }
                            ),
                            flush=True,
                        )
                if certificate is not None:
                    break
            if certificate is not None:
                break
        if certificate is not None:
            break

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": "Theorem 4.11, arXiv:2603.14204v1",
        "parameters": {
            "k_max": args.k_max,
            "ell_max": args.ell_max,
            "n_max": args.n_max,
            "m_max": args.m_max,
            "degree_max": args.degree_max,
            "tail": args.tail,
        },
        "checked": checked,
        "closest": closest,
        "counterexample": certificate,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"event": "done", **report}), flush=True)


if __name__ == "__main__":
    main()
