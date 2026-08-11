#!/usr/bin/env python3
"""Kill-test low-order rooted-tree invariants for STP2 product closure.

For a rooted n-vertex tree let I be its independence polynomial and let E be
the polynomial after deleting the root.  Necessarily

    I[0] = E[0] = 1,
    I[1] = n, E[1] = n - 1,
    I[2] = C(n - 1, 2),
    I - E = x R,

where R is the independence polynomial after deleting the closed root
neighborhood.  In particular, if q = R[1], then

    E[2] = I[2] - q

and R[2] is the number of independent pairs in a forest on q vertices.

This script exhausts short integer pairs satisfying those identities plus the
abstract shape package used in Formal/STP2Closure.lean.  It then tests whether
the product pair preserves the guarded ladder/STP2 inequalities.  A returned
witness is an obstruction to this *proposed invariant*, not a tree.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from dataclasses import dataclass


def conv(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def coeff(a: tuple[int, ...], k: int) -> int:
    return a[k] if 0 <= k < len(a) else 0


def is_lc(a: tuple[int, ...]) -> bool:
    return all(a[k] * a[k] >= a[k - 1] * coeff(a, k + 1)
               for k in range(1, len(a)))


def is_stp2(i_poly: tuple[int, ...], e_poly: tuple[int, ...]) -> bool:
    limit = max(len(i_poly), len(e_poly)) + 1
    return all(coeff(i_poly, k) * coeff(e_poly, k)
               >= coeff(i_poly, k - 1) * coeff(e_poly, k + 1)
               for k in range(1, limit))


def first_stp2_failure(
    i_poly: tuple[int, ...], e_poly: tuple[int, ...]
) -> tuple[int, int] | None:
    limit = max(len(i_poly), len(e_poly)) + 1
    for k in range(1, limit):
        gap = (coeff(i_poly, k) * coeff(e_poly, k)
               - coeff(i_poly, k - 1) * coeff(e_poly, k + 1))
        if gap < 0:
            return k, gap
    return None


@dataclass(frozen=True)
class Pair:
    n: int
    q: int
    i_poly: tuple[int, ...]
    e_poly: tuple[int, ...]
    r_poly: tuple[int, ...]


def forest_pair_bounds(q: int) -> tuple[int, int]:
    """Bounds on independent pairs in an arbitrary forest on q vertices."""
    if q <= 1:
        return 0, 0
    total_pairs = math.comb(q, 2)
    return max(0, total_pairs - (q - 1)), total_pairs


def generate_pairs(max_n: int, max_last: int, degree: int) -> list[Pair]:
    if degree not in (2, 3, 4):
        raise ValueError("degree must be 2, 3, or 4")

    pairs: list[Pair] = []
    for n in range(3, max_n + 1):
        i2 = math.comb(n - 1, 2)
        for q in range(n - 1):
            e2 = i2 - q
            if e2 < 0:
                continue

            r2_lo, r2_hi = forest_pair_bounds(q)
            if degree == 2:
                r_tails = [()]
            elif degree == 3:
                r_tails = [(r2,) for r2 in range(r2_lo, r2_hi + 1)]
            else:
                r_tails = [
                    (r2, r3)
                    for r2 in range(r2_lo, r2_hi + 1)
                    for r3 in range(max_last + 1)
                ]

            for r_tail in r_tails:
                r_poly = (1, q) + r_tail
                if not is_lc(r_poly):
                    continue

                # I[k+1] - E[k+1] = R[k].  Enumerate E beyond order 2;
                # the low-order rooted-tree identities then determine I.
                tail_len = degree - 2
                for e_tail in itertools.product(range(max_last + 1), repeat=tail_len):
                    e_poly = (1, n - 1, e2) + e_tail
                    i_tail = tuple(
                        e_tail[j] + r_poly[j + 2]
                        for j in range(tail_len)
                    )
                    i_poly = (1, n, i2) + i_tail

                    # Contiguous positive support at the truncated degrees.
                    if any(x <= 0 for x in i_poly) or any(x <= 0 for x in e_poly):
                        continue
                    if not (is_lc(i_poly) and is_lc(e_poly)):
                        continue
                    if not is_stp2(i_poly, e_poly):
                        continue
                    if any(e > i for e, i in zip(e_poly, i_poly)):
                        continue
                    pairs.append(Pair(n, q, i_poly, e_poly, r_poly))
    return pairs


def random_lc_tail(
    prefix: tuple[int, ...], degree: int, rng: random.Random
) -> tuple[int, ...] | None:
    """Extend an integer LC prefix to exactly ``degree`` with positive terms."""
    out = list(prefix)
    while len(out) <= degree:
        upper = out[-1] * out[-1] // out[-2]
        if upper < 1:
            return None
        # Mix boundary-biased and log-uniform choices.  Extremal failures often
        # live near an LC boundary, while log sampling explores long tails.
        if rng.random() < 0.35:
            value = max(1, upper - rng.randrange(min(upper, 8)))
        else:
            value = max(1, int(math.exp(rng.random() * math.log(upper + 1))))
        out.append(value)
    return tuple(out)


def random_pair(max_n: int, rng: random.Random) -> Pair | None:
    n = rng.randint(4, max_n)
    q = rng.randint(1, n - 2)
    i2 = math.comb(n - 1, 2)
    e2 = i2 - q

    min_deg_e = max((n - 1 + 1) // 2, n - 1 - q)
    min_deg_r = (q + 1) // 2
    deg_e = rng.randint(max(2, min_deg_e), n - 1)
    deg_r = rng.randint(max(1, min_deg_r), q)

    e_poly = random_lc_tail((1, n - 1, e2), deg_e, rng)
    r2_lo, r2_hi = forest_pair_bounds(q)
    if deg_r == 1:
        if r2_lo > 0:
            return None
        r_poly = (1, q)
    else:
        r2 = rng.randint(r2_lo, r2_hi)
        r_poly = random_lc_tail((1, q, r2), deg_r, rng)
    if e_poly is None or r_poly is None:
        return None

    if any(value > math.comb(n - 1, k) for k, value in enumerate(e_poly)):
        return None
    if any(value > math.comb(q, k) for k, value in enumerate(r_poly)):
        return None

    degree_i = max(deg_e, deg_r + 1)
    i_poly = tuple(
        coeff(e_poly, k) + coeff(r_poly, k - 1)
        for k in range(degree_i + 1)
    )
    if not is_lc(i_poly) or not is_stp2(i_poly, e_poly):
        return None
    if any(value > math.comb(n, k) for k, value in enumerate(i_poly)):
        return None
    return Pair(n, q, i_poly, e_poly, r_poly)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=14)
    parser.add_argument("--max-last", type=int, default=120)
    parser.add_argument("--degree", type=int, choices=(2, 3, 4), default=3)
    parser.add_argument("--distinct", action="store_true",
                        help="test all ordered pairs, not only self-products")
    parser.add_argument("--samples", type=int, default=0,
                        help="sample this many ordered products (0 means exhaustive)")
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--random-trials", type=int, default=0,
                        help="generate and test this many random integer pairs")
    args = parser.parse_args()

    if args.random_trials:
        rng = random.Random(args.seed)
        pairs = []
        tested = 0
        previous: list[Pair] = []
        while tested < args.random_trials:
            pair = random_pair(args.max_n, rng)
            tested += 1
            if pair is None:
                continue
            partners = [pair]
            if previous:
                partners.extend(rng.choice(previous) for _ in range(min(4, len(previous))))
            for other in partners:
                prod_i = conv(pair.i_poly, other.i_poly)
                prod_e = conv(pair.e_poly, other.e_poly)
                failure = first_stp2_failure(prod_i, prod_e)
                if failure is None:
                    continue
                k, gap = failure
                print(json.dumps({
                    "found": True,
                    "random_trials": tested,
                    "failure_k": k,
                    "failure_gap": gap,
                    "left": pair.__dict__,
                    "right": other.__dict__,
                    "product_i": prod_i,
                    "product_e": prod_e,
                }, default=list, indent=2))
                return
            previous.append(pair)
            if len(previous) > 5000:
                previous[rng.randrange(len(previous))] = previous[-1]
                previous.pop()
        print(json.dumps({
            "found": False,
            "random_trials": tested,
            "valid_pairs": len(previous),
        }, indent=2))
        return

    pairs = generate_pairs(args.max_n, args.max_last, args.degree)
    print(json.dumps({"candidate_pairs": len(pairs), "degree": args.degree}))

    if args.samples:
        rng = random.Random(args.seed)
        products = ((rng.choice(pairs), rng.choice(pairs)) for _ in range(args.samples))
    elif args.distinct:
        products = itertools.product(pairs, repeat=2)
    else:
        products = ((p, p) for p in pairs)
    tested = 0
    for left, right in products:
        tested += 1
        prod_i = conv(left.i_poly, right.i_poly)
        prod_e = conv(left.e_poly, right.e_poly)
        failure = first_stp2_failure(prod_i, prod_e)
        if failure is None:
            continue
        k, gap = failure
        print(json.dumps({
            "found": True,
            "tested_products": tested,
            "failure_k": k,
            "failure_gap": gap,
            "left": left.__dict__,
            "right": right.__dict__,
            "product_i": prod_i,
            "product_e": prod_e,
        }, default=list, indent=2))
        return

    print(json.dumps({"found": False, "tested_products": tested}, indent=2))


if __name__ == "__main__":
    main()
