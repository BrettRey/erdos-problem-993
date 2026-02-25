#!/usr/bin/env python3
"""Search fixed-lambda rooted-component multiset collisions on anchored aggregates.

This is a component-level probe for the missing injectivity bridge:
  (d, mu1, mu2, rho, sigma) -> N
at a fixed lambda.

Definitions for rooted components C_j:
  F_j(x) = I(C_j; x)
  G_j(x) = I(C_j - root_j; x)
  a_j    = lambda * F'_j(lambda) / F_j(lambda)
  b_j    = lambda^2 * F''_j(lambda) / F_j(lambda)
  r_j    = G_j(lambda) / F_j(lambda)
  g_j    = lambda * G'_j(lambda) / G_j(lambda)
  d_j    = deg(F_j)
  n_j    = |V(C_j)|

Aggregate key for a multiset:
  d      = sum d_j
  mu1    = sum a_j
  mu2    = sum b_j + 2 * sum_{i<j} a_i a_j
  rho    = lambda * product r_j
  sigma  = rho * (1 + sum g_j)
  N      = sum n_j

The script reports the first collision with same key but different N.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations_with_replacement
from typing import Any

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from indpoly import independence_poly
from trees import trees_geng_raw


def parse_fraction(s: str) -> Fraction:
    if "/" in s:
        a, b = s.split("/", 1)
        return Fraction(int(a), int(b))
    return Fraction(int(s), 1)


def poly_eval(poly: list[int], x: Fraction) -> Fraction:
    out = Fraction(0, 1)
    for c in reversed(poly):
        out = out * x + c
    return out


def deriv_eval(poly: list[int], order: int, x: Fraction) -> Fraction:
    out = Fraction(0, 1)
    for i, c in enumerate(poly):
        if i < order or c == 0:
            continue
        ff = 1
        for t in range(order):
            ff *= (i - t)
        out += c * ff * (x ** (i - order))
    return out


def remove_vertex(adj: list[list[int]], root: int) -> list[list[int]]:
    n = len(adj)
    keep = [v for v in range(n) if v != root]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for w in adj[v]:
            if w == root:
                continue
            out[vv].append(idx[w])
    return out


def frac_pair(x: Fraction) -> list[int]:
    return [x.numerator, x.denominator]


@dataclass(frozen=True)
class Component:
    n: int
    d: int
    f_poly: tuple[int, ...]
    g_poly: tuple[int, ...]
    a: Fraction
    b: Fraction
    r: Fraction
    g: Fraction


def build_components(max_comp_n: int, lam: Fraction) -> list[Component]:
    # De-duplicate rooted components by (F,G,n); this keeps stats exact.
    seen: dict[tuple[tuple[int, ...], tuple[int, ...], int], Component] = {}

    for n in range(1, max_comp_n + 1):
        for nn, adj, _ in trees_geng_raw(n):
            assert nn == n
            f_poly = independence_poly(n, adj)
            f_lam = poly_eval(f_poly, lam)
            if f_lam == 0:
                continue

            for root in range(n):
                g_adj = remove_vertex(adj, root)
                g_poly = independence_poly(n - 1, g_adj) if n > 1 else [1]
                g_lam = poly_eval(g_poly, lam)
                if g_lam == 0:
                    continue

                a = lam * deriv_eval(f_poly, 1, lam) / f_lam
                b = lam * lam * deriv_eval(f_poly, 2, lam) / f_lam
                r = g_lam / f_lam
                g = lam * deriv_eval(g_poly, 1, lam) / g_lam

                key = (tuple(f_poly), tuple(g_poly), n)
                if key not in seen:
                    seen[key] = Component(
                        n=n,
                        d=len(f_poly) - 1,
                        f_poly=tuple(f_poly),
                        g_poly=tuple(g_poly),
                        a=a,
                        b=b,
                        r=r,
                        g=g,
                    )

    return list(seen.values())


def aggregate_key(
    comps: list[Component], idxs: tuple[int, ...], lam: Fraction
) -> tuple[int, ...]:
    d = 0
    mu1 = Fraction(0, 1)
    t_sum = Fraction(0, 1)  # sum (b - a^2)
    r_prod = Fraction(1, 1)
    g_sum = Fraction(0, 1)
    for i in idxs:
        c = comps[i]
        d += c.d
        mu1 += c.a
        t_sum += c.b - c.a * c.a
        r_prod *= c.r
        g_sum += c.g
    mu2 = mu1 * mu1 + t_sum
    rho = lam * r_prod
    sigma = rho * (1 + g_sum)
    return (
        d,
        mu1.numerator,
        mu1.denominator,
        mu2.numerator,
        mu2.denominator,
        rho.numerator,
        rho.denominator,
        sigma.numerator,
        sigma.denominator,
    )


def aggregate_N(comps: list[Component], idxs: tuple[int, ...]) -> int:
    return sum(comps[i].n for i in idxs)


def encode_multiset(comps: list[Component], idxs: tuple[int, ...]) -> list[dict[str, Any]]:
    out = []
    for i in idxs:
        c = comps[i]
        out.append(
            {
                "n": c.n,
                "d": c.d,
                "F": list(c.f_poly),
                "G": list(c.g_poly),
                "a": frac_pair(c.a),
                "b": frac_pair(c.b),
                "r": frac_pair(c.r),
                "g": frac_pair(c.g),
            }
        )
    return out


def scan(
    lam: Fraction,
    max_comp_n: int,
    multiset_size: int,
    progress_every: int,
) -> dict[str, Any]:
    comps = build_components(max_comp_n=max_comp_n, lam=lam)
    m = len(comps)

    seen: dict[tuple[int, ...], tuple[tuple[int, ...], int]] = {}
    tested = 0
    collisions = 0
    witness: dict[str, Any] | None = None

    for idxs in combinations_with_replacement(range(m), multiset_size):
        tested += 1
        if progress_every and tested % progress_every == 0:
            print(f"tested={tested} seen={len(seen)} collisions={collisions}", flush=True)

        key = aggregate_key(comps, idxs, lam)
        nsum = aggregate_N(comps, idxs)
        prev = seen.get(key)
        if prev is None:
            seen[key] = (idxs, nsum)
            continue

        collisions += 1
        prev_idxs, prev_n = prev
        if prev_n != nsum:
            witness = {
                "lambda": frac_pair(lam),
                "max_comp_n": max_comp_n,
                "multiset_size": multiset_size,
                "key": list(key),
                "A_N": prev_n,
                "B_N": nsum,
                "A_multiset": encode_multiset(comps, prev_idxs),
                "B_multiset": encode_multiset(comps, idxs),
            }
            break

    return {
        "scan": "component_aggregate_collision_fixed_lambda",
        "lambda": frac_pair(lam),
        "max_comp_n": max_comp_n,
        "multiset_size": multiset_size,
        "component_library_size": m,
        "tested_multisets": tested,
        "unique_keys": len(seen),
        "collisions": collisions,
        "split_found": witness is not None,
        "split": witness,
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Fixed-lambda rooted-component aggregate collision scan."
    )
    ap.add_argument("--lambda", dest="lam", required=True, help="lambda as a/b or int")
    ap.add_argument("--max-comp-n", type=int, default=8)
    ap.add_argument("--multiset-size", type=int, default=2)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    lam = parse_fraction(args.lam)
    payload = scan(
        lam=lam,
        max_comp_n=args.max_comp_n,
        multiset_size=args.multiset_size,
        progress_every=args.progress_every,
    )

    print(
        f"done tested={payload['tested_multisets']} unique={payload['unique_keys']} "
        f"collisions={payload['collisions']} split_found={payload['split_found']}",
        flush=True,
    )
    if payload["split_found"]:
        print(
            f"N split: {payload['split']['A_N']} vs {payload['split']['B_N']}",
            flush=True,
        )

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()

