#!/usr/bin/env python3
"""Exhaustive canonical K* split scan with explicit min-u tie-break.

Canonical triplet rule:
  among all admissible degree-2 bridge triplets (leaf,support,u),
  pick lexicographically by (u, support, leaf) minimum.

Admissible triplet:
  - leaf has degree 1
  - support = unique neighbor of leaf
  - deg(support) = 2
  - u = other neighbor of support

Gate:
  - is_dleaf_le_1(n, adj) is True
  - at least one admissible triplet exists

For canonical triplet (leaf,support,u):
  B = T \\ {leaf,support}, rooted at u
  P = dp0[u], Q = dp1[u], I(T) = independence_poly(T)
  m = leftmost mode index of I(T)
  lambda = i_{m-1}/i_m
  K* = (d,m,lambda,mu1,mu2,rho,sigma), where
    d = deg(P)
    mu1 = lambda*P'(lambda)/P(lambda)
    mu2 = lambda^2*P''(lambda)/P(lambda)
    rho = Q(lambda)/P(lambda)
    sigma = lambda*Q'(lambda)/P(lambda)

Split condition:
  same K* but different N=[x]P (equiv. i1-3).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from typing import Any

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly
from trees import trees_geng_raw


def coeff(poly: list[int], k: int) -> int:
    return poly[k] if 0 <= k < len(poly) else 0


def frac_pair(x: Fraction) -> list[int]:
    return [x.numerator, x.denominator]


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


def rooted_dp(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[list[int]]]:
    n = len(adj)
    children = [[] for _ in range(n)]
    parent = [-1] * n
    order: list[int] = []
    stack = [root]
    parent[root] = root
    while stack:
        v = stack.pop()
        order.append(v)
        for w in adj[v]:
            if parent[w] != -1:
                continue
            parent[w] = v
            children[v].append(w)
            stack.append(w)

    dp0: list[list[int]] = [[] for _ in range(n)]
    dp1: list[list[int]] = [[] for _ in range(n)]

    for v in reversed(order):
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
            continue
        prod0 = [1]
        prod1 = [1]
        for c in children[v]:
            prod0 = _polymul(prod0, _polyadd(dp0[c], dp1[c]))
            prod1 = _polymul(prod1, dp0[c])
        dp0[v] = prod0
        dp1[v] = [0] + prod1
    return dp0, dp1


def admissible_triplets(adj: list[list[int]]) -> list[tuple[int, int, int]]:
    n = len(adj)
    out: list[tuple[int, int, int]] = []
    for leaf in range(n):
        if len(adj[leaf]) != 1:
            continue
        support = adj[leaf][0]
        if len(adj[support]) != 2:
            continue
        u = adj[support][0] if adj[support][1] == leaf else adj[support][1]
        out.append((leaf, support, u))
    return out


def choose_triplet_min_u(adj: list[list[int]]) -> tuple[int, int, int] | None:
    trips = admissible_triplets(adj)
    if not trips:
        return None
    return min(trips, key=lambda t: (t[2], t[1], t[0]))


def build_record_from_triplet(n: int, adj: list[list[int]], g6: str, trip: tuple[int, int, int]) -> dict[str, Any] | None:
    leaf, support, u = trip
    i_poly = independence_poly(n, adj)
    m = mode_index_leftmost(i_poly)
    i_m = coeff(i_poly, m)
    if i_m == 0:
        return None
    lam = Fraction(coeff(i_poly, m - 1), i_m)

    b_adj = remove_vertices(adj, {leaf, support})
    keep = [v for v in range(n) if v not in {leaf, support}]
    idx = {old: new for new, old in enumerate(keep)}
    if u not in idx:
        return None
    u_in_b = idx[u]
    dp0, dp1 = rooted_dp(b_adj, u_in_b)
    p_poly = dp0[u_in_b]
    q_poly = dp1[u_in_b]

    p_lam = poly_eval(p_poly, lam)
    if p_lam == 0:
        return None
    mu1 = lam * deriv_eval(p_poly, 1, lam) / p_lam
    mu2 = lam * lam * deriv_eval(p_poly, 2, lam) / p_lam
    rho = poly_eval(q_poly, lam) / p_lam
    sigma = lam * deriv_eval(q_poly, 1, lam) / p_lam

    return {
        "n": n,
        "g6": g6,
        "leaf": leaf,
        "support": support,
        "u": u,
        "triplet": [leaf, support, u],
        "triplet_count": len(admissible_triplets(adj)),
        "d": len(p_poly) - 1,
        "m": m,
        "lambda": frac_pair(lam),
        "mu1": frac_pair(mu1),
        "mu2": frac_pair(mu2),
        "rho": frac_pair(rho),
        "sigma": frac_pair(sigma),
        "i1": coeff(i_poly, 1),
        "N": coeff(p_poly, 1),
        "P": p_poly,
        "Q": q_poly,
        "I": i_poly,
    }


def key_from_record(rec: dict[str, Any]) -> tuple[int, ...]:
    d = rec["d"]
    m = rec["m"]
    lam_n, lam_d = rec["lambda"]
    mu1_n, mu1_d = rec["mu1"]
    mu2_n, mu2_d = rec["mu2"]
    rho_n, rho_d = rec["rho"]
    sig_n, sig_d = rec["sigma"]
    return (
        d,
        m,
        lam_n,
        lam_d,
        mu1_n,
        mu1_d,
        mu2_n,
        mu2_d,
        rho_n,
        rho_d,
        sig_n,
        sig_d,
    )


def rebuild_record(n: int, g6: str) -> dict[str, Any]:
    nn, adj = parse_graph6(g6.encode("ascii"))
    if nn != n:
        raise RuntimeError(f"n mismatch for {g6}: expected {n}, got {nn}")
    trip = choose_triplet_min_u(adj)
    if trip is None:
        raise RuntimeError(f"no admissible triplet for {g6}")
    rec = build_record_from_triplet(n, adj, g6, trip)
    if rec is None:
        raise RuntimeError(f"could not build record for {g6}")
    return rec


def scan(
    min_n: int,
    max_n: int,
    progress_every: int,
    within_n_only: bool,
    m_min: int,
) -> dict[str, Any]:
    started = time.time()
    seen: dict[tuple[int, ...], tuple[int, str, int, int]] = {}

    checked_total = 0
    skipped_dleaf = 0
    skipped_no_triplet = 0
    skipped_bad = 0
    skipped_m = 0
    collisions = 0
    split: dict[str, Any] | None = None
    per_n: list[dict[str, Any]] = []

    for n in range(min_n, max_n + 1):
        if within_n_only:
            seen.clear()
        total_n = 0
        checked_n = 0
        t0 = time.time()
        for nn, adj, raw in trees_geng_raw(n):
            total_n += 1
            if progress_every > 0 and total_n % progress_every == 0:
                print(
                    f"progress n={n} total={total_n} checked={checked_n} "
                    f"unique={len(seen)} collisions={collisions}",
                    flush=True,
                )
            if not is_dleaf_le_1(nn, adj):
                skipped_dleaf += 1
                continue
            trip = choose_triplet_min_u(adj)
            if trip is None:
                skipped_no_triplet += 1
                continue
            g6 = raw.decode("ascii").strip()
            rec = build_record_from_triplet(nn, adj, g6, trip)
            if rec is None:
                skipped_bad += 1
                continue
            if rec["m"] < m_min:
                skipped_m += 1
                continue
            key = key_from_record(rec)
            prev = seen.get(key)
            if prev is None:
                seen[key] = (nn, g6, rec["i1"], rec["N"])
            else:
                collisions += 1
                pn, pg6, pi1, pN = prev
                if pi1 != rec["i1"] or pN != rec["N"]:
                    split = {"A": rebuild_record(pn, pg6), "B": rec}
                    break
            checked_total += 1
            checked_n += 1

        per_n.append(
            {
                "n": n,
                "total_trees": total_n,
                "checked": checked_n,
                "unique_keys": len(seen),
                "collisions": collisions,
                "elapsed_sec": time.time() - t0,
            }
        )
        print(
            f"n={n:2d} total={total_n:8d} checked={checked_n:7d} "
            f"unique={len(seen):7d} collisions={collisions:7d} "
            f"time={time.time()-t0:.2f}s",
            flush=True,
        )
        if split is not None:
            break

    return {
        "scan": "canonical_kstar_split_scan_minu",
        "triplet_rule": "min_u_then_support_then_leaf",
        "key": "(d,m,lambda,mu1,mu2,rho,sigma)",
        "min_n": min_n,
        "max_n": per_n[-1]["n"] if per_n else min_n,
        "within_n_only": within_n_only,
        "m_min": m_min,
        "checked_total": checked_total,
        "unique_keys": len(seen),
        "collisions": collisions,
        "skipped_dleaf": skipped_dleaf,
        "skipped_no_triplet": skipped_no_triplet,
        "skipped_bad": skipped_bad,
        "skipped_m": skipped_m,
        "split_found": split is not None,
        "split": split,
        "per_n": per_n,
        "elapsed_sec": time.time() - started,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=22)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument(
        "--within-n-only",
        action="store_true",
        help="Reset key table per n (intra-layer collisions only).",
    )
    ap.add_argument(
        "--m-min",
        type=int,
        default=0,
        help="Ignore instances with mode index m < this value.",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    payload = scan(
        min_n=args.min_n,
        max_n=args.max_n,
        progress_every=args.progress_every,
        within_n_only=args.within_n_only,
        m_min=args.m_min,
    )
    print(
        f"done checked={payload['checked_total']} unique={payload['unique_keys']} "
        f"collisions={payload['collisions']} split_found={payload['split_found']} "
        f"elapsed={payload['elapsed_sec']:.2f}s",
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
