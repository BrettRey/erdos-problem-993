#!/usr/bin/env python3
"""Targeted scan for same-K* / different-N splits in a star-of-components family.

Construct trees of the form:
  leaf(0) - support(1) - u(2)
and attach a multiset of rooted component trees to u.

For each constructed tree passing canonical gates via bridge_decomposition(..., require_dleaf=True),
compute
  K* = (d, m, lambda_hat, mu1, mu2, rho, sigma)
where
  d      = deg(P)
  m      = leftmost mode index of I(T)
  lambda = i_{m-1}/i_m
  mu1    = lambda P'(lambda)/P(lambda)
  mu2    = lambda^2 P''(lambda)/P(lambda)
  rho    = Q(lambda)/P(lambda)
  sigma  = lambda Q'(lambda)/P(lambda)

Report first key collision with different i1/N.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations_with_replacement
from typing import Any

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from attack4_common import bridge_decomposition
from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import independence_poly
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


def component_is_dleaf_compatible(adj: list[list[int]], root: int) -> bool:
    """Whether attaching a parent to `root` preserves d_leaf<=1 locally."""
    n = len(adj)
    ext = [list(nei) for nei in adj] + [[]]
    parent = n
    ext[root].append(parent)
    ext[parent].append(root)
    ext = [sorted(nei) for nei in ext]
    return is_dleaf_le_1(n + 1, ext)


def maybe_g6(adj: list[list[int]]) -> str:
    try:
        import networkx as nx  # type: ignore
    except Exception:
        return ""
    g = nx.Graph()
    n = len(adj)
    g.add_nodes_from(range(n))
    for v in range(n):
        for w in adj[v]:
            if v < w:
                g.add_edge(v, w)
    return nx.to_graph6_bytes(g, header=False).decode("ascii").strip()


@dataclass(frozen=True)
class RootedComponent:
    n: int
    root: int
    adj: tuple[tuple[int, ...], ...]
    f_poly: tuple[int, ...]
    g_poly: tuple[int, ...]
    deg_root: int
    signature: tuple[Any, ...]


def build_rooted_library(
    min_comp_n: int,
    max_comp_n: int,
    root_min_deg: int,
    require_component_dleaf: bool,
) -> list[RootedComponent]:
    seen: dict[tuple[Any, ...], RootedComponent] = {}
    for n in range(min_comp_n, max_comp_n + 1):
        for nn, adj, _ in trees_geng_raw(n):
            assert nn == n
            f_poly = independence_poly(n, adj)
            for root in range(n):
                if len(adj[root]) < root_min_deg:
                    continue
                if require_component_dleaf and not component_is_dleaf_compatible(adj, root):
                    continue
                g_adj = remove_vertex(adj, root)
                g_poly = independence_poly(n - 1, g_adj) if n > 1 else [1]
                sig = (n, tuple(f_poly), tuple(g_poly))
                if sig in seen:
                    continue
                comp = RootedComponent(
                    n=n,
                    root=root,
                    adj=tuple(tuple(nei) for nei in adj),
                    f_poly=tuple(f_poly),
                    g_poly=tuple(g_poly),
                    deg_root=len(adj[root]),
                    signature=sig,
                )
                seen[sig] = comp
    return list(seen.values())


def build_tree_from_components(comps: list[RootedComponent], idxs: tuple[int, ...]) -> list[list[int]]:
    # Base vertices: 0-1-2 where (0,1,2) are intended (leaf,support,u)
    total = 3 + sum(comps[i].n for i in idxs)
    adj = [[] for _ in range(total)]

    def add_edge(a: int, b: int) -> None:
        adj[a].append(b)
        adj[b].append(a)

    add_edge(0, 1)
    add_edge(1, 2)

    off = 3
    for i in idxs:
        c = comps[i]
        # copy component internal edges
        for v in range(c.n):
            for w in c.adj[v]:
                if v < w:
                    add_edge(off + v, off + w)
        # attach root to u=2
        add_edge(2, off + c.root)
        off += c.n

    return [sorted(nei) for nei in adj]


def kstar_key_and_record(decomp: Any) -> tuple[tuple[int, ...], dict[str, Any]] | None:
    p = decomp.p_poly
    q = decomp.q_poly
    i_poly = decomp.poly_t
    m = decomp.m_t

    i_m = coeff(i_poly, m)
    if i_m == 0:
        return None
    lam = Fraction(coeff(i_poly, m - 1), i_m)
    p_lam = poly_eval(p, lam)
    if p_lam == 0:
        return None

    mu1 = lam * deriv_eval(p, 1, lam) / p_lam
    mu2 = lam * lam * deriv_eval(p, 2, lam) / p_lam
    rho = poly_eval(q, lam) / p_lam
    sigma = lam * deriv_eval(q, 1, lam) / p_lam

    key = (
        len(p) - 1,
        m,
        lam.numerator,
        lam.denominator,
        mu1.numerator,
        mu1.denominator,
        mu2.numerator,
        mu2.denominator,
        rho.numerator,
        rho.denominator,
        sigma.numerator,
        sigma.denominator,
    )
    rec = {
        "n": decomp.n,
        "g6": decomp.g6,
        "leaf": decomp.leaf,
        "support": decomp.support,
        "u": decomp.u,
        "d": len(p) - 1,
        "m": m,
        "lambda_hat": frac_pair(lam),
        "mu1": frac_pair(mu1),
        "mu2": frac_pair(mu2),
        "rho": frac_pair(rho),
        "sigma": frac_pair(sigma),
        "i1": coeff(i_poly, 1),
        "N": coeff(p, 1),
        "P": p,
        "Q": q,
        "I": i_poly,
    }
    return key, rec


def scan(
    min_comp_n: int,
    max_comp_n: int,
    root_min_deg: int,
    require_component_dleaf: bool,
    multiset_size: int,
    max_tested: int,
    random_samples: int,
    seed: int,
    progress_every: int,
) -> dict[str, Any]:
    started = time.time()

    comps = build_rooted_library(
        min_comp_n=min_comp_n,
        max_comp_n=max_comp_n,
        root_min_deg=root_min_deg,
        require_component_dleaf=require_component_dleaf,
    )
    m = len(comps)

    seen: dict[tuple[int, ...], dict[str, Any]] = {}
    tested = 0
    passed_gate = 0
    collisions = 0
    split: dict[str, Any] | None = None

    if random_samples > 0:
        rng = random.Random(seed)
        def iter_multisets():
            for _ in range(random_samples):
                draw = [rng.randrange(m) for _ in range(multiset_size)]
                draw.sort()
                yield tuple(draw)
    else:
        def iter_multisets():
            for idxs in combinations_with_replacement(range(m), multiset_size):
                yield idxs

    for idxs in iter_multisets():
        tested += 1
        if max_tested > 0 and tested > max_tested:
            break
        if progress_every > 0 and tested % progress_every == 0:
            print(
                f"tested={tested} passed_gate={passed_gate} unique={len(seen)} collisions={collisions}",
                flush=True,
            )

        adj = build_tree_from_components(comps, idxs)
        n = len(adj)
        g6 = maybe_g6(adj)
        decomp = bridge_decomposition(n, adj, g6, require_dleaf=True)
        if decomp is None:
            continue
        passed_gate += 1

        out = kstar_key_and_record(decomp)
        if out is None:
            continue
        key, rec = out
        rec["component_indices"] = list(idxs)
        rec["component_signatures"] = [list(comps[i].signature) for i in idxs]

        prev = seen.get(key)
        if prev is None:
            seen[key] = rec
            continue
        collisions += 1
        if prev["N"] != rec["N"] or prev["i1"] != rec["i1"]:
            split = {"A": prev, "B": rec}
            break

    return {
        "scan": "star_component_kstar_split_search",
        "min_comp_n": min_comp_n,
        "max_comp_n": max_comp_n,
        "root_min_deg": root_min_deg,
        "require_component_dleaf": require_component_dleaf,
        "multiset_size": multiset_size,
        "max_tested": max_tested,
        "random_samples": random_samples,
        "seed": seed,
        "component_library_size": m,
        "tested_multisets": tested,
        "passed_gate": passed_gate,
        "unique_keys": len(seen),
        "collisions": collisions,
        "split_found": split is not None,
        "split": split,
        "elapsed_sec": time.time() - started,
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Targeted scan for same-K* / different-N splits in star-component family."
    )
    ap.add_argument("--min-comp-n", type=int, default=2)
    ap.add_argument("--max-comp-n", type=int, default=6)
    ap.add_argument("--root-min-deg", type=int, default=1)
    ap.add_argument(
        "--require-component-dleaf",
        action="store_true",
        help="Filter rooted components to those satisfying d_leaf<=1 after root-parent attachment.",
    )
    ap.add_argument("--multiset-size", type=int, default=4)
    ap.add_argument("--max-tested", type=int, default=0, help="0 = no cap")
    ap.add_argument(
        "--random-samples",
        type=int,
        default=0,
        help="If >0, sample this many random multisets (ignores lexicographic full scan).",
    )
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    payload = scan(
        min_comp_n=args.min_comp_n,
        max_comp_n=args.max_comp_n,
        root_min_deg=args.root_min_deg,
        require_component_dleaf=args.require_component_dleaf,
        multiset_size=args.multiset_size,
        max_tested=args.max_tested,
        random_samples=args.random_samples,
        seed=args.seed,
        progress_every=args.progress_every,
    )

    print(
        f"done tested={payload['tested_multisets']} passed_gate={payload['passed_gate']} "
        f"unique={payload['unique_keys']} collisions={payload['collisions']} "
        f"split_found={payload['split_found']} elapsed={payload['elapsed_sec']:.2f}s",
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
