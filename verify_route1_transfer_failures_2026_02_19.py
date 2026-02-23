#!/usr/bin/env python3
"""Extract and analyze Route-1 transfer exact-violation witnesses.

This script reproduces the canonical leaf choice used in
`conjecture_a_route1_transfer_scan.py` and outputs all witnesses with
`exact_excess = D - 1/(1+lambda) > tol`, where

  D = mu_B(lambda) - mu_P(lambda),   lambda = i_{m-1}(T)/i_m(T).

It also records the decomposition terms:
  D = p_u * (1 - S),   S = sum_c delta_c,
  p_u = lambda R/(1 + lambda R),   R = prod_c (1 - p_c),
and child-level statistics that help characterize extremal failures.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from dataclasses import asdict, dataclass
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def choose_min_support_leaf(adj: list[list[int]]) -> tuple[int, int]:
    deg = [len(nb) for nb in adj]
    leaves = [v for v, d in enumerate(deg) if d == 1]
    parent = {l: adj[l][0] for l in leaves}
    # minimum support degree, tie by largest leaf id
    leaf = max(leaves, key=lambda l: (-deg[parent[l]], -l))
    return leaf, parent[leaf]


def remove_vertices(adj: list[list[int]], remove_set: set[int]) -> tuple[list[list[int]], dict[int, int]]:
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out, idx


def rooted_dp(adj: list[list[int]], root: int):
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in range(n)]
    parent[root] = root
    queue = [root]
    for v in queue:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v
                children[v].append(w)
                queue.append(w)

    order: list[int] = []
    stack: list[tuple[int, bool]] = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    dp0: list[list[int]] = [[] for _ in range(n)]
    dp1: list[list[int]] = [[] for _ in range(n)]
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
            continue
        p0 = [1]
        for c in children[v]:
            p0 = _polymul(p0, _polyadd(dp0[c], dp1[c]))
        dp0[v] = p0

        p1 = [1]
        for c in children[v]:
            p1 = _polymul(p1, dp0[c])
        dp1[v] = [0] + p1

    return dp0, dp1, children, parent, order


def eval_poly(poly: list[int], lam: float) -> float:
    val = 0.0
    p = 1.0
    for ck in poly:
        val += ck * p
        p *= lam
    return val


def mean_at_lambda(poly: list[int], lam: float) -> float:
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z if z else 0.0


@dataclass
class ChildStats:
    child: int
    deg: int
    n_descendants: int
    n_grandchildren: int
    n_leaf_grandchildren: int
    p_c: float
    delta_c: float
    mu_Ic: float
    mu_qc: float
    mu_dc: float
    sum_delta_gc: float


def subtree_size(children: list[list[int]], root: int) -> int:
    total = 0
    stack = [root]
    while stack:
        v = stack.pop()
        total += 1
        stack.extend(children[v])
    return total


def analyze_single_tree(nn: int, adj: list[list[int]], tol: float = 1e-12) -> dict[str, Any] | None:
    if not is_dleaf_le_1(nn, adj):
        return None

    poly_t = independence_poly(nn, adj)
    m = mode_index_leftmost(poly_t)
    if m == 0 or m - 1 >= len(poly_t) or poly_t[m - 1] == 0 or poly_t[m] == 0:
        return None
    lam = poly_t[m - 1] / poly_t[m]

    leaf, support = choose_min_support_leaf(adj)
    if len(adj[support]) != 2:
        return None
    u = adj[support][0] if adj[support][1] == leaf else adj[support][1]

    b_adj, idx = remove_vertices(adj, {leaf, support})
    u_in_b = idx[u]
    dp0, dp1, children_b, _, _ = rooted_dp(b_adj, u_in_b)

    p_poly = dp0[u_in_b]
    q_poly = dp1[u_in_b]
    b_poly = _polyadd(p_poly, q_poly)

    mu_p = mean_at_lambda(p_poly, lam)
    mu_b = mean_at_lambda(b_poly, lam)
    d_val = mu_b - mu_p

    threshold = 1.0 / (1.0 + lam)
    exact_excess = d_val - threshold
    if exact_excess <= tol:
        return None

    z_p = eval_poly(p_poly, lam)
    z_q = eval_poly(q_poly, lam)
    z_b = z_p + z_q
    p_u = z_q / z_b if z_b else 0.0

    children_stats: list[ChildStats] = []
    sum_delta = 0.0
    sum_p = 0.0
    r_prod = 1.0

    for c in children_b[u_in_b]:
        ic_poly = _polyadd(dp0[c], dp1[c])
        z_ic = eval_poly(ic_poly, lam)
        z_dc = eval_poly(dp1[c], lam)
        p_c = z_dc / z_ic if z_ic else 0.0
        mu_ic = mean_at_lambda(ic_poly, lam)
        mu_qc = mean_at_lambda(dp0[c], lam)
        mu_dc = mean_at_lambda(dp1[c], lam)
        delta_c = mu_ic - mu_qc

        # delta_c = p_c * (1 - sum_{gc} delta_gc)
        gc_sum = 0.0
        leaf_gc = 0
        for gc in children_b[c]:
            igc_poly = _polyadd(dp0[gc], dp1[gc])
            mu_igc = mean_at_lambda(igc_poly, lam)
            mu_qgc = mean_at_lambda(dp0[gc], lam)
            d_gc = mu_igc - mu_qgc
            gc_sum += d_gc
            if not children_b[gc]:
                leaf_gc += 1

        children_stats.append(
            ChildStats(
                child=c,
                deg=len(b_adj[c]),
                n_descendants=subtree_size(children_b, c),
                n_grandchildren=len(children_b[c]),
                n_leaf_grandchildren=leaf_gc,
                p_c=p_c,
                delta_c=delta_c,
                mu_Ic=mu_ic,
                mu_qc=mu_qc,
                mu_dc=mu_dc,
                sum_delta_gc=gc_sum,
            )
        )
        sum_delta += delta_c
        sum_p += p_c
        r_prod *= (1.0 - p_c)

    one_minus_sum_delta = 1.0 - sum_delta
    d_from_decomp = p_u * one_minus_sum_delta

    # algebraic numerator for exact_excess sign
    num = lam * r_prod * (lam - (1.0 + lam) * sum_delta) - 1.0
    den = (1.0 + lam * r_prod) * (1.0 + lam)
    exact_excess_from_num = num / den

    # Primary witness record
    out = {
        "n": nn,
        "g6": None,  # set by caller
        "mode": m,
        "lambda": lam,
        "leaf": leaf,
        "support": support,
        "u_in_T": u,
        "deg_u_in_T": len(adj[u]),
        "deg_u_in_B": len(b_adj[u_in_b]),
        "num_children_u": len(children_b[u_in_b]),
        "mu_p": mu_p,
        "mu_b": mu_b,
        "D": d_val,
        "threshold_1_over_1plam": threshold,
        "exact_excess": exact_excess,
        "p_u": p_u,
        "sum_p_children": sum_p,
        "sum_delta_children": sum_delta,
        "one_minus_sum_delta": one_minus_sum_delta,
        "R_prod_1_minus_pc": r_prod,
        "D_from_decomp": d_from_decomp,
        "exact_excess_from_num": exact_excess_from_num,
        "num_exact_sign": num,
        "den_exact_sign": den,
        "sum_delta_negative_children": sum(1 for ch in children_stats if ch.delta_c < 0.0),
        "sum_p_minus_sum_delta": sum_p - sum_delta,
        "children": [asdict(ch) for ch in children_stats],
    }
    return out


def run_extract(
    min_n: int,
    max_n: int,
    geng: str,
    res: int | None,
    mod: int | None,
    out: str,
    tol: float,
) -> None:
    if (res is None) != (mod is None):
        raise ValueError("Specify both --res and --mod together, or neither.")
    if mod is not None and not (0 <= res < mod):
        raise ValueError("Require 0 <= res < mod.")

    params = {
        "min_n": min_n,
        "max_n": max_n,
        "res": res,
        "mod": mod,
        "tol": tol,
    }

    seen = 0
    considered = 0
    checked = 0
    no_deg2_support = 0
    witnesses: list[dict[str, Any]] = []

    for n in range(min_n, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        if mod is not None:
            cmd.append(f"{res}/{mod}")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        for raw in proc.stdout:
            seen += 1
            nn, adj = parse_graph6(raw)
            if not is_dleaf_le_1(nn, adj):
                continue
            considered += 1

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                no_deg2_support += 1
                continue
            checked += 1

            rec = analyze_single_tree(nn, adj, tol=tol)
            if rec is None:
                continue
            rec["g6"] = raw.decode("ascii").strip()
            witnesses.append(rec)

        proc.wait()

    witnesses.sort(key=lambda r: r["exact_excess"], reverse=True)
    summary = {
        "seen": seen,
        "considered": considered,
        "checked": checked,
        "no_deg2_support": no_deg2_support,
        "exact_violation_count": len(witnesses),
        "max_exact_excess": witnesses[0]["exact_excess"] if witnesses else None,
        "max_exact_excess_witness": witnesses[0] if witnesses else None,
    }

    payload = {
        "params": params,
        "summary": summary,
        "witnesses": witnesses,
    }

    out_dir = os.path.dirname(out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"Wrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Extract Route-1 exact transfer violations with decomposition.")
    ap.add_argument("--min-n", type=int, required=True)
    ap.add_argument("--max-n", type=int, required=True)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--res", type=int, default=None)
    ap.add_argument("--mod", type=int, default=None)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    run_extract(
        min_n=args.min_n,
        max_n=args.max_n,
        geng=args.geng,
        res=args.res,
        mod=args.mod,
        out=args.out,
        tol=args.tol,
    )


if __name__ == "__main__":
    main()
