#!/usr/bin/env python3
"""Scan the Route-2 -> Route-1 transfer gap on canonical degree-2 leaves.

For each d_leaf<=1 tree T:
  - choose canonical leaf l by minimum support degree (tie: largest leaf id)
  - require deg(s)=2, with support s and other neighbor u
  - let B = T-{l,s}
  - let P = dp_B[u][0]
  - m = mode(I(T)), lambda = i_{m-1}(T)/i_m(T)

Define transfer difference:
  D := mu_B(lambda) - mu_P(lambda)

Two useful benchmark excesses:
  half_excess  := D - 1/2
  exact_excess := D - (1 - lambda/(1+lambda))

If exact_excess <= 0, then Route-2 exact threshold
  mu_B >= m - 1 - lambda/(1+lambda)
would imply Route-1 target
  mu_P >= m - 2.

This scanner quantifies how far exact_excess is from <= 0.
Supports split runs via nauty residue classes (r/m).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
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


def rooted_dp_u(adj: list[list[int]], root: int) -> tuple[list[int], list[int]]:
    """Return (dp[root][0], dp[root][1]) for tree rooted at root."""
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

    return dp0[root], dp1[root]


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


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "no_deg2_support": 0,
        "half_violation": 0,
        "exact_violation": 0,
        "max_D": None,
        "max_D_witness": None,
        "max_half_excess": None,
        "max_half_excess_witness": None,
        "max_exact_excess": None,
        "max_exact_excess_witness": None,
        "max_p_u": None,
        "max_p_u_witness": None,
        "wall_s": 0.0,
    }


def maybe_update_max(
    stats: dict[str, Any],
    key: str,
    witness_key: str,
    value: float,
    witness: dict[str, Any],
) -> None:
    cur = stats[key]
    if cur is None or value > cur:
        stats[key] = value
        stats[witness_key] = witness


def merge_stats(dst: dict[str, Any], src: dict[str, Any]) -> None:
    for k in ["seen", "considered", "checked", "no_deg2_support", "half_violation", "exact_violation"]:
        dst[k] += int(src.get(k, 0))

    for key, witness_key in [
        ("max_D", "max_D_witness"),
        ("max_half_excess", "max_half_excess_witness"),
        ("max_exact_excess", "max_exact_excess_witness"),
        ("max_p_u", "max_p_u_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_max(dst, key, witness_key, float(src[key]), src[witness_key])


def aggregate(per_n: dict[str, dict[str, Any]]) -> dict[str, Any]:
    out = fresh_stats()
    for key in sorted(per_n, key=lambda s: int(s)):
        merge_stats(out, per_n[key])
    out.pop("wall_s", None)
    return out


def write_payload(out_path: str, params: dict[str, Any], per_n: dict[str, dict[str, Any]]) -> None:
    payload = {
        "params": params,
        "summary": aggregate(per_n),
        "per_n": {k: per_n[k] for k in sorted(per_n, key=lambda s: int(s))},
    }
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def fmt_opt(v: float | None) -> str:
    return "None" if v is None else f"{v:.6g}"


def main() -> None:
    ap = argparse.ArgumentParser(description="Route-1 transfer-gap scan (checkpointed).")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument("--res", type=int, default=None, help="Optional split residue r in r/m.")
    ap.add_argument("--mod", type=int, default=None, help="Optional split modulus m in r/m.")
    ap.add_argument("--out", default="results/whnc_route1_transfer_scan_n23.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    if (args.res is None) != (args.mod is None):
        raise ValueError("Specify both --res and --mod together, or neither.")
    if args.mod is not None and not (0 <= args.res < args.mod):
        raise ValueError("Require 0 <= res < mod.")

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "all_trees": args.all_trees,
        "tol": args.tol,
        "res": args.res,
        "mod": args.mod,
    }

    per_n: dict[str, dict[str, Any]] = {}
    if (not args.no_resume) and args.out and os.path.exists(args.out):
        with open(args.out, "r", encoding="utf-8") as f:
            old = json.load(f)
        per_n = old.get("per_n", {})
        done = ", ".join(sorted(per_n, key=lambda s: int(s)))
        print(f"Resuming from {args.out}; completed n: {done}", flush=True)

    scope = "all trees" if args.all_trees else "d_leaf<=1"
    split_desc = f", split={args.res}/{args.mod}" if args.mod is not None else ""
    print(
        f"Route-1 transfer scan on {scope}, n={args.min_n}..{args.max_n}{split_desc}",
        flush=True,
    )
    print("Canonical leaf: minimum support degree (tie: largest leaf id)", flush=True)
    print(f"Output: {args.out}", flush=True)
    print("-" * 96, flush=True)

    t_all = time.time()
    for n in range(args.min_n, args.max_n + 1):
        n_key = str(n)
        if n_key in per_n:
            print(f"n={n:2d}: already complete, skipping", flush=True)
            continue

        t0 = time.time()
        stats = fresh_stats()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        if args.mod is not None:
            cmd.append(f"{args.res}/{args.mod}")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        for raw in proc.stdout:
            stats["seen"] += 1
            nn, adj = parse_graph6(raw)

            if (not args.all_trees) and (not is_dleaf_le_1(nn, adj)):
                continue
            stats["considered"] += 1

            poly_t = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_t)
            if m == 0 or m - 1 >= len(poly_t) or poly_t[m - 1] == 0 or poly_t[m] == 0:
                continue
            lam = poly_t[m - 1] / poly_t[m]

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                stats["no_deg2_support"] += 1
                continue

            stats["checked"] += 1
            u = adj[support][0] if adj[support][1] == leaf else adj[support][1]
            b_adj, idx = remove_vertices(adj, {leaf, support})
            u_in_b = idx[u]

            p_poly, q_poly = rooted_dp_u(b_adj, u_in_b)
            b_poly = independence_poly(len(b_adj), b_adj)

            mu_p = mean_at_lambda(p_poly, lam)
            mu_b = mean_at_lambda(b_poly, lam)
            d_val = mu_b - mu_p

            a = lam / (1.0 + lam)
            half_excess = d_val - 0.5
            exact_excess = d_val - (1.0 - a)

            z_p = eval_poly(p_poly, lam)
            z_q = eval_poly(q_poly, lam)
            p_u = z_q / (z_p + z_q) if (z_p + z_q) else 0.0

            g6 = raw.decode("ascii").strip()
            witness = {
                "n": nn,
                "g6": g6,
                "mode": m,
                "lambda_mode": lam,
                "leaf": leaf,
                "support": support,
                "mu_p": mu_p,
                "mu_b": mu_b,
                "D": d_val,
                "a": a,
                "half_excess": half_excess,
                "exact_excess": exact_excess,
                "p_u": p_u,
            }

            if half_excess > args.tol:
                stats["half_violation"] += 1
            if exact_excess > args.tol:
                stats["exact_violation"] += 1

            maybe_update_max(stats, "max_D", "max_D_witness", d_val, witness)
            maybe_update_max(stats, "max_half_excess", "max_half_excess_witness", half_excess, witness)
            maybe_update_max(stats, "max_exact_excess", "max_exact_excess_witness", exact_excess, witness)
            maybe_update_max(stats, "max_p_u", "max_p_u_witness", p_u, witness)

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"checked={stats['checked']:8d} half_v={stats['half_violation']:5d} "
            f"exact_v={stats['exact_violation']:5d} maxD={fmt_opt(stats['max_D'])} "
            f"maxHalf={fmt_opt(stats['max_half_excess'])} "
            f"maxExact={fmt_opt(stats['max_exact_excess'])} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    wall = time.time() - t_all
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"checked={summary['checked']:,} half_v={summary['half_violation']:,} "
        f"exact_v={summary['exact_violation']:,} wall={wall:.1f}s",
        flush=True,
    )
    print(f"max_D={summary['max_D']}", flush=True)
    print(f"max_half_excess={summary['max_half_excess']}", flush=True)
    print(f"max_exact_excess={summary['max_exact_excess']}", flush=True)
    print(f"max_p_u={summary['max_p_u']}", flush=True)
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
