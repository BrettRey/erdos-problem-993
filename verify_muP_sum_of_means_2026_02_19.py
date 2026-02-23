#!/usr/bin/env python3
"""Unified verifier for Path-1 direct mu_P inequality via sum-of-means.

For each d_leaf<=1 tree T (unless --all-trees):
  - choose canonical leaf l by min support degree (tie: largest leaf id),
  - require deg(s)=2 for support s of l,
  - let u be the other neighbor of s,
  - set B = T - {l, s}, root B at u,
  - define P = dp_B[u][0], and child subtrees T_c for c child of u in B,
  - m = mode(I(T)), lambda = i_{m-1}(T) / i_m(T), a = lambda/(1+lambda).

The script checks, per tree:
  1) sum-of-means identity:
       mu_P(lambda) = sum_c mu_{T_c}(lambda)
  2) exact chain identity:
       mu_P - (m-2)
       = [mu_B - (m-1-a)] - [D - (1-a)],  where D := mu_B - mu_P
  3) target inequality:
       mu_P(lambda) >= m - 2

It records extremal slack values and witness trees, and supports per-n checkpoints.
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


def rooted_dp(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[list[int]], list[list[int]], list[int]]:
    """Return (dp0, dp1, children, order_post) for tree rooted at root."""
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

    return dp0, dp1, children, order


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
        "muP_fail": 0,
        "sum_identity_fail": 0,
        "chain_identity_fail": 0,
        "partition_size_fail": 0,
        "transfer_cap_violation": 0,
        "checked_k1": 0,
        "checked_k2plus": 0,
        "min_muP_gap": None,
        "min_muP_gap_witness": None,
        "min_muP_gap_k1": None,
        "min_muP_gap_k1_witness": None,
        "min_muP_gap_k2plus": None,
        "min_muP_gap_k2plus_witness": None,
        "min_exact_slack_B": None,
        "min_exact_slack_B_witness": None,
        "max_exact_excess_D": None,
        "max_exact_excess_D_witness": None,
        "min_chain_gap": None,
        "min_chain_gap_witness": None,
        "max_abs_sum_identity_err": None,
        "max_abs_sum_identity_err_witness": None,
        "max_abs_chain_identity_err": None,
        "max_abs_chain_identity_err_witness": None,
        "wall_s": 0.0,
    }


def maybe_update_min(
    stats: dict[str, Any],
    key: str,
    witness_key: str,
    value: float,
    witness: dict[str, Any],
) -> None:
    cur = stats[key]
    if cur is None or value < cur:
        stats[key] = value
        stats[witness_key] = witness


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
    for k in [
        "seen",
        "considered",
        "checked",
        "no_deg2_support",
        "muP_fail",
        "sum_identity_fail",
        "chain_identity_fail",
        "partition_size_fail",
        "transfer_cap_violation",
        "checked_k1",
        "checked_k2plus",
    ]:
        dst[k] += int(src.get(k, 0))

    for key, witness_key in [
        ("min_muP_gap", "min_muP_gap_witness"),
        ("min_muP_gap_k1", "min_muP_gap_k1_witness"),
        ("min_muP_gap_k2plus", "min_muP_gap_k2plus_witness"),
        ("min_exact_slack_B", "min_exact_slack_B_witness"),
        ("min_chain_gap", "min_chain_gap_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])

    for key, witness_key in [
        ("max_exact_excess_D", "max_exact_excess_D_witness"),
        ("max_abs_sum_identity_err", "max_abs_sum_identity_err_witness"),
        ("max_abs_chain_identity_err", "max_abs_chain_identity_err_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_max(dst, key, witness_key, float(src[key]), src[witness_key])


def aggregate(per_n: dict[str, dict[str, Any]]) -> dict[str, Any]:
    out = fresh_stats()
    for n_key in sorted(per_n, key=lambda s: int(s)):
        merge_stats(out, per_n[n_key])
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
    ap = argparse.ArgumentParser(description="Verify direct mu_P bound via sum-of-means and transfer chain.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--identity-tol", type=float, default=1e-9)
    ap.add_argument("--transfer-cap", type=float, default=0.006)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument("--res", type=int, default=None, help="Optional split residue r in r/m.")
    ap.add_argument("--mod", type=int, default=None, help="Optional split modulus m in r/m.")
    ap.add_argument("--out", default="results/verify_muP_sum_of_means_n23.json")
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
        "identity_tol": args.identity_tol,
        "transfer_cap": args.transfer_cap,
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
        f"verify_muP_sum_of_means on {scope}, n={args.min_n}..{args.max_n}{split_desc}",
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
            if m == 0 or m - 1 >= len(poly_t) or poly_t[m - 1] <= 0 or poly_t[m] <= 0:
                continue
            lam = poly_t[m - 1] / poly_t[m]
            a = lam / (1.0 + lam)

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                stats["no_deg2_support"] += 1
                continue

            u = adj[support][0] if adj[support][1] == leaf else adj[support][1]
            b_adj, idx = remove_vertices(adj, {leaf, support})
            u_in_b = idx[u]

            dp0, dp1, children, order = rooted_dp(b_adj, u_in_b)
            p_poly = dp0[u_in_b]
            q_poly = dp1[u_in_b]
            b_poly = _polyadd(p_poly, q_poly)

            mu_p = mean_at_lambda(p_poly, lam)
            mu_b = mean_at_lambda(b_poly, lam)
            d_val = mu_b - mu_p

            exact_slack_b = mu_b - (m - 1 - a)
            exact_excess_d = d_val - (1 - a)
            mu_p_gap = mu_p - (m - 2)
            chain_gap = exact_slack_b - exact_excess_d
            chain_err = mu_p_gap - chain_gap

            # sum_c mu_{T_c}(lambda)
            sum_mu_children = 0.0
            for c in children[u_in_b]:
                i_c = _polyadd(dp0[c], dp1[c])
                sum_mu_children += mean_at_lambda(i_c, lam)
            sum_err = mu_p - sum_mu_children

            # partition-size identity: sum_c |T_c| = |B| - 1 = n - 3
            subtree_size = [0] * len(b_adj)
            for v in order:
                sz = 1
                for c in children[v]:
                    sz += subtree_size[c]
                subtree_size[v] = sz
            part_size = sum(subtree_size[c] for c in children[u_in_b])
            size_ok = (part_size == len(b_adj) - 1 == nn - 3)

            k_children = len(children[u_in_b])

            stats["checked"] += 1
            if k_children == 1:
                stats["checked_k1"] += 1
            else:
                stats["checked_k2plus"] += 1

            if mu_p_gap < -args.tol:
                stats["muP_fail"] += 1
            if abs(sum_err) > args.identity_tol:
                stats["sum_identity_fail"] += 1
            if abs(chain_err) > args.identity_tol:
                stats["chain_identity_fail"] += 1
            if not size_ok:
                stats["partition_size_fail"] += 1
            if exact_excess_d > args.transfer_cap + args.tol:
                stats["transfer_cap_violation"] += 1

            g6 = raw.decode("ascii").strip()
            base = {
                "n": nn,
                "g6": g6,
                "mode": m,
                "lambda_mode": lam,
                "leaf": leaf,
                "support": support,
                "u_in_B": u_in_b,
                "k_children": k_children,
                "mu_p": mu_p,
                "mu_b": mu_b,
                "D": d_val,
                "a": a,
                "mu_p_gap": mu_p_gap,
                "exact_slack_B": exact_slack_b,
                "exact_excess_D": exact_excess_d,
                "chain_gap": chain_gap,
                "sum_mu_children": sum_mu_children,
                "sum_identity_err": sum_err,
                "chain_identity_err": chain_err,
                "partition_size_sum": part_size,
                "partition_size_target": nn - 3,
            }

            maybe_update_min(stats, "min_muP_gap", "min_muP_gap_witness", mu_p_gap, base)
            if k_children == 1:
                maybe_update_min(stats, "min_muP_gap_k1", "min_muP_gap_k1_witness", mu_p_gap, base)
            else:
                maybe_update_min(
                    stats,
                    "min_muP_gap_k2plus",
                    "min_muP_gap_k2plus_witness",
                    mu_p_gap,
                    base,
                )

            maybe_update_min(
                stats,
                "min_exact_slack_B",
                "min_exact_slack_B_witness",
                exact_slack_b,
                base,
            )
            maybe_update_max(
                stats,
                "max_exact_excess_D",
                "max_exact_excess_D_witness",
                exact_excess_d,
                base,
            )
            maybe_update_min(stats, "min_chain_gap", "min_chain_gap_witness", chain_gap, base)
            maybe_update_max(
                stats,
                "max_abs_sum_identity_err",
                "max_abs_sum_identity_err_witness",
                abs(sum_err),
                base,
            )
            maybe_update_max(
                stats,
                "max_abs_chain_identity_err",
                "max_abs_chain_identity_err_witness",
                abs(chain_err),
                base,
            )

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"checked={stats['checked']:8d} muP_fail={stats['muP_fail']:4d} "
            f"sum_fail={stats['sum_identity_fail']:4d} chain_fail={stats['chain_identity_fail']:4d} "
            f"min_gap={fmt_opt(stats['min_muP_gap'])} "
            f"min_exact={fmt_opt(stats['min_exact_slack_B'])} "
            f"max_excess={fmt_opt(stats['max_exact_excess_D'])} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    wall = time.time() - t_all
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} checked={summary['checked']:,} "
        f"muP_fail={summary['muP_fail']:,} sum_fail={summary['sum_identity_fail']:,} "
        f"chain_fail={summary['chain_identity_fail']:,} part_fail={summary['partition_size_fail']:,} "
        f"cap_v={summary['transfer_cap_violation']:,} wall={wall:.1f}s",
        flush=True,
    )
    print(
        "Extrema: "
        f"min_muP_gap={summary['min_muP_gap']}, "
        f"min_exact_slack_B={summary['min_exact_slack_B']}, "
        f"max_exact_excess_D={summary['max_exact_excess_D']}, "
        f"min_chain_gap={summary['min_chain_gap']}",
        flush=True,
    )
    print(
        "Identity errors: "
        f"max_abs_sum_identity_err={summary['max_abs_sum_identity_err']}, "
        f"max_abs_chain_identity_err={summary['max_abs_chain_identity_err']}",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
