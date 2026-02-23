#!/usr/bin/env python3
"""Scan pendant-bonus mean bounds at the degree-2 leaf bridge.

For each eligible tree T (default: d_leaf<=1), with m = leftmost mode of I(T)
and lambda_m = i_{m-1}/i_m, choose either:
  - every leaf l whose support s has degree 2, or
  - one canonical such leaf per tree (smallest leaf index among deg2 supports).

Define B = T - {l, s}.  This scanner checks two bounds:

1) Route-2 target (weaker):
     mu_B(lambda_m) >= m - 3/2

2) Pendant-bonus sufficient bound (stronger):
     mu_B(lambda_m) >= m - 1 - lambda_m/(1+lambda_m)

It writes checkpointed per-n JSON so long runs can be resumed.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def remove_vertices(adj: list[list[int]], remove_set: set[int]) -> list[list[int]]:
    keep = [v for v in range(len(adj)) if v not in remove_set]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
    return out


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
        "checked_trees": 0,
        "checked_leaves": 0,
        "fail_half": 0,
        "fail_exact": 0,
        "min_half_slack": None,
        "min_half_witness": None,
        "min_exact_slack": None,
        "min_exact_witness": None,
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


def merge_stats(dst: dict[str, Any], src: dict[str, Any]) -> None:
    for k in ["seen", "considered", "checked_trees", "checked_leaves", "fail_half", "fail_exact"]:
        dst[k] += int(src.get(k, 0))
    for key, witness_key in [
        ("min_half_slack", "min_half_witness"),
        ("min_exact_slack", "min_exact_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])


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


def choose_canonical_deg2_leaf(adj: list[list[int]]) -> int | None:
    deg = [len(nb) for nb in adj]
    cand = []
    for v, d in enumerate(deg):
        if d != 1:
            continue
        s = adj[v][0]
        if deg[s] == 2:
            cand.append(v)
    if not cand:
        return None
    return min(cand)


def fmt_opt(val: float | None) -> str:
    return "None" if val is None else f"{val:.6g}"


def main() -> None:
    ap = argparse.ArgumentParser(description="Pendant-bonus route-2 mean scan (checkpointed).")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument(
        "--canonical-leaf-only",
        action="store_true",
        help="Use one canonical degree-2-support leaf per tree instead of all such leaves.",
    )
    ap.add_argument("--out", default="results/whnc_pendant_bonus_scan_n23.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "tol": args.tol,
        "all_trees": args.all_trees,
        "canonical_leaf_only": args.canonical_leaf_only,
    }

    per_n: dict[str, dict[str, Any]] = {}
    if (not args.no_resume) and args.out and os.path.exists(args.out):
        with open(args.out, "r", encoding="utf-8") as f:
            old = json.load(f)
        per_n = old.get("per_n", {})
        done = ", ".join(sorted(per_n, key=lambda s: int(s)))
        print(f"Resuming from {args.out}; completed n: {done}", flush=True)

    scope = "all trees" if args.all_trees else "d_leaf<=1"
    leaf_mode = "canonical leaf" if args.canonical_leaf_only else "all deg2-support leaves"
    print(
        f"Pendant-bonus scan on {scope}, {leaf_mode}, n={args.min_n}..{args.max_n}",
        flush=True,
    )
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
        proc = subprocess.Popen(
            [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            stats["seen"] += 1

            if (not args.all_trees) and (not is_dleaf_le_1(nn, adj)):
                continue
            stats["considered"] += 1

            poly = independence_poly(nn, adj)
            m = mode_index_leftmost(poly)
            if m == 0 or poly[m - 1] <= 0 or poly[m] <= 0:
                continue
            lam = poly[m - 1] / poly[m]
            stats["checked_trees"] += 1

            deg = [len(nb) for nb in adj]
            if args.canonical_leaf_only:
                leaf = choose_canonical_deg2_leaf(adj)
                if leaf is None:
                    continue
                leaves = [leaf]
            else:
                leaves = []
                for v, d in enumerate(deg):
                    if d != 1:
                        continue
                    s = adj[v][0]
                    if deg[s] == 2:
                        leaves.append(v)
                if not leaves:
                    continue

            g6 = line.decode("ascii").strip()
            deg_sig = dict(sorted(Counter(deg).items()))
            for leaf in leaves:
                support = adj[leaf][0]
                b_adj = remove_vertices(adj, {leaf, support})
                b_poly = independence_poly(len(b_adj), b_adj)
                mu_b = mean_at_lambda(b_poly, lam)

                half_slack = mu_b - (m - 1.5)
                exact_threshold = m - 1 - (lam / (1 + lam))
                exact_slack = mu_b - exact_threshold

                stats["checked_leaves"] += 1
                if half_slack < -args.tol:
                    stats["fail_half"] += 1
                if exact_slack < -args.tol:
                    stats["fail_exact"] += 1

                base_witness = {
                    "n": nn,
                    "g6": g6,
                    "mode": m,
                    "lambda_mode": lam,
                    "leaf": leaf,
                    "support": support,
                    "mu_b": mu_b,
                    "degree_signature": deg_sig,
                }
                maybe_update_min(
                    stats,
                    "min_half_slack",
                    "min_half_witness",
                    half_slack,
                    {**base_witness, "half_slack": half_slack},
                )
                maybe_update_min(
                    stats,
                    "min_exact_slack",
                    "min_exact_witness",
                    exact_slack,
                    {**base_witness, "exact_slack": exact_slack},
                )

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"trees={stats['checked_trees']:8d} leaves={stats['checked_leaves']:9d} "
            f"fail_half={stats['fail_half']:4d} fail_exact={stats['fail_exact']:4d} "
            f"min_half={fmt_opt(stats['min_half_slack'])} "
            f"min_exact={fmt_opt(stats['min_exact_slack'])} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"trees={summary['checked_trees']:,} leaves={summary['checked_leaves']:,} "
        f"fail_half={summary['fail_half']:,} fail_exact={summary['fail_exact']:,} "
        f"wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(
        f"Global minima: half={summary['min_half_slack']}, exact={summary['min_exact_slack']}",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
