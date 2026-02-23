#!/usr/bin/env python3
"""Scan route-2 compensation profile at the degree-2 leaf bridge.

For each eligible tree T (default: d_leaf<=1), let m be the leftmost mode index
of I(T) and lambda_m = i_{m-1}/i_m. Choose either:
  - one canonical leaf l whose support s has degree 2 (default), or
  - every such leaf (with --all-deg2-leaves).

Set B = T-{l,s}. Define tau = b_{m-2}/b_{m-1} from I(B)=sum b_k x^k.

At each checked leaf, this scanner records:
  - route-2 target slack at lambda_m:
      mu_B(lambda_m) - (m - 3/2)
  - stronger pendant-bonus slack:
      mu_B(lambda_m) - (m - 1 - lambda_m/(1+lambda_m))
  - tau-deficit:
      (m - 3/2) - mu_B(tau)
  - compensation lift from tau to lambda_m:
      gain = mu_B(lambda_m) - mu_B(tau), gap = lambda_m - tau,
      avg_slope = gain/gap (when tau-deficit > 0 and gap > 0).

Output is checkpointed by n so long runs can resume.
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


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked_trees": 0,
        "checked_leaves": 0,
        "fail_half": 0,
        "fail_exact": 0,
        "deficient_tau": 0,
        "deficient_gap_nonpos": 0,
        "min_route2_slack": None,
        "min_route2_witness": None,
        "min_exact_slack": None,
        "min_exact_witness": None,
        "max_deficit_tau": None,
        "max_deficit_tau_witness": None,
        "min_gap_def": None,
        "min_gap_def_witness": None,
        "min_gain_def": None,
        "min_gain_def_witness": None,
        "min_avg_slope_def": None,
        "min_avg_slope_def_witness": None,
        "min_tie_gap": None,
        "min_tie_gap_witness": None,
        "wall_s": 0.0,
    }


def merge_stats(dst: dict[str, Any], src: dict[str, Any]) -> None:
    for k in [
        "seen",
        "considered",
        "checked_trees",
        "checked_leaves",
        "fail_half",
        "fail_exact",
        "deficient_tau",
        "deficient_gap_nonpos",
    ]:
        dst[k] += int(src.get(k, 0))
    for key, witness_key in [
        ("min_route2_slack", "min_route2_witness"),
        ("min_exact_slack", "min_exact_witness"),
        ("min_gap_def", "min_gap_def_witness"),
        ("min_gain_def", "min_gain_def_witness"),
        ("min_avg_slope_def", "min_avg_slope_def_witness"),
        ("min_tie_gap", "min_tie_gap_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])
    for key, witness_key in [("max_deficit_tau", "max_deficit_tau_witness")]:
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


def fmt_opt(val: float | None) -> str:
    return "None" if val is None else f"{val:.6g}"


def main() -> None:
    ap = argparse.ArgumentParser(description="Route-2 tau-deficit compensation scanner (checkpointed).")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument(
        "--all-deg2-leaves",
        action="store_true",
        help="Check every degree-2-support leaf instead of one canonical leaf per tree.",
    )
    ap.add_argument("--out", default="results/whnc_route2_compensation_scan_n23_canonical.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "tol": args.tol,
        "all_trees": args.all_trees,
        "all_deg2_leaves": args.all_deg2_leaves,
    }

    per_n: dict[str, dict[str, Any]] = {}
    if (not args.no_resume) and args.out and os.path.exists(args.out):
        with open(args.out, "r", encoding="utf-8") as f:
            old = json.load(f)
        per_n = old.get("per_n", {})
        done = ", ".join(sorted(per_n, key=lambda s: int(s)))
        print(f"Resuming from {args.out}; completed n: {done}", flush=True)

    scope = "all trees" if args.all_trees else "d_leaf<=1"
    leaf_mode = "all deg2-support leaves" if args.all_deg2_leaves else "canonical leaf"
    print(
        f"Route-2 compensation scan on {scope}, {leaf_mode}, n={args.min_n}..{args.max_n}",
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

            deg = [len(nb) for nb in adj]
            if args.all_deg2_leaves:
                leaves = []
                for v, d in enumerate(deg):
                    if d != 1:
                        continue
                    s = adj[v][0]
                    if deg[s] == 2:
                        leaves.append(v)
                if not leaves:
                    continue
            else:
                leaf = choose_canonical_deg2_leaf(adj)
                if leaf is None:
                    continue
                leaves = [leaf]

            poly = independence_poly(nn, adj)
            m = mode_index_leftmost(poly)
            if m < 2 or poly[m - 1] <= 0 or poly[m] <= 0:
                continue

            lam = poly[m - 1] / poly[m]
            stats["checked_trees"] += 1
            g6 = line.decode("ascii").strip()
            deg_sig = dict(sorted(Counter(deg).items()))

            for leaf in leaves:
                support = adj[leaf][0]
                b_adj = remove_vertices(adj, {leaf, support})
                b_poly = independence_poly(len(b_adj), b_adj)
                if m - 1 >= len(b_poly) or b_poly[m - 2] <= 0 or b_poly[m - 1] <= 0:
                    continue

                tau = b_poly[m - 2] / b_poly[m - 1]
                mu_tau = mean_at_lambda(b_poly, tau)
                mu_lam = mean_at_lambda(b_poly, lam)

                half_threshold = m - 1.5
                exact_threshold = m - 1 - (lam / (1 + lam))
                route2_slack = mu_lam - half_threshold
                exact_slack = mu_lam - exact_threshold

                tie_gap = lam - tau
                deficit_tau = half_threshold - mu_tau
                gain = mu_lam - mu_tau

                stats["checked_leaves"] += 1
                if route2_slack < -args.tol:
                    stats["fail_half"] += 1
                if exact_slack < -args.tol:
                    stats["fail_exact"] += 1

                base_witness = {
                    "n": nn,
                    "g6": g6,
                    "mode": m,
                    "lambda_mode": lam,
                    "tau_b": tau,
                    "leaf": leaf,
                    "support": support,
                    "mu_tau": mu_tau,
                    "mu_b": mu_lam,
                    "degree_signature": deg_sig,
                }
                maybe_update_min(
                    stats,
                    "min_route2_slack",
                    "min_route2_witness",
                    route2_slack,
                    {**base_witness, "route2_slack": route2_slack},
                )
                maybe_update_min(
                    stats,
                    "min_exact_slack",
                    "min_exact_witness",
                    exact_slack,
                    {**base_witness, "exact_slack": exact_slack},
                )
                maybe_update_min(
                    stats,
                    "min_tie_gap",
                    "min_tie_gap_witness",
                    tie_gap,
                    {**base_witness, "tie_gap": tie_gap},
                )

                if deficit_tau > args.tol:
                    stats["deficient_tau"] += 1
                    def_witness = {
                        **base_witness,
                        "deficit_tau": deficit_tau,
                        "tie_gap": tie_gap,
                        "gain": gain,
                    }
                    maybe_update_max(
                        stats,
                        "max_deficit_tau",
                        "max_deficit_tau_witness",
                        deficit_tau,
                        def_witness,
                    )
                    if tie_gap <= args.tol:
                        stats["deficient_gap_nonpos"] += 1
                    else:
                        avg_slope = gain / tie_gap
                        maybe_update_min(
                            stats,
                            "min_gap_def",
                            "min_gap_def_witness",
                            tie_gap,
                            {**def_witness, "avg_slope": avg_slope},
                        )
                        maybe_update_min(
                            stats,
                            "min_gain_def",
                            "min_gain_def_witness",
                            gain,
                            {**def_witness, "avg_slope": avg_slope},
                        )
                        maybe_update_min(
                            stats,
                            "min_avg_slope_def",
                            "min_avg_slope_def_witness",
                            avg_slope,
                            {**def_witness, "avg_slope": avg_slope},
                        )

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"trees={stats['checked_trees']:8d} leaves={stats['checked_leaves']:8d} "
            f"fail_half={stats['fail_half']:4d} fail_exact={stats['fail_exact']:4d} "
            f"def_tau={stats['deficient_tau']:8d} def_gap_nonpos={stats['deficient_gap_nonpos']:4d} "
            f"min_half={fmt_opt(stats['min_route2_slack'])} "
            f"max_def_tau={fmt_opt(stats['max_deficit_tau'])} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"trees={summary['checked_trees']:,} leaves={summary['checked_leaves']:,} "
        f"fail_half={summary['fail_half']:,} fail_exact={summary['fail_exact']:,} "
        f"def_tau={summary['deficient_tau']:,} def_gap_nonpos={summary['deficient_gap_nonpos']:,} "
        f"wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(
        "Extrema: "
        f"min_half={summary['min_route2_slack']}, "
        f"min_exact={summary['min_exact_slack']}, "
        f"max_def_tau={summary['max_deficit_tau']}, "
        f"min_gap_def={summary['min_gap_def']}, "
        f"min_gain_def={summary['min_gain_def']}, "
        f"min_avg_slope_def={summary['min_avg_slope_def']}",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
