#!/usr/bin/env python3
"""Verify slope-based closure criteria for Route-2 compensation.

For each eligible bridge leaf (l,s) with deg(s)=2 in a tree T:
  B = T - {l,s}
  m = mode(I(T)) (leftmost)
  lam = lambda_m(T) = i_{m-1}(T) / i_m(T)
  tau = b_{m-2} / b_{m-1} from I(B)=sum b_k x^k

Route-2 target:
  mu_B(lam) >= m - 3/2

This checker records:
  deficit_tau = (m-3/2) - mu_B(tau)
  gain        = mu_B(lam) - mu_B(tau)
  gap         = lam - tau

and slope-based margins (deficit cases only):
  M_tau = (Var_B(tau)/tau) * gap - deficit_tau
  M_lam = (Var_B(lam)/lam) * gap - deficit_tau

If mu_B is concave on [tau,lam], then gain >= (Var_B(lam)/lam)*gap,
so M_lam >= 0 is a sufficient closure condition.
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


def mean_and_var_at_lambda(poly: list[int], lam: float) -> tuple[float, float]:
    z = 0.0
    mu_num = 0.0
    mu2_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        mu2_num += k * k * w
        p *= lam
    if z == 0.0:
        return 0.0, 0.0
    mu = mu_num / z
    var = (mu2_num / z) - (mu * mu)
    return mu, var


def concavity_profile(
    poly: list[int],
    left: float,
    right: float,
    grid: int,
) -> tuple[bool, float]:
    """Check discrete concavity of mu_B over [left,right] on a uniform grid.

    Returns:
      (is_concave, max_second_diff)
    where max_second_diff > 0 indicates a concavity violation.
    """
    if grid < 2 or right <= left:
        return True, 0.0

    mus = []
    for i in range(grid + 1):
        t = i / grid
        lam = left + (right - left) * t
        mu, _ = mean_and_var_at_lambda(poly, lam)
        mus.append(mu)

    max_d2 = 0.0
    for i in range(1, grid):
        d2 = mus[i + 1] - 2.0 * mus[i] + mus[i - 1]
        if d2 > max_d2:
            max_d2 = d2
    return max_d2 <= 0.0, max_d2


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
        "route2_fail": 0,
        "exact_fail": 0,
        "deficit_tau_cases": 0,
        "deficit_nonpos_gap": 0,
        "min_route2_slack": None,
        "min_route2_witness": None,
        "min_exact_slack": None,
        "min_exact_witness": None,
        "min_gain_minus_deficit": None,
        "min_gain_minus_deficit_witness": None,
        "min_tau_endpoint_margin": None,
        "min_tau_endpoint_margin_witness": None,
        "min_lam_endpoint_margin": None,
        "min_lam_endpoint_margin_witness": None,
        "min_gap_deficit": None,
        "min_gap_deficit_witness": None,
        "max_deficit_tau": None,
        "max_deficit_tau_witness": None,
        "concavity_checked": 0,
        "concavity_fail": 0,
        "max_concavity_d2": None,
        "max_concavity_d2_witness": None,
        "wall_s": 0.0,
    }


def merge_stats(dst: dict[str, Any], src: dict[str, Any]) -> None:
    for k in [
        "seen",
        "considered",
        "checked_trees",
        "checked_leaves",
        "route2_fail",
        "exact_fail",
        "deficit_tau_cases",
        "deficit_nonpos_gap",
        "concavity_checked",
        "concavity_fail",
    ]:
        dst[k] += int(src.get(k, 0))

    for key, witness_key in [
        ("min_route2_slack", "min_route2_witness"),
        ("min_exact_slack", "min_exact_witness"),
        ("min_gain_minus_deficit", "min_gain_minus_deficit_witness"),
        ("min_tau_endpoint_margin", "min_tau_endpoint_margin_witness"),
        ("min_lam_endpoint_margin", "min_lam_endpoint_margin_witness"),
        ("min_gap_deficit", "min_gap_deficit_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])

    for key, witness_key in [
        ("max_deficit_tau", "max_deficit_tau_witness"),
        ("max_concavity_d2", "max_concavity_d2_witness"),
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


def fmt_opt(val: float | None) -> str:
    return "None" if val is None else f"{val:.6g}"


def main() -> None:
    ap = argparse.ArgumentParser(description="Route-2 slope-compensation verifier (checkpointed).")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument("--all-deg2-leaves", action="store_true")
    ap.add_argument(
        "--concavity-grid",
        type=int,
        default=0,
        help="If >0, check discrete concavity of mu_B on [tau,lambda] using this grid size.",
    )
    ap.add_argument(
        "--concavity-all",
        action="store_true",
        help="Check concavity for all checked leaves (default: only deficit+positive-gap cases).",
    )
    ap.add_argument("--mod", type=int, default=0, help="Process only trees with seen %% mod == res.")
    ap.add_argument("--res", type=int, default=0, help="Residue for --mod split.")
    ap.add_argument("--out", default="results/route2_slope_compensation_check.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    if args.mod < 0:
        raise ValueError("--mod must be >= 0")
    if args.mod == 0 and args.res != 0:
        raise ValueError("--res requires --mod > 0")
    if args.mod > 0 and not (0 <= args.res < args.mod):
        raise ValueError("--res must satisfy 0 <= res < mod")

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "tol": args.tol,
        "all_trees": args.all_trees,
        "all_deg2_leaves": args.all_deg2_leaves,
        "concavity_grid": args.concavity_grid,
        "concavity_all": args.concavity_all,
        "mod": args.mod,
        "res": args.res,
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
    split = f"split res={args.res}/{args.mod}" if args.mod > 0 else "no split"
    print(
        f"Route-2 slope-compensation check on {scope}, {leaf_mode}, n={args.min_n}..{args.max_n} ({split})",
        flush=True,
    )
    print(f"Output: {args.out}", flush=True)
    print("-" * 108, flush=True)

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
            if args.mod > 0 and (stats["seen"] - 1) % args.mod != args.res:
                continue

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

            poly_t = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_t)
            if m < 2 or poly_t[m - 1] <= 0 or poly_t[m] <= 0:
                continue

            lam = poly_t[m - 1] / poly_t[m]
            stats["checked_trees"] += 1
            g6 = line.decode("ascii").strip()
            deg_sig = dict(sorted(Counter(deg).items()))

            for leaf in leaves:
                support = adj[leaf][0]
                b_adj = remove_vertices(adj, {leaf, support})
                poly_b = independence_poly(len(b_adj), b_adj)
                if m - 1 >= len(poly_b) or poly_b[m - 2] <= 0 or poly_b[m - 1] <= 0:
                    continue

                tau = poly_b[m - 2] / poly_b[m - 1]
                mu_tau, var_tau = mean_and_var_at_lambda(poly_b, tau)
                mu_lam, var_lam = mean_and_var_at_lambda(poly_b, lam)

                threshold_half = m - 1.5
                threshold_exact = m - 1 - (lam / (1 + lam))
                route2_slack = mu_lam - threshold_half
                exact_slack = mu_lam - threshold_exact

                gap = lam - tau
                deficit_tau = threshold_half - mu_tau
                gain = mu_lam - mu_tau
                gain_minus_deficit = gain - max(deficit_tau, 0.0)

                stats["checked_leaves"] += 1
                if route2_slack < -args.tol:
                    stats["route2_fail"] += 1
                if exact_slack < -args.tol:
                    stats["exact_fail"] += 1

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

                if deficit_tau > args.tol:
                    stats["deficit_tau_cases"] += 1
                    if gap <= args.tol:
                        stats["deficit_nonpos_gap"] += 1
                    else:
                        slope_tau = var_tau / tau if tau > 0 else 0.0
                        slope_lam = var_lam / lam if lam > 0 else 0.0
                        tau_margin = slope_tau * gap - deficit_tau
                        lam_margin = slope_lam * gap - deficit_tau

                        def_witness = {
                            **base_witness,
                            "deficit_tau": deficit_tau,
                            "gain": gain,
                            "gap": gap,
                            "var_tau": var_tau,
                            "var_lam": var_lam,
                            "tau_margin": tau_margin,
                            "lam_margin": lam_margin,
                        }

                        maybe_update_min(
                            stats,
                            "min_gain_minus_deficit",
                            "min_gain_minus_deficit_witness",
                            gain_minus_deficit,
                            {**def_witness, "gain_minus_deficit": gain_minus_deficit},
                        )
                        maybe_update_min(
                            stats,
                            "min_tau_endpoint_margin",
                            "min_tau_endpoint_margin_witness",
                            tau_margin,
                            def_witness,
                        )
                        maybe_update_min(
                            stats,
                            "min_lam_endpoint_margin",
                            "min_lam_endpoint_margin_witness",
                            lam_margin,
                            def_witness,
                        )
                        maybe_update_min(
                            stats,
                            "min_gap_deficit",
                            "min_gap_deficit_witness",
                            gap,
                            def_witness,
                        )
                        maybe_update_max(
                            stats,
                            "max_deficit_tau",
                            "max_deficit_tau_witness",
                            deficit_tau,
                            def_witness,
                        )

                do_concavity = args.concavity_grid > 0 and (
                    args.concavity_all or (deficit_tau > args.tol and gap > args.tol)
                )
                if do_concavity:
                    ok_conc, max_d2 = concavity_profile(poly_b, tau, lam, args.concavity_grid)
                    stats["concavity_checked"] += 1
                    witness = {
                        **base_witness,
                        "deficit_tau": deficit_tau,
                        "gap": gap,
                        "max_d2": max_d2,
                    }
                    maybe_update_max(
                        stats,
                        "max_concavity_d2",
                        "max_concavity_d2_witness",
                        max_d2,
                        witness,
                    )
                    if (not ok_conc) and max_d2 > args.tol:
                        stats["concavity_fail"] += 1

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"trees={stats['checked_trees']:8d} leaves={stats['checked_leaves']:8d} "
            f"route2_fail={stats['route2_fail']:4d} exact_fail={stats['exact_fail']:4d} "
            f"def_tau={stats['deficit_tau_cases']:8d} nonpos_gap={stats['deficit_nonpos_gap']:4d} "
            f"min_r2={fmt_opt(stats['min_route2_slack'])} "
            f"min_lam_margin={fmt_opt(stats['min_lam_endpoint_margin'])} "
            f"conc_fail={stats['concavity_fail']:4d} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    print("-" * 108, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"trees={summary['checked_trees']:,} leaves={summary['checked_leaves']:,} "
        f"route2_fail={summary['route2_fail']:,} exact_fail={summary['exact_fail']:,} "
        f"def_tau={summary['deficit_tau_cases']:,} nonpos_gap={summary['deficit_nonpos_gap']:,} "
        f"conc_checked={summary['concavity_checked']:,} conc_fail={summary['concavity_fail']:,} "
        f"wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(
        "Extrema: "
        f"min_r2={summary['min_route2_slack']}, "
        f"min_exact={summary['min_exact_slack']}, "
        f"min_gain_minus_def={summary['min_gain_minus_deficit']}, "
        f"min_tau_margin={summary['min_tau_endpoint_margin']}, "
        f"min_lam_margin={summary['min_lam_endpoint_margin']}, "
        f"max_def_tau={summary['max_deficit_tau']}, "
        f"max_d2={summary['max_concavity_d2']}",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
