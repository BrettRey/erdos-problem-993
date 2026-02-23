#!/usr/bin/env python3
"""Scan Route-1 local invariant: Phi_{m-1}(P; lambda_m(T)) >= 0.

Setup (degree-2 bridge):
  - choose a leaf l whose support s has minimum degree
  - require deg(s)=2 (observed always in d_leaf<=1 trees)
  - let u be the other neighbor of s
  - B = T - {l, s}
  - P = dp_B[u][0]  (u excluded)

At m = mode(I(T)) and lambda_m = i_{m-1}(T)/i_m(T), this script tracks:
  1) Phi_{m-1}(P; lambda_m)
  2) lambda gap: lambda_m - p_{m-2}/p_{m-1}
  3) ULC status of P
  4) mode(P) - (m-1) and degree(P) - (m-1)
  5) tie-mean gap at tau_P = p_{m-2}/p_{m-1}

Checkpointed per n with resume support.
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
from indpoly import _polyadd, _polymul, independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


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


def rooted_dp_u_excluded(adj: list[list[int]], root: int) -> list[int]:
    """Return dp[root][0] for a tree rooted at root."""
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

    return dp0[root]


def phi_q(poly: list[int], q: int, lam: float) -> float:
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    mu = mu_num / z
    return z * (mu - (q - 1))


def mean_at_lambda(poly: list[int], lam: float) -> float:
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z


def is_ulc(poly: list[int]) -> bool:
    # k * a_k^2 >= (k+1) * a_{k-1} * a_{k+1}, for k=1..d-1
    for k in range(1, len(poly) - 1):
        if k * poly[k] * poly[k] < (k + 1) * poly[k - 1] * poly[k + 1]:
            return False
    return True


def choose_min_support_leaf(adj: list[list[int]]) -> tuple[int, int]:
    deg = [len(nb) for nb in adj]
    leaves = [v for v, d in enumerate(deg) if d == 1]
    parent = {l: adj[l][0] for l in leaves}
    leaf = max(leaves, key=lambda l: (-deg[parent[l]], -l))
    return leaf, parent[leaf]


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "no_deg2_support": 0,
        "phiP_fail": 0,
        "lambda_gap_fail": 0,
        "modeP_fail": 0,
        "ulc_fail": 0,
        "tau_skipped": 0,
        "modeP_shift_dist": {},
        "degP_shift_dist": {},
        "min_phiP": None,
        "min_phiP_witness": None,
        "min_mu_gap": None,
        "min_mu_gap_witness": None,
        "min_lambda_gap": None,
        "min_lambda_gap_witness": None,
        "min_tau_mu_gap": None,
        "min_tau_mu_gap_witness": None,
        "small_deg_shift_cases": [],
        "wall_s": 0.0,
    }


def inc_dist(dist: dict[str, int], key: int) -> None:
    sk = str(key)
    dist[sk] = dist.get(sk, 0) + 1


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
    for k in [
        "seen",
        "considered",
        "checked",
        "no_deg2_support",
        "phiP_fail",
        "lambda_gap_fail",
        "modeP_fail",
        "ulc_fail",
        "tau_skipped",
    ]:
        dst[k] += int(src.get(k, 0))

    for key in ["modeP_shift_dist", "degP_shift_dist"]:
        for kk, vv in src.get(key, {}).items():
            dst[key][kk] = dst[key].get(kk, 0) + int(vv)

    for key, witness_key in [
        ("min_phiP", "min_phiP_witness"),
        ("min_mu_gap", "min_mu_gap_witness"),
        ("min_lambda_gap", "min_lambda_gap_witness"),
        ("min_tau_mu_gap", "min_tau_mu_gap_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])

    # Keep at most first 200 representative small-shift cases.
    cases = dst["small_deg_shift_cases"] + src.get("small_deg_shift_cases", [])
    dst["small_deg_shift_cases"] = cases[:200]


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
    ap = argparse.ArgumentParser(description="Checkpointed Route-1 Phi(P) scan.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument(
        "--res",
        type=int,
        default=None,
        help="Optional nauty split residue r in r/m.",
    )
    ap.add_argument(
        "--mod",
        type=int,
        default=None,
        help="Optional nauty split modulus m in r/m.",
    )
    ap.add_argument("--out", default="results/whnc_phiP_scan_n23.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    if (args.res is None) != (args.mod is None):
        raise ValueError("Specify both --res and --mod together, or neither.")
    if args.mod is not None and not (0 <= args.res < args.mod):
        raise ValueError("Require 0 <= res < mod.")

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "tol": args.tol,
        "all_trees": args.all_trees,
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
        f"Route-1 Phi(P) scan on {scope}, n={args.min_n}..{args.max_n}{split_desc}",
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
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        if args.mod is not None:
            cmd.append(f"{args.res}/{args.mod}")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            stats["seen"] += 1

            if (not args.all_trees) and (not is_dleaf_le_1(nn, adj)):
                continue
            stats["considered"] += 1

            poly_t = independence_poly(nn, adj)
            m = mode_index_leftmost(poly_t)
            if m == 0 or poly_t[m - 1] <= 0 or poly_t[m] <= 0:
                continue
            lam = poly_t[m - 1] / poly_t[m]

            leaf, support = choose_min_support_leaf(adj)
            if len(adj[support]) != 2:
                stats["no_deg2_support"] += 1
                continue

            u = adj[support][0] if adj[support][1] == leaf else adj[support][1]
            b_adj, idx = remove_vertices(adj, {leaf, support})
            u_in_b = idx[u]
            p_poly = rooted_dp_u_excluded(b_adj, u_in_b)

            stats["checked"] += 1
            g6 = line.decode("ascii").strip()
            deg_sig = dict(sorted(Counter(len(nb) for nb in adj).items()))

            # Route-1 target
            phi_p = phi_q(p_poly, m - 1, lam)
            mu_p = mean_at_lambda(p_poly, lam)
            mu_gap = mu_p - (m - 2)

            if phi_p < -args.tol:
                stats["phiP_fail"] += 1
            if not is_ulc(p_poly):
                stats["ulc_fail"] += 1

            mode_p = mode_index_leftmost(p_poly)
            mode_shift = mode_p - (m - 1)
            inc_dist(stats["modeP_shift_dist"], mode_shift)
            if mode_shift < 0:
                stats["modeP_fail"] += 1

            deg_shift = (len(p_poly) - 1) - (m - 1)
            inc_dist(stats["degP_shift_dist"], deg_shift)

            base = {
                "n": nn,
                "g6": g6,
                "mode": m,
                "lambda_mode": lam,
                "leaf": leaf,
                "support": support,
                "mode_p": mode_p,
                "mode_p_shift": mode_shift,
                "deg_p": len(p_poly) - 1,
                "deg_p_shift": deg_shift,
                "phi_p": phi_p,
                "mu_p": mu_p,
                "mu_gap": mu_gap,
                "degree_signature": deg_sig,
            }

            maybe_update_min(stats, "min_phiP", "min_phiP_witness", phi_p, base)
            maybe_update_min(stats, "min_mu_gap", "min_mu_gap_witness", mu_gap, base)

            if deg_shift <= 1 and len(stats["small_deg_shift_cases"]) < 200:
                stats["small_deg_shift_cases"].append(base)

            # Tie comparison for P at r = m-1
            if m - 2 < 0 or m - 1 >= len(p_poly) or p_poly[m - 1] == 0:
                stats["tau_skipped"] += 1
                continue

            tau = p_poly[m - 2] / p_poly[m - 1]
            lam_gap = lam - tau
            tau_mu_gap = mean_at_lambda(p_poly, tau) - (m - 2)
            tie_w = {
                **base,
                "tau_p": tau,
                "lambda_gap": lam_gap,
                "tau_mu_gap": tau_mu_gap,
            }

            if lam_gap < -args.tol:
                stats["lambda_gap_fail"] += 1
            maybe_update_min(stats, "min_lambda_gap", "min_lambda_gap_witness", lam_gap, tie_w)
            maybe_update_min(stats, "min_tau_mu_gap", "min_tau_mu_gap_witness", tau_mu_gap, tie_w)

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"checked={stats['checked']:8d} phi_fail={stats['phiP_fail']:4d} "
            f"lam_fail={stats['lambda_gap_fail']:4d} ulc_fail={stats['ulc_fail']:4d} "
            f"min_phi={fmt_opt(stats['min_phiP'])} min_mu={fmt_opt(stats['min_mu_gap'])} "
            f"min_lgap={fmt_opt(stats['min_lambda_gap'])} min_taugap={fmt_opt(stats['min_tau_mu_gap'])} "
            f"({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    wall = time.time() - t_all
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"checked={summary['checked']:,} phi_fail={summary['phiP_fail']:,} "
        f"lam_fail={summary['lambda_gap_fail']:,} mode_fail={summary['modeP_fail']:,} "
        f"ulc_fail={summary['ulc_fail']:,} wall={wall:.1f}s",
        flush=True,
    )
    print(f"modeP_shift_dist={summary['modeP_shift_dist']}", flush=True)
    print(f"degP_shift_dist={summary['degP_shift_dist']}", flush=True)
    print(f"min_phiP={summary['min_phiP']}", flush=True)
    print(f"min_mu_gap={summary['min_mu_gap']}", flush=True)
    print(f"min_lambda_gap={summary['min_lambda_gap']}", flush=True)
    print(f"min_tau_mu_gap={summary['min_tau_mu_gap']}", flush=True)
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
