#!/usr/bin/env python3
"""Attack 3 scan: exact mean-lift mu_T - mu_P at mode tie fugacity.

Setup per decomposition:
  - T is a d_leaf<=1 tree.
  - Pick leaf l whose support s has deg(s)=2, and u the other neighbor of s.
  - B = T - {l,s}, rooted at u.
  - P = dp_B[u][0], Q = dp_B[u][1], so I(B)=P+Q.
  - I(T) = (1+2x)P + (1+x)Q.
  - m = mode(I(T)), lambda = i_{m-1}(T)/i_m(T).

Main derived identities verified numerically:
  lift := mu_T(lambda) - mu_P(lambda)

  (Z form)
    lift =
      [lambda*(2 Z_P + Z_Q) + (1+lambda) Z_Q (mu_Q - mu_P)] / Z_T
    where Z_T = (1+2lambda)Z_P + (1+lambda)Z_Q.

  (p_u, D form)
    p_u := Z_Q / (Z_P + Z_Q), D := mu_B - mu_P,
    lift = [lambda*(2 - p_u) + (1+lambda)D] / [1 + 2lambda - lambda p_u].

  (A+B decomposition)
    I(A)=(1+x)P+Q, I(T)=I(A)+xI(B),
    mu_T = (Z_A/Z_T)mu_A + (lambda Z_B/Z_T)(1+mu_B).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def choose_canonical_leaf(adj: list[list[int]]) -> tuple[int, int]:
    """Choose canonical leaf by min support degree, tie by largest leaf id."""
    deg = [len(nb) for nb in adj]
    leaves = [v for v, d in enumerate(deg) if d == 1]
    parent = {l: adj[l][0] for l in leaves}
    leaf = max(leaves, key=lambda l: (-deg[parent[l]], -l))
    return leaf, parent[leaf]


def all_deg2_support_leaves(adj: list[list[int]]) -> list[tuple[int, int]]:
    deg = [len(nb) for nb in adj]
    pairs: list[tuple[int, int]] = []
    for leaf in range(len(adj)):
        if deg[leaf] != 1:
            continue
        support = adj[leaf][0]
        if deg[support] == 2:
            pairs.append((leaf, support))
    return pairs


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
    z = 0.0
    p = 1.0
    for ck in poly:
        z += ck * p
        p *= lam
    return z


def mean_at_lambda(poly: list[int], lam: float) -> float:
    z = 0.0
    mu_num = 0.0
    p = 1.0
    for k, ck in enumerate(poly):
        w = ck * p
        z += w
        mu_num += k * w
        p *= lam
    return mu_num / z if z else float("nan")


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked_trees": 0,
        "checked_leaves": 0,
        "skip_invalid_mode_tie": 0,
        "no_deg2_support_trees": 0,
        "lift_ge_1_fail": 0,
        "muT_below_m_minus_1_fail": 0,
        "muT_above_m_fail": 0,
        "formula_z_fail": 0,
        "formula_pu_fail": 0,
        "formula_a_fail": 0,
        "mode_p_lt_m_minus_1_fail": 0,
        "mode_p_lt_floor_muP_fail": 0,
        "muP_lt_m_minus_1_fail": 0,
        "max_lift": None,
        "max_lift_witness": None,
        "min_lift": None,
        "min_lift_witness": None,
        "min_one_minus_lift": None,
        "min_one_minus_lift_witness": None,
        "max_D": None,
        "max_D_witness": None,
        "max_p_u": None,
        "max_p_u_witness": None,
        "max_muQ_minus_muP": None,
        "max_muQ_minus_muP_witness": None,
        "max_base_term": None,
        "max_base_term_witness": None,
        "max_transfer_term": None,
        "max_transfer_term_witness": None,
        "max_abs_formula_z_err": None,
        "max_abs_formula_z_err_witness": None,
        "max_abs_formula_pu_err": None,
        "max_abs_formula_pu_err_witness": None,
        "max_abs_formula_a_err": None,
        "max_abs_formula_a_err_witness": None,
        "min_muT_minus_m_minus_1": None,
        "min_muT_minus_m_minus_1_witness": None,
        "min_m_minus_muT": None,
        "min_m_minus_muT_witness": None,
        "min_modeP_minus_m_minus_1": None,
        "min_modeP_minus_m_minus_1_witness": None,
        "min_modeP_minus_floor_muP": None,
        "min_modeP_minus_floor_muP_witness": None,
        "min_muP_minus_m_minus_1": None,
        "min_muP_minus_m_minus_1_witness": None,
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
        "checked_trees",
        "checked_leaves",
        "skip_invalid_mode_tie",
        "no_deg2_support_trees",
        "lift_ge_1_fail",
        "muT_below_m_minus_1_fail",
        "muT_above_m_fail",
        "formula_z_fail",
        "formula_pu_fail",
        "formula_a_fail",
        "mode_p_lt_m_minus_1_fail",
        "mode_p_lt_floor_muP_fail",
        "muP_lt_m_minus_1_fail",
    ]:
        dst[k] += int(src.get(k, 0))

    for key, witness_key in [
        ("min_lift", "min_lift_witness"),
        ("min_one_minus_lift", "min_one_minus_lift_witness"),
        ("min_muT_minus_m_minus_1", "min_muT_minus_m_minus_1_witness"),
        ("min_m_minus_muT", "min_m_minus_muT_witness"),
        ("min_modeP_minus_m_minus_1", "min_modeP_minus_m_minus_1_witness"),
        ("min_modeP_minus_floor_muP", "min_modeP_minus_floor_muP_witness"),
        ("min_muP_minus_m_minus_1", "min_muP_minus_m_minus_1_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])

    for key, witness_key in [
        ("max_lift", "max_lift_witness"),
        ("max_D", "max_D_witness"),
        ("max_p_u", "max_p_u_witness"),
        ("max_muQ_minus_muP", "max_muQ_minus_muP_witness"),
        ("max_base_term", "max_base_term_witness"),
        ("max_transfer_term", "max_transfer_term_witness"),
        ("max_abs_formula_z_err", "max_abs_formula_z_err_witness"),
        ("max_abs_formula_pu_err", "max_abs_formula_pu_err_witness"),
        ("max_abs_formula_a_err", "max_abs_formula_a_err_witness"),
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
    ap = argparse.ArgumentParser(description="Attack3 scan: exact mu_T-mu_P lift identities and extrema.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--identity-tol", type=float, default=1e-9)
    ap.add_argument("--floor-eps", type=float, default=1e-12)
    ap.add_argument("--leaf-mode", choices=["canonical", "all_deg2"], default="canonical")
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument("--res", type=int, default=None, help="Optional split residue r in r/m.")
    ap.add_argument("--mod", type=int, default=None, help="Optional split modulus m in r/m.")
    ap.add_argument("--out", default="results/attack3_mean_lift_scan_n23_canonical_2026_02_19.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    if (args.res is None) != (args.mod is None):
        raise ValueError("Specify both --res and --mod together, or neither.")
    if args.mod is not None and not (0 <= args.res < args.mod):
        raise ValueError("Require 0 <= res < mod.")

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "leaf_mode": args.leaf_mode,
        "all_trees": args.all_trees,
        "tol": args.tol,
        "identity_tol": args.identity_tol,
        "floor_eps": args.floor_eps,
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
        f"attack3 mean-lift scan on {scope}, n={args.min_n}..{args.max_n}, "
        f"leaf_mode={args.leaf_mode}{split_desc}",
        flush=True,
    )
    print(f"Output: {args.out}", flush=True)
    print("-" * 100, flush=True)

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
            if m == 0 or m >= len(poly_t) or poly_t[m - 1] <= 0 or poly_t[m] <= 0:
                stats["skip_invalid_mode_tie"] += 1
                continue

            lam = poly_t[m - 1] / poly_t[m]
            z_t = eval_poly(poly_t, lam)
            mu_t = mean_at_lambda(poly_t, lam)
            g6 = raw.decode("ascii").strip()

            if args.leaf_mode == "canonical":
                leaf, support = choose_canonical_leaf(adj)
                leaf_support_pairs = [(leaf, support)]
            else:
                leaf_support_pairs = all_deg2_support_leaves(adj)

            if not leaf_support_pairs:
                stats["no_deg2_support_trees"] += 1
                continue

            valid_in_tree = 0
            for leaf, support in leaf_support_pairs:
                if len(adj[support]) != 2:
                    continue
                u = adj[support][0] if adj[support][1] == leaf else adj[support][1]
                b_adj, idx = remove_vertices(adj, {leaf, support})
                u_in_b = idx[u]

                p_poly, q_poly = rooted_dp_u(b_adj, u_in_b)
                b_poly = _polyadd(p_poly, q_poly)

                z_p = eval_poly(p_poly, lam)
                z_q = eval_poly(q_poly, lam)
                z_b = z_p + z_q
                if z_b <= 0.0 or z_t <= 0.0:
                    continue

                mu_p = mean_at_lambda(p_poly, lam)
                mu_q = mean_at_lambda(q_poly, lam)
                mu_b = mean_at_lambda(b_poly, lam)

                p_u = z_q / z_b
                d_val = mu_b - mu_p
                dq = mu_q - mu_p

                lift = mu_t - mu_p
                base_term = lam * (2.0 * z_p + z_q) / z_t
                transfer_term = (1.0 + lam) * z_q * dq / z_t
                lift_formula_z = base_term + transfer_term

                denom = 1.0 + 2.0 * lam - lam * p_u
                lift_formula_pu = (lam * (2.0 - p_u) + (1.0 + lam) * d_val) / denom

                # A-decomposition formula:
                # mu_T = (Z_A/Z_T)mu_A + (lambda Z_B/Z_T)(1+mu_B)
                z_a = (1.0 + lam) * z_p + z_q
                mu_a_num = lam * z_p + (1.0 + lam) * mu_p * z_p + mu_q * z_q
                mu_a = mu_a_num / z_a
                mu_t_formula_a = (z_a / z_t) * mu_a + (lam * z_b / z_t) * (1.0 + mu_b)

                err_z = abs(lift - lift_formula_z)
                err_pu = abs(lift - lift_formula_pu)
                err_a = abs(mu_t - mu_t_formula_a)

                mode_p = mode_index_leftmost(p_poly)
                floor_mu_p = math.floor(mu_p + args.floor_eps)
                mode_gap_m = mode_p - (m - 1)
                mode_gap_floor = mode_p - floor_mu_p
                mu_p_gap_m1 = mu_p - (m - 1)
                mu_t_gap_lower = mu_t - (m - 1)
                mu_t_gap_upper = m - mu_t
                one_minus_lift = 1.0 - lift

                witness = {
                    "n": nn,
                    "g6": g6,
                    "mode_T": m,
                    "lambda_mode": lam,
                    "leaf": leaf,
                    "support": support,
                    "u": u,
                    "mode_P": mode_p,
                    "floor_mu_P": floor_mu_p,
                    "mu_T": mu_t,
                    "mu_P": mu_p,
                    "mu_Q": mu_q,
                    "mu_B": mu_b,
                    "lift": lift,
                    "one_minus_lift": one_minus_lift,
                    "D": d_val,
                    "p_u": p_u,
                    "mu_Q_minus_mu_P": dq,
                    "base_term": base_term,
                    "transfer_term": transfer_term,
                    "z_P": z_p,
                    "z_Q": z_q,
                    "z_B": z_b,
                    "z_T": z_t,
                    "formula_z_err": err_z,
                    "formula_pu_err": err_pu,
                    "formula_a_err": err_a,
                    "modeP_minus_m_minus_1": mode_gap_m,
                    "modeP_minus_floor_muP": mode_gap_floor,
                    "muP_minus_m_minus_1": mu_p_gap_m1,
                    "muT_minus_m_minus_1": mu_t_gap_lower,
                    "m_minus_muT": mu_t_gap_upper,
                }

                valid_in_tree += 1
                stats["checked_leaves"] += 1

                if lift >= 1.0 + args.tol:
                    stats["lift_ge_1_fail"] += 1
                if mu_t < (m - 1) - args.tol:
                    stats["muT_below_m_minus_1_fail"] += 1
                if mu_t > m + args.tol:
                    stats["muT_above_m_fail"] += 1
                if err_z > args.identity_tol:
                    stats["formula_z_fail"] += 1
                if err_pu > args.identity_tol:
                    stats["formula_pu_fail"] += 1
                if err_a > args.identity_tol:
                    stats["formula_a_fail"] += 1
                if mode_gap_m < 0:
                    stats["mode_p_lt_m_minus_1_fail"] += 1
                if mode_gap_floor < 0:
                    stats["mode_p_lt_floor_muP_fail"] += 1
                if mu_p_gap_m1 < -args.tol:
                    stats["muP_lt_m_minus_1_fail"] += 1

                maybe_update_max(stats, "max_lift", "max_lift_witness", lift, witness)
                maybe_update_min(stats, "min_lift", "min_lift_witness", lift, witness)
                maybe_update_min(stats, "min_one_minus_lift", "min_one_minus_lift_witness", one_minus_lift, witness)
                maybe_update_max(stats, "max_D", "max_D_witness", d_val, witness)
                maybe_update_max(stats, "max_p_u", "max_p_u_witness", p_u, witness)
                maybe_update_max(stats, "max_muQ_minus_muP", "max_muQ_minus_muP_witness", dq, witness)
                maybe_update_max(stats, "max_base_term", "max_base_term_witness", base_term, witness)
                maybe_update_max(stats, "max_transfer_term", "max_transfer_term_witness", transfer_term, witness)
                maybe_update_max(stats, "max_abs_formula_z_err", "max_abs_formula_z_err_witness", err_z, witness)
                maybe_update_max(stats, "max_abs_formula_pu_err", "max_abs_formula_pu_err_witness", err_pu, witness)
                maybe_update_max(stats, "max_abs_formula_a_err", "max_abs_formula_a_err_witness", err_a, witness)
                maybe_update_min(stats, "min_muT_minus_m_minus_1", "min_muT_minus_m_minus_1_witness", mu_t_gap_lower, witness)
                maybe_update_min(stats, "min_m_minus_muT", "min_m_minus_muT_witness", mu_t_gap_upper, witness)
                maybe_update_min(stats, "min_modeP_minus_m_minus_1", "min_modeP_minus_m_minus_1_witness", float(mode_gap_m), witness)
                maybe_update_min(stats, "min_modeP_minus_floor_muP", "min_modeP_minus_floor_muP_witness", float(mode_gap_floor), witness)
                maybe_update_min(stats, "min_muP_minus_m_minus_1", "min_muP_minus_m_minus_1_witness", mu_p_gap_m1, witness)

            if valid_in_tree == 0:
                stats["no_deg2_support_trees"] += 1
            else:
                stats["checked_trees"] += 1

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"checked_trees={stats['checked_trees']:8d} checked_leaves={stats['checked_leaves']:9d} "
            f"max_lift={fmt_opt(stats['max_lift'])} min(1-lift)={fmt_opt(stats['min_one_minus_lift'])} "
            f"lift>=1_fail={stats['lift_ge_1_fail']:3d} modeP<m-1={stats['mode_p_lt_m_minus_1_fail']:3d} "
            f"modeP<floor(muP)={stats['mode_p_lt_floor_muP_fail']:3d} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    total_s = time.time() - t_all
    summary = aggregate(per_n)
    print("-" * 100, flush=True)
    print(f"Completed in {total_s:.1f}s", flush=True)
    print(
        f"TOTAL: seen={summary['seen']:,} considered={summary['considered']:,} "
        f"checked_trees={summary['checked_trees']:,} checked_leaves={summary['checked_leaves']:,}",
        flush=True,
    )
    print(
        f"max_lift={fmt_opt(summary['max_lift'])}, min(1-lift)={fmt_opt(summary['min_one_minus_lift'])}, "
        f"lift>=1_fail={summary['lift_ge_1_fail']}",
        flush=True,
    )
    print(
        f"max_D={fmt_opt(summary['max_D'])}, max_p_u={fmt_opt(summary['max_p_u'])}, "
        f"max(muQ-muP)={fmt_opt(summary['max_muQ_minus_muP'])}",
        flush=True,
    )
    print(
        f"modeP<m-1_fail={summary['mode_p_lt_m_minus_1_fail']}, "
        f"modeP<floor(muP)_fail={summary['mode_p_lt_floor_muP_fail']}, "
        f"muP<m-1_fail={summary['muP_lt_m_minus_1_fail']}",
        flush=True,
    )
    print(f"results written to: {args.out}", flush=True)


if __name__ == "__main__":
    main()
