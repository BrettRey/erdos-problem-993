#!/usr/bin/env python3
"""Profile focused mode-tie leaf bridge under min-support-degree leaf choice.

For each tree T (default: d_leaf<=1), let m be the leftmost mode of I(T) and
lambda_m = i_{m-1}/i_m. Choose a leaf l whose support has minimum degree
(tie-break: smallest leaf index), and define:

  A = T - l
  B = T - {l, s}  where s is the support of l

Record:
- mode(A) - m
- mode(B) - (m - 1)
- Phi_m(A; lambda_m) and Phi_{m-1}(B; lambda_m)

The script checkpoints per-n results to JSON after each n so interrupted runs
can be resumed without rerunning completed ranges.
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


def choose_leaf_min_parent_degree(adj: list[list[int]]) -> tuple[int, int]:
    deg = [len(nb) for nb in adj]
    leaves = [v for v, d in enumerate(deg) if d == 1]
    parent = {l: adj[l][0] for l in leaves}
    # Same tie policy as conjecture_a_mode_tie_leaf_bridge_scan.py:
    # max with key (-parent_deg, -leaf) picks smallest parent degree and leaf.
    leaf = max(leaves, key=lambda l: (-deg[parent[l]], -l))
    return leaf, parent[leaf]


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "fail_phi": 0,
        "chosen_support_deg2": 0,
        "chosen_support_deg_dist": {},
        "mode_a_shift_dist": {},
        "mode_b_shift_dist": {},
        "min_phi_a": None,
        "min_phi_b": None,
        "min_phi_a_witness": None,
        "min_phi_b_witness": None,
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
    for k in ["seen", "considered", "checked", "fail_phi", "chosen_support_deg2"]:
        dst[k] += int(src.get(k, 0))
    for k in ["chosen_support_deg_dist", "mode_a_shift_dist", "mode_b_shift_dist"]:
        for key, val in src.get(k, {}).items():
            dst[k][key] = dst[k].get(key, 0) + int(val)

    if src.get("min_phi_a") is not None:
        maybe_update_min(
            dst,
            "min_phi_a",
            "min_phi_a_witness",
            float(src["min_phi_a"]),
            src["min_phi_a_witness"],
        )
    if src.get("min_phi_b") is not None:
        maybe_update_min(
            dst,
            "min_phi_b",
            "min_phi_b_witness",
            float(src["min_phi_b"]),
            src["min_phi_b_witness"],
        )


def aggregate(per_n: dict[str, dict[str, Any]]) -> dict[str, Any]:
    out = fresh_stats()
    for key in sorted(per_n, key=lambda s: int(s)):
        merge_stats(out, per_n[key])
    out.pop("wall_s", None)
    return out


def write_payload(
    out_path: str,
    params: dict[str, Any],
    per_n: dict[str, dict[str, Any]],
) -> None:
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


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Profile min-support-degree leaf bridge mode shifts (checkpointed)."
    )
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument("--out", default="results/whnc_mode_tie_minleaf_profile_n23.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "tol": args.tol,
        "all_trees": args.all_trees,
    }

    per_n: dict[str, dict[str, Any]] = {}
    if (not args.no_resume) and args.out and os.path.exists(args.out):
        with open(args.out, "r", encoding="utf-8") as f:
            old = json.load(f)
        per_n = old.get("per_n", {})
        print(
            f"Resuming from {args.out}; completed n: "
            f"{', '.join(sorted(per_n, key=lambda s: int(s)))}",
            flush=True,
        )

    scope = "all trees" if args.all_trees else "d_leaf<=1"
    print(
        f"Min-parent-degree leaf profile on {scope} for n={args.min_n}..{args.max_n}",
        flush=True,
    )
    print("Checkpoint output:", args.out, flush=True)
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
            leaf, support = choose_leaf_min_parent_degree(adj)
            deg = [len(nb) for nb in adj]
            sdeg = deg[support]

            stats["checked"] += 1
            if sdeg == 2:
                stats["chosen_support_deg2"] += 1
            inc_dist(stats["chosen_support_deg_dist"], sdeg)

            a_adj = remove_vertices(adj, {leaf})
            b_adj = remove_vertices(adj, {leaf, support})
            a_poly = independence_poly(len(a_adj), a_adj)
            b_poly = independence_poly(len(b_adj), b_adj)

            mode_a = mode_index_leftmost(a_poly)
            mode_b = mode_index_leftmost(b_poly) if b_poly else 0
            inc_dist(stats["mode_a_shift_dist"], mode_a - m)
            inc_dist(stats["mode_b_shift_dist"], mode_b - (m - 1))

            phi_a = phi_q(a_poly, m, lam)
            phi_b = phi_q(b_poly, m - 1, lam)
            if phi_a < -args.tol or phi_b < -args.tol:
                stats["fail_phi"] += 1

            deg_sig = dict(sorted(Counter(deg).items()))
            base_witness = {
                "n": nn,
                "g6": line.decode("ascii").strip(),
                "mode": m,
                "lambda_mode": lam,
                "leaf": leaf,
                "support": support,
                "support_degree": sdeg,
                "mode_a": mode_a,
                "mode_b": mode_b,
                "degree_signature": deg_sig,
            }
            maybe_update_min(
                stats,
                "min_phi_a",
                "min_phi_a_witness",
                phi_a,
                {**base_witness, "phi_a": phi_a, "phi_b": phi_b},
            )
            maybe_update_min(
                stats,
                "min_phi_b",
                "min_phi_b_witness",
                phi_b,
                {**base_witness, "phi_a": phi_a, "phi_b": phi_b},
            )

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats

        if args.out:
            write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"checked={stats['checked']:8d} fail_phi={stats['fail_phi']:4d} "
            f"deg2={stats['chosen_support_deg2']:8d} "
            f"min_phi_a={stats['min_phi_a']:.6f} min_phi_b={stats['min_phi_b']:.6f} "
            f"({stats['wall_s']:.1f}s)",
            flush=True,
        )

    total = aggregate(per_n)
    wall = time.time() - t_all
    print("-" * 96, flush=True)
    print(
        f"TOTAL checked={total['checked']:,} fail_phi={total['fail_phi']:,} "
        f"deg2={total['chosen_support_deg2']:,} wall={wall:.1f}s",
        flush=True,
    )
    print(f"mode_a_shift_dist={total['mode_a_shift_dist']}", flush=True)
    print(f"mode_b_shift_dist={total['mode_b_shift_dist']}", flush=True)
    print(f"min_phi_a={total['min_phi_a']}", flush=True)
    print(f"min_phi_b={total['min_phi_b']}", flush=True)
    if args.out:
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
