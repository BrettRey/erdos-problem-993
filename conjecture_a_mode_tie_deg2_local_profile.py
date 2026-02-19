#!/usr/bin/env python3
"""Profile degree-2-support local bridge invariants for focused mode-tie.

For each eligible tree T (default: d_leaf<=1), let m be the leftmost mode of
I(T) at lambda=1 and define lambda_m = i_{m-1}/i_m.

For each leaf l whose support s has degree 2, define:
  A = T - l
  B = T - {l, s}

and record the following local quantities:
  - Phi_m(A; lambda_m)
  - Phi_{m-1}(B; lambda_m)
  - mode(B) - (m-1)
  - tie gap: lambda_m(T) - lambda_{m-1}(B)
    where lambda_{m-1}(B) := b_{m-2}/b_{m-1}

The script checkpoints after each n and supports resume.
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


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked_trees": 0,
        "checked_deg2_leaves": 0,
        "fail_phi": 0,
        "fail_mode_b": 0,
        "fail_tie_gap": 0,
        "tie_skipped": 0,
        "mode_b_shift_dist": {},
        "min_phi_a": None,
        "min_phi_a_witness": None,
        "min_phi_b": None,
        "min_phi_b_witness": None,
        "min_tie_gap": None,
        "min_tie_gap_witness": None,
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
        "checked_trees",
        "checked_deg2_leaves",
        "fail_phi",
        "fail_mode_b",
        "fail_tie_gap",
        "tie_skipped",
    ]:
        dst[k] += int(src.get(k, 0))

    for key, val in src.get("mode_b_shift_dist", {}).items():
        dst["mode_b_shift_dist"][key] = dst["mode_b_shift_dist"].get(key, 0) + int(val)

    for key, witness_key in [
        ("min_phi_a", "min_phi_a_witness"),
        ("min_phi_b", "min_phi_b_witness"),
        ("min_tie_gap", "min_tie_gap_witness"),
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


def fmt_opt(val: float | None) -> str:
    return "None" if val is None else f"{val:.6g}"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Checkpointed profile of degree-2 local bridge invariants."
    )
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--all-trees", action="store_true")
    ap.add_argument("--out", default="results/whnc_mode_tie_deg2_local_profile_n23.json")
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
        done = ", ".join(sorted(per_n, key=lambda s: int(s)))
        print(f"Resuming from {args.out}; completed n: {done}", flush=True)

    scope = "all trees" if args.all_trees else "d_leaf<=1"
    print(
        f"Degree-2 local bridge profile on {scope}, n={args.min_n}..{args.max_n}",
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
            stats["checked_trees"] += 1

            deg = [len(nb) for nb in adj]
            g6 = line.decode("ascii").strip()
            for leaf, d in enumerate(deg):
                if d != 1:
                    continue
                support = adj[leaf][0]
                if deg[support] != 2:
                    continue

                stats["checked_deg2_leaves"] += 1

                a_adj = remove_vertices(adj, {leaf})
                b_adj = remove_vertices(adj, {leaf, support})
                a_poly = independence_poly(len(a_adj), a_adj)
                b_poly = independence_poly(len(b_adj), b_adj)

                phi_a = phi_q(a_poly, m, lam)
                phi_b = phi_q(b_poly, m - 1, lam)

                mode_b = mode_index_leftmost(b_poly) if b_poly else 0
                mode_b_shift = mode_b - (m - 1)
                inc_dist(stats["mode_b_shift_dist"], mode_b_shift)

                base_witness = {
                    "n": nn,
                    "g6": g6,
                    "mode": m,
                    "lambda_mode": lam,
                    "leaf": leaf,
                    "support": support,
                    "phi_a": phi_a,
                    "phi_b": phi_b,
                    "mode_b": mode_b,
                    "mode_b_shift": mode_b_shift,
                }

                maybe_update_min(
                    stats,
                    "min_phi_a",
                    "min_phi_a_witness",
                    phi_a,
                    base_witness,
                )
                maybe_update_min(
                    stats,
                    "min_phi_b",
                    "min_phi_b_witness",
                    phi_b,
                    base_witness,
                )

                if phi_a < -args.tol or phi_b < -args.tol:
                    stats["fail_phi"] += 1

                if mode_b < m - 1:
                    stats["fail_mode_b"] += 1

                if m - 2 < 0 or m - 1 >= len(b_poly) or b_poly[m - 1] == 0:
                    stats["tie_skipped"] += 1
                    continue

                lam_b = b_poly[m - 2] / b_poly[m - 1]
                tie_gap = lam - lam_b
                tie_witness = dict(base_witness)
                tie_witness["lambda_tie_b"] = lam_b
                tie_witness["tie_gap"] = tie_gap
                maybe_update_min(
                    stats,
                    "min_tie_gap",
                    "min_tie_gap_witness",
                    tie_gap,
                    tie_witness,
                )

                if tie_gap < -args.tol:
                    stats["fail_tie_gap"] += 1

        proc.wait()
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"trees={stats['checked_trees']:8d} deg2={stats['checked_deg2_leaves']:9d} "
            f"phi_fail={stats['fail_phi']:4d} modeB_fail={stats['fail_mode_b']:4d} "
            f"tie_fail={stats['fail_tie_gap']:4d} "
            f"min_phi_a={fmt_opt(stats['min_phi_a'])} "
            f"min_phi_b={fmt_opt(stats['min_phi_b'])} "
            f"min_tie_gap={fmt_opt(stats['min_tie_gap'])} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"trees={summary['checked_trees']:,} deg2={summary['checked_deg2_leaves']:,} "
        f"phi_fail={summary['fail_phi']:,} modeB_fail={summary['fail_mode_b']:,} "
        f"tie_fail={summary['fail_tie_gap']:,} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(
        f"Global minima: phi_a={summary['min_phi_a']}, phi_b={summary['min_phi_b']}, "
        f"tie_gap={summary['min_tie_gap']}",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
