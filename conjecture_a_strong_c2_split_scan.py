#!/usr/bin/env python3
"""Profile the STRONG C2 determinant split on canonical degree-2 leaves.

For each d_leaf<=1 tree T, choose the canonical leaf used by the successful
bridge heuristic: a leaf whose support has minimum degree (tie: smallest leaf id).

With m = mode(I(T)) at lambda=1, leaf l, support s (deg(s)=2), and
B = T-{l,s}, define:

  mismatch := p_{m-2} b_{m-1} - p_{m-1} b_{m-2}
  lc_surplus := b_{m-1}^2 - b_m b_{m-2}
  combined := lc_surplus + mismatch

where P = I(B-u) (u is the non-leaf neighbor of s), and b_k, p_k are
coefficients of I(B), P respectively.

The STRONG C2 cross-tree inequality is equivalent to combined >= 0.
This scanner quantifies how often mismatch < 0 and how much lc_surplus
compensates it.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly


def maybe_update_min(stats: dict[str, Any], key: str, witness_key: str, value: float, witness: dict[str, Any]) -> None:
    cur = stats[key]
    if cur is None or value < cur:
        stats[key] = value
        stats[witness_key] = witness


def maybe_update_max(stats: dict[str, Any], key: str, witness_key: str, value: float, witness: dict[str, Any]) -> None:
    cur = stats[key]
    if cur is None or value > cur:
        stats[key] = value
        stats[witness_key] = witness


def fresh_stats() -> dict[str, Any]:
    return {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "skip_no_mode_tie": 0,
        "skip_non_deg2_min_support": 0,
        "skip_short_B": 0,
        "mismatch_neg": 0,
        "mismatch_pos": 0,
        "mismatch_zero": 0,
        "combined_neg": 0,
        "mode_b_shift_dist": {},
        "min_mismatch": None,
        "min_mismatch_witness": None,
        "min_lc_surplus": None,
        "min_lc_surplus_witness": None,
        "min_combined": None,
        "min_combined_witness": None,
        "max_neg_over_lc": None,
        "max_neg_over_lc_witness": None,
        "max_neg_over_rise": None,
        "max_neg_over_rise_witness": None,
        "max_neg_over_drop": None,
        "max_neg_over_drop_witness": None,
        "bound_rise_fail": 0,  # (-mismatch) > b_{m-1}(b_{m-1}-b_{m-2})
        "bound_drop_fail": 0,  # (-mismatch) > b_{m-2}(b_{m-1}-b_m), only when drop>=0
        "wall_s": 0.0,
    }


def inc_dist(dist: dict[str, int], key: int) -> None:
    s = str(key)
    dist[s] = dist.get(s, 0) + 1


def merge_stats(dst: dict[str, Any], src: dict[str, Any]) -> None:
    for k in [
        "seen",
        "considered",
        "checked",
        "skip_no_mode_tie",
        "skip_non_deg2_min_support",
        "skip_short_B",
        "mismatch_neg",
        "mismatch_pos",
        "mismatch_zero",
        "combined_neg",
        "bound_rise_fail",
        "bound_drop_fail",
    ]:
        dst[k] += int(src.get(k, 0))

    for k, v in src.get("mode_b_shift_dist", {}).items():
        dst["mode_b_shift_dist"][k] = dst["mode_b_shift_dist"].get(k, 0) + int(v)

    for key, witness_key in [
        ("min_mismatch", "min_mismatch_witness"),
        ("min_lc_surplus", "min_lc_surplus_witness"),
        ("min_combined", "min_combined_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])

    for key, witness_key in [
        ("max_neg_over_lc", "max_neg_over_lc_witness"),
        ("max_neg_over_rise", "max_neg_over_rise_witness"),
        ("max_neg_over_drop", "max_neg_over_drop_witness"),
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
    ap = argparse.ArgumentParser(description="Scan STRONG C2 split: lc_surplus + mismatch.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=22)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument(
        "--mod",
        type=int,
        default=1,
        help="Split geng enumeration into residue classes res/mod (recommended 8 for n=23).",
    )
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--out", default="results/whnc_strong_c2_split_n22.json")
    ap.add_argument("--no-resume", action="store_true")
    args = ap.parse_args()

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "mod": args.mod,
        "tol": args.tol,
    }

    per_n: dict[str, dict[str, Any]] = {}
    if (not args.no_resume) and args.out and os.path.exists(args.out):
        with open(args.out, "r", encoding="utf-8") as f:
            old = json.load(f)
        per_n = old.get("per_n", {})
        done = ", ".join(sorted(per_n, key=lambda s: int(s)))
        print(f"Resuming from {args.out}; completed n: {done}", flush=True)

    print(
        f"STRONG C2 split scan on d_leaf<=1, n={args.min_n}..{args.max_n}",
        flush=True,
    )
    print("Canonical leaf: minimum support degree (tie: smallest leaf id)", flush=True)
    print("Output:", args.out, flush=True)
    print("-" * 96, flush=True)

    t_all = time.time()
    for n in range(args.min_n, args.max_n + 1):
        n_key = str(n)
        if n_key in per_n:
            print(f"n={n:2d}: already complete, skipping", flush=True)
            continue

        t0 = time.time()
        stats = fresh_stats()

        for res in range(args.mod):
            cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
            if args.mod > 1:
                cmd.append(f"{res}/{args.mod}")

            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert proc.stdout is not None

            for raw in proc.stdout:
                stats["seen"] += 1
                nn, adj = parse_graph6(raw)

                if not is_dleaf_le_1(nn, adj):
                    continue
                stats["considered"] += 1

                poly_T = independence_poly(nn, adj)
                m = mode_index_leftmost(poly_T)
                if m == 0 or m - 1 >= len(poly_T) or poly_T[m - 1] == 0 or poly_T[m] == 0:
                    stats["skip_no_mode_tie"] += 1
                    continue

                deg = [len(nb) for nb in adj]
                leaves = [v for v in range(nn) if deg[v] == 1]
                min_parent_deg = min(deg[adj[l][0]] for l in leaves)
                leaf = min(l for l in leaves if deg[adj[l][0]] == min_parent_deg)
                support = adj[leaf][0]

                if deg[support] != 2:
                    stats["skip_non_deg2_min_support"] += 1
                    continue

                stats["checked"] += 1
                g6 = raw.decode("ascii").strip()

                u = [x for x in adj[support] if x != leaf][0]
                b_adj = remove_vertices(adj, {leaf, support})
                b_poly = independence_poly(len(b_adj), b_adj)

                if m - 2 < 0 or m >= len(b_poly) or m - 1 >= len(b_poly) or b_poly[m - 1] == 0:
                    stats["skip_short_B"] += 1
                    continue

                keep = [v for v in range(nn) if v not in {leaf, support}]
                idx = {v: i for i, v in enumerate(keep)}
                u_in_B = idx[u]
                P, _Q = compute_hub_polys(b_adj, u_in_B)

                p_m2 = P[m - 2] if m - 2 < len(P) else 0
                p_m1 = P[m - 1] if m - 1 < len(P) else 0
                b_m2 = b_poly[m - 2]
                b_m1 = b_poly[m - 1]
                b_m = b_poly[m]

                mismatch = p_m2 * b_m1 - p_m1 * b_m2
                lc_surplus = b_m1 * b_m1 - b_m * b_m2
                combined = lc_surplus + mismatch
                mode_b = mode_index_leftmost(b_poly) if b_poly else 0
                mode_b_shift = mode_b - (m - 1)
                inc_dist(stats["mode_b_shift_dist"], mode_b_shift)

                rise_term = b_m1 * (b_m1 - b_m2)
                drop_term = b_m2 * (b_m1 - b_m)

                w = {
                    "n": nn,
                    "g6": g6,
                    "mode_T": m,
                    "mode_B": mode_b,
                    "mode_b_shift": mode_b_shift,
                    "leaf": leaf,
                    "support": support,
                    "mismatch": mismatch,
                    "lc_surplus": lc_surplus,
                    "combined": combined,
                    "p_m2": p_m2,
                    "p_m1": p_m1,
                    "b_m2": b_m2,
                    "b_m1": b_m1,
                    "b_m": b_m,
                    "rise_term": rise_term,
                    "drop_term": drop_term,
                }

                maybe_update_min(stats, "min_mismatch", "min_mismatch_witness", mismatch, w)
                maybe_update_min(stats, "min_lc_surplus", "min_lc_surplus_witness", lc_surplus, w)
                maybe_update_min(stats, "min_combined", "min_combined_witness", combined, w)

                if mismatch < -args.tol:
                    stats["mismatch_neg"] += 1
                    neg = -mismatch
                    if lc_surplus > 0:
                        maybe_update_max(stats, "max_neg_over_lc", "max_neg_over_lc_witness", neg / lc_surplus, w)
                    if rise_term > 0:
                        maybe_update_max(stats, "max_neg_over_rise", "max_neg_over_rise_witness", neg / rise_term, w)
                    if drop_term > 0:
                        maybe_update_max(stats, "max_neg_over_drop", "max_neg_over_drop_witness", neg / drop_term, w)

                    if neg > rise_term:
                        stats["bound_rise_fail"] += 1
                    if drop_term >= 0 and neg > drop_term:
                        stats["bound_drop_fail"] += 1
                elif mismatch > args.tol:
                    stats["mismatch_pos"] += 1
                else:
                    stats["mismatch_zero"] += 1

                if combined < -args.tol:
                    stats["combined_neg"] += 1

            stderr = proc.stderr.read().decode("utf-8", errors="replace")
            ret = proc.wait()
            if ret != 0:
                raise RuntimeError(
                    f"geng failed at n={n}, res={res}/{args.mod}, returncode={ret}, stderr={stderr.strip()}"
                )
        stats["wall_s"] = time.time() - t0
        per_n[n_key] = stats
        write_payload(args.out, params, per_n)

        print(
            f"n={n:2d}: seen={stats['seen']:9d} considered={stats['considered']:8d} "
            f"checked={stats['checked']:8d} mismatch-={stats['mismatch_neg']:5d} "
            f"comb-={stats['combined_neg']:3d} min_comb={fmt_opt(stats['min_combined'])} "
            f"max(-mis/lc)={fmt_opt(stats['max_neg_over_lc'])} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"checked={summary['checked']:,} mismatch-={summary['mismatch_neg']:,} "
        f"comb-={summary['combined_neg']:,} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(
        f"Global mins/maxes: min_mismatch={summary['min_mismatch']}, "
        f"min_lc={summary['min_lc_surplus']}, min_comb={summary['min_combined']}, "
        f"max(-mis/lc)={summary['max_neg_over_lc']}, "
        f"max(-mis/rise)={summary['max_neg_over_rise']}",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
