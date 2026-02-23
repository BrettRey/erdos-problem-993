#!/usr/bin/env python3
"""Profile obstruction terms for the STRONG C2 split on canonical degree-2 leaves.

For each d_leaf<=1 tree T, choose the canonical leaf:
  - minimum support degree, tie by smallest leaf id

Require deg(s)=2 for support s, with leaf l and other neighbor u.
Set B = T-{l,s}, P = dp_B[u][0], Q = dp_B[u][1], and let m = leftmost mode(I(T)).

At index m-1, define:
  b0=b_{m-2}, b1=b_{m-1}, b2=b_m
  p0=p_{m-2}, p1=p_{m-1}
  q0=q_{m-2}, q1=q_{m-1}

Core quantities:
  mismatch := p0*b1 - p1*b0
  lc_surplus := b1^2 - b2*b0
  combined := lc_surplus + mismatch
  neg := -mismatch (only when mismatch<0)
  rise := b1*(b1-b0)

This scanner tracks where mismatch<0 occurs and how it relates to:
  - q-drop (q1<q0),
  - B-mode shift,
  - rise compensation condition neg <= rise.

It is checkpointed by n, with optional nauty residue splitting via --mod.
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
        "rise_fail": 0,
        "mismatch_neg_shift0": 0,
        "mismatch_neg_shift1": 0,
        "mismatch_neg_shift_other": 0,
        "q_drop_all": 0,
        "q_drop_shift0": 0,
        "q_drop_shift1": 0,
        "q_drop_shift_other": 0,
        "q_drop_mismatch_neg": 0,
        "q_drop_mismatch_neg_shift0": 0,
        "q_drop_mismatch_neg_shift1": 0,
        "q_drop_mismatch_neg_shift_other": 0,
        "q_drop_ratio_undef": 0,
        "mismatch_neg_drop_neg": 0,  # b1 < b2 among mismatch-neg
        "min_delta_q": None,
        "min_delta_q_witness": None,
        "min_delta_q_mismatch_neg": None,
        "min_delta_q_mismatch_neg_witness": None,
        "min_rise_margin_neg": None,  # min of rise - neg over mismatch-neg
        "min_rise_margin_neg_witness": None,
        "max_neg_over_rise": None,
        "max_neg_over_rise_witness": None,
        "min_ratio_gap_neg_qdrop": None,  # (p1/b1) - ((-dq)/db)
        "min_ratio_gap_neg_qdrop_witness": None,
        "max_transfer_ratio_neg_qdrop": None,  # (-dq)/db
        "max_transfer_ratio_neg_qdrop_witness": None,
        "min_need_ratio_neg_qdrop": None,  # p1/b1
        "min_need_ratio_neg_qdrop_witness": None,
        "wall_s": 0.0,
    }


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
        "rise_fail",
        "mismatch_neg_shift0",
        "mismatch_neg_shift1",
        "mismatch_neg_shift_other",
        "q_drop_all",
        "q_drop_shift0",
        "q_drop_shift1",
        "q_drop_shift_other",
        "q_drop_mismatch_neg",
        "q_drop_mismatch_neg_shift0",
        "q_drop_mismatch_neg_shift1",
        "q_drop_mismatch_neg_shift_other",
        "q_drop_ratio_undef",
        "mismatch_neg_drop_neg",
    ]:
        dst[k] += int(src.get(k, 0))

    for key, witness_key in [
        ("min_delta_q", "min_delta_q_witness"),
        ("min_delta_q_mismatch_neg", "min_delta_q_mismatch_neg_witness"),
        ("min_rise_margin_neg", "min_rise_margin_neg_witness"),
        ("min_ratio_gap_neg_qdrop", "min_ratio_gap_neg_qdrop_witness"),
        ("min_need_ratio_neg_qdrop", "min_need_ratio_neg_qdrop_witness"),
    ]:
        if src.get(key) is not None:
            maybe_update_min(dst, key, witness_key, float(src[key]), src[witness_key])

    for key, witness_key in [
        ("max_neg_over_rise", "max_neg_over_rise_witness"),
        ("max_transfer_ratio_neg_qdrop", "max_transfer_ratio_neg_qdrop_witness"),
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
    ap = argparse.ArgumentParser(description="Obstruction profile for STRONG C2 split.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=22)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument(
        "--mod",
        type=int,
        default=1,
        help="Split geng enumeration into residue classes res/mod.",
    )
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--out", default="results/whnc_strong_c2_obstruction_n22.json")
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
        f"STRONG C2 obstruction scan on d_leaf<=1, n={args.min_n}..{args.max_n}",
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

                poly_t = independence_poly(nn, adj)
                m = mode_index_leftmost(poly_t)
                if m == 0 or m - 1 >= len(poly_t) or poly_t[m - 1] == 0 or poly_t[m] == 0:
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
                u_in_b = idx[u]
                p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

                p0 = p_poly[m - 2] if m - 2 < len(p_poly) else 0
                p1 = p_poly[m - 1] if m - 1 < len(p_poly) else 0
                q0 = q_poly[m - 2] if m - 2 < len(q_poly) else 0
                q1 = q_poly[m - 1] if m - 1 < len(q_poly) else 0
                b0 = b_poly[m - 2]
                b1 = b_poly[m - 1]
                b2 = b_poly[m]

                mismatch = p0 * b1 - p1 * b0
                lc_surplus = b1 * b1 - b2 * b0
                combined = lc_surplus + mismatch
                neg = -mismatch
                rise = b1 * (b1 - b0)
                rise_margin = rise - neg

                mode_b = mode_index_leftmost(b_poly) if b_poly else 0
                shift = mode_b - (m - 1)

                delta_q = q1 - q0
                delta_b = b1 - b0
                q_drop = delta_q < -args.tol

                w = {
                    "n": nn,
                    "g6": g6,
                    "mode_T": m,
                    "mode_B": mode_b,
                    "mode_b_shift": shift,
                    "leaf": leaf,
                    "support": support,
                    "p_m2": p0,
                    "p_m1": p1,
                    "q_m2": q0,
                    "q_m1": q1,
                    "b_m2": b0,
                    "b_m1": b1,
                    "b_m": b2,
                    "mismatch": mismatch,
                    "lc_surplus": lc_surplus,
                    "combined": combined,
                    "neg": neg,
                    "rise": rise,
                    "rise_margin": rise_margin,
                    "delta_q": delta_q,
                    "delta_b": delta_b,
                }

                maybe_update_min(stats, "min_delta_q", "min_delta_q_witness", float(delta_q), w)

                if q_drop:
                    stats["q_drop_all"] += 1
                    if shift == 0:
                        stats["q_drop_shift0"] += 1
                    elif shift == 1:
                        stats["q_drop_shift1"] += 1
                    else:
                        stats["q_drop_shift_other"] += 1

                if mismatch < -args.tol:
                    stats["mismatch_neg"] += 1
                    if shift == 0:
                        stats["mismatch_neg_shift0"] += 1
                    elif shift == 1:
                        stats["mismatch_neg_shift1"] += 1
                    else:
                        stats["mismatch_neg_shift_other"] += 1

                    maybe_update_min(
                        stats,
                        "min_delta_q_mismatch_neg",
                        "min_delta_q_mismatch_neg_witness",
                        float(delta_q),
                        w,
                    )
                    maybe_update_min(
                        stats,
                        "min_rise_margin_neg",
                        "min_rise_margin_neg_witness",
                        float(rise_margin),
                        w,
                    )

                    if q_drop:
                        stats["q_drop_mismatch_neg"] += 1
                        if shift == 0:
                            stats["q_drop_mismatch_neg_shift0"] += 1
                        elif shift == 1:
                            stats["q_drop_mismatch_neg_shift1"] += 1
                        else:
                            stats["q_drop_mismatch_neg_shift_other"] += 1

                        if delta_b > args.tol and b1 > 0:
                            transfer_ratio = (-delta_q) / delta_b
                            need_ratio = p1 / b1
                            ratio_gap = need_ratio - transfer_ratio
                            w_ratio = {
                                **w,
                                "transfer_ratio": transfer_ratio,
                                "need_ratio": need_ratio,
                                "ratio_gap": ratio_gap,
                            }
                            maybe_update_min(
                                stats,
                                "min_ratio_gap_neg_qdrop",
                                "min_ratio_gap_neg_qdrop_witness",
                                float(ratio_gap),
                                w_ratio,
                            )
                            maybe_update_max(
                                stats,
                                "max_transfer_ratio_neg_qdrop",
                                "max_transfer_ratio_neg_qdrop_witness",
                                float(transfer_ratio),
                                w_ratio,
                            )
                            maybe_update_min(
                                stats,
                                "min_need_ratio_neg_qdrop",
                                "min_need_ratio_neg_qdrop_witness",
                                float(need_ratio),
                                w_ratio,
                            )
                        else:
                            stats["q_drop_ratio_undef"] += 1

                    if rise > 0:
                        maybe_update_max(
                            stats,
                            "max_neg_over_rise",
                            "max_neg_over_rise_witness",
                            float(neg / rise),
                            w,
                        )
                    if rise_margin < -args.tol:
                        stats["rise_fail"] += 1
                    if b1 < b2:
                        stats["mismatch_neg_drop_neg"] += 1
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
            f"qdrop_neg={stats['q_drop_mismatch_neg']:4d} rise_fail={stats['rise_fail']:3d} "
            f"min_rise_margin={fmt_opt(stats['min_rise_margin_neg'])} ({stats['wall_s']:.1f}s)",
            flush=True,
        )

    summary = aggregate(per_n)
    print("-" * 96, flush=True)
    print(
        f"TOTAL seen={summary['seen']:,} considered={summary['considered']:,} "
        f"checked={summary['checked']:,} mismatch-={summary['mismatch_neg']:,} "
        f"qdrop_neg={summary['q_drop_mismatch_neg']:,} rise_fail={summary['rise_fail']:,} "
        f"comb-={summary['combined_neg']:,} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(
        "Extrema: "
        f"min_delta_q={summary['min_delta_q']}, "
        f"min_delta_q_neg={summary['min_delta_q_mismatch_neg']}, "
        f"min_rise_margin_neg={summary['min_rise_margin_neg']}, "
        f"min_ratio_gap_neg_qdrop={summary['min_ratio_gap_neg_qdrop']}, "
        f"max_neg_over_rise={summary['max_neg_over_rise']}",
        flush=True,
    )
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
