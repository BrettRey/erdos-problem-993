#!/usr/bin/env python3
"""Verify the STRONG C2 rise-compensation algebraic split on canonical leaves.

For each d_leaf<=1 tree T:
  - choose canonical leaf (minimum support degree, tie by leaf id)
  - require deg(support)=2
  - set B = T-{leaf, support}, and P,Q from root-u decomposition in B

At m = leftmost mode index of I(T), define:
  p0 = p_{m-2}, p1 = p_{m-1}
  q0 = q_{m-2}, q1 = q_{m-1}
  b0 = b_{m-2}, b1 = b_{m-1}, b2 = b_m

Core terms:
  mismatch   = p0*b1 - p1*b0
  lc_surplus = b1*b1 - b2*b0
  combined   = lc_surplus + mismatch
  neg        = -mismatch
  rise       = b1*(b1-b0)

Identity checked exactly:
  rise - neg = p1*(b1-b0) + b1*(q1-q0)

Hard regime (for rise-compensation):
  mismatch < 0 and q1 < q0
  In this regime, rise-neg >= 0 is equivalent to
    (- (q1-q0) / (b1-b0)) <= p1/b1.
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


def maybe_min(stats: dict[str, Any], key: str, wit_key: str, value: float, witness: dict[str, Any]) -> None:
    cur = stats.get(key)
    if cur is None or value < cur:
        stats[key] = value
        stats[wit_key] = witness


def maybe_max(stats: dict[str, Any], key: str, wit_key: str, value: float, witness: dict[str, Any]) -> None:
    cur = stats.get(key)
    if cur is None or value > cur:
        stats[key] = value
        stats[wit_key] = witness


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify rise-compensation split and hard-regime ratio condition.")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--mod", type=int, default=1, help="Residue modulus for geng split.")
    ap.add_argument("--res", type=int, default=-1, help="Single residue in [0,mod-1]; -1 means all residues.")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument(
        "--out",
        default="results/verify_strong_c2_rise_identity_2026_02_19.json",
    )
    args = ap.parse_args()

    if args.mod <= 0:
        raise ValueError("--mod must be >= 1")
    if args.res >= args.mod:
        raise ValueError("--res must be < --mod")

    residues = list(range(args.mod)) if args.res < 0 else [args.res]

    stats: dict[str, Any] = {
        "seen": 0,
        "considered": 0,
        "checked": 0,
        "skip_no_mode_tie": 0,
        "skip_non_deg2_min_support": 0,
        "skip_short_B": 0,
        "identity_fail": 0,
        "combined_neg": 0,
        "mismatch_neg": 0,
        "rise_fail": 0,
        "easy_regime": 0,
        "hard_regime": 0,
        "hard_ratio_fail": 0,
        "hard_ratio_undef": 0,
        "shift0_mismatch_neg": 0,
        "shift1_mismatch_neg": 0,
        "shift_other_mismatch_neg": 0,
        "min_combined": None,
        "min_combined_witness": None,
        "min_rise_minus_neg": None,
        "min_rise_minus_neg_witness": None,
        "min_hard_ratio_gap": None,
        "min_hard_ratio_gap_witness": None,
        "max_hard_transfer_ratio": None,
        "max_hard_transfer_ratio_witness": None,
        "min_hard_need_ratio": None,
        "min_hard_need_ratio_witness": None,
        "mismatch_neg_witnesses": [],
    }

    params = {
        "min_n": args.min_n,
        "max_n": args.max_n,
        "mod": args.mod,
        "res": args.res,
        "tol": args.tol,
        "canonical_leaf": "min_support_degree_then_leaf_id",
    }

    t_all = time.time()
    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        n_seen = n_cons = n_checked = n_neg = 0

        for res in residues:
            cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
            if args.mod > 1:
                cmd.append(f"{res}/{args.mod}")

            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert proc.stdout is not None

            for raw in proc.stdout:
                stats["seen"] += 1
                n_seen += 1

                nn, adj = parse_graph6(raw)
                if not is_dleaf_le_1(nn, adj):
                    continue

                stats["considered"] += 1
                n_cons += 1

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
                n_checked += 1

                b_adj = remove_vertices(adj, {leaf, support})
                b_poly = independence_poly(len(b_adj), b_adj)
                if m - 2 < 0 or m >= len(b_poly) or m - 1 >= len(b_poly) or b_poly[m - 1] == 0:
                    stats["skip_short_B"] += 1
                    continue

                u = [x for x in adj[support] if x != leaf][0]
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
                db = b1 - b0
                dq = q1 - q0
                rise = b1 * db
                rise_minus_neg = rise - neg

                rhs = p1 * db + b1 * dq
                if rise_minus_neg != rhs:
                    stats["identity_fail"] += 1

                mode_b = mode_index_leftmost(b_poly)
                shift = mode_b - (m - 1)

                witness = {
                    "n": nn,
                    "g6": raw.decode("ascii").strip(),
                    "m": m,
                    "mode_B": mode_b,
                    "shift": shift,
                    "leaf": leaf,
                    "support": support,
                    "p0": p0,
                    "p1": p1,
                    "q0": q0,
                    "q1": q1,
                    "b0": b0,
                    "b1": b1,
                    "b2": b2,
                    "mismatch": mismatch,
                    "lc_surplus": lc_surplus,
                    "combined": combined,
                    "neg": neg,
                    "rise": rise,
                    "db": db,
                    "dq": dq,
                    "rise_minus_neg": rise_minus_neg,
                }

                maybe_min(stats, "min_combined", "min_combined_witness", float(combined), witness)

                if mismatch < -args.tol:
                    stats["mismatch_neg"] += 1
                    n_neg += 1
                    stats["mismatch_neg_witnesses"].append(witness)

                    if shift == 0:
                        stats["shift0_mismatch_neg"] += 1
                    elif shift == 1:
                        stats["shift1_mismatch_neg"] += 1
                    else:
                        stats["shift_other_mismatch_neg"] += 1

                    maybe_min(
                        stats,
                        "min_rise_minus_neg",
                        "min_rise_minus_neg_witness",
                        float(rise_minus_neg),
                        witness,
                    )

                    if rise_minus_neg < -args.tol:
                        stats["rise_fail"] += 1

                    if dq >= -args.tol:
                        stats["easy_regime"] += 1
                    else:
                        stats["hard_regime"] += 1
                        if db <= args.tol or b1 <= 0:
                            stats["hard_ratio_undef"] += 1
                        else:
                            transfer_ratio = (-dq) / db
                            need_ratio = p1 / b1
                            ratio_gap = need_ratio - transfer_ratio
                            w_ratio = {
                                **witness,
                                "transfer_ratio": transfer_ratio,
                                "need_ratio": need_ratio,
                                "ratio_gap": ratio_gap,
                            }
                            maybe_min(
                                stats,
                                "min_hard_ratio_gap",
                                "min_hard_ratio_gap_witness",
                                float(ratio_gap),
                                w_ratio,
                            )
                            maybe_max(
                                stats,
                                "max_hard_transfer_ratio",
                                "max_hard_transfer_ratio_witness",
                                float(transfer_ratio),
                                w_ratio,
                            )
                            maybe_min(
                                stats,
                                "min_hard_need_ratio",
                                "min_hard_need_ratio_witness",
                                float(need_ratio),
                                w_ratio,
                            )
                            if ratio_gap < -args.tol:
                                stats["hard_ratio_fail"] += 1

                if combined < -args.tol:
                    stats["combined_neg"] += 1

            stderr = proc.stderr.read().decode("utf-8", errors="replace")
            ret = proc.wait()
            if ret != 0:
                raise RuntimeError(
                    f"geng failed at n={n}, res={res}/{args.mod}, returncode={ret}, stderr={stderr.strip()}"
                )

        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_cons:8d} checked={n_checked:8d} "
            f"mismatch_neg={n_neg:4d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    payload = {
        "params": params,
        "summary": {
            k: stats[k]
            for k in stats
            if k != "mismatch_neg_witnesses"
        },
        "mismatch_neg_witnesses": stats["mismatch_neg_witnesses"],
        "wall_s": time.time() - t_all,
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print("-" * 96)
    print(
        f"TOTAL seen={stats['seen']:,} considered={stats['considered']:,} checked={stats['checked']:,} "
        f"mismatch_neg={stats['mismatch_neg']:,} hard={stats['hard_regime']:,} "
        f"rise_fail={stats['rise_fail']:,} combined_neg={stats['combined_neg']:,} "
        f"identity_fail={stats['identity_fail']:,} wall={time.time()-t_all:.1f}s"
    )
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
