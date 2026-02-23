#!/usr/bin/env python3
"""Direct coefficient attack for mode(P) >= m-1.

Setup (canonical degree-2 bridge choice):
  - T is a d_leaf <= 1 tree.
  - Choose a leaf l whose support s has minimum support degree;
    require deg(s)=2 and let u be the other neighbor of s.
  - B = T - {l,s}, rooted at u.
  - P = dp_B[u][0], Q = dp_B[u][1].
  - I(T) = (1+2x)P + (1+x)Q.

At mode index m of I(T), define:
  g = p_{m-1} - p_{m-2}     (target: g >= 0).

Mode inequalities and decompositions:
  I1 := a_m - a_{m-1}
      = (p_m - p_{m-2}) + g + (q_m - q_{m-2})
      = p_term1 + q_term1.

  I2 := a_m - a_{m+1}
      = (2p_{m-1} - p_m - p_{m+1}) + (q_{m-1} - q_{m+1})
      = p_term2 + q_term2.

This script records:
  1) exact values of these terms per tree (aggregated),
  2) compensation patterns (whether Q-term rescues a negative P-term, etc.),
  3) tight witnesses for g, I1, I2, and compensation ratios.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, remove_vertices
from graph6 import parse_graph6
from indpoly import _polyadd, independence_poly, is_log_concave, is_unimodal


def coeff(poly: list[int], k: int) -> int:
    if k < 0 or k >= len(poly):
        return 0
    return poly[k]


def mode_index_leftmost(poly: list[int]) -> int:
    return max(range(len(poly)), key=lambda i: poly[i])


def choose_canonical_leaf(adj: list[list[int]]) -> int | None:
    """Choose canonical leaf as in prior bridge diagnostics.

    Among all leaves, choose one with minimum support degree; tie by leaf id.
    For d_leaf <= 1 this minimum should be 2. If it is not 2, return None.
    """
    deg = [len(nb) for nb in adj]
    leaves = [v for v, d in enumerate(deg) if d == 1]
    if not leaves:
        return None
    min_parent_deg = min(deg[adj[l][0]] for l in leaves)
    leaf = min(l for l in leaves if deg[adj[l][0]] == min_parent_deg)
    if deg[adj[leaf][0]] != 2:
        return None
    return leaf


def poly_from_pq(p_poly: list[int], q_poly: list[int]) -> list[int]:
    """Build I(T) from I(T)=(1+2x)P+(1+x)Q."""
    twox_p = [0] + [2 * c for c in p_poly]
    x_q = [0] + q_poly
    return _polyadd(_polyadd(p_poly, twox_p), _polyadd(q_poly, x_q))


def update_min_witness(
    container: dict[str, Any],
    key: str,
    value: float,
    witness: dict[str, Any],
) -> None:
    cur = container.get(key)
    if cur is None or value < cur:
        container[key] = value
        container[f"{key}_witness"] = witness


def update_max_witness(
    container: dict[str, Any],
    key: str,
    value: float,
    witness: dict[str, Any],
) -> None:
    cur = container.get(key)
    if cur is None or value > cur:
        container[key] = value
        container[f"{key}_witness"] = witness


def witness_payload(
    n: int,
    g6: str,
    m: int,
    leaf: int,
    support: int,
    u: int,
    p_m2: int,
    p_m1: int,
    p_m: int,
    p_m1p: int,
    q_m2: int,
    q_m1: int,
    q_m: int,
    q_m1p: int,
    a_m1: int,
    a_m: int,
    a_m1p: int,
    g: int,
    i1: int,
    i2: int,
    p_term1: int,
    q_term1: int,
    p_term2: int,
    q_term2: int,
    half_margin: int,
) -> dict[str, Any]:
    return {
        "n": n,
        "g6": g6,
        "m": m,
        "leaf": leaf,
        "support": support,
        "u": u,
        "p_m_minus_2": p_m2,
        "p_m_minus_1": p_m1,
        "p_m": p_m,
        "p_m_plus_1": p_m1p,
        "q_m_minus_2": q_m2,
        "q_m_minus_1": q_m1,
        "q_m": q_m,
        "q_m_plus_1": q_m1p,
        "a_m_minus_1": a_m1,
        "a_m": a_m,
        "a_m_plus_1": a_m1p,
        "g_target": g,
        "I1_a_m_minus_a_m_minus_1": i1,
        "I2_a_m_minus_a_m_plus_1": i2,
        "p_term1": p_term1,
        "q_term1": q_term1,
        "p_term2": p_term2,
        "q_term2": q_term2,
        "half_margin_2p_m_minus_1_minus_a_m_minus_1": half_margin,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument(
        "--out",
        default="results/attack2_coeff_scan_n23.json",
        help="JSON output path",
    )
    ap.add_argument(
        "--verify-identity-every",
        type=int,
        default=0,
        help=(
            "If >0, every k-th checked tree verifies I(T) directly via "
            "independence_poly(T)."
        ),
    )
    ap.add_argument("--tol", type=float, default=1e-12)
    args = ap.parse_args()

    summary: dict[str, Any] = {
        "seen": 0,
        "considered_dleaf_le_1": 0,
        "checked": 0,
        "no_deg2_support_canonical": 0,
        "mode_p_fail": 0,
        "mode_p_equal": 0,
        "i1_fail": 0,
        "i2_fail": 0,
        "i1_identity_fail": 0,
        "i2_identity_fail": 0,
        "p_ge_pprime_fail": 0,
        "q_shift_fail": 0,
        "p_log_concave_fail": 0,
        "q_log_concave_fail": 0,
        "p_unimodal_fail": 0,
        "half_condition_true": 0,
        "half_condition_false": 0,
        "p_term1_negative": 0,
        "q_term1_negative": 0,
        "p_term2_negative": 0,
        "q_term2_negative": 0,
        "i1_q_compensates_negative_p": 0,
        "i1_p_compensates_negative_q": 0,
        "i2_q_compensates_negative_p": 0,
        "i2_p_compensates_negative_q": 0,
        "sampled_identity_checks": 0,
        "sampled_identity_fail": 0,
    }

    minima: dict[str, Any] = {
        "min_g": None,
        "min_i1": None,
        "min_i2": None,
        "min_half_margin": None,
        "min_q_comp_ratio_i1": None,
        "min_q_comp_ratio_i2": None,
    }
    maxima: dict[str, Any] = {
        "max_required_q_from_i1": None,
        "max_required_q_from_i2": None,
        "max_negative_p_term1": None,
        "max_negative_p_term2": None,
    }

    per_n: dict[str, dict[str, int]] = {}
    t0 = time.time()

    for n in range(args.min_n, args.max_n + 1):
        n_stats = {
            "seen": 0,
            "considered": 0,
            "checked": 0,
            "mode_p_fail": 0,
            "mode_p_equal": 0,
            "i1_fail": 0,
            "i2_fail": 0,
        }
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        assert proc.stdout is not None

        for raw in proc.stdout:
            g6 = raw.decode("ascii").strip()
            nn, adj = parse_graph6(raw)
            summary["seen"] += 1
            n_stats["seen"] += 1
            if not is_dleaf_le_1(nn, adj):
                continue

            summary["considered_dleaf_le_1"] += 1
            n_stats["considered"] += 1
            leaf = choose_canonical_leaf(adj)
            if leaf is None:
                summary["no_deg2_support_canonical"] += 1
                continue

            support = adj[leaf][0]
            u = next(v for v in adj[support] if v != leaf)
            b_adj = remove_vertices(adj, {leaf, support})
            keep = [v for v in range(nn) if v not in {leaf, support}]
            idx = {v: i for i, v in enumerate(keep)}
            u_in_b = idx[u]

            p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)
            t_poly = poly_from_pq(p_poly, q_poly)
            m = mode_index_leftmost(t_poly)
            if m == 0:
                continue

            summary["checked"] += 1
            n_stats["checked"] += 1

            # Optional expensive identity check against direct independence DP.
            if (
                args.verify_identity_every > 0
                and summary["checked"] % args.verify_identity_every == 0
            ):
                summary["sampled_identity_checks"] += 1
                poly_true = independence_poly(nn, adj)
                if poly_true != t_poly:
                    summary["sampled_identity_fail"] += 1

            p_m2 = coeff(p_poly, m - 2)
            p_m1 = coeff(p_poly, m - 1)
            p_m = coeff(p_poly, m)
            p_m1p = coeff(p_poly, m + 1)

            q_m2 = coeff(q_poly, m - 2)
            q_m1 = coeff(q_poly, m - 1)
            q_m = coeff(q_poly, m)
            q_m1p = coeff(q_poly, m + 1)

            a_m2 = coeff(t_poly, m - 2)
            a_m1 = coeff(t_poly, m - 1)
            a_m = coeff(t_poly, m)
            a_m1p = coeff(t_poly, m + 1)

            g = p_m1 - p_m2
            i1 = a_m - a_m1
            i2 = a_m - a_m1p
            p_term1 = (p_m - p_m2) + g
            q_term1 = q_m - q_m2
            p_term2 = 2 * p_m1 - p_m - p_m1p
            q_term2 = q_m1 - q_m1p
            half_margin = 2 * p_m1 - a_m1

            witness = witness_payload(
                n=n,
                g6=g6,
                m=m,
                leaf=leaf,
                support=support,
                u=u,
                p_m2=p_m2,
                p_m1=p_m1,
                p_m=p_m,
                p_m1p=p_m1p,
                q_m2=q_m2,
                q_m1=q_m1,
                q_m=q_m,
                q_m1p=q_m1p,
                a_m1=a_m1,
                a_m=a_m,
                a_m1p=a_m1p,
                g=g,
                i1=i1,
                i2=i2,
                p_term1=p_term1,
                q_term1=q_term1,
                p_term2=p_term2,
                q_term2=q_term2,
                half_margin=half_margin,
            )

            # Identity checks for the mode inequalities at m.
            if i1 != p_term1 + q_term1:
                summary["i1_identity_fail"] += 1
            if i2 != p_term2 + q_term2:
                summary["i2_identity_fail"] += 1

            if i1 < 0:
                summary["i1_fail"] += 1
                n_stats["i1_fail"] += 1
            if i2 < 0:
                summary["i2_fail"] += 1
                n_stats["i2_fail"] += 1

            if g < 0:
                summary["mode_p_fail"] += 1
                n_stats["mode_p_fail"] += 1
            elif g == 0:
                summary["mode_p_equal"] += 1
                n_stats["mode_p_equal"] += 1

            pprime = q_poly[1:] if len(q_poly) >= 2 else []
            # q_k = pprime_{k-1}
            for k in range(len(q_poly)):
                if coeff(q_poly, k) != coeff(pprime, k - 1):
                    summary["q_shift_fail"] += 1
                    break

            max_len = max(len(p_poly), len(pprime))
            for k in range(max_len):
                if coeff(p_poly, k) < coeff(pprime, k):
                    summary["p_ge_pprime_fail"] += 1
                    break

            if not is_log_concave(p_poly):
                summary["p_log_concave_fail"] += 1
            if not is_log_concave(q_poly):
                summary["q_log_concave_fail"] += 1
            if not is_unimodal(p_poly):
                summary["p_unimodal_fail"] += 1

            if half_margin >= 0:
                summary["half_condition_true"] += 1
            else:
                summary["half_condition_false"] += 1

            if p_term1 < 0:
                summary["p_term1_negative"] += 1
            if q_term1 < 0:
                summary["q_term1_negative"] += 1
            if p_term2 < 0:
                summary["p_term2_negative"] += 1
            if q_term2 < 0:
                summary["q_term2_negative"] += 1

            if p_term1 < 0 and q_term1 >= -p_term1:
                summary["i1_q_compensates_negative_p"] += 1
                ratio = q_term1 / (-p_term1) if p_term1 != 0 else float("inf")
                update_min_witness(
                    minima,
                    "min_q_comp_ratio_i1",
                    ratio,
                    witness,
                )
            if q_term1 < 0 and p_term1 >= -q_term1:
                summary["i1_p_compensates_negative_q"] += 1

            if p_term2 < 0 and q_term2 >= -p_term2:
                summary["i2_q_compensates_negative_p"] += 1
                ratio = q_term2 / (-p_term2) if p_term2 != 0 else float("inf")
                update_min_witness(
                    minima,
                    "min_q_comp_ratio_i2",
                    ratio,
                    witness,
                )
            if q_term2 < 0 and p_term2 >= -q_term2:
                summary["i2_p_compensates_negative_q"] += 1

            req_q_i1 = max(0, -(g + (p_m - p_m2)))
            req_q_i2 = max(0, -p_term2)
            update_max_witness(maxima, "max_required_q_from_i1", float(req_q_i1), witness)
            update_max_witness(maxima, "max_required_q_from_i2", float(req_q_i2), witness)
            update_max_witness(maxima, "max_negative_p_term1", float(-p_term1), witness)
            update_max_witness(maxima, "max_negative_p_term2", float(-p_term2), witness)

            update_min_witness(minima, "min_g", float(g), witness)
            update_min_witness(minima, "min_i1", float(i1), witness)
            update_min_witness(minima, "min_i2", float(i2), witness)
            update_min_witness(minima, "min_half_margin", float(half_margin), witness)

        proc.wait()
        per_n[str(n)] = n_stats
        print(
            (
                f"n={n:2d}: seen={n_stats['seen']:8d} "
                f"considered={n_stats['considered']:8d} "
                f"checked={n_stats['checked']:8d} "
                f"modeP_fail={n_stats['mode_p_fail']:6d} "
                f"modeP_eq={n_stats['mode_p_equal']:6d}"
            ),
            flush=True,
        )

    wall = time.time() - t0
    summary["wall_seconds"] = wall

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "tol": args.tol,
            "geng": args.geng,
            "verify_identity_every": args.verify_identity_every,
        },
        "summary": summary,
        "minima": minima,
        "maxima": maxima,
        "per_n": per_n,
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print("-" * 100, flush=True)
    print(f"Checked: {summary['checked']}", flush=True)
    print(f"mode(P) >= m-1 failures: {summary['mode_p_fail']}", flush=True)
    print(f"min g = p_(m-1)-p_(m-2): {minima['min_g']}", flush=True)
    print(f"I1 failures: {summary['i1_fail']}", flush=True)
    print(f"I2 failures: {summary['i2_fail']}", flush=True)
    print(f"half condition true count: {summary['half_condition_true']}", flush=True)
    print(f"Output written to: {args.out}", flush=True)
    print(f"Wall seconds: {wall:.2f}", flush=True)


if __name__ == "__main__":
    main()

