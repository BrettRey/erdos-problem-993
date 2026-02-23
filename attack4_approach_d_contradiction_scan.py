#!/usr/bin/env python3
"""Approach D diagnostics: contradiction from assuming mode(P) <= m-2.

Real-tree part (d_leaf<=1, n-range):
- Detect cases with mode(P) <= m-2.
- For each such case, record the exact compensation inequality
    q_m - q_{m-2} >= (p_{m-1}-p_m) + 2(p_{m-2}-p_{m-1}).

Synthetic part (bounded exhaustive search):
- Test whether extracted local constraints are logically inconsistent.
- If feasible tuples exist, then the local inequalities alone do not imply a
  contradiction; extra structure is needed.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from typing import Any

from attack4_common import bridge_decomposition, get_coeff, iter_tree_g6
from diagnose_bridge_decomposition import mode_index_leftmost


def update_min(stats: dict[str, Any], key: str, value: int, witness: dict[str, Any]) -> None:
    cur = stats.get(key)
    if cur is None or value < cur:
        stats[key] = value
        stats[key + "_witness"] = witness


def run_real_scan(min_n: int, max_n: int, geng: str, all_trees: bool) -> dict[str, Any]:
    stats: dict[str, Any] = {
        "checked": 0,
        "assumption_modeP_le_mminus2": 0,
        "assumption_and_igap_fail": 0,
        "assumption_min_required_qgain": None,
        "assumption_min_qgain_minus_required": None,
        "assumption_min_gap_modeP_minus_mminus1": None,
    }

    for n, adj, g6 in iter_tree_g6(min_n, max_n, geng):
        decomp = bridge_decomposition(n=n, adj=adj, g6=g6, require_dleaf=not all_trees)
        if decomp is None:
            continue

        stats["checked"] += 1

        m = decomp.m_t
        mode_p = mode_index_leftmost(decomp.p_poly)
        gap_mode = mode_p - (m - 1)
        update_min(
            stats,
            "assumption_min_gap_modeP_minus_mminus1",
            gap_mode,
            {
                "n": n,
                "g6": g6,
                "m": m,
                "mode_p": mode_p,
            },
        )

        if mode_p > m - 2:
            continue

        stats["assumption_modeP_le_mminus2"] += 1

        p2 = get_coeff(decomp.p_poly, m - 2)
        p1 = get_coeff(decomp.p_poly, m - 1)
        p0 = get_coeff(decomp.p_poly, m)

        q2 = get_coeff(decomp.q_poly, m - 2)
        q0 = get_coeff(decomp.q_poly, m)

        # Exact compensation inequality under the assumption.
        required = (p1 - p0) + 2 * (p2 - p1)
        qgain = q0 - q2
        margin = qgain - required

        update_min(
            stats,
            "assumption_min_required_qgain",
            required,
            {
                "n": n,
                "g6": g6,
                "m": m,
                "p_m2": p2,
                "p_m1": p1,
                "p_m": p0,
                "required": required,
            },
        )
        update_min(
            stats,
            "assumption_min_qgain_minus_required",
            margin,
            {
                "n": n,
                "g6": g6,
                "m": m,
                "p_m2": p2,
                "p_m1": p1,
                "p_m": p0,
                "q_m2": q2,
                "q_m": q0,
                "required": required,
                "qgain": qgain,
                "margin": margin,
            },
        )

        # Should never fail because m is a mode of I(T), but keep explicit check.
        i_m = get_coeff(decomp.poly_t, m)
        i_m1 = get_coeff(decomp.poly_t, m - 1)
        if i_m < i_m1:
            stats["assumption_and_igap_fail"] += 1

    return stats


def synthetic_minimal_feasible(
    max_coeff: int,
) -> tuple[int, dict[str, int] | None, dict[str, int] | None]:
    """Model A: direct local variables (p2,p1,p0,q2,q0)."""
    found = 0
    first: dict[str, int] | None = None
    first_nontrivial: dict[str, int] | None = None

    for p2 in range(max_coeff + 1):
        for p1 in range(p2 + 1):
            for p0 in range(p1 + 1):
                required = (p1 - p0) + 2 * (p2 - p1)
                for q2 in range(max_coeff + 1):
                    for q0 in range(max_coeff + 1):
                        qgain = q0 - q2
                        if qgain < required:
                            continue
                        found += 1
                        if first is None:
                            first = {
                                "p_m2": p2,
                                "p_m1": p1,
                                "p_m": p0,
                                "q_m2": q2,
                                "q_m": q0,
                                "required": required,
                                "qgain": qgain,
                                "margin": qgain - required,
                            }
                        if first_nontrivial is None and (
                            required > 0
                            or p2 > 0
                            or p1 > 0
                            or p0 > 0
                            or q2 > 0
                            or q0 > 0
                        ):
                            first_nontrivial = {
                                "p_m2": p2,
                                "p_m1": p1,
                                "p_m": p0,
                                "q_m2": q2,
                                "q_m": q0,
                                "required": required,
                                "qgain": qgain,
                                "margin": qgain - required,
                            }
    return found, first, first_nontrivial


def synthetic_structured_feasible(max_coeff: int) -> dict[str, Any]:
    """Model B/C: local G+E structure with optional local LC/unimodality checks."""
    stats: dict[str, Any] = {
        "model_B_found": 0,
        "model_B_first": None,
        "model_B_first_nontrivial": None,
        "model_C_found": 0,
        "model_C_first": None,
        "model_C_first_nontrivial": None,
    }

    # Variables:
    # g0=g_{m-3}, g1=g_{m-2}, g2=g_{m-1}, g3=g_m
    # e1=E_{m-2}, e2=E_{m-1}, e3=E_m
    # p2=g1+e1, p1=g2+e2, p0=g3+e3
    for g0 in range(max_coeff + 1):
        for g1 in range(max_coeff + 1):
            for g2 in range(max_coeff + 1):
                # Local log-concavity for G at m-2: g1^2 >= g0*g2
                if g1 * g1 < g0 * g2:
                    continue
                for g3 in range(max_coeff + 1):
                    # Local log-concavity for G at m-1: g2^2 >= g1*g3
                    if g2 * g2 < g1 * g3:
                        continue

                    for e1 in range(max_coeff + 1):
                        p2 = g1 + e1
                        for e2 in range(max_coeff + 1):
                            p1 = g2 + e2
                            if p1 > p2:
                                continue
                            for e3 in range(max_coeff + 1):
                                p0 = g3 + e3
                                if p0 > p1:
                                    continue

                                # Local log-concavity for P at m-1.
                                if p1 * p1 < p2 * p0:
                                    continue

                                # Compensation inequality from I_m >= I_{m-1}.
                                required = (p1 - p0) + 2 * (p2 - p1)
                                qgain = g2 - g0
                                if qgain < required:
                                    continue

                                stats["model_B_found"] += 1
                                if stats["model_B_first"] is None:
                                    stats["model_B_first"] = {
                                        "g_m3": g0,
                                        "g_m2": g1,
                                        "g_m1": g2,
                                        "g_m": g3,
                                        "e_m2": e1,
                                        "e_m1": e2,
                                        "e_m": e3,
                                        "p_m2": p2,
                                        "p_m1": p1,
                                        "p_m": p0,
                                        "required": required,
                                        "qgain": qgain,
                                        "margin": qgain - required,
                                    }
                                if stats["model_B_first_nontrivial"] is None and (
                                    required > 0
                                    or g0 > 0
                                    or g1 > 0
                                    or g2 > 0
                                    or g3 > 0
                                    or e1 > 0
                                    or e2 > 0
                                    or e3 > 0
                                ):
                                    stats["model_B_first_nontrivial"] = {
                                        "g_m3": g0,
                                        "g_m2": g1,
                                        "g_m1": g2,
                                        "g_m": g3,
                                        "e_m2": e1,
                                        "e_m1": e2,
                                        "e_m": e3,
                                        "p_m2": p2,
                                        "p_m1": p1,
                                        "p_m": p0,
                                        "required": required,
                                        "qgain": qgain,
                                        "margin": qgain - required,
                                    }

                                # Model C adds local unimodality-at-m-1 for Q:
                                # q_{m-1} >= q_{m-2}, q_{m-1} >= q_m, i.e. g1>=g0 and g1>=g2.
                                if g1 >= g0 and g1 >= g2:
                                    stats["model_C_found"] += 1
                                    if stats["model_C_first"] is None:
                                        stats["model_C_first"] = {
                                            "g_m3": g0,
                                            "g_m2": g1,
                                            "g_m1": g2,
                                            "g_m": g3,
                                            "e_m2": e1,
                                            "e_m1": e2,
                                            "e_m": e3,
                                            "p_m2": p2,
                                            "p_m1": p1,
                                            "p_m": p0,
                                            "required": required,
                                            "qgain": qgain,
                                            "margin": qgain - required,
                                        }
                                    if stats["model_C_first_nontrivial"] is None and (
                                        required > 0
                                        or g0 > 0
                                        or g1 > 0
                                        or g2 > 0
                                        or g3 > 0
                                        or e1 > 0
                                        or e2 > 0
                                        or e3 > 0
                                    ):
                                        stats["model_C_first_nontrivial"] = {
                                            "g_m3": g0,
                                            "g_m2": g1,
                                            "g_m1": g2,
                                            "g_m": g3,
                                            "e_m2": e1,
                                            "e_m1": e2,
                                            "e_m": e3,
                                            "p_m2": p2,
                                            "p_m1": p1,
                                            "p_m": p0,
                                            "required": required,
                                            "qgain": qgain,
                                            "margin": qgain - required,
                                        }

    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument(
        "--all-trees",
        action="store_true",
        help="Disable d_leaf<=1 filter for the real-tree scan.",
    )
    ap.add_argument(
        "--synthetic-max-coeff",
        type=int,
        default=6,
        help="Upper bound for synthetic coefficient search.",
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    t0 = time.time()

    real_stats = run_real_scan(
        min_n=args.min_n,
        max_n=args.max_n,
        geng=args.geng,
        all_trees=args.all_trees,
    )

    model_a_count, model_a_first, model_a_first_nontrivial = synthetic_minimal_feasible(
        args.synthetic_max_coeff
    )
    structured = synthetic_structured_feasible(args.synthetic_max_coeff)

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "geng": args.geng,
            "all_trees": args.all_trees,
            "synthetic_max_coeff": args.synthetic_max_coeff,
        },
        "real_scan": real_stats,
        "synthetic": {
            "model_A_found": model_a_count,
            "model_A_first": model_a_first,
            "model_A_first_nontrivial": model_a_first_nontrivial,
            **structured,
        },
        "wall_s": time.time() - t0,
    }

    print(json.dumps(payload, indent=2, sort_keys=True))

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
