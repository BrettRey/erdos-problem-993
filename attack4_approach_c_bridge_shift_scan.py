#!/usr/bin/env python3
"""Approach C diagnostics: bridge-shift inequalities around index m.

Given I(T) = (1+2x)P + (1+x)Q and m = mode(I(T)), define

  Delta_I = I_m - I_{m-1}
  Delta_P = (p_m - p_{m-1}) + 2(p_{m-1} - p_{m-2})
  Delta_Q = q_m - q_{m-2}

Then Delta_I = Delta_P + Delta_Q and Delta_I >= 0.

This script tests candidate inequalities that could force p_{m-1} >= p_{m-2}
(and thus mode(P) >= m-1 under unimodality/log-concavity heuristics).
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


def update_max(stats: dict[str, Any], key: str, value: int, witness: dict[str, Any]) -> None:
    cur = stats.get(key)
    if cur is None or value > cur:
        stats[key] = value
        stats[key + "_witness"] = witness


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    stats: dict[str, Any] = {
        "checked": 0,
        "identity_fail": 0,
        "modeP_ge_mminus1_fail": 0,
        # sign profiles
        "deltaP_neg": 0,
        "deltaP_nonneg": 0,
        "deltaQ_pos": 0,
        "deltaQ_nonpos": 0,
        # candidate claims
        "claim_deltaP_nonneg_fail": 0,
        "claim_deltaQ_nonpos_fail": 0,
        "claim_deltaQ_le_pm1_minus_pm_fail": 0,
        "claim_deltaQ_le_pm1_minus_pm_plus_pm1_minus_pm2_fail": 0,
        # extremal margins
        "min_deltaI": None,
        "min_deltaP": None,
        "max_deltaQ": None,
        "min_deltaQ": None,
        "min_gap_deltaQ_bound1": None,
        "min_gap_deltaQ_bound2": None,
        "wall_s": 0.0,
    }

    t0 = time.time()

    for n, adj, g6 in iter_tree_g6(args.min_n, args.max_n, args.geng):
        decomp = bridge_decomposition(n=n, adj=adj, g6=g6, require_dleaf=True)
        if decomp is None:
            continue

        stats["checked"] += 1

        m = decomp.m_t
        p2 = get_coeff(decomp.p_poly, m - 2)
        p1 = get_coeff(decomp.p_poly, m - 1)
        p0 = get_coeff(decomp.p_poly, m)

        q2 = get_coeff(decomp.q_poly, m - 2)
        q0 = get_coeff(decomp.q_poly, m)

        i_m = get_coeff(decomp.poly_t, m)
        i_m1 = get_coeff(decomp.poly_t, m - 1)

        delta_i = i_m - i_m1
        delta_p = (p0 - p1) + 2 * (p1 - p2)
        delta_q = q0 - q2

        if delta_i != delta_p + delta_q:
            stats["identity_fail"] += 1

        mode_p = mode_index_leftmost(decomp.p_poly)
        if mode_p < m - 1:
            stats["modeP_ge_mminus1_fail"] += 1

        if delta_p < 0:
            stats["deltaP_neg"] += 1
            stats["claim_deltaP_nonneg_fail"] += 1
        else:
            stats["deltaP_nonneg"] += 1

        if delta_q > 0:
            stats["deltaQ_pos"] += 1
            stats["claim_deltaQ_nonpos_fail"] += 1
        else:
            stats["deltaQ_nonpos"] += 1

        # Candidate upper bounds on Delta_Q that would force p_{m-1} >= p_{m-2}.
        bound1 = p1 - p0
        bound2 = (p1 - p0) + (p1 - p2)

        gap1 = bound1 - delta_q
        gap2 = bound2 - delta_q

        if gap1 < 0:
            stats["claim_deltaQ_le_pm1_minus_pm_fail"] += 1
        if gap2 < 0:
            stats["claim_deltaQ_le_pm1_minus_pm_plus_pm1_minus_pm2_fail"] += 1

        witness = {
            "n": n,
            "g6": g6,
            "m": m,
            "delta_i": delta_i,
            "delta_p": delta_p,
            "delta_q": delta_q,
            "p_m2": p2,
            "p_m1": p1,
            "p_m": p0,
            "q_m2": q2,
            "q_m": q0,
        }

        update_min(stats, "min_deltaI", delta_i, witness)
        update_min(stats, "min_deltaP", delta_p, witness)
        update_max(stats, "max_deltaQ", delta_q, witness)
        update_min(stats, "min_deltaQ", delta_q, witness)
        update_min(stats, "min_gap_deltaQ_bound1", gap1, witness)
        update_min(stats, "min_gap_deltaQ_bound2", gap2, witness)

    stats["wall_s"] = time.time() - t0

    payload = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "geng": args.geng,
        },
        "stats": stats,
    }

    print(json.dumps(payload, indent=2, sort_keys=True))

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
