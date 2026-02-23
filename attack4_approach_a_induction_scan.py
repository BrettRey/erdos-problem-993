#!/usr/bin/env python3
"""Approach A diagnostics: induction on n via B = T - {l,s}.

Key question: can mode(P_T) >= m_T - 1 be inherited from B?

For T in d_leaf<=1:
1) Build canonical degree-2 bridge decomposition for T.
2) Form B = T-{l,s}.
3) Check closure (is B also d_leaf<=1?).
4) On closure + non-base cases, test transfer relations:
   - m_B vs m_T
   - mode(P_T) vs m_B - 1
   - mode(P_T) vs mode(P_B) when B has its own decomposition.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter
from typing import Any

from attack4_common import bridge_decomposition, iter_tree_g6
from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices


def update_min(stats: dict[str, Any], key: str, value: int, witness: dict[str, Any]) -> None:
    cur = stats.get(key)
    if cur is None or value < cur:
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
        "checked_T": 0,
        "bound_modeP_ge_mminus1_fail": 0,
        "B_is_dleaf": 0,
        "B_not_dleaf": 0,
        "B_n_distribution": {},
        "closure_nonbase": 0,   # B size >= 4 and d_leaf<=1
        "closure_base": 0,      # B size <= 3 and d_leaf<=1
        "B_decomp_missing_in_closure_nonbase": 0,
        "B_has_deg2_support_any": 0,
        "B_no_deg2_support_any": 0,
        "B_any_modeP_ge_mminus1_fail": 0,
        "B_any_existential_choice_fail": 0,
        "B_any_no_deg2_choice": 0,
        "claim_mB_ge_mT_minus1_fail": 0,
        "claim_modePT_ge_mB_minus1_fail": 0,
        "claim_modePT_ge_modePB_fail": 0,
        "delta_mB_minus_mT_distribution": {},
        "min_gap_modePT_minus_mTminus1": None,
        "min_gap_modePT_minus_mBminus1": None,
        "min_gap_modePT_minus_modePB": None,
        "min_best_gap_modePBexist_minus_mBminus1": None,
        "wall_s": 0.0,
    }

    t0 = time.time()

    for n, adj, g6 in iter_tree_g6(args.min_n, args.max_n, args.geng):
        decomp_t = bridge_decomposition(n=n, adj=adj, g6=g6, require_dleaf=True)
        if decomp_t is None:
            continue

        stats["checked_T"] += 1
        m_t = decomp_t.m_t
        mode_p_t = mode_index_leftmost(decomp_t.p_poly)

        gap_t = mode_p_t - (m_t - 1)
        update_min(
            stats,
            "min_gap_modePT_minus_mTminus1",
            gap_t,
            {
                "n": n,
                "g6": g6,
                "m_t": m_t,
                "mode_p_t": mode_p_t,
                "leaf": decomp_t.leaf,
                "support": decomp_t.support,
                "u": decomp_t.u,
            },
        )
        if gap_t < 0:
            stats["bound_modeP_ge_mminus1_fail"] += 1

        b_n = len(decomp_t.b_adj)
        stats["B_n_distribution"][str(b_n)] = stats["B_n_distribution"].get(str(b_n), 0) + 1

        if is_dleaf_le_1(b_n, decomp_t.b_adj):
            stats["B_is_dleaf"] += 1
            if b_n <= 3:
                stats["closure_base"] += 1
                continue
            stats["closure_nonbase"] += 1
        else:
            stats["B_not_dleaf"] += 1

        # Analyze B in broader class (no d_leaf filter).
        b_decomp_any = bridge_decomposition(
            n=b_n,
            adj=decomp_t.b_adj,
            g6=f"{g6}|B",
            require_dleaf=False,
        )
        if b_decomp_any is None:
            stats["B_no_deg2_support_any"] += 1
        else:
            stats["B_has_deg2_support_any"] += 1
            m_b_any = b_decomp_any.m_t
            mode_p_b_any = mode_index_leftmost(b_decomp_any.p_poly)
            if mode_p_b_any < m_b_any - 1:
                stats["B_any_modeP_ge_mminus1_fail"] += 1

            # Existential variant: maximize gap over all degree-2 support leaves in B.
            deg_b = [len(nb) for nb in decomp_t.b_adj]
            leaves_b = [v for v in range(b_n) if deg_b[v] == 1]
            m_b_from_poly = mode_index_leftmost(decomp_t.b_poly)
            best_gap = None

            for leaf_b in leaves_b:
                support_b = decomp_t.b_adj[leaf_b][0]
                if deg_b[support_b] != 2:
                    continue

                u_b = (
                    decomp_t.b_adj[support_b][0]
                    if decomp_t.b_adj[support_b][1] == leaf_b
                    else decomp_t.b_adj[support_b][1]
                )
                bb_adj = remove_vertices(decomp_t.b_adj, {leaf_b, support_b})
                keep_b = [v for v in range(b_n) if v not in {leaf_b, support_b}]
                idx_b = {old: new for new, old in enumerate(keep_b)}
                u_in_bb = idx_b[u_b]
                p_choice, _ = compute_hub_polys(bb_adj, u_in_bb)
                gap_choice = mode_index_leftmost(p_choice) - (m_b_from_poly - 1)

                if best_gap is None or gap_choice > best_gap:
                    best_gap = gap_choice

            if best_gap is None:
                stats["B_any_no_deg2_choice"] += 1
            else:
                update_min(
                    stats,
                    "min_best_gap_modePBexist_minus_mBminus1",
                    best_gap,
                    {
                        "n": n,
                        "g6": g6,
                        "m_b": m_b_from_poly,
                        "best_gap": best_gap,
                    },
                )
                if best_gap < 0:
                    stats["B_any_existential_choice_fail"] += 1

        # For actual induction transfer checks, require closure and non-base.
        if not is_dleaf_le_1(b_n, decomp_t.b_adj) or b_n <= 3:
            continue

        b_decomp = bridge_decomposition(
            n=b_n,
            adj=decomp_t.b_adj,
            g6=f"{g6}|B",
            require_dleaf=True,
        )
        if b_decomp is None:
            stats["B_decomp_missing_in_closure_nonbase"] += 1
            continue

        m_b = b_decomp.m_t
        mode_p_b = mode_index_leftmost(b_decomp.p_poly)

        delta = m_b - m_t
        stats["delta_mB_minus_mT_distribution"][str(delta)] = (
            stats["delta_mB_minus_mT_distribution"].get(str(delta), 0) + 1
        )

        if m_b < m_t - 1:
            stats["claim_mB_ge_mT_minus1_fail"] += 1

        gap_mb = mode_p_t - (m_b - 1)
        update_min(
            stats,
            "min_gap_modePT_minus_mBminus1",
            gap_mb,
            {
                "n": n,
                "g6": g6,
                "m_t": m_t,
                "m_b": m_b,
                "mode_p_t": mode_p_t,
                "leaf": decomp_t.leaf,
                "support": decomp_t.support,
                "u": decomp_t.u,
            },
        )
        if gap_mb < 0:
            stats["claim_modePT_ge_mB_minus1_fail"] += 1

        gap_pb = mode_p_t - mode_p_b
        update_min(
            stats,
            "min_gap_modePT_minus_modePB",
            gap_pb,
            {
                "n": n,
                "g6": g6,
                "m_t": m_t,
                "m_b": m_b,
                "mode_p_t": mode_p_t,
                "mode_p_b": mode_p_b,
                "leaf": decomp_t.leaf,
                "support": decomp_t.support,
                "u": decomp_t.u,
            },
        )
        if gap_pb < 0:
            stats["claim_modePT_ge_modePB_fail"] += 1

    stats["wall_s"] = time.time() - t0

    closure_total = stats["B_is_dleaf"]
    if closure_total > 0:
        stats["closure_rate"] = closure_total / stats["checked_T"]
    else:
        stats["closure_rate"] = 0.0

    if stats["closure_nonbase"] > 0:
        stats["closure_nonbase_success_rate"] = (
            stats["closure_nonbase"] - stats["B_decomp_missing_in_closure_nonbase"]
        ) / stats["closure_nonbase"]
    else:
        stats["closure_nonbase_success_rate"] = 0.0

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
