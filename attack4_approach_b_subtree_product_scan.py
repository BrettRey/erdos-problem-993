#!/usr/bin/env python3
"""Approach B diagnostics: subtree decomposition claims for P = prod_i I(T_{c_i}).

Tested intermediate claims:
1) Superadditivity candidate: mode(P) >= sum_i mode(I(T_{c_i})).
2) Subtree-sum target: sum_i mode(I(T_{c_i})) >= m(T) - 1.
3) Size-only proxy: sum_i floor(|T_{c_i}| / 3) >= m(T) - 1.
4) Per-child empirical baseline: mode(I(T_{c_i})) >= floor(|T_{c_i}| / 3).
"""

from __future__ import annotations

import argparse
import json
import os
import time
from typing import Any

from attack4_common import bridge_decomposition, iter_tree_g6
from diagnose_bridge_decomposition import mode_index_leftmost


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
        "checked": 0,
        "superadditivity_fail": 0,
        "sum_child_modes_ge_mminus1_fail": 0,
        "size_proxy_floor_v3_ge_mminus1_fail": 0,
        "per_child_floor_v3_fail": 0,
        "num_children_distribution": {},
        "min_gap_modeP_minus_sumChildModes": None,
        "min_gap_sumChildModes_minus_mminus1": None,
        "min_gap_sumFloorV3_minus_mminus1": None,
        "min_gap_child_mode_minus_floorv3": None,
        "wall_s": 0.0,
    }

    t0 = time.time()

    for n, adj, g6 in iter_tree_g6(args.min_n, args.max_n, args.geng):
        decomp = bridge_decomposition(n=n, adj=adj, g6=g6, require_dleaf=True)
        if decomp is None:
            continue

        stats["checked"] += 1
        m = decomp.m_t

        child_modes = [mode_index_leftmost(poly) for poly in decomp.f_list]
        mode_p = mode_index_leftmost(decomp.p_poly)
        sum_child_modes = sum(child_modes)

        num_children = len(decomp.children_of_u)
        stats["num_children_distribution"][str(num_children)] = (
            stats["num_children_distribution"].get(str(num_children), 0) + 1
        )

        # Claim 1: mode(P) >= sum_i mode(f_i)
        gap_super = mode_p - sum_child_modes
        update_min(
            stats,
            "min_gap_modeP_minus_sumChildModes",
            gap_super,
            {
                "n": n,
                "g6": g6,
                "m": m,
                "mode_p": mode_p,
                "child_modes": child_modes,
                "child_sizes": decomp.child_sizes,
            },
        )
        if gap_super < 0:
            stats["superadditivity_fail"] += 1

        # Claim 2: sum_i mode(f_i) >= m-1
        gap_sum_modes = sum_child_modes - (m - 1)
        update_min(
            stats,
            "min_gap_sumChildModes_minus_mminus1",
            gap_sum_modes,
            {
                "n": n,
                "g6": g6,
                "m": m,
                "sum_child_modes": sum_child_modes,
                "child_modes": child_modes,
                "child_sizes": decomp.child_sizes,
            },
        )
        if gap_sum_modes < 0:
            stats["sum_child_modes_ge_mminus1_fail"] += 1

        # Claim 3: size-only proxy sum floor(v_i/3) >= m-1
        sum_floor_v3 = sum(v // 3 for v in decomp.child_sizes)
        gap_floor = sum_floor_v3 - (m - 1)
        update_min(
            stats,
            "min_gap_sumFloorV3_minus_mminus1",
            gap_floor,
            {
                "n": n,
                "g6": g6,
                "m": m,
                "sum_floor_v3": sum_floor_v3,
                "child_sizes": decomp.child_sizes,
                "child_modes": child_modes,
            },
        )
        if gap_floor < 0:
            stats["size_proxy_floor_v3_ge_mminus1_fail"] += 1

        # Claim 4: per-child mode >= floor(v/3)
        for idx, (v, mu) in enumerate(zip(decomp.child_sizes, child_modes)):
            gap_child = mu - (v // 3)
            update_min(
                stats,
                "min_gap_child_mode_minus_floorv3",
                gap_child,
                {
                    "n": n,
                    "g6": g6,
                    "m": m,
                    "child_index": idx,
                    "child_size": v,
                    "child_mode": mu,
                },
            )
            if gap_child < 0:
                stats["per_child_floor_v3_fail"] += 1

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
