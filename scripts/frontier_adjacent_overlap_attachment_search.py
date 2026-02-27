#!/usr/bin/env python3
"""Frontier odd/even overlap lifting search via bounded rooted gadget attachments.

Target pair defaults to the known same-(m,lambda) adjacent pair:
  A g6: O??????_A?C?E?@_WG@j?   (N=13)
  B g6: P????A?OD?E?E?B??o?E?OO? (N=14)

Search strategy:
  - Build rooted gadget library from d_leaf<=1 trees up to gadget-max-n.
  - Deduplicate rooted gadgets by (size, F, G), where
      F = I(component), G = I(component - root).
  - Enumerate multiset attachment patterns of up to max-gadgets gadgets,
    total added size <= max-added-size, and even total added size.
  - Attach the SAME pattern to both base trees at canonical root implicitly via
    canonical DP multiplication:
      P' = P0 * Fmul, Q' = Q0 * Gmul, I'=(1+2x)P'+(1+x)Q'
  - Keep records with m >= m_min and compare by key (deltaN, m, lambda).
  - Report whether any common key also has equal rho (witness).
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from fractions import Fraction
from typing import Any

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import mode_index_leftmost
from indpoly import _polyadd, _polymul
from scripts.canonical_kstar_split_scan_minu import rebuild_record, rooted_dp
from trees import trees_geng_raw


def poly_eval(poly: list[int], x: Fraction) -> Fraction:
    out = Fraction(0, 1)
    for c in reversed(poly):
        out = out * x + c
    return out


def frac_pair(fr: Fraction) -> list[int]:
    return [int(fr.numerator), int(fr.denominator)]


def make_i_poly(P: list[int], Q: list[int]) -> list[int]:
    return _polyadd(_polymul([1, 2], P), _polymul([1, 1], Q))


def rooted_fg(adj: list[list[int]], root: int) -> tuple[list[int], list[int]]:
    dp0, dp1 = rooted_dp(adj, root)
    G = dp0[root]
    F = _polyadd(dp0[root], dp1[root])
    return F, G


def build_gadget_library(max_n: int) -> list[dict[str, Any]]:
    # Dedup by (size, F, G)
    seen: dict[tuple[Any, ...], dict[str, Any]] = {}
    for n in range(2, max_n + 1):  # exclude trivial size-1 gadget
        for nn, adj, raw in trees_geng_raw(n):
            if not is_dleaf_le_1(nn, adj):
                continue
            for root in range(nn):
                F, G = rooted_fg(adj, root)
                key = (nn, tuple(F), tuple(G))
                if key in seen:
                    continue
                seen[key] = {
                    "size": nn,
                    "F": F,
                    "G": G,
                    "root": int(root),
                    "g6": raw.decode("ascii").strip(),
                }
    # Stable order for reproducibility
    gadgets = list(seen.values())
    gadgets.sort(key=lambda g: (g["size"], tuple(g["F"]), tuple(g["G"])))
    return gadgets


def build_patterns(
    gadgets: list[dict[str, Any]],
    max_gadgets: int,
    max_added_size: int,
) -> list[dict[str, Any]]:
    sizes = [g["size"] for g in gadgets]
    patterns: list[dict[str, Any]] = []
    for k in range(0, max_gadgets + 1):
        for combo in itertools.combinations_with_replacement(range(len(gadgets)), k):
            added = sum(sizes[i] for i in combo)
            if added > max_added_size:
                continue
            if added % 2 != 0:
                continue
            fmul = [1]
            gmul = [1]
            for idx in combo:
                fmul = _polymul(fmul, gadgets[idx]["F"])
                gmul = _polymul(gmul, gadgets[idx]["G"])
            patterns.append(
                {
                    "combo": list(combo),
                    "k": int(k),
                    "added_size": int(added),
                    "Fmul": fmul,
                    "Gmul": gmul,
                }
            )
    return patterns


def eval_pattern_on_base(
    base: dict[str, Any],
    pat: dict[str, Any],
    m_min: int,
) -> dict[str, Any] | None:
    P = _polymul(base["P"], pat["Fmul"])
    Q = _polymul(base["Q"], pat["Gmul"])
    I = make_i_poly(P, Q)
    m = mode_index_leftmost(I)
    if m < m_min:
        return None
    if m <= 0 or I[m] == 0:
        return None
    lam = Fraction(I[m - 1], I[m])
    p_lam = poly_eval(P, lam)
    if p_lam == 0:
        return None
    rho = poly_eval(Q, lam) / p_lam
    deltaN = int(P[1] - base["N"])
    return {
        "deltaN": deltaN,
        "m": int(m),
        "lambda": lam,
        "rho": rho,
        "N": int(P[1]),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pair-a-g6", default="O??????_A?C?E?@_WG@j?")
    ap.add_argument("--pair-b-g6", default="P????A?OD?E?E?B??o?E?OO?")
    ap.add_argument("--gadget-max-n", type=int, default=8)
    ap.add_argument("--max-gadgets", type=int, default=4)
    ap.add_argument("--max-added-size", type=int, default=24)
    ap.add_argument("--m-min", type=int, default=4)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    base_a = rebuild_record(16, args.pair_a_g6)
    base_b = rebuild_record(17, args.pair_b_g6)

    gadgets = build_gadget_library(args.gadget_max_n)
    patterns = build_patterns(gadgets, args.max_gadgets, args.max_added_size)

    # key = (deltaN,m,lam_num,lam_den) -> list of rho values (with pattern refs)
    map_a: dict[tuple[int, int, int, int], list[dict[str, Any]]] = {}
    map_b: dict[tuple[int, int, int, int], list[dict[str, Any]]] = {}

    checked = 0
    checked_a = 0
    checked_b = 0

    for idx, pat in enumerate(patterns):
        ra = eval_pattern_on_base(base_a, pat, args.m_min)
        rb = eval_pattern_on_base(base_b, pat, args.m_min)
        if ra is not None:
            checked += 1
            checked_a += 1
            ka = (ra["deltaN"], ra["m"], ra["lambda"].numerator, ra["lambda"].denominator)
            map_a.setdefault(ka, []).append(
                {
                    "rho": frac_pair(ra["rho"]),
                    "pattern_index": idx,
                    "added_size": pat["added_size"],
                    "k": pat["k"],
                    "N": ra["N"],
                }
            )
        if rb is not None:
            checked += 1
            checked_b += 1
            kb = (rb["deltaN"], rb["m"], rb["lambda"].numerator, rb["lambda"].denominator)
            map_b.setdefault(kb, []).append(
                {
                    "rho": frac_pair(rb["rho"]),
                    "pattern_index": idx,
                    "added_size": pat["added_size"],
                    "k": pat["k"],
                    "N": rb["N"],
                }
            )

    common_keys = sorted(set(map_a.keys()) & set(map_b.keys()))
    common_keys_fmt = [
        {
            "deltaN": k[0],
            "m": k[1],
            "lambda": [k[2], k[3]],
        }
        for k in common_keys
    ]

    witness = None
    closest = None

    for k in common_keys:
        ras = map_a[k]
        rbs = map_b[k]
        set_a = {tuple(x["rho"]) for x in ras}
        set_b = {tuple(x["rho"]) for x in rbs}
        inter = set_a & set_b
        if inter and witness is None:
            rho_pair = next(iter(inter))
            a_ex = next(x for x in ras if tuple(x["rho"]) == rho_pair)
            b_ex = next(x for x in rbs if tuple(x["rho"]) == rho_pair)
            witness = {
                "key": {"deltaN": k[0], "m": k[1], "lambda": [k[2], k[3]], "rho": list(rho_pair)},
                "A": a_ex,
                "B": b_ex,
            }

        # Closest rho gap for this shared (deltaN,m,lambda)
        min_gap = None
        min_pair = None
        for xa in ras:
            fa = Fraction(xa["rho"][0], xa["rho"][1])
            for xb in rbs:
                fb = Fraction(xb["rho"][0], xb["rho"][1])
                g = abs(fa - fb)
                if min_gap is None or g < min_gap:
                    min_gap = g
                    min_pair = (xa, xb)
        if min_gap is not None:
            rec = {
                "key": {"deltaN": k[0], "m": k[1], "lambda": [k[2], k[3]]},
                "rho_gap": [min_gap.numerator, min_gap.denominator],
                "A": min_pair[0],
                "B": min_pair[1],
            }
            if closest is None:
                closest = rec
            else:
                g0 = Fraction(closest["rho_gap"][0], closest["rho_gap"][1])
                if min_gap < g0:
                    closest = rec

    payload = {
        "scan": "frontier_adjacent_overlap_attachment_search",
        "pair": {
            "A": {"g6": args.pair_a_g6, "N0": int(base_a["N"]), "m0": int(base_a["m"]), "lambda0": base_a["lambda"]},
            "B": {"g6": args.pair_b_g6, "N0": int(base_b["N"]), "m0": int(base_b["m"]), "lambda0": base_b["lambda"]},
        },
        "params": {
            "gadget_max_n": args.gadget_max_n,
            "max_gadgets": args.max_gadgets,
            "max_added_size": args.max_added_size,
            "m_min": args.m_min,
            "even_added_size_only": True,
        },
        "totals": {
            "gadgets_count": len(gadgets),
            "patterns_count": len(patterns),
            "checked": checked,
            "checked_A": checked_a,
            "checked_B": checked_b,
            "common_size_m_lambda_count": len(common_keys),
            "collision_count": 1 if witness is not None else 0,
            "split_found": witness is not None,
        },
        "common_size_m_lambda_keys": common_keys_fmt,
        "first_split": witness,
        "closest_overlap_by_rho_gap": closest,
    }

    print(
        "done "
        f"gadgets={len(gadgets)} patterns={len(patterns)} checked={checked} "
        f"common_keys={len(common_keys)} split_found={witness is not None}",
        flush=True,
    )

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
