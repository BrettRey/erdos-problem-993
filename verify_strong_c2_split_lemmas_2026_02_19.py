#!/usr/bin/env python3
"""Post-process rise-identity scan output into lemma-style checks.

This script does not re-enumerate trees. It consumes
`verify_strong_c2_rise_identity_2026_02_19.py` output and verifies the key
algebraic decompositions and regime-specific inequalities on recorded
mismatch-negative witnesses.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any


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
    ap = argparse.ArgumentParser(description="Verify lemma-level inequalities on mismatch-negative witnesses.")
    ap.add_argument(
        "--input",
        default="results/verify_strong_c2_rise_identity_2026_02_19.json",
        help="Output JSON from verify_strong_c2_rise_identity_2026_02_19.py",
    )
    ap.add_argument(
        "--out",
        default="results/verify_strong_c2_split_lemmas_2026_02_19.json",
    )
    ap.add_argument("--tol", type=float, default=1e-12)
    args = ap.parse_args()

    with open(args.input, "r", encoding="utf-8") as f:
        payload = json.load(f)

    rows = payload.get("mismatch_neg_witnesses", [])

    stats: dict[str, Any] = {
        "mismatch_neg": len(rows),
        # Exact decomposition checks
        "decomp_fail": 0,  # combined != (rise-neg) + b0*(b1-b2)
        "ratio_form_fail": 0,  # in dq<0 regime, rise-neg ratio form fails
        # Regime counts
        "hard": 0,
        "easy": 0,
        "hard_shift0": 0,
        "hard_shift1": 0,
        "hard_shift_other": 0,
        # Candidate structural inequalities
        "shift1_qdrop": 0,
        "shift1_comp_bound_fail": 0,
        "shift1_comp_bound_min_margin": None,
        "shift1_comp_bound_min_margin_witness": None,
        # Min margins
        "min_decomp_margin": None,
        "min_decomp_margin_witness": None,
        "min_ratio_gap_hard": None,
        "min_ratio_gap_hard_witness": None,
        "max_transfer_ratio_hard": None,
        "max_transfer_ratio_hard_witness": None,
        "min_need_ratio_hard": None,
        "min_need_ratio_hard_witness": None,
    }

    for r in rows:
        b0 = int(r["b0"])
        b1 = int(r["b1"])
        b2 = int(r["b2"])
        p1 = int(r["p1"])
        q0 = int(r["q0"])
        q1 = int(r["q1"])
        dq = int(r["dq"])
        db = int(r["db"])
        shift = int(r["shift"])

        rise = int(r["rise"])
        neg = int(r["neg"])
        combined = int(r["combined"])

        decomp_rhs = (rise - neg) + b0 * (b1 - b2)
        if combined != decomp_rhs:
            stats["decomp_fail"] += 1

        decomp_margin = combined - decomp_rhs
        maybe_min(stats, "min_decomp_margin", "min_decomp_margin_witness", float(decomp_margin), r)

        if dq < 0:
            stats["hard"] += 1
            if shift == 0:
                stats["hard_shift0"] += 1
            elif shift == 1:
                stats["hard_shift1"] += 1
            else:
                stats["hard_shift_other"] += 1

            if db <= 0 or b1 <= 0:
                stats["ratio_form_fail"] += 1
            else:
                transfer_ratio = (-dq) / db
                need_ratio = p1 / b1
                ratio_gap = need_ratio - transfer_ratio

                w_ratio = {
                    **r,
                    "transfer_ratio": transfer_ratio,
                    "need_ratio": need_ratio,
                    "ratio_gap": ratio_gap,
                }
                maybe_min(stats, "min_ratio_gap_hard", "min_ratio_gap_hard_witness", float(ratio_gap), w_ratio)
                maybe_max(
                    stats,
                    "max_transfer_ratio_hard",
                    "max_transfer_ratio_hard_witness",
                    float(transfer_ratio),
                    w_ratio,
                )
                maybe_min(stats, "min_need_ratio_hard", "min_need_ratio_hard_witness", float(need_ratio), w_ratio)
                if ratio_gap < -args.tol:
                    stats["ratio_form_fail"] += 1
        else:
            stats["easy"] += 1

        # Shift-1 candidate bound:
        #   p1*(b1-b0) >= b0*(b2-b1)
        # This plus dq>=0 implies combined>=0 in shift=1 via
        #   combined = [p1*(b1-b0)-b0*(b2-b1)] + b1*dq.
        if shift == 1:
            if dq < 0:
                stats["shift1_qdrop"] += 1
            margin = p1 * (b1 - b0) - b0 * (b2 - b1)
            maybe_min(
                stats,
                "shift1_comp_bound_min_margin",
                "shift1_comp_bound_min_margin_witness",
                float(margin),
                r,
            )
            if margin < -args.tol:
                stats["shift1_comp_bound_fail"] += 1

    out = {
        "input": args.input,
        "summary": stats,
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print(json.dumps(out["summary"], indent=2))
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
