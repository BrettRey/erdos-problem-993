#!/usr/bin/env python3
"""Calibrate Path-1 direct chain from existing scan artifacts.

Checks the inequality
  mu_P - (m-2) = exact_slack_B - exact_excess_D
using precomputed global extrema:
  exact_slack_B := mu_B - (m - 1 - lambda/(1+lambda))
  exact_excess_D := D - (1 - lambda/(1+lambda)), D := mu_B - mu_P.

If:
  min exact_slack_B >= sigma
  max exact_excess_D <= kappa
then:
  mu_P - (m-2) >= sigma - kappa.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any


def load_json(path: str) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def min_from_per_n(payload: dict[str, Any], value_key: str, witness_key: str) -> tuple[float | None, Any]:
    best = None
    best_w = None
    per_n = payload.get("per_n", {})
    for n_key in sorted(per_n, key=lambda s: int(s)):
        stats = per_n[n_key]
        v = stats.get(value_key)
        if v is None:
            continue
        if best is None or v < best:
            best = float(v)
            best_w = stats.get(witness_key)
    return best, best_w


def max_from_per_n(payload: dict[str, Any], value_key: str, witness_key: str) -> tuple[float | None, Any]:
    best = None
    best_w = None
    per_n = payload.get("per_n", {})
    for n_key in sorted(per_n, key=lambda s: int(s)):
        stats = per_n[n_key]
        v = stats.get(value_key)
        if v is None:
            continue
        if best is None or v > best:
            best = float(v)
            best_w = stats.get(witness_key)
    return best, best_w


def main() -> None:
    ap = argparse.ArgumentParser(description="Verify Path-1 transfer-chain constants from existing result files.")
    ap.add_argument(
        "--transfer-json",
        default="results/whnc_route1_transfer_scan_n23_merged.json",
        help="Route-1 transfer scan artifact (contains max_exact_excess).",
    )
    ap.add_argument(
        "--pendant-json",
        nargs="+",
        default=[
            "results/whnc_pendant_bonus_scan_n20_canonical.json",
            "results/whnc_pendant_bonus_scan_n23_canonical_tail.json",
        ],
        help="One or more pendant-bonus artifacts to combine (contains min_exact_slack).",
    )
    ap.add_argument(
        "--phi-json",
        default="results/whnc_phiP_scan_n23_merged.json",
        help="Phi(P) artifact (for direct min mu_P gap cross-check).",
    )
    ap.add_argument("--transfer-cap", type=float, default=0.006)
    ap.add_argument("--out", default="results/path1_direct_chain_from_existing_results_2026_02_19.json")
    args = ap.parse_args()

    transfer = load_json(args.transfer_json)
    transfer_max_excess, transfer_max_witness = max_from_per_n(
        transfer,
        "max_exact_excess",
        "max_exact_excess_witness",
    )
    if transfer_max_excess is None:
        s = transfer.get("summary", {})
        transfer_max_excess = s.get("max_exact_excess")
        transfer_max_witness = s.get("max_exact_excess_witness")
    if transfer_max_excess is None:
        raise RuntimeError("Could not extract max_exact_excess from transfer artifact.")

    pendant_checked_total = 0
    pendant_seen_total = 0
    pendant_considered_total = 0
    pendant_min_exact = None
    pendant_min_witness = None

    for path in args.pendant_json:
        payload = load_json(path)
        v, w = min_from_per_n(payload, "min_exact_slack", "min_exact_witness")
        if v is None:
            s = payload.get("summary", {})
            v = s.get("min_exact_slack")
            w = s.get("min_exact_witness")
        if v is None:
            raise RuntimeError(f"Could not extract min_exact_slack from {path}")
        if pendant_min_exact is None or v < pendant_min_exact:
            pendant_min_exact = v
            pendant_min_witness = w

        s = payload.get("summary", {})
        pendant_checked_total += int(s.get("checked_leaves", 0))
        pendant_seen_total += int(s.get("seen", 0))
        pendant_considered_total += int(s.get("considered", 0))

    if pendant_min_exact is None:
        raise RuntimeError("No pendant min_exact_slack extracted.")

    phi = load_json(args.phi_json)
    phi_min_gap, phi_min_gap_witness = min_from_per_n(phi, "min_mu_gap", "min_mu_gap_witness")
    if phi_min_gap is None:
        s = phi.get("summary", {})
        phi_min_gap = s.get("min_mu_gap")
        phi_min_gap_witness = s.get("min_mu_gap_witness")

    lower_bound_chain = pendant_min_exact - transfer_max_excess
    cap_gap = args.transfer_cap - transfer_max_excess
    cap_holds = transfer_max_excess <= args.transfer_cap

    payload_out = {
        "formula": {
            "identity": "mu_P-(m-2) = exact_slack_B - exact_excess_D",
            "exact_slack_B": "mu_B - (m - 1 - lambda/(1+lambda))",
            "exact_excess_D": "D - (1 - lambda/(1+lambda))",
            "D": "mu_B - mu_P",
        },
        "inputs": {
            "transfer_json": args.transfer_json,
            "pendant_json": args.pendant_json,
            "phi_json": args.phi_json,
            "transfer_cap": args.transfer_cap,
        },
        "extrema": {
            "min_exact_slack_B": pendant_min_exact,
            "min_exact_slack_B_witness": pendant_min_witness,
            "max_exact_excess_D": transfer_max_excess,
            "max_exact_excess_D_witness": transfer_max_witness,
            "min_muP_gap_direct_phi_scan": phi_min_gap,
            "min_muP_gap_direct_phi_scan_witness": phi_min_gap_witness,
        },
        "derived": {
            "lower_bound_from_chain": lower_bound_chain,
            "cap_holds": cap_holds,
            "cap_margin": cap_gap,
            "implies_muP_ge_m_minus_2": lower_bound_chain >= 0.0,
        },
        "coverage": {
            "pendant_checked_leaves_total": pendant_checked_total,
            "pendant_seen_total": pendant_seen_total,
            "pendant_considered_total": pendant_considered_total,
            "transfer_checked": int(transfer.get("summary", {}).get("checked", 0)),
        },
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload_out, f, indent=2)

    print("Path-1 chain calibration from existing artifacts")
    print("-" * 96)
    print(f"min_exact_slack_B = {pendant_min_exact}")
    print(f"max_exact_excess_D = {transfer_max_excess}")
    print(f"lower_bound_from_chain = {lower_bound_chain}")
    print(f"transfer_cap = {args.transfer_cap} (holds={cap_holds}, margin={cap_gap})")
    if phi_min_gap is not None:
        print(f"direct min mu_P-(m-2) from phi scan = {phi_min_gap}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
