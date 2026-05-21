"""Run the sampled fixed-r Route-2 certificate suite.

This is a convenience driver for the currently sampled lanes

    r in {4,8,12,16,20,24,32,40,60,80}.

For each lane it checks:

1. exact finite Route-2 for a <= 199;
2. hub-off symbolic margins/reserve for a >= 200;
3. hub-on mode domination for a >= 200;
4. hub-on Route-2 perturbation for a >= 200.

The driver is intentionally a certificate runner, not a proof generator.  It
records which callable certificates passed and writes a compact JSON summary.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from fixed_r_finite_route2_check import DEFAULT_R_VALUES, check_lanes
from fixed_r_huboff_certificate import check_r as check_huboff
from fixed_r_hubon_mode_certificate import check_r as check_hubon_mode
from fixed_r_hubon_route2_perturbation import check_r as check_hubon_route2


def denominators_for_r(r: int) -> tuple[int, int]:
    """Return (margin_denom, reserve_denom) used by the sampled certificates."""
    if r == 8:
        return 10, 4
    return 1000, 1000


def run_suite(r_values: list[int], finite_max: int, threshold: int) -> dict:
    finite = check_lanes(r_values, finite_max)
    lanes = []
    all_ok = (
        finite["total_failures"] == 0
        and finite["total_stronger_failures"] == 0
    )

    for r in r_values:
        margin_denom, reserve_denom = denominators_for_r(r)
        print(f"\n=== sampled fixed-r lane r={r} ===")
        huboff_ok = check_huboff(
            r,
            threshold=threshold,
            margin_denom=margin_denom,
            reserve_denom=reserve_denom,
            exact_max=0,
        )
        mode_ok = check_hubon_mode(
            r,
            threshold=threshold,
            margin_denom=margin_denom,
        )
        perturb_ok = check_hubon_route2(
            r,
            threshold=threshold,
            reserve_denom=reserve_denom,
        )
        lane_ok = huboff_ok and mode_ok and perturb_ok
        all_ok = all_ok and lane_ok
        lanes.append(
            {
                "r": r,
                "margin_denom": margin_denom,
                "reserve_denom": reserve_denom,
                "huboff_ok": huboff_ok,
                "hubon_mode_ok": mode_ok,
                "hubon_route2_perturbation_ok": perturb_ok,
                "asymptotic_ok": lane_ok,
            }
        )

    return {
        "params": {
            "r_values": r_values,
            "finite_max": finite_max,
            "threshold": threshold,
        },
        "all_ok": all_ok,
        "finite": finite,
        "asymptotic_lanes": lanes,
    }


def parse_r_values(raw: str) -> list[int]:
    if raw.strip().lower() == "selected":
        return DEFAULT_R_VALUES[:]
    return [int(part) for part in raw.split(",") if part.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description="Run sampled fixed-r Route-2 certificates.")
    ap.add_argument("--r-values", default="selected")
    ap.add_argument("--finite-max", type=int, default=199)
    ap.add_argument("--threshold", type=int, default=200)
    ap.add_argument(
        "--out",
        default="results/fixed_r_sampled_certificate_suite.json",
    )
    args = ap.parse_args()

    result = run_suite(parse_r_values(args.r_values), args.finite_max, args.threshold)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))
    print(f"\nwrote {out}")
    print(f"all_ok={result['all_ok']}")
    assert result["all_ok"]


if __name__ == "__main__":
    main()
