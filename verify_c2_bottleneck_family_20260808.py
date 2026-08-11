#!/usr/bin/env python3
"""Exact audit of the asymptotically tight C2 base-invariant family.

The algebraic claim certified here is only the top diagonal inequality for
arm vectors (1,1) and (2h+1,2h+1).  Full partial synchronization is checked
through a finite parameter bound and is reported as computational evidence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from indpoly import is_log_concave
from verify_c2_bounded_pendant_core_20260808 import (
    adjacent_and_contracted,
    rooted_arm_state,
)
from verify_connector_partial_sync_route_20260808 import shift
from verify_two_hub_ratio_order_obstructions_20260808 import (
    partial_synchronization_failure,
)


def state(arms: tuple[int, ...]) -> tuple[tuple[int, ...], list[int], list[int]]:
    return (arms, *rooted_arm_state(arms))


def verify(maximum_long_arm: int) -> dict[str, object]:
    left = state((1, 1))
    minimum_top_margin: int | None = None
    minimum_top_margin_parameter: int | None = None

    for long_arm in range(1, maximum_long_arm + 1):
        adjacent, contracted = adjacent_and_contracted(
            left,
            state((long_arm, long_arm)),
        )
        assert is_log_concave(adjacent)
        assert is_log_concave(contracted)
        assert partial_synchronization_failure(contracted, adjacent) is None
        assert partial_synchronization_failure(
            shift(contracted),
            adjacent,
        ) is None

        if long_arm % 2 == 1:
            h = (long_arm - 1) // 2
            top_degree = 2 * h + 4
            assert len(adjacent) - 1 == top_degree
            assert len(contracted) - 1 == top_degree
            assert adjacent[top_degree] == 1
            assert contracted[top_degree] == 1
            assert adjacent[top_degree - 1] == 2 * h * h + 5 * h + 6
            assert contracted[top_degree - 1] == h * h + 3 * h + 4
            top_margin = (
                2 * contracted[top_degree - 1]
                - adjacent[top_degree - 1]
            )
            assert top_margin == h + 2
            if minimum_top_margin is None or top_margin < minimum_top_margin:
                minimum_top_margin = top_margin
                minimum_top_margin_parameter = long_arm

    return {
        "family": {
            "left_pendant_arms": [1, 1],
            "right_pendant_arms": "[m,m]",
            "parameter_range_checked": [1, maximum_long_arm],
        },
        "theorem": {
            "scope": "odd m=2h+1, top-degree diagonal only, all h>=0",
            "top_degree": "N=2h+4",
            "B_N_minus_1": "2h^2+5h+6",
            "C_N_minus_1": "h^2+3h+4",
            "partial_sync_margin": "2*C_(N-1)-B_(N-1)=h+2>0",
        },
        "computational_evidence": {
            "scope": (
                "Exact full log-concavity and both full two-index base "
                "partial-synchronization relations only through the stated "
                "finite m bound."
            ),
            "parameters_checked": maximum_long_arm,
            "failures": 0,
            "minimum_odd_top_margin": minimum_top_margin,
            "minimum_odd_top_margin_at_m": minimum_top_margin_parameter,
        },
        "not_claimed": (
            "The finite audit is not a proof of full partial synchronization "
            "for all m, of universal C2 log-concavity, or of Erdos #993."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-long-arm", type=int, default=500)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/c2_bottleneck_family_20260808.json"),
    )
    args = parser.parse_args()
    report = verify(args.max_long_arm)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"event": "verified", "out": str(args.out), **report}))


if __name__ == "__main__":
    main()
