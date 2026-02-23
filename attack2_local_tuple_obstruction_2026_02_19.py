#!/usr/bin/env python3
"""Search abstract local tuples showing limits of direct coefficient inequalities.

Variables (near mode index m):
  p2 = p_{m-2}, p1 = p_{m-1}, p0 = p_m, pp = p_{m+1}
  q2 = q_{m-2}, q1 = q_{m-1}, q0 = q_m, qp = q_{m+1}

Goal:
  Find tuples with p1 < p2 (target failure) that still satisfy:
    - I1: (p0 + p1 - 2p2) + (q0 - q2) >= 0
    - I2: (2p1 - p0 - pp) + (q1 - qp) >= 0
    - local log-concavity of P and Q at m-1,m
    - local dominance bounds q_{k+1} <= p_k.

This does not prove realizability by trees; it isolates what the two inequalities
can and cannot force without additional structure.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from typing import Any


@dataclass
class LocalWitness:
    p_m_minus_2: int
    p_m_minus_1: int
    p_m: int
    p_m_plus_1: int
    q_m_minus_2: int
    q_m_minus_1: int
    q_m: int
    q_m_plus_1: int
    g_target: int
    i1_mode_left: int
    i2_mode_right: int
    p_term1: int
    q_term1: int
    p_term2: int
    q_term2: int


def is_lc_trip(a: int, b: int, c: int) -> bool:
    return b * b >= a * c


def search_first(max_val: int, require_positive_tail: bool, enforce_q_shape: bool) -> LocalWitness | None:
    for p2 in range(1, max_val + 1):
        for p1 in range(0, p2):
            for p0 in range(0, p1 + 1):
                pp_min = 1 if require_positive_tail else 0
                for pp in range(pp_min, p0 + 1):
                    # Local shape compatible with mode(P)<=m-2 assumption.
                    if not (p2 >= p1 >= p0 >= pp):
                        continue
                    if not is_lc_trip(p2, p1, p0):
                        continue
                    if not is_lc_trip(p1, p0, pp):
                        continue
                    # Target failure.
                    if p1 - p2 >= 0:
                        continue

                    q2_min = 1 if require_positive_tail else 0
                    q1_min = 1 if require_positive_tail else 0
                    q0_min = 1 if require_positive_tail else 0
                    qp_min = 1 if require_positive_tail else 0

                    for q2 in range(q2_min, p2 + 1):
                        for q1 in range(q1_min, p2 + 1):
                            for q0 in range(q0_min, p1 + 1):
                                for qp in range(qp_min, p0 + 1):
                                    # Dominance: q_{k+1} <= p_k.
                                    if q0 > p1 or q1 > p2 or qp > p0:
                                        continue
                                    # Local LC for Q.
                                    if not is_lc_trip(q2, q1, q0):
                                        continue
                                    if not is_lc_trip(q1, q0, qp):
                                        continue
                                    if enforce_q_shape:
                                        # Optional "rise-to-m then fall" local shape.
                                        if not (q0 >= q1 >= q2):
                                            continue
                                        if not (q0 >= qp):
                                            continue

                                    p_term1 = p0 + p1 - 2 * p2
                                    q_term1 = q0 - q2
                                    p_term2 = 2 * p1 - p0 - pp
                                    q_term2 = q1 - qp
                                    i1 = p_term1 + q_term1
                                    i2 = p_term2 + q_term2
                                    if i1 < 0 or i2 < 0:
                                        continue

                                    return LocalWitness(
                                        p_m_minus_2=p2,
                                        p_m_minus_1=p1,
                                        p_m=p0,
                                        p_m_plus_1=pp,
                                        q_m_minus_2=q2,
                                        q_m_minus_1=q1,
                                        q_m=q0,
                                        q_m_plus_1=qp,
                                        g_target=p1 - p2,
                                        i1_mode_left=i1,
                                        i2_mode_right=i2,
                                        p_term1=p_term1,
                                        q_term1=q_term1,
                                        p_term2=p_term2,
                                        q_term2=q_term2,
                                    )
    return None


def run_regime(
    max_val: int,
    require_positive_tail: bool,
    enforce_q_shape: bool,
) -> dict[str, Any]:
    witness = search_first(
        max_val=max_val,
        require_positive_tail=require_positive_tail,
        enforce_q_shape=enforce_q_shape,
    )
    return {
        "max_val": max_val,
        "require_positive_tail": require_positive_tail,
        "enforce_q_shape": enforce_q_shape,
        "found": witness is not None,
        "witness": asdict(witness) if witness is not None else None,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-val", type=int, default=40)
    ap.add_argument(
        "--out",
        default="results/attack2_local_tuple_obstruction.json",
        help="JSON output path",
    )
    args = ap.parse_args()

    regimes = [
        ("minimal_constraints", False, False),
        ("with_q_shape", False, True),
        ("with_q_shape_and_positive_tail", True, True),
    ]

    payload: dict[str, Any] = {"params": {"max_val": args.max_val}, "regimes": {}}
    for name, positive_tail, q_shape in regimes:
        payload["regimes"][name] = run_regime(
            max_val=args.max_val,
            require_positive_tail=positive_tail,
            enforce_q_shape=q_shape,
        )

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"Output written to: {args.out}", flush=True)
    for name in regimes:
        label = name[0]
        reg = payload["regimes"][label]
        print(f"{label}: found={reg['found']}", flush=True)
        if reg["found"]:
            print(f"  witness={reg['witness']}", flush=True)


if __name__ == "__main__":
    main()

