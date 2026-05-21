"""Probe eventual mode shifts for fixed-r spider lanes S(2^a,r).

For fixed ``r``,

    I_{a,r}(x) = P_r(x)(1+2x)^a + x P_{r-1}(x)(1+x)^a,

where P_r is the independence polynomial of a path on r vertices.  The second
term is exponentially small near the relevant saddle, so the mode should have
an eventually constant shift by residue class:

    3 m_{a,r} - 2a = D_{r, a mod 3}.

This script records those observed shifts and flags non-stabilization in the
chosen window.  It is a proof-target generator, not a proof.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from route2_spider_lane_scan import mode_left, path_polys, spider_poly_from_counts


def shifts_for_r(r: int, a_start: int, a_end: int, paths: list[list[int]]) -> dict:
    rows = []
    by_residue: dict[int, set[int]] = {0: set(), 1: set(), 2: set()}
    for a in range(a_start, a_end + 1):
        counts = {2: a}
        if r == 2:
            counts[2] += 1
        else:
            counts[r] = 1
        poly = spider_poly_from_counts(counts, paths)
        m = mode_left(poly)
        shift = 3 * m - 2 * a
        residue = a % 3
        by_residue[residue].add(shift)
        rows.append({"a": a, "m": m, "residue": residue, "shift": shift})

    stable = {q: len(values) == 1 for q, values in by_residue.items()}
    shift_by_residue = {
        q: next(iter(values)) if len(values) == 1 else sorted(values)
        for q, values in by_residue.items()
    }
    return {
        "r": r,
        "stable": all(stable.values()),
        "stable_by_residue": stable,
        "shift_by_residue": shift_by_residue,
        "last_rows": rows[-9:],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Fixed-r S(2^a,r) mode-shift probe.")
    ap.add_argument("--r-max", type=int, default=80)
    ap.add_argument("--a-start", type=int, default=200)
    ap.add_argument("--a-end", type=int, default=260)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    paths = path_polys(max(2, args.r_max))
    records = [
        shifts_for_r(r, args.a_start, args.a_end, paths)
        for r in range(2, args.r_max + 1)
    ]
    unstable = [rec for rec in records if not rec["stable"]]

    print(f"checked r=2..{args.r_max}, a={args.a_start}..{args.a_end}")
    print(f"unstable records in window: {len(unstable)}")
    if unstable:
        for rec in unstable[:10]:
            print(f"  r={rec['r']}: {rec['shift_by_residue']}")

    print("\nselected shifts D_q = 3m - 2a")
    selected = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 40, 80]
    for rec in records:
        if rec["r"] not in selected:
            continue
        shifts = rec["shift_by_residue"]
        print(
            f"r={rec['r']:2d}: "
            f"q0={shifts[0]!s:>3} q1={shifts[1]!s:>3} q2={shifts[2]!s:>3}"
        )

    if args.out:
        payload = {
            "params": {
                "r_max": args.r_max,
                "a_start": args.a_start,
                "a_end": args.a_end,
            },
            "records": records,
        }
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
