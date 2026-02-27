#!/usr/bin/env python3
"""Validate and summarize adjacent rho split scan output.

Reads JSON from scripts/adjacent_rho_split_scan_minu.py and prints a compact
status report. If an adjacent split is present, validates:
  - same (m, lambda, rho) projection key,
  - |delta N| = 1.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def frac_pair(x: Any) -> tuple[int, int]:
    """Parse fraction-like values from supported JSON shapes."""
    if isinstance(x, dict) and "num" in x and "den" in x:
        return (int(x["num"]), int(x["den"]))
    if isinstance(x, (list, tuple)) and len(x) == 2:
        return (int(x[0]), int(x[1]))
    raise ValueError(f"unsupported fraction shape: {x!r}")


def projection_from_record(rec: dict[str, Any]) -> tuple[int, tuple[int, int], tuple[int, int]]:
    """Extract (m, lambda, rho) from a rebuilt record."""
    m = int(rec["m"])
    lam = frac_pair(rec["lambda"])
    rho = frac_pair(rec["rho"])
    return (m, lam, rho)


def projection_from_compact(rec: dict[str, Any]) -> tuple[int, tuple[int, int], tuple[int, int]]:
    """Extract (m, lambda, rho) from compact witness form."""
    inv = rec["full_invariants"]
    m = int(inv["m"])
    lam = frac_pair(inv["lambda"])
    rho = frac_pair(inv["rho"])
    return (m, lam, rho)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--in",
        dest="in_path",
        default="results/adjacent_rho_split_scan_minu_mge4_n27_exact.json",
        help="Input JSON from adjacent_rho_split_scan_minu.py",
    )
    args = ap.parse_args()

    path = Path(args.in_path)
    if not path.exists():
        raise SystemExit(f"missing file: {path}")

    data = json.loads(path.read_text(encoding="utf-8"))
    proj = data.get("projection_result", {})

    checked = int(data.get("checked_total", 0))
    full_unique = int(data.get("full_unique_keys", 0))
    full_collisions = int(data.get("full_collisions", 0))
    p_unique = int(proj.get("unique_keys", 0))
    p_collisions = int(proj.get("collisions", 0))
    split_found = bool(proj.get("split_found", False))
    adj_found = bool(proj.get("adjacent_split_found", False))
    done = bool(data.get("done", False))

    print(f"file={path}")
    print(f"done={done}")
    print(
        "totals "
        f"checked_total={checked} "
        f"full_unique_keys={full_unique} "
        f"full_collisions={full_collisions}"
    )
    print(
        "projection "
        f"unique_keys={p_unique} "
        f"collisions={p_collisions} "
        f"split_found={split_found} "
        f"adjacent_split_found={adj_found}"
    )

    if not adj_found:
        return

    fs = proj.get("first_adjacent_split")
    if not isinstance(fs, dict):
        raise SystemExit("adjacent_split_found=true but first_adjacent_split missing")

    a = fs["A"]
    b = fs["B"]
    delta = abs(int(a["N"]) - int(b["N"]))
    print(
        "adjacent_witness "
        f"A_N={int(a['N'])} B_N={int(b['N'])} |delta_N|={delta} "
        f"A_g6={a['g6']} B_g6={b['g6']}"
    )
    if delta != 1:
        raise SystemExit(f"invalid witness: expected |delta_N|=1, got {delta}")

    if "rebuilt" in fs:
        ra = fs["rebuilt"]["A"]
        rb = fs["rebuilt"]["B"]
        pa = projection_from_record(ra)
        pb = projection_from_record(rb)
    else:
        pa = projection_from_compact(a)
        pb = projection_from_compact(b)

    if pa != pb:
        raise SystemExit(
            "invalid witness: projection mismatch "
            f"A={pa} B={pb}"
        )

    print(
        "projection_key "
        f"m={pa[0]} lambda={pa[1][0]}/{pa[1][1]} rho={pa[2][0]}/{pa[2][1]}"
    )
    print("witness_validation=ok")


if __name__ == "__main__":
    main()
