#!/usr/bin/env python3
"""Summarize m-lambda-rho related JSON artifacts.

Supports:
  - canonical_projection_battery_minu*.json
  - adjacent_rho_split_scan_minu*.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _fmt_frac(x: Any) -> str:
    if isinstance(x, dict) and "num" in x and "den" in x:
        return f"{int(x['num'])}/{int(x['den'])}"
    if isinstance(x, (list, tuple)) and len(x) == 2:
        return f"{int(x[0])}/{int(x[1])}"
    return str(x)


def summarize_projection_battery(path: Path, data: dict[str, Any]) -> None:
    print(f"[projection-battery] {path}")
    print(
        "  "
        f"checked_total={data.get('checked_total')} "
        f"full_unique_keys={data.get('full_unique_keys')} "
        f"full_collisions={data.get('full_collisions')}"
    )
    for pr in data.get("projection_results", []):
        print(
            "  "
            f"{pr.get('name')}: split_found={pr.get('split_found')} "
            f"unique_keys={pr.get('unique_keys')} collisions={pr.get('collisions')}"
        )
        c1c2 = pr.get("c1c2_analysis")
        if isinstance(c1c2, dict):
            print(
                "    "
                f"c1_true={c1c2.get('c1_true')} "
                f"c2_true={c1c2.get('c2_true')} "
                f"c1c2_true={c1c2.get('c1c2_true')} "
                f"keys_total={c1c2.get('keys_total')}"
            )


def summarize_adjacent_scan(path: Path, data: dict[str, Any]) -> None:
    proj = data.get("projection_result", {})
    print(f"[adjacent-scan] {path}")
    print(
        "  "
        f"done={data.get('done')} "
        f"checked_total={data.get('checked_total')} "
        f"full_unique_keys={data.get('full_unique_keys')} "
        f"full_collisions={data.get('full_collisions')}"
    )
    print(
        "  "
        f"projection={proj.get('name')} "
        f"split_found={proj.get('split_found')} "
        f"adjacent_split_found={proj.get('adjacent_split_found')} "
        f"unique_keys={proj.get('unique_keys')} "
        f"collisions={proj.get('collisions')}"
    )
    fs = proj.get("first_adjacent_split")
    if not isinstance(fs, dict):
        return
    a = fs["A"]
    b = fs["B"]
    print(
        "  witness "
        f"A.N={a.get('N')} B.N={b.get('N')} "
        f"A.g6={a.get('g6')} B.g6={b.get('g6')}"
    )
    pkey = fs.get("projection_key")
    if isinstance(pkey, dict):
        print(
            "  key "
            f"m={pkey.get('m')} "
            f"lambda={_fmt_frac(pkey.get('lambda'))} "
            f"rho={_fmt_frac(pkey.get('rho'))}"
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "paths",
        nargs="+",
        help="One or more JSON artifact paths to summarize.",
    )
    args = ap.parse_args()

    for p in args.paths:
        path = Path(p)
        if not path.exists():
            print(f"[missing] {path}")
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        if "projection_result" in data and "scan" in data and "adjacent_rho" in str(data["scan"]):
            summarize_adjacent_scan(path, data)
        elif "projection_results" in data:
            summarize_projection_battery(path, data)
        else:
            print(f"[unknown-schema] {path}")


if __name__ == "__main__":
    main()
