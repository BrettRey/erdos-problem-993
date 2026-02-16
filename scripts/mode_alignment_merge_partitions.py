#!/usr/bin/env python3
"""Merge partitioned mode-alignment JSON artifacts."""

from __future__ import annotations

import argparse
import glob
import json
from typing import Any


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pattern",
        required=True,
        help="Glob for partition JSON files, e.g. results/mode_alignment_n21_res*_mod8.json",
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    paths = sorted(glob.glob(args.pattern))
    if not paths:
        raise SystemExit(f"no files matched: {args.pattern}")

    files_read = 0
    n_value: int | None = None
    mod_value: int | None = None
    seen_res: set[int] = set()

    total_trees = 0
    total_vertex_cases = 0
    total_failures = 0
    total_equalities = 0
    global_max_gap = -(10**9)
    global_max_abs_mode_diff = 0
    sample_failures: list[dict[str, Any]] = []
    sample_equalities: list[dict[str, Any]] = []

    for path in paths:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        files_read += 1

        n = int(data["n"])
        res = int(data["res"])
        mod = int(data["mod"])
        if n_value is None:
            n_value = n
        elif n_value != n:
            raise SystemExit(f"mixed n values: saw {n_value} and {n} ({path})")
        if mod_value is None:
            mod_value = mod
        elif mod_value != mod:
            raise SystemExit(f"mixed mod values: saw {mod_value} and {mod} ({path})")
        if res in seen_res:
            raise SystemExit(f"duplicate res partition: {res} ({path})")
        seen_res.add(res)

        total_trees += int(data["trees"])
        total_vertex_cases += int(data["vertex_cases"])
        total_failures += int(data["failures"])
        total_equalities += int(data["equality_count"])
        global_max_gap = max(global_max_gap, int(data["max_mode_gap"]))
        global_max_abs_mode_diff = max(
            global_max_abs_mode_diff, int(data["max_abs_mode_diff"])
        )

        for ex in data.get("sample_failures", []):
            if len(sample_failures) < 10:
                sample_failures.append(ex)
        for ex in data.get("sample_equalities", []):
            if len(sample_equalities) < 10:
                sample_equalities.append(ex)

    merged = {
        "claim": "mode(I(T-w)) <= d(I(T))",
        "n": n_value,
        "mod": mod_value,
        "partitions_found": sorted(seen_res),
        "partitions_expected": list(range(mod_value)) if mod_value is not None else [],
        "is_complete": (mod_value is not None and seen_res == set(range(mod_value))),
        "files_read": files_read,
        "total_trees": total_trees,
        "total_vertex_cases": total_vertex_cases,
        "total_failures": total_failures,
        "total_equality_count": total_equalities,
        "max_mode_gap": global_max_gap,
        "max_abs_mode_diff": global_max_abs_mode_diff,
        "sample_failures": sample_failures,
        "sample_equalities": sample_equalities,
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(merged, f, indent=2)

    status = "PASS" if total_failures == 0 else "FAIL"
    completeness = "complete" if merged["is_complete"] else "incomplete"
    print(
        f"{status} ({completeness}): n={n_value}, files={files_read}, "
        f"trees={total_trees:,}, vertex_cases={total_vertex_cases:,}, "
        f"failures={total_failures:,}, max_gap={global_max_gap}, "
        f"max|diff|={global_max_abs_mode_diff}"
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

