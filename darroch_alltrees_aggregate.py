#!/usr/bin/env python3
"""Aggregate chunk JSON outputs from darroch_alltrees_chunk_scan.py."""

from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path
from typing import Any


def main() -> None:
    ap = argparse.ArgumentParser(description="Aggregate chunked Darroch scan outputs.")
    ap.add_argument("--glob", required=True, help="Glob for chunk JSON files.")
    ap.add_argument("--out", required=True, help="Output aggregate JSON path.")
    ap.add_argument(
        "--expect-mod",
        type=int,
        default=0,
        help="If >0, require exactly this many distinct chunk residues.",
    )
    args = ap.parse_args()

    files = sorted(glob.glob(args.glob))
    if not files:
        raise FileNotFoundError(f"No files matched glob: {args.glob}")

    chunks: dict[int, dict[str, Any]] = {}
    n_value: int | None = None
    mod_value: int | None = None

    total_trees = 0
    mode_ceil_fail = 0
    darroch_fail = 0
    min_mode_ceil_margin: int | None = None
    min_set_margin: int | None = None
    wall_s = 0.0
    first_mode_ceil_fail: dict[str, Any] | None = None
    first_darroch_fail: dict[str, Any] | None = None

    for path in files:
        payload = json.loads(Path(path).read_text())
        s = payload["summary"]
        n = int(s["n"])
        res = int(s["chunk"]["res"])
        mod = int(s["chunk"]["mod"])

        if n_value is None:
            n_value = n
        elif n != n_value:
            raise ValueError(f"Inconsistent n across files: {n} vs {n_value} ({path})")

        if mod_value is None:
            mod_value = mod
        elif mod != mod_value:
            raise ValueError(
                f"Inconsistent mod across files: {mod} vs {mod_value} ({path})"
            )

        if res in chunks:
            raise ValueError(f"Duplicate residue {res} across files; found again in {path}")

        chunks[res] = {
            "path": path,
            "total_trees": int(s["total_trees"]),
            "mode_ceil_fail": int(s["mode_ceil_fail"]),
            "darroch_fail": int(s["darroch_fail"]),
            "min_mode_ceil_margin": s["min_mode_ceil_margin"],
            "min_set_margin": s["min_set_margin"],
            "wall_s": float(s["wall_s"]),
        }

        total_trees += int(s["total_trees"])
        mode_ceil_fail += int(s["mode_ceil_fail"])
        darroch_fail += int(s["darroch_fail"])
        wall_s += float(s["wall_s"])

        mm = s["min_mode_ceil_margin"]
        if mm is not None and (min_mode_ceil_margin is None or mm < min_mode_ceil_margin):
            min_mode_ceil_margin = mm

        sm = s["min_set_margin"]
        if sm is not None and (min_set_margin is None or sm < min_set_margin):
            min_set_margin = sm

        if first_mode_ceil_fail is None and s["first_mode_ceil_fail"] is not None:
            first_mode_ceil_fail = {
                "chunk_res": res,
                "chunk_mod": mod,
                "path": path,
                "witness": s["first_mode_ceil_fail"],
            }
        if first_darroch_fail is None and s["first_darroch_fail"] is not None:
            first_darroch_fail = {
                "chunk_res": res,
                "chunk_mod": mod,
                "path": path,
                "witness": s["first_darroch_fail"],
            }

    assert n_value is not None
    assert mod_value is not None

    missing_residues = [r for r in range(mod_value) if r not in chunks]
    if args.expect_mod > 0 and mod_value != args.expect_mod:
        raise ValueError(f"Expected mod={args.expect_mod}, but files have mod={mod_value}")

    if args.expect_mod > 0 and missing_residues:
        raise ValueError(f"Missing residues for full partition: {missing_residues}")

    out_payload = {
        "params": {
            "glob": args.glob,
            "expect_mod": args.expect_mod,
        },
        "summary": {
            "n": n_value,
            "mod": mod_value,
            "chunks_found": len(chunks),
            "missing_residues": missing_residues,
            "total_trees": total_trees,
            "mode_ceil_fail": mode_ceil_fail,
            "darroch_fail": darroch_fail,
            "min_mode_ceil_margin": min_mode_ceil_margin,
            "min_set_margin": min_set_margin,
            "first_mode_ceil_fail": first_mode_ceil_fail,
            "first_darroch_fail": first_darroch_fail,
            "sum_chunk_wall_s": wall_s,
        },
        "chunks": {str(r): chunks[r] for r in sorted(chunks)},
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out_payload, f, indent=2)

    print(f"Wrote {args.out}")
    print(
        f"n={n_value}, mod={mod_value}, chunks={len(chunks)}, "
        f"total_trees={total_trees:,}, mode<=ceil fails={mode_ceil_fail:,}, "
        f"darroch fails={darroch_fail:,}"
    )


if __name__ == "__main__":
    main()
