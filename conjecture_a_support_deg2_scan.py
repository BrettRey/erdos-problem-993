#!/usr/bin/env python3
"""Resumable structural scan for support-degree patterns in d_leaf<=1 trees.

For each n in [min_n, max_n], enumerate all n-vertex trees via geng, filter to
d_leaf<=1, and record:
  - whether each tree has a support vertex of degree 2
  - whether the minimum support degree is 2
  - histogram of minimum support degree

The output JSON is checkpointed after every n so interrupted runs can resume.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from graph6 import parse_graph6


def load_existing(path: str) -> dict[str, Any]:
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_payload(path: str, payload: dict[str, Any]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    os.replace(tmp, path)


def merge_hist(dst: Counter[int], src: Counter[int]) -> None:
    for k, v in src.items():
        dst[k] += v


def recompute_summary(per_n: dict[str, dict[str, Any]]) -> dict[str, Any]:
    total_dleaf = 0
    trees_with_support = 0
    has_deg2_support = 0
    min_support_deg_eq2 = 0
    min_support_deg_gt2 = 0
    support_deg_hist = Counter()
    first_no_deg2 = None
    seen_total = 0

    for n_str in sorted(per_n.keys(), key=lambda x: int(x)):
        row = per_n[n_str]
        seen_total += int(row.get("seen", 0))
        total_dleaf += int(row.get("dleaf", 0))
        trees_with_support += int(row.get("trees_with_support", 0))
        has_deg2_support += int(row.get("has_deg2_support", 0))
        min_support_deg_eq2 += int(row.get("min_support_deg_eq2", 0))
        min_support_deg_gt2 += int(row.get("min_support_deg_gt2", 0))
        merge_hist(
            support_deg_hist,
            Counter({int(k): int(v) for k, v in row.get("support_deg_hist", {}).items()}),
        )
        if first_no_deg2 is None and row.get("first_no_deg2") is not None:
            first_no_deg2 = row["first_no_deg2"]

    return {
        "seen_total": seen_total,
        "total_dleaf": total_dleaf,
        "trees_with_support": trees_with_support,
        "has_deg2_support": has_deg2_support,
        "min_support_deg_eq2": min_support_deg_eq2,
        "min_support_deg_gt2": min_support_deg_gt2,
        "support_deg_hist": {str(k): support_deg_hist[k] for k in sorted(support_deg_hist)},
        "first_no_deg2": first_no_deg2,
    }


def scan_n(n: int, geng: str, mod: int) -> dict[str, Any]:
    seen = 0
    dleaf = 0
    trees_with_support = 0
    has_deg2_support = 0
    min_support_deg_eq2 = 0
    min_support_deg_gt2 = 0
    support_deg_hist: Counter[int] = Counter()
    first_no_deg2 = None
    part_stats: list[dict[str, Any]] = []

    for res in range(mod):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        if mod > 1:
            cmd.append(f"{res}/{mod}")

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert proc.stdout is not None

        part_seen = 0
        for raw in proc.stdout:
            part_seen += 1
            seen += 1
            nn, adj = parse_graph6(raw)
            if not is_dleaf_le_1(nn, adj):
                continue

            dleaf += 1
            deg = [len(nb) for nb in adj]
            leaves = [v for v, dv in enumerate(deg) if dv == 1]
            supports = sorted({adj[l][0] for l in leaves})

            if not supports:
                continue

            trees_with_support += 1
            sdeg = [deg[s] for s in supports]
            ms = min(sdeg)
            support_deg_hist[ms] += 1

            if any(x == 2 for x in sdeg):
                has_deg2_support += 1
            elif first_no_deg2 is None:
                first_no_deg2 = {
                    "n": nn,
                    "g6": raw.decode("ascii").strip(),
                    "degree_signature": dict(sorted(Counter(deg).items())),
                    "support_degrees": sdeg,
                }

            if ms == 2:
                min_support_deg_eq2 += 1
            else:
                min_support_deg_gt2 += 1

        stderr = proc.stderr.read().decode("utf-8", errors="replace")
        ret = proc.wait()
        part_stats.append({"res": res, "seen": part_seen, "returncode": ret})
        if ret != 0:
            raise RuntimeError(
                f"geng failed at n={n}, res={res}/{mod}, returncode={ret}, stderr={stderr.strip()}"
            )

    return {
        "seen": seen,
        "dleaf": dleaf,
        "trees_with_support": trees_with_support,
        "has_deg2_support": has_deg2_support,
        "min_support_deg_eq2": min_support_deg_eq2,
        "min_support_deg_gt2": min_support_deg_gt2,
        "support_deg_hist": {str(k): support_deg_hist[k] for k in sorted(support_deg_hist)},
        "first_no_deg2": first_no_deg2,
        "parts": part_stats,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Resumable support-degree scan for d_leaf<=1 trees")
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", default="results/conjA_support_deg2_structure_n23.json")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument(
        "--mod",
        type=int,
        default=1,
        help="Split geng enumeration into residue classes res/mod (recommended: 8 for n=23).",
    )
    ap.add_argument(
        "--recompute-n",
        type=int,
        nargs="*",
        default=[],
        help="Force recomputation for listed n values even when present in output JSON.",
    )
    args = ap.parse_args()

    t0 = time.time()
    payload: dict[str, Any] = {
        "params": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "geng": args.geng,
            "resume": args.resume,
            "mod": args.mod,
        },
        "per_n": {},
        "summary": {},
    }

    if args.resume:
        existing = load_existing(args.out)
        if existing:
            payload["per_n"] = existing.get("per_n", {})
            payload["params"] = existing.get("params", payload["params"])

    for n in range(args.min_n, args.max_n + 1):
        key = str(n)
        if key in payload["per_n"] and n not in set(args.recompute_n):
            print(f"n={n}: already present, skipping", flush=True)
            continue

        tn = time.time()
        row = scan_n(n, args.geng, args.mod)
        row["wall_s"] = time.time() - tn
        payload["per_n"][key] = row

        payload["summary"] = recompute_summary(payload["per_n"])
        payload["summary"]["wall_s"] = time.time() - t0
        save_payload(args.out, payload)

        print(
            f"n={n}: dleaf={row['dleaf']:,}, has_deg2={row['has_deg2_support']:,}, "
            f"min_deg2={row['min_support_deg_eq2']:,}, min_gt2={row['min_support_deg_gt2']:,} "
            f"({row['wall_s']:.1f}s)",
            flush=True,
        )

    payload["summary"] = recompute_summary(payload["per_n"])
    payload["summary"]["wall_s"] = time.time() - t0
    save_payload(args.out, payload)
    print(f"Wrote {args.out}", flush=True)
    print(json.dumps(payload["summary"], indent=2), flush=True)


if __name__ == "__main__":
    main()
