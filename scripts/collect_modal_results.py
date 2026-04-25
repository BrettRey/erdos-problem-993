#!/usr/bin/env python3
"""Collect and merge Modal dict-backed exhaustive run outputs.

Supports the two dict schemas used in this repo:
  - unimodality search dicts from ``search_modal_exhaustive*.py``
  - LC / near-miss analysis dicts from ``analyze_modal_lc_nm*.py``

Examples:
  python3 scripts/collect_modal_results.py status --kind unimodality --n 28 --workers 1024
  python3 scripts/collect_modal_results.py collect --kind unimodality --n 28 --workers 1024 \
    --out results/analysis_n28_modal_unimodality.json
  python3 scripts/collect_modal_results.py collect --kind lc_nm --n 28 --workers 1024 \
    --top-k 200 --lc-top-k 200 --out results/analysis_n28_modal_lc_nm.json
"""

from __future__ import annotations

import argparse
import heapq
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from search import A000055_COUNTS


def dict_name(kind: str, n: int) -> str:
    if kind == "unimodality":
        return f"erdos-993-n{n}-unimodality-results"
    if kind == "lc_nm":
        return f"erdos-993-n{n}-lc-nm-results"
    raise ValueError(f"unknown kind: {kind}")


def load_modal_dict(name: str, limit: int) -> dict[str, Any]:
    cmd = ["modal", "dict", "items", name, str(limit), "--json"]
    raw = subprocess.check_output(cmd, text=True)
    rows = json.loads(raw)
    out: dict[str, Any] = {}
    for row in rows:
        out[row["Key"]] = row["Value"]
    return out


def summarize_unimodality(rows: dict[str, Any], n: int, workers: int, expected: int | None) -> dict[str, Any]:
    total_trees = 0
    counterexamples: list[dict[str, Any]] = []
    for value in rows.values():
        total_trees += value["count"]
        cx = value.get("counterexample")
        if cx is not None:
            counterexamples.append(cx)
    return {
        "n": n,
        "workers": workers,
        "completed": len(rows),
        "total_trees_so_far": total_trees,
        "expected": expected,
        "expected_match_so_far": (total_trees == expected) if expected is not None and len(rows) == workers else None,
        "counterexamples_found": len(counterexamples),
        "counterexample_examples": counterexamples[:3],
        "dict_name": dict_name("unimodality", n),
        "complete": len(rows) == workers,
    }


def collect_unimodality(rows: dict[str, Any], n: int, workers: int, expected: int | None) -> dict[str, Any]:
    status = summarize_unimodality(rows, n, workers, expected)
    return {
        "n": n,
        "workers": workers,
        "total_trees": status["total_trees_so_far"],
        "expected": expected,
        "match": (status["total_trees_so_far"] == expected) if expected is not None else None,
        "counterexamples": status["counterexamples_found"],
        "counterexample_examples": status["counterexample_examples"],
        "platform": "Modal",
        "dict_name": status["dict_name"],
        "date": time.strftime("%Y-%m-%d"),
        "complete": status["complete"],
        "completed_partitions": status["completed"],
    }


def resolve_expected(n: int, cli_expected: int | None) -> int | None:
    """Return the OEIS A000055 expected count, checking any CLI override."""
    oeis_expected = A000055_COUNTS.get(n)
    if cli_expected is not None and oeis_expected is not None and cli_expected != oeis_expected:
        raise ValueError(
            f"--expected={cli_expected} disagrees with OEIS A000055 count "
            f"for n={n}: {oeis_expected}"
        )
    return oeis_expected if oeis_expected is not None else cli_expected


def collect_lc_nm(
    rows: dict[str, Any],
    n: int,
    workers: int,
    top_k: int,
    lc_top_k: int,
) -> dict[str, Any]:
    agg_count = 0
    agg_non_unimodal = 0
    agg_lc_fail = 0

    global_worst_lc_ratio = -1.0
    global_worst_lc_item: dict[str, Any] | None = None

    merged_nm: list[tuple[float, int, dict[str, Any]]] = []
    merged_lc: list[tuple[float, int, dict[str, Any]]] = []
    all_lc_failures: list[dict[str, Any]] = []
    serial = 0
    collect_all_lc = False

    for value in rows.values():
        agg_count += value["count"]
        agg_non_unimodal += value["non_unimodal_count"]
        agg_lc_fail += value["lc_fail_count"]

        if value["worst_lc_ratio"] > global_worst_lc_ratio:
            global_worst_lc_ratio = value["worst_lc_ratio"]
            global_worst_lc_item = value["worst_lc_item"]

        for item in value["top_nm"]:
            serial += 1
            ratio = item["nm_ratio"]
            if len(merged_nm) < top_k:
                heapq.heappush(merged_nm, (ratio, serial, item))
            elif ratio > merged_nm[0][0]:
                heapq.heapreplace(merged_nm, (ratio, serial, item))

        for item in value["top_lc"]:
            serial += 1
            ratio = item["lc_ratio"]
            if len(merged_lc) < lc_top_k:
                heapq.heappush(merged_lc, (ratio, serial, item))
            elif ratio > merged_lc[0][0]:
                heapq.heapreplace(merged_lc, (ratio, serial, item))

        if value.get("all_lc_failures") is not None:
            collect_all_lc = True
            all_lc_failures.extend(value["all_lc_failures"])

    top_nm = [d for _, _, d in sorted(merged_nm, key=lambda x: -x[0])]
    top_lc = [d for _, _, d in sorted(merged_lc, key=lambda x: -x[0])]

    return {
        "n": n,
        "workers": workers,
        "total_trees": agg_count,
        "non_unimodal_count": agg_non_unimodal,
        "lc_failure_count": agg_lc_fail,
        "worst_lc_ratio": global_worst_lc_ratio,
        "worst_lc_item": global_worst_lc_item,
        "top_near_misses": top_nm,
        "top_lc_failures": top_lc,
        "all_lc_failures": all_lc_failures if collect_all_lc else None,
        "platform": "Modal",
        "dict_name": dict_name("lc_nm", n),
        "date": time.strftime("%Y-%m-%d"),
        "complete": len(rows) == workers,
        "completed_partitions": len(rows),
    }


def format_status(summary: dict[str, Any], kind: str) -> str:
    lines = [
        f"kind={kind}",
        f"n={summary['n']}",
        f"workers={summary['workers']}",
        f"completed_partitions={summary['completed'] if 'completed' in summary else summary['completed_partitions']}",
        f"complete={summary['complete']}",
    ]
    if kind == "unimodality":
        lines.extend(
            [
                f"total_trees_so_far={summary['total_trees_so_far']}",
                f"counterexamples_found={summary['counterexamples_found']}",
                f"dict_name={summary['dict_name']}",
            ]
        )
    else:
        lines.extend(
            [
                f"total_trees={summary['total_trees']}",
                f"non_unimodal_count={summary['non_unimodal_count']}",
                f"lc_failure_count={summary['lc_failure_count']}",
                f"best_nm={summary['top_near_misses'][0]['nm_ratio'] if summary['top_near_misses'] else None}",
                f"dict_name={summary['dict_name']}",
            ]
        )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Collect Modal dict-backed result sets")
    sub = parser.add_subparsers(dest="command", required=True)

    for command in ("status", "collect"):
        p = sub.add_parser(command)
        p.add_argument("--kind", choices=("unimodality", "lc_nm"), required=True)
        p.add_argument("--n", type=int, required=True)
        p.add_argument("--workers", type=int, default=1024)
        p.add_argument("--expected", type=int, default=None)
        p.add_argument("--top-k", type=int, default=200)
        p.add_argument("--lc-top-k", type=int, default=200)
        p.add_argument("--dict-name", default="")
        p.add_argument("--out", default="")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    expected = resolve_expected(args.n, args.expected) if args.kind == "unimodality" else args.expected
    name = args.dict_name or dict_name(args.kind, args.n)
    rows = load_modal_dict(name, max(args.workers, 2000))

    if args.kind == "unimodality":
        status = summarize_unimodality(rows, args.n, args.workers, expected)
        if args.command == "status":
            print(format_status(status, args.kind))
            return 0
        summary = collect_unimodality(rows, args.n, args.workers, expected)
    else:
        summary = collect_lc_nm(rows, args.n, args.workers, args.top_k, args.lc_top_k)
        if args.command == "status":
            print(format_status(summary, args.kind))
            return 0

    if args.out:
        out_path = Path(args.out)
        out_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"saved={out_path}")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
