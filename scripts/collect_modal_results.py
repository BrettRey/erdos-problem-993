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
  python3 scripts/collect_modal_results.py collect --kind unimodality --n 29 --workers 1024 \
    --raw-in results/raw_n29_modal_rows.json --out /tmp/replayed_n29.json
"""

from __future__ import annotations

import argparse
import hashlib
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


RAW_EXPORT_SCHEMA = "erdos993-modal-dict-rows-v1"


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


def _canonical_json_bytes(value: Any) -> bytes:
    """Return a stable JSON encoding suitable for artifact hashes."""
    return json.dumps(
        value,
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _partition_sort_key(key: str) -> tuple[int, int, str]:
    try:
        residue_s, workers_s = key.split("/", 1)
        return (int(workers_s), int(residue_s), key)
    except (ValueError, TypeError):
        return (sys.maxsize, sys.maxsize, key)


def expected_partition_keys(workers: int) -> list[str]:
    if workers < 1:
        raise ValueError("workers must be positive")
    return [f"{residue}/{workers}" for residue in range(workers)]


def audit_partition_keys(rows: dict[str, Any], workers: int) -> dict[str, Any]:
    """Check the exact residue-key set, including stale or foreign keys."""
    expected = expected_partition_keys(workers)
    expected_set = set(expected)
    actual = sorted(rows, key=_partition_sort_key)
    actual_set = set(actual)
    missing = [key for key in expected if key not in actual_set]
    unexpected = [key for key in actual if key not in expected_set]
    present_expected = workers - len(missing)
    return {
        "expected_count": workers,
        "present_expected_count": present_expected,
        "actual_count": len(actual),
        "missing_keys": missing,
        "unexpected_keys": unexpected,
        "exact_match": not missing and not unexpected,
        "expected_keys_sha256": _sha256_json(expected),
        "actual_keys_sha256": _sha256_json(actual),
    }


def expected_rows_only(rows: dict[str, Any], workers: int) -> dict[str, Any]:
    """Exclude stale keys from status totals while retaining partial results."""
    return {key: rows[key] for key in expected_partition_keys(workers) if key in rows}


def require_exact_partition_keys(rows: dict[str, Any], workers: int) -> dict[str, Any]:
    audit = audit_partition_keys(rows, workers)
    if not audit["exact_match"]:
        missing = ", ".join(audit["missing_keys"][:5]) or "none"
        unexpected = ", ".join(audit["unexpected_keys"][:5]) or "none"
        raise ValueError(
            "Modal partition keys do not exactly match the expected residue set: "
            f"missing={len(audit['missing_keys'])} ({missing}); "
            f"unexpected={len(audit['unexpected_keys'])} ({unexpected})"
        )
    return audit


def sorted_raw_rows(rows: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"key": key, "value": rows[key]}
        for key in sorted(rows, key=_partition_sort_key)
    ]


def build_raw_export(
    rows: dict[str, Any],
    *,
    kind: str,
    n: int,
    workers: int,
    expected: int | None,
    name: str,
) -> dict[str, Any]:
    """Build a deterministic, self-hashing export of the raw Modal rows."""
    raw_rows = sorted_raw_rows(rows)
    payload = {
        "schema": RAW_EXPORT_SCHEMA,
        "kind": kind,
        "n": n,
        "workers": workers,
        "expected_tree_count": expected,
        "dict_name": name,
        "row_count": len(raw_rows),
        "partition_key_audit": audit_partition_keys(rows, workers),
        "rows_sha256": _sha256_json(raw_rows),
        "rows": raw_rows,
    }
    payload["content_sha256"] = _sha256_json(payload)
    return payload


def load_raw_export(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load and validate a deterministic raw-row export without Modal."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != RAW_EXPORT_SCHEMA:
        raise ValueError(f"unsupported raw export schema: {payload.get('schema')!r}")
    claimed_content_hash = payload.get("content_sha256")
    content_without_hash = {
        key: value for key, value in payload.items() if key != "content_sha256"
    }
    if claimed_content_hash != _sha256_json(content_without_hash):
        raise ValueError("raw export content_sha256 mismatch")
    for field, expected_type in (
        ("kind", str),
        ("n", int),
        ("workers", int),
        ("dict_name", str),
    ):
        if not isinstance(payload.get(field), expected_type):
            raise ValueError(f"raw export {field} has the wrong type")
    raw_rows = payload.get("rows")
    if not isinstance(raw_rows, list):
        raise ValueError("raw export rows must be a list")
    if payload.get("row_count") != len(raw_rows):
        raise ValueError("raw export row_count does not match its rows")
    if payload.get("rows_sha256") != _sha256_json(raw_rows):
        raise ValueError("raw export rows_sha256 mismatch")

    rows: dict[str, Any] = {}
    for row in raw_rows:
        if not isinstance(row, dict) or set(row) != {"key", "value"}:
            raise ValueError("each raw export row must contain exactly key and value")
        key = row["key"]
        if not isinstance(key, str):
            raise ValueError("raw export row key must be a string")
        if key in rows:
            raise ValueError(f"duplicate raw export key: {key}")
        rows[key] = row["value"]

    workers = payload.get("workers")
    if not isinstance(workers, int):
        raise ValueError("raw export workers must be an integer")
    calculated_audit = audit_partition_keys(rows, workers)
    if payload.get("partition_key_audit") != calculated_audit:
        raise ValueError("raw export partition_key_audit mismatch")
    if raw_rows != sorted_raw_rows(rows):
        raise ValueError("raw export rows are not in canonical partition order")
    return rows, payload


def summarize_unimodality(rows: dict[str, Any], n: int, workers: int, expected: int | None) -> dict[str, Any]:
    partition_audit = audit_partition_keys(rows, workers)
    rows = expected_rows_only(rows, workers)
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
        "complete": partition_audit["exact_match"],
        "partition_key_audit": partition_audit,
        "rows_sha256": _sha256_json(sorted_raw_rows(rows)),
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
        "partition_key_audit": status["partition_key_audit"],
        "rows_sha256": status["rows_sha256"],
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
    partition_audit = audit_partition_keys(rows, workers)
    rows = expected_rows_only(rows, workers)
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
        "complete": partition_audit["exact_match"],
        "completed_partitions": len(rows),
        "partition_key_audit": partition_audit,
        "rows_sha256": _sha256_json(sorted_raw_rows(rows)),
    }


def format_status(summary: dict[str, Any], kind: str) -> str:
    lines = [
        f"kind={kind}",
        f"n={summary['n']}",
        f"workers={summary['workers']}",
        f"completed_partitions={summary['completed'] if 'completed' in summary else summary['completed_partitions']}",
        f"complete={summary['complete']}",
        f"missing_partition_keys={len(summary['partition_key_audit']['missing_keys'])}",
        f"unexpected_partition_keys={len(summary['partition_key_audit']['unexpected_keys'])}",
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
        p.add_argument(
            "--raw-in",
            type=Path,
            help="Validate and replay a prior raw-row export instead of contacting Modal",
        )
        if command == "collect":
            p.add_argument(
                "--raw-out",
                type=Path,
                help="Write sorted raw Modal rows with hashes for offline replay",
            )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    expected = resolve_expected(args.n, args.expected) if args.kind == "unimodality" else args.expected
    name = args.dict_name or dict_name(args.kind, args.n)
    if args.raw_in:
        try:
            rows, raw_payload = load_raw_export(args.raw_in)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            parser.error(f"cannot validate --raw-in: {exc}")
        for field, requested in (
            ("kind", args.kind),
            ("n", args.n),
            ("workers", args.workers),
        ):
            if raw_payload.get(field) != requested:
                parser.error(
                    f"--raw-in metadata {field}={raw_payload.get(field)!r} "
                    f"does not match requested {requested!r}"
                )
        if args.dict_name and raw_payload.get("dict_name") != args.dict_name:
            parser.error("--raw-in dict_name does not match --dict-name")
        if raw_payload.get("expected_tree_count") != expected:
            parser.error(
                "--raw-in expected_tree_count does not match the resolved expected count"
            )
        name = raw_payload["dict_name"]
    else:
        rows = load_modal_dict(name, max(args.workers, 2000))

    if args.command == "collect":
        try:
            require_exact_partition_keys(rows, args.workers)
        except ValueError as exc:
            parser.error(str(exc))

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

    summary["dict_name"] = name

    if args.command == "collect" and args.raw_out:
        raw_export = build_raw_export(
            rows,
            kind=args.kind,
            n=args.n,
            workers=args.workers,
            expected=expected,
            name=name,
        )
        args.raw_out.parent.mkdir(parents=True, exist_ok=True)
        args.raw_out.write_text(
            json.dumps(raw_export, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        summary["raw_export"] = {
            "path": str(args.raw_out),
            "rows_sha256": raw_export["rows_sha256"],
            "content_sha256": raw_export["content_sha256"],
            "row_count": raw_export["row_count"],
        }
        print(f"saved_raw={args.raw_out}")

    if args.out:
        out_path = Path(args.out)
        out_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"saved={out_path}")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
