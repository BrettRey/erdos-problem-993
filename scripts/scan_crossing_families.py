#!/usr/bin/env python3
"""Deterministic family scan for post-descent crossing pressure."""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indpoly import independence_poly, is_log_concave, is_unimodal, log_concavity_ratio, near_miss_ratio  # noqa: E402
from scripts.analyze_prufer_corpus import graph6_from_adj, mode_interval  # noqa: E402
from targeted import make_broom  # noqa: E402


def parse_int_values(spec: str) -> list[int]:
    values: set[int] = set()
    for raw_part in spec.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_s, end_s = part.split("-", 1)
            start = int(start_s)
            end = int(end_s)
            if start > end:
                raise ValueError(f"empty descending range {part!r}")
            values.update(range(start, end + 1))
        else:
            values.add(int(part))
    return sorted(values)


def make_multiarm_star(s: int, arms: tuple[int, ...]) -> tuple[int, list[list[int]]]:
    """Build M(s; arms): hub with s leaves and path arms of given lengths."""
    if s < 0:
        raise ValueError("s must be nonnegative")
    n = 1 + s + sum(arms)
    adj: list[list[int]] = [[] for _ in range(n)]
    hub = 0
    v = 1
    for _ in range(s):
        adj[hub].append(v)
        adj[v].append(hub)
        v += 1
    for arm_len in arms:
        prev = hub
        for _ in range(arm_len):
            adj[prev].append(v)
            adj[v].append(prev)
            prev = v
            v += 1
    return n, adj


def first_descent(poly: list[int]) -> int | None:
    for k in range(1, len(poly)):
        if poly[k] < poly[k - 1]:
            return k
    return None


def post_descent_upward_count(poly: list[int]) -> int:
    descent = first_descent(poly)
    if descent is None:
        return 0
    return sum(1 for k in range(descent, len(poly) - 1) if poly[k + 1] > poly[k])


def evaluate_row(
    *,
    family: str,
    label: str,
    params: dict[str, Any],
    adj: list[list[int]],
    include_graph6: bool,
) -> dict[str, Any]:
    n = len(adj)
    poly = independence_poly(n, adj)
    pressure, pos = near_miss_ratio(poly)
    mode_first, mode_last = mode_interval(poly)
    lc_ratio, lc_pos = log_concavity_ratio(poly)
    upward_count = post_descent_upward_count(poly)
    return {
        "family": family,
        "label": label,
        "params": params,
        "n": n,
        "alpha": len(poly) - 1,
        "mode_first": mode_first,
        "mode_last": mode_last,
        "first_descent": first_descent(poly),
        "crossing_pressure": pressure,
        "crossing_pos": pos,
        "crossing_reserve": 1.0 - pressure,
        "post_descent_upward_count": upward_count,
        "unimodal": is_unimodal(poly),
        "log_concave": is_log_concave(poly),
        "lc_ratio": lc_ratio,
        "lc_pos": lc_pos,
        "graph6": graph6_from_adj(adj) if include_graph6 else None,
    }


def arm_tuples(values: list[int], max_arms: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for size in range(2, max_arms + 1):
        out.extend(itertools.combinations_with_replacement(values, size))
    return out


def summarize_numeric(rows: list[dict[str, Any]], key: str) -> dict[str, float | int | None]:
    values = sorted(row[key] for row in rows if row[key] is not None)
    if not values:
        return {"min": None, "max": None, "mean": None, "median": None}
    return {
        "min": values[0],
        "max": values[-1],
        "mean": sum(values) / len(values),
        "median": values[len(values) // 2],
    }


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "family",
        "label",
        "params",
        "n",
        "alpha",
        "mode_first",
        "first_descent",
        "crossing_pressure",
        "crossing_pos",
        "crossing_reserve",
        "post_descent_upward_count",
        "unimodal",
        "log_concave",
        "lc_ratio",
        "lc_pos",
        "graph6",
    ]
    return {key: row[key] for key in keys}


def build_rows(
    *,
    min_n: int,
    max_n: int,
    step_n: int,
    broom_max_path: int,
    arm_values: list[int],
    max_arms: int,
    include_graph6: bool,
) -> list[dict[str, Any]]:
    rows = []
    arms_grid = arm_tuples(arm_values, max_arms)
    for n in range(min_n, max_n + 1, step_n):
        for path_len in range(2, min(n, broom_max_path + 1)):
            star_size = n - path_len
            if star_size < 1:
                continue
            _, adj = make_broom(path_len, star_size)
            rows.append(
                evaluate_row(
                    family="broom",
                    label=f"broom({path_len},{star_size})",
                    params={"path_len": path_len, "star_size": star_size},
                    adj=adj,
                    include_graph6=include_graph6,
                )
            )
        for arms in arms_grid:
            s = n - 1 - sum(arms)
            if s < 1:
                continue
            _, adj = make_multiarm_star(s, arms)
            rows.append(
                evaluate_row(
                    family="multiarm_star",
                    label=f"M({s};{','.join(map(str, arms))})",
                    params={"s": s, "arms": list(arms)},
                    adj=adj,
                    include_graph6=include_graph6,
                )
            )
    return rows


def build_summary(
    rows: list[dict[str, Any]],
    *,
    min_n: int,
    max_n: int,
    step_n: int,
    broom_max_path: int,
    arm_values: list[int],
    max_arms: int,
    top: int,
) -> dict[str, Any]:
    return {
        "source": {
            "kind": "crossing_family_scan",
            "min_n": min_n,
            "max_n": max_n,
            "step_n": step_n,
            "broom_max_path": broom_max_path,
            "arm_values": arm_values,
            "max_arms": max_arms,
        },
        "processed": len(rows),
        "counts": {
            "non_unimodal": sum(1 for row in rows if not row["unimodal"]),
            "non_log_concave": sum(1 for row in rows if not row["log_concave"]),
            "post_descent_upward_rows": sum(
                1 for row in rows if row["post_descent_upward_count"] > 0
            ),
            "broom_rows": sum(1 for row in rows if row["family"] == "broom"),
            "multiarm_star_rows": sum(1 for row in rows if row["family"] == "multiarm_star"),
        },
        "numeric": {
            "crossing_pressure": summarize_numeric(rows, "crossing_pressure"),
            "crossing_reserve": summarize_numeric(rows, "crossing_reserve"),
            "lc_ratio": summarize_numeric(rows, "lc_ratio"),
        },
        "top_by_crossing_pressure": [
            compact(row)
            for row in sorted(rows, key=lambda row: row["crossing_pressure"], reverse=True)[:top]
        ],
        "top_by_lc_ratio": [
            compact(row) for row in sorted(rows, key=lambda row: row["lc_ratio"], reverse=True)[:top]
        ],
        "post_descent_upward_rows": [
            compact(row) for row in rows if row["post_descent_upward_count"] > 0
        ][:top],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-n", type=int, default=50)
    parser.add_argument("--max-n", type=int, default=500)
    parser.add_argument("--step-n", type=int, default=10)
    parser.add_argument("--broom-max-path", type=int, default=80)
    parser.add_argument("--arm-values", default="2-12")
    parser.add_argument("--max-arms", type=int, default=4)
    parser.add_argument("--include-graph6", action="store_true")
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/crossing_family_scan_2026-07-03.json"),
    )
    args = parser.parse_args()
    if args.min_n < 2 or args.max_n < args.min_n or args.step_n < 1:
        parser.error("invalid n range")
    if args.broom_max_path < 2:
        parser.error("--broom-max-path must be at least 2")
    if args.max_arms < 2:
        parser.error("--max-arms must be at least 2")
    arm_values = parse_int_values(args.arm_values)

    rows = build_rows(
        min_n=args.min_n,
        max_n=args.max_n,
        step_n=args.step_n,
        broom_max_path=args.broom_max_path,
        arm_values=arm_values,
        max_arms=args.max_arms,
        include_graph6=args.include_graph6,
    )
    summary = build_summary(
        rows,
        min_n=args.min_n,
        max_n=args.max_n,
        step_n=args.step_n,
        broom_max_path=args.broom_max_path,
        arm_values=arm_values,
        max_arms=args.max_arms,
        top=args.top,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Scanned {summary['processed']} family rows; "
        f"non-unimodal={summary['counts']['non_unimodal']}, "
        f"post-descent-upward={summary['counts']['post_descent_upward_rows']}, "
        f"best-pressure={summary['numeric']['crossing_pressure']['max']}, "
        f"wrote {args.out}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
