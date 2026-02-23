#!/usr/bin/env python3
"""Summarize local branch-shape structure of Route-1 transfer violations.

Input files are the JSON outputs from
`verify_route1_transfer_failures_2026_02_19.py`.
For each witness, we compute:
  - rooted profile of B = T-{l,s} at u,
  - branch multiset under the unique/selected child of u,
  - compact canonical signatures for those branch subtrees.
"""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from typing import Any

from graph6 import parse_graph6
from verify_route1_transfer_failures_2026_02_19 import remove_vertices


def rooted_children(adj: list[list[int]], root: int) -> tuple[list[list[int]], list[int]]:
    parent = [-1] * len(adj)
    parent[root] = root
    children = [[] for _ in adj]
    queue = [root]
    for v in queue:
        for w in adj[v]:
            if parent[w] == -1:
                parent[w] = v
                children[v].append(w)
                queue.append(w)
    return children, parent


def encode_shape(children: list[list[int]], v: int) -> tuple:
    """AHU-style rooted tree shape encoding."""
    return tuple(sorted(encode_shape(children, c) for c in children[v]))


def subtree_size(children: list[list[int]], root: int) -> int:
    total = 0
    stack = [root]
    while stack:
        v = stack.pop()
        total += 1
        stack.extend(children[v])
    return total


def subtree_height(children: list[list[int]], root: int) -> int:
    max_depth = 0
    stack = [(root, 0)]
    while stack:
        v, d = stack.pop()
        if not children[v]:
            max_depth = max(max_depth, d)
        for c in children[v]:
            stack.append((c, d + 1))
    return max_depth


def depth_hist(children: list[list[int]], root: int) -> dict[str, int]:
    depth = {root: 0}
    queue = [root]
    for v in queue:
        for c in children[v]:
            depth[c] = depth[v] + 1
            queue.append(c)
    out: dict[str, int] = {}
    for d in depth.values():
        key = str(d)
        out[key] = out.get(key, 0) + 1
    return out


def summarize_witness(w: dict[str, Any]) -> dict[str, Any]:
    nn, adj = parse_graph6(w["g6"].encode())
    b_adj, idx = remove_vertices(adj, {w["leaf"], w["support"]})
    u = idx[w["u_in_T"]]
    children, _ = rooted_children(b_adj, u)

    root_children = children[u]
    selected_child = root_children[0] if root_children else None
    if selected_child is None:
        branch_rows: list[dict[str, Any]] = []
    else:
        branch_rows = []
        for gc in children[selected_child]:
            shape = encode_shape(children, gc)
            row = {
                "gc": gc,
                "size": subtree_size(children, gc),
                "height": subtree_height(children, gc),
                "shape": str(shape),
            }
            branch_rows.append(row)

    branch_rows.sort(key=lambda r: (r["size"], r["height"], r["shape"]))
    branch_counter = Counter((r["size"], r["height"], r["shape"]) for r in branch_rows)
    branch_multiset = [
        {
            "size": size,
            "height": height,
            "shape": shape,
            "count": count,
        }
        for (size, height, shape), count in sorted(branch_counter.items())
    ]

    return {
        "g6": w["g6"],
        "n": w["n"],
        "mode": w["mode"],
        "lambda": w["lambda"],
        "D": w["D"],
        "threshold_1_over_1plam": w["threshold_1_over_1plam"],
        "exact_excess": w["exact_excess"],
        "deg_u_in_T": w["deg_u_in_T"],
        "deg_u_in_B": w["deg_u_in_B"],
        "num_children_u": w["num_children_u"],
        "p_u": w["p_u"],
        "sum_delta_children": w["sum_delta_children"],
        "one_minus_sum_delta": w["one_minus_sum_delta"],
        "depth_hist_B": depth_hist(children, u),
        "selected_child_deg_in_B": len(b_adj[selected_child]) if selected_child is not None else 0,
        "selected_child_branch_count": len(branch_rows),
        "selected_child_branches": branch_rows,
        "selected_child_branch_multiset": branch_multiset,
    }


def load_witnesses(paths: list[str]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for p in paths:
        with open(p, "r", encoding="utf-8") as f:
            payload = json.load(f)
        rows = payload.get("witnesses", [])
        out.extend(rows)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize branch-shape signatures of Route-1 exact failures.")
    ap.add_argument("--in", dest="inputs", nargs="+", required=True, help="Input JSON files from failure extractor.")
    ap.add_argument("--out", required=True, help="Output JSON report path.")
    args = ap.parse_args()

    witnesses = load_witnesses(args.inputs)
    witnesses.sort(key=lambda w: w["exact_excess"], reverse=True)
    summaries = [summarize_witness(w) for w in witnesses]

    payload = {
        "inputs": args.inputs,
        "count": len(summaries),
        "summaries": summaries,
    }

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"summaries={len(summaries)}")
    for s in summaries:
        print(
            f"n={s['n']} exact_excess={s['exact_excess']:.12f} "
            f"deg_u_B={s['deg_u_in_B']} child_branches={s['selected_child_branch_count']} "
            f"g6={s['g6']}"
        )
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
