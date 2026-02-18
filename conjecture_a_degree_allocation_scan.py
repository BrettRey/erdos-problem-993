#!/usr/bin/env python3
"""Test degree-only fractional allocation for WHNC.

For heavy set H = {h : P(h) > 1/3}, consider edge weights x_{h,u} = 1/deg(h)
on heavy-neighbor edges (h,u), u in N(H).

This would prove WHNC if c(u) := sum_{h~u, h in H} 1/deg(h) <= 1 for all u in N(H),
since each heavy vertex receives unit total weight.

This script checks c(u) <= 1 exhaustively.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Degree-only allocation scan for heavy-neighborhood overlap."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument(
        "--all-trees",
        action="store_true",
        help="Do not filter by d_leaf<=1; scan all trees.",
    )
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_seen = 0
    total_considered = 0
    fail_trees = 0
    worst_c = -1.0
    worst = None

    t_all = time.time()
    print("Degree-only allocation scan")
    print("Checking c(u) = sum_{h~u,h in H} 1/deg(h) <= 1 for u in N(H)")
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_fail = 0
        n_worst = -1.0

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1

            if (not args.all_trees) and (not is_dleaf_le_1(n0, adj)):
                continue
            n_considered += 1
            total_considered += 1

            probs = hard_core_probs(n0, adj)
            h_nodes = [v for v, pv in enumerate(probs) if pv > 1.0 / 3.0 + args.tol]
            if not h_nodes:
                continue
            h_set = set(h_nodes)
            deg_h = {h: len(adj[h]) for h in h_nodes}

            bad_tree = False
            for u in {x for h in h_nodes for x in adj[h] if x not in h_set}:
                c = 0.0
                heavy_neighbors = []
                heavy_degrees = []
                for h in adj[u]:
                    if h in h_set:
                        c += 1.0 / deg_h[h]
                        heavy_neighbors.append(h)
                        heavy_degrees.append(deg_h[h])

                if c > n_worst:
                    n_worst = c
                if c > worst_c:
                    worst_c = c
                    worst = {
                        "n": n0,
                        "g6": line.decode("ascii").strip(),
                        "u": u,
                        "heavy_neighbors": heavy_neighbors,
                        "heavy_degrees": heavy_degrees,
                        "c_u": c,
                    }

                if c > 1.0 + args.tol:
                    bad_tree = True

            if bad_tree:
                n_fail += 1
                fail_trees += 1

        proc.wait()
        worst_str = "N/A" if n_worst < 0 else f"{n_worst:.6f}"
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"fail_trees={n_fail:7d} worst_c={worst_str:>8} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"fail_trees={fail_trees:,} worst_c={worst_c:.12f} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Worst witness: {worst}", flush=True)

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        payload = {
            "params": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "all_trees": args.all_trees,
                "tol": args.tol,
            },
            "summary": {
                "seen": total_seen,
                "considered": total_considered,
                "fail_trees": fail_trees,
                "worst_c": worst_c,
                "worst": worst,
            },
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
