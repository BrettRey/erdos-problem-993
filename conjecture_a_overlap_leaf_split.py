#!/usr/bin/env python3
"""Split local overlap-control failures by leaf-heavy adjacency.

For u in N(H), local overlap-control checks:
  local_excess(u) = sum_{h in N(u) cap H} (P(h) - 1/3)
  deficit(u)      = 1/3 - P(u)

This script counts failures local_excess(u) > deficit(u) and splits by whether
u has at least one heavy leaf neighbor.
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
        description="Leaf-split diagnostic for local overlap-control failures."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    one_third = 1.0 / 3.0
    total_seen = 0
    total_considered = 0
    fail_total = 0
    fail_with_leafheavy = 0
    fail_without_leafheavy = 0

    first_with = None
    first_without = None

    t_all = time.time()
    print("Local overlap leaf-split scan")
    print("Checking local_excess(u) > deficit(u) over u in N(H)")
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_fail = 0

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1
            if not is_dleaf_le_1(n0, adj):
                continue

            n_considered += 1
            total_considered += 1
            probs = hard_core_probs(n0, adj)
            h_set = {v for v, pv in enumerate(probs) if pv > one_third + args.tol}
            if not h_set:
                continue

            u_nodes = {u for h in h_set for u in adj[h] if u not in h_set}
            for u in u_nodes:
                hs = [h for h in adj[u] if h in h_set]
                if not hs:
                    continue
                deficit = one_third - probs[u]
                excess = sum(probs[h] - one_third for h in hs)
                if excess > deficit + args.tol:
                    n_fail += 1
                    fail_total += 1
                    has_leafheavy = any(len(adj[h]) == 1 for h in hs)
                    ratio = excess / deficit if deficit > args.tol else float("inf")

                    if has_leafheavy:
                        fail_with_leafheavy += 1
                        if first_with is None:
                            first_with = {
                                "n": n0,
                                "g6": line.decode("ascii").strip(),
                                "u": u,
                                "heavy_neighbors": hs,
                                "heavy_degrees": [len(adj[h]) for h in hs],
                                "ratio": ratio,
                            }
                    else:
                        fail_without_leafheavy += 1
                        if first_without is None:
                            first_without = {
                                "n": n0,
                                "g6": line.decode("ascii").strip(),
                                "u": u,
                                "heavy_neighbors": hs,
                                "heavy_degrees": [len(adj[h]) for h in hs],
                                "ratio": ratio,
                            }

        proc.wait()
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"local_fail={n_fail:7d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"fail_total={fail_total:,} with_leafheavy={fail_with_leafheavy:,} "
        f"without_leafheavy={fail_without_leafheavy:,} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"First with_leafheavy: {first_with}", flush=True)
    print(f"First without_leafheavy: {first_without}", flush=True)

    if args.out:
        out_dir = os.path.dirname(args.out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        payload = {
            "params": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "tol": args.tol,
            },
            "summary": {
                "seen": total_seen,
                "considered": total_considered,
                "fail_total": fail_total,
                "fail_with_leafheavy": fail_with_leafheavy,
                "fail_without_leafheavy": fail_without_leafheavy,
                "first_with_leafheavy": first_with,
                "first_without_leafheavy": first_without,
            },
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
