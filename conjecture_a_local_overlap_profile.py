#!/usr/bin/env python3
"""Profile local overlap behavior by heavy-neighbor type.

For d_leaf <= 1 trees, with
  H = {v : P(v) > 1/3},
  U = N(H) \\ H,
define for each u in U:
  deficit(u) = 1/3 - P(u),
  excess(u)  = sum_{h in N(u) cap H} (P(h) - 1/3),
  slack(u)   = deficit(u) - excess(u).

This script aggregates u-vertices by key:
  (deg_H(u), leaf_H(u)),
where
  deg_H(u)  = |N(u) cap H|,
  leaf_H(u) = number of heavy leaf neighbors of u.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import defaultdict

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Local-overlap profile by (deg_H(u), heavy-leaf count)."
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
    total_u = 0
    total_fail = 0
    per_n: dict[str, dict] = {}

    # key=(deg_H, leaf_H)
    stats: dict[tuple[int, int], dict[str, float | int]] = defaultdict(
        lambda: {
            "count": 0,
            "fail_count": 0,
            "min_slack": float("inf"),
            "max_ratio": 0.0,
        }
    )

    t_all = time.time()
    print("Local-overlap profile scan", flush=True)
    print("Aggregating by key (deg_H(u), heavy_leaf_count(u))", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_u = 0
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
                n_u += 1
                total_u += 1

                deg_h = len(hs)
                leaf_h = sum(1 for h in hs if len(adj[h]) == 1)
                key = (deg_h, leaf_h)

                deficit = one_third - probs[u]
                excess = sum(probs[h] - one_third for h in hs)
                slack = deficit - excess

                st = stats[key]
                st["count"] += 1
                if slack < st["min_slack"]:
                    st["min_slack"] = slack
                if deficit > args.tol:
                    ratio = excess / deficit
                    if ratio > st["max_ratio"]:
                        st["max_ratio"] = ratio
                else:
                    # Should not happen in this regime, but keep diagnostics robust.
                    st["max_ratio"] = float("inf")

                if slack < -args.tol:
                    st["fail_count"] += 1
                    n_fail += 1
                    total_fail += 1

        proc.wait()
        per_n[str(n)] = {
            "seen": n_seen,
            "considered": n_considered,
            "u_count": n_u,
            "fail_count": n_fail,
            "wall_s": time.time() - t0,
        }
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"u={n_u:9d} fail={n_fail:7d} ({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"u={total_u:,} fail={total_fail:,} wall={time.time() - t_all:.1f}s",
        flush=True,
    )

    stats_out = {}
    for key in sorted(stats):
        deg_h, leaf_h = key
        st = stats[key]
        count = int(st["count"])
        fail_count = int(st["fail_count"])
        fail_pct = 100.0 * fail_count / count if count else 0.0
        stats_out[f"degH={deg_h},leafH={leaf_h}"] = {
            "deg_H": deg_h,
            "leaf_H": leaf_h,
            "count": count,
            "fail_count": fail_count,
            "fail_pct": fail_pct,
            "min_slack": float(st["min_slack"]),
            "max_ratio": float(st["max_ratio"]),
        }

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
                "u_count": total_u,
                "fail_count": total_fail,
                "wall_s": time.time() - t_all,
            },
            "profile": stats_out,
            "per_n": per_n,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
