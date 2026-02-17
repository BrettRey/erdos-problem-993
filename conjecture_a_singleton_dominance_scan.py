#!/usr/bin/env python3
"""Scan stronger singleton-dominance inequalities for WHNC.

Setup:
  H = {h : P(h) > 1/3}, U = N(H),
  demand(h) = P(h) - 1/3, supply(u) = 1/3 - P(u),
  slack(S) = supply(N(S)) - demand(S) for non-empty S subseteq H.

Checks:
  (A) global singleton argmin:
      min_{non-empty S subseteq H} slack(S) == min_{h in H} slack({h})
  (B) stronger local singleton dominance:
      slack(S) >= min_{h in S} slack({h}) for every non-empty S subseteq H
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time

from conjecture_a_hall_subset_scan import hard_core_probs, is_dleaf_le_1, parse_graph6


def scan_tree(
    adj: list[list[int]],
    probs: list[float],
    tol: float,
) -> dict:
    one_third = 1.0 / 3.0
    h_nodes = [v for v, pv in enumerate(probs) if pv > one_third + tol]
    if not h_nodes:
        return {
            "m": 0,
            "global_ok": True,
            "local_ok": True,
            "min_global_gap": 0.0,
            "min_local_gap": 0.0,
            "argmin_size": 0,
            "argmin_mask": 0,
            "global_min_slack": 0.0,
            "global_min_singleton": 0.0,
            "local_bad_mask": 0,
        }

    h_set = set(h_nodes)
    u_nodes = sorted({u for h in h_nodes for u in adj[h] if u not in h_set})
    u_idx = {u: i for i, u in enumerate(u_nodes)}

    dem = [probs[h] - one_third for h in h_nodes]
    sup = [one_third - probs[u] for u in u_nodes]

    h_masks: list[int] = []
    for h in h_nodes:
        m = 0
        for u in adj[h]:
            j = u_idx.get(u)
            if j is not None:
                m |= 1 << j
        h_masks.append(m)

    u_pow = 1 << len(u_nodes)
    sup_sum = [0.0] * u_pow
    for s in range(1, u_pow):
        lsb = s & -s
        b = lsb.bit_length() - 1
        sup_sum[s] = sup_sum[s ^ lsb] + sup[b]

    m = len(h_nodes)
    h_pow = 1 << m
    dem_sum = [0.0] * h_pow
    nbr_mask = [0] * h_pow
    min_singleton_in = [float("inf")] * h_pow

    singleton_slack = [0.0] * m
    min_singleton = float("inf")
    for i in range(m):
        sl = sup_sum[h_masks[i]] - dem[i]
        singleton_slack[i] = sl
        if sl < min_singleton:
            min_singleton = sl

    min_slack = float("inf")
    argmin_size = 0
    argmin_mask = 0

    min_global_gap = float("inf")
    min_local_gap = float("inf")
    local_bad_mask = 0
    local_ok = True

    for s in range(1, h_pow):
        lsb = s & -s
        b = lsb.bit_length() - 1
        prev = s ^ lsb

        dem_sum[s] = dem_sum[prev] + dem[b]
        nbr_mask[s] = nbr_mask[prev] | h_masks[b]
        min_singleton_in[s] = min(min_singleton_in[prev], singleton_slack[b])

        sl = sup_sum[nbr_mask[s]] - dem_sum[s]
        if sl < min_slack:
            min_slack = sl
            argmin_size = s.bit_count()
            argmin_mask = s

        ggap = sl - min_singleton
        if ggap < min_global_gap:
            min_global_gap = ggap

        lgap = sl - min_singleton_in[s]
        if lgap < min_local_gap:
            min_local_gap = lgap
            if lgap < -tol:
                local_bad_mask = s
        if lgap < -tol:
            local_ok = False

    global_ok = min_global_gap >= -tol
    return {
        "m": m,
        "global_ok": global_ok,
        "local_ok": local_ok,
        "min_global_gap": min_global_gap,
        "min_local_gap": min_local_gap,
        "argmin_size": argmin_size,
        "argmin_mask": argmin_mask,
        "global_min_slack": min_slack,
        "global_min_singleton": min_singleton,
        "local_bad_mask": local_bad_mask,
    }


def decode_subset(mask: int, h_nodes: list[int]) -> list[int]:
    out: list[int] = []
    rem = mask
    while rem:
        lsb = rem & -rem
        i = lsb.bit_length() - 1
        rem ^= lsb
        out.append(h_nodes[i])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan singleton-dominance properties over heavy subsets."
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
    global_fail = 0
    local_fail = 0

    global_worst = None
    local_worst = None
    argmin_size_dist: dict[int, int] = {}

    t_all = time.time()
    print("Singleton-dominance scan")
    print("Checks:")
    print("  (A) min_{S!=empty} slack(S) == min_{h in H} slack({h})")
    print("  (B) slack(S) >= min_{h in S} slack({h}) for all non-empty S")
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_seen = 0
        n_considered = 0
        n_global_fail = 0
        n_local_fail = 0

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
            res = scan_tree(adj, probs, args.tol)

            argmin_size = res["argmin_size"]
            argmin_size_dist[argmin_size] = argmin_size_dist.get(argmin_size, 0) + 1

            if not res["global_ok"]:
                n_global_fail += 1
                global_fail += 1
            if not res["local_ok"]:
                n_local_fail += 1
                local_fail += 1

            g6 = line.decode("ascii").strip()
            h_nodes = [v for v, pv in enumerate(probs) if pv > 1.0 / 3.0 + args.tol]

            if (global_worst is None) or (res["min_global_gap"] < global_worst["min_global_gap"]):
                global_worst = {
                    "n": n0,
                    "g6": g6,
                    "H": h_nodes,
                    "argmin_size": res["argmin_size"],
                    "argmin_subset": decode_subset(res["argmin_mask"], h_nodes),
                    "min_slack": res["global_min_slack"],
                    "min_singleton": res["global_min_singleton"],
                    "min_global_gap": res["min_global_gap"],
                }

            if (local_worst is None) or (res["min_local_gap"] < local_worst["min_local_gap"]):
                local_worst = {
                    "n": n0,
                    "g6": g6,
                    "H": h_nodes,
                    "argmin_size": res["argmin_size"],
                    "argmin_subset": decode_subset(res["argmin_mask"], h_nodes),
                    "local_bad_subset": decode_subset(res["local_bad_mask"], h_nodes),
                    "min_slack": res["global_min_slack"],
                    "min_singleton": res["global_min_singleton"],
                    "min_local_gap": res["min_local_gap"],
                }

        proc.wait()
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_considered:8d} "
            f"global_fail={n_global_fail:7d} local_fail={n_local_fail:7d} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_considered:,} "
        f"global_fail={global_fail:,} local_fail={local_fail:,} "
        f"wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    if global_worst is not None:
        print("Worst global gap witness:", global_worst, flush=True)
    if local_worst is not None:
        print("Worst local gap witness:", local_worst, flush=True)
    print(f"argmin_size_dist={argmin_size_dist}", flush=True)

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
                "global_fail": global_fail,
                "local_fail": local_fail,
                "global_worst": global_worst,
                "local_worst": local_worst,
                "argmin_size_dist": argmin_size_dist,
            },
        }
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
