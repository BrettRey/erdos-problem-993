#!/usr/bin/env python3
"""Scan weighted Hall slacks over non-empty subsets S subseteq H.

Setup (d_leaf <= 1 tree):
  H = {h : P(h) > 1/3}
  U = N(H)
  demand(h) = P(h) - 1/3
  supply(u) = 1/3 - P(u)

Weighted Hall slack for non-empty S subseteq H:
  slack(S) = supply(N(S)) - demand(S).

This script brute-forces all non-empty subsets of H (bitmask DP) and records:
  - minimum non-empty slack
  - argmin subset size
  - whether argmin can be chosen singleton
  - failures (negative slack)
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter


def parse_graph6(line: bytes) -> tuple[int, list[list[int]]]:
    s = line.decode("ascii").strip()
    data = [ord(c) - 63 for c in s]
    n = data[0]
    idx = 1
    adj = [[] for _ in range(n)]
    bit = 5
    word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5
                idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj


def is_dleaf_le_1(n: int, adj: list[list[int]]) -> bool:
    leaves = {v for v in range(n) if len(adj[v]) == 1}
    for v in range(n):
        leaf_children = sum(1 for u in adj[v] if u in leaves)
        if leaf_children > 1:
            return False
    return True


def hard_core_probs(n: int, adj: list[list[int]]) -> list[float]:
    if n == 1:
        return [0.5]

    parent = [-1] * n
    children = [[] for _ in range(n)]
    order: list[int] = []
    queue = [0]
    seen = [False] * n
    seen[0] = True
    qi = 0
    while qi < len(queue):
        v = queue[qi]
        qi += 1
        order.append(v)
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    r_up = [0.0] * n
    for v in reversed(order):
        prod = 1.0
        for c in children[v]:
            prod *= 1.0 + r_up[c]
        r_up[v] = 1.0 / prod

    r_down = [0.0] * n
    for v in order:
        for c in children[v]:
            prod = 1.0
            if parent[v] != -1:
                prod *= 1.0 + r_down[v]
            for s in children[v]:
                if s != c:
                    prod *= 1.0 + r_up[s]
            r_down[c] = 1.0 / prod

    p = [0.0] * n
    for v in range(n):
        prod = 1.0
        if parent[v] != -1:
            prod *= 1.0 + r_down[v]
        for c in children[v]:
            prod *= 1.0 + r_up[c]
        r_full = 1.0 / prod
        p[v] = r_full / (1.0 + r_full)
    return p


def min_nonempty_hall_slack(
    adj: list[list[int]],
    probs: list[float],
    tol: float,
) -> tuple[float, int, float]:
    one_third = 1.0 / 3.0
    h_nodes = [v for v, pv in enumerate(probs) if pv > one_third + tol]
    if not h_nodes:
        return 0.0, 0, 0.0

    h_set = set(h_nodes)
    u_nodes = sorted({u for h in h_nodes for u in adj[h] if u not in h_set})
    u_idx = {u: i for i, u in enumerate(u_nodes)}

    dem = [probs[h] - one_third for h in h_nodes]
    sup = [one_third - probs[u] for u in u_nodes]

    h_masks = []
    for h in h_nodes:
        m = 0
        for u in adj[h]:
            if u in u_idx:
                m |= 1 << u_idx[u]
        h_masks.append(m)

    u_pow = 1 << len(u_nodes)
    sup_sum = [0.0] * u_pow
    for m in range(1, u_pow):
        lsb = m & -m
        b = lsb.bit_length() - 1
        sup_sum[m] = sup_sum[m ^ lsb] + sup[b]

    h_pow = 1 << len(h_nodes)
    dem_sum = [0.0] * h_pow
    nbr_mask = [0] * h_pow

    min_slack = float("inf")
    argmin_size = 0
    min_singleton = float("inf")

    for i, d in enumerate(dem):
        singleton_mask = 1 << i
        m = h_masks[i]
        sl = sup_sum[m] - d
        if sl < min_singleton:
            min_singleton = sl

    for s in range(1, h_pow):
        lsb = s & -s
        b = lsb.bit_length() - 1
        prev = s ^ lsb
        dem_sum[s] = dem_sum[prev] + dem[b]
        nbr_mask[s] = nbr_mask[prev] | h_masks[b]
        sl = sup_sum[nbr_mask[s]] - dem_sum[s]
        if sl < min_slack:
            min_slack = sl
            argmin_size = s.bit_count()

    return min_slack, argmin_size, min_singleton


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exhaustive non-empty weighted Hall subset scan on H."
    )
    parser.add_argument("--min-n", type=int, default=4)
    parser.add_argument("--max-n", type=int, default=20)
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
    total_dleaf = 0
    hall_fail = 0
    singleton_best = 0
    singleton_not_best = 0
    argmin_size_dist: Counter[int] = Counter()
    global_min_slack = float("inf")
    global_min_singleton = float("inf")
    worst = None
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("WHNC weighted Hall subset scan")
    print("Checking min_{nonempty S subseteq H} supply(N(S)) - demand(S)", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_dleaf = 0
        n_seen = 0
        n_fail = 0
        n_singleton_best = 0
        n_singleton_not_best = 0

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj = parse_graph6(line)
            n_seen += 1
            total_seen += 1
            if (not args.all_trees) and (not is_dleaf_le_1(n0, adj)):
                continue

            n_dleaf += 1
            total_dleaf += 1
            probs = hard_core_probs(n0, adj)
            min_slack, argmin_size, min_singleton = min_nonempty_hall_slack(
                adj, probs, args.tol
            )

            argmin_size_dist[argmin_size] += 1
            if min_slack < -1e-10:
                n_fail += 1
                hall_fail += 1
            if min_singleton <= min_slack + 1e-12:
                n_singleton_best += 1
                singleton_best += 1
            else:
                n_singleton_not_best += 1
                singleton_not_best += 1

            if min_slack < global_min_slack:
                global_min_slack = min_slack
                global_min_singleton = min_singleton
                worst = {
                    "n": n0,
                    "g6": line.decode("ascii").strip(),
                    "min_slack": min_slack,
                    "argmin_size": argmin_size,
                    "min_singleton": min_singleton,
                }

        proc.wait()
        rate = (
            n_singleton_best / n_dleaf if n_dleaf else 0.0
        )
        print(
            f"n={n:2d}: seen={n_seen:9d} considered={n_dleaf:8d} "
            f"hall_fail={n_fail:6d} "
            f"singleton_best={n_singleton_best:8d} ({100*rate:6.2f}%) "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )
        per_n[str(n)] = {
            "seen": n_seen,
            "dleaf": n_dleaf,
            "hall_fail": n_fail,
            "singleton_best": n_singleton_best,
            "singleton_not_best": n_singleton_not_best,
            "singleton_best_pct": 100.0 * rate,
            "elapsed_s": time.time() - t0,
        }

    print("-" * 80, flush=True)
    print(
        f"TOTAL seen={total_seen:,} considered={total_dleaf:,} "
        f"hall_fail={hall_fail} "
        f"singleton_best={singleton_best:,} singleton_not_best={singleton_not_best:,} "
        f"wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Global minimum non-empty slack: {global_min_slack:.12f}", flush=True)
    print(f"Global minimum singleton slack: {global_min_singleton:.12f}", flush=True)
    print("Argmin size distribution:", flush=True)
    for k in sorted(argmin_size_dist):
        v = argmin_size_dist[k]
        pct = 100.0 * v / total_dleaf if total_dleaf else 0.0
        print(f"  size={k:2d}: {v:10d} ({pct:6.2f}%)", flush=True)
    if worst is not None:
        print("Worst example:", worst, flush=True)

    if args.out:
        payload = {
            "config": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "geng": args.geng,
                "tol": args.tol,
                "all_trees": args.all_trees,
            },
            "summary": {
                "total_seen": total_seen,
                "total_dleaf": total_dleaf,
                "hall_fail": hall_fail,
                "singleton_best": singleton_best,
                "singleton_not_best": singleton_not_best,
                "global_min_slack": global_min_slack,
                "global_min_singleton": global_min_singleton,
                "worst": worst,
                "wall_s": time.time() - t_all,
            },
            "argmin_size_dist": dict(argmin_size_dist),
            "per_n": per_n,
        }
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Saved {args.out}", flush=True)


if __name__ == "__main__":
    main()
