#!/usr/bin/env python3
"""Check local overlap control for WHNC in d_leaf <= 1 trees.

For H = {v : P(v) > 1/3} and u in N(H), define:
  deficit(u) = 1/3 - P(u)
  local_excess(u) = sum_{h in N(u) cap H} (P(h) - 1/3)

This script checks whether local_excess(u) <= deficit(u) always holds.
If true for every u in N(H), WHNC follows immediately by summing over u.

It also tracks overlap stress metrics requested in the research note:
  - ratio_excess = local_excess / deficit
  - ratio_count = |N(u) cap H| / deficit
"""

from __future__ import annotations

import argparse
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
    q_idx = 0

    while q_idx < len(queue):
        v = queue[q_idx]
        q_idx += 1
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

    probs = [0.0] * n
    for v in range(n):
        prod = 1.0
        if parent[v] != -1:
            prod *= 1.0 + r_down[v]
        for c in children[v]:
            prod *= 1.0 + r_up[c]
        r_full = 1.0 / prod
        probs[v] = r_full / (1.0 + r_full)
    return probs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Local overlap control test for WHNC on d_leaf <= 1 trees."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    one_third = 1.0 / 3.0

    total_dleaf = 0
    total_u = 0
    local_fail = 0
    worst_ratio_excess = -1.0
    worst_ratio_count = -1.0
    worst_example = None
    overlap_deg_dist: Counter[int] = Counter()

    t_all = time.time()
    print("WHNC local-overlap control scan")
    print("Checking local_excess(u) <= deficit(u) for u in N(H)", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_dleaf = 0
        n_u = 0
        n_fail = 0
        n_worst_ratio = -1.0

        assert proc.stdout is not None
        for line in proc.stdout:
            nn, adj = parse_graph6(line)
            if not is_dleaf_le_1(nn, adj):
                continue
            n_dleaf += 1
            total_dleaf += 1

            probs = hard_core_probs(nn, adj)
            heavy = {v for v, p in enumerate(probs) if p > one_third + args.tol}
            if not heavy:
                continue

            nh = set()
            for h in heavy:
                nh.update(adj[h])
            nh -= heavy

            for u in nh:
                hu = [h for h in adj[u] if h in heavy]
                if not hu:
                    continue
                n_u += 1
                total_u += 1
                overlap_deg_dist[len(hu)] += 1

                deficit = one_third - probs[u]
                excess = sum(probs[h] - one_third for h in hu)

                ratio_excess = excess / deficit if deficit > args.tol else float("inf")
                ratio_count = len(hu) / deficit if deficit > args.tol else float("inf")

                if ratio_excess > n_worst_ratio:
                    n_worst_ratio = ratio_excess
                if ratio_excess > worst_ratio_excess:
                    worst_ratio_excess = ratio_excess
                    worst_ratio_count = ratio_count
                    worst_example = {
                        "n": nn,
                        "g6": line.decode("ascii").strip(),
                        "u": u,
                        "heavy_neighbors": hu,
                        "P_u": probs[u],
                        "deficit": deficit,
                        "excess": excess,
                        "ratio_excess": ratio_excess,
                        "ratio_count": ratio_count,
                    }

                if excess > deficit + args.tol:
                    n_fail += 1
                    local_fail += 1

        proc.wait()
        ratio_str = "N/A" if n_worst_ratio < 0 else f"{n_worst_ratio:.6f}"
        print(
            f"n={n:2d}: d_leaf<=1={n_dleaf:8d}  u_in_N(H)={n_u:8d}  "
            f"fail={n_fail:6d}  worst_excess/deficit={ratio_str:>8}  "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL d_leaf<=1 trees: {total_dleaf:,}  total_u: {total_u:,}  "
        f"local_fail={local_fail}  wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Worst excess/deficit ratio: {worst_ratio_excess:.12f}", flush=True)
    print(f"Worst count/deficit ratio: {worst_ratio_count:.12f}", flush=True)
    print("Overlap degree distribution |N(u) cap H|:", flush=True)
    for k in sorted(overlap_deg_dist):
        cnt = overlap_deg_dist[k]
        pct = 100.0 * cnt / total_u if total_u else 0.0
        print(f"  {k}: {cnt:12d} ({pct:6.3f}%)", flush=True)

    if worst_example is not None:
        print("Worst local example:", flush=True)
        print(
            f"  n={worst_example['n']} g6={worst_example['g6']} u={worst_example['u']}",
            flush=True,
        )
        print(
            f"  heavy_neighbors={worst_example['heavy_neighbors']} "
            f"P_u={worst_example['P_u']:.12f}",
            flush=True,
        )
        print(
            f"  deficit={worst_example['deficit']:.12f} "
            f"excess={worst_example['excess']:.12f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
