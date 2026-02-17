#!/usr/bin/env python3
"""Check the equal-share heavy allocation inequality.

For each heavy vertex h in H={v:P(v)>1/3}, test:

  P(h) - 1/3 <= sum_{u~h} (1/3 - P(u)) / deg_H(u),

where deg_H(u)=|N(u) cap H|.

This is a concrete candidate from the cavity-allocation approach.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time


def parse_graph6(line: bytes) -> tuple[int, list[list[int]], str]:
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
    return n, adj, s


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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check equal-share allocation inequality on d_leaf<=1 trees."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    one_third = 1.0 / 3.0
    total_dleaf = 0
    total_fail = 0
    worst_ratio = -1.0
    worst = None
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("Equal-share allocation scan")
    print("Testing ex(h) <= sum_{u~h} deficit(u)/deg_H(u)", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        n_dleaf = 0
        n_fail = 0
        n_worst = -1.0

        assert proc.stdout is not None
        for line in proc.stdout:
            n0, adj, g6 = parse_graph6(line)
            if not is_dleaf_le_1(n0, adj):
                continue
            n_dleaf += 1
            total_dleaf += 1

            p = hard_core_probs(n0, adj)
            h = [v for v, pv in enumerate(p) if pv > one_third + args.tol]
            if not h:
                continue

            h_set = set(h)
            deg_h = [sum(1 for x in adj[u] if x in h_set) for u in range(n0)]

            for hv in h:
                ex = p[hv] - one_third
                rhs = 0.0
                for u in adj[hv]:
                    if deg_h[u] > 0:
                        rhs += (one_third - p[u]) / deg_h[u]
                ratio = ex / rhs if rhs > 1e-15 else float("inf")

                if ratio > n_worst:
                    n_worst = ratio
                if ratio > worst_ratio:
                    worst_ratio = ratio
                    worst = {
                        "n": n0,
                        "g6": g6,
                        "h": hv,
                        "ratio": ratio,
                        "excess": ex,
                        "rhs": rhs,
                    }
                if ex > rhs + args.tol:
                    n_fail += 1
                    total_fail += 1

        proc.wait()
        per_n[str(n)] = {
            "dleaf": n_dleaf,
            "fail": n_fail,
            "worst_ratio": n_worst,
            "elapsed_s": time.time() - t0,
        }
        wstr = "N/A" if n_worst < 0 else f"{n_worst:.6f}"
        print(
            f"n={n:2d}: d_leaf<=1={n_dleaf:8d} fail={n_fail:8d} "
            f"worst_ratio={wstr:>8} ({time.time()-t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL d_leaf<=1={total_dleaf:,} equal_share_fail={total_fail} "
        f"worst_ratio={worst_ratio:.12f} wall={time.time()-t_all:.1f}s",
        flush=True,
    )
    if worst is not None:
        print("Worst example:", worst, flush=True)

    if args.out:
        payload = {
            "config": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "geng": args.geng,
                "tol": args.tol,
            },
            "summary": {
                "total_dleaf": total_dleaf,
                "equal_share_fail": total_fail,
                "worst_ratio": worst_ratio,
                "worst": worst,
                "wall_s": time.time() - t_all,
            },
            "per_n": per_n,
        }
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Saved {args.out}", flush=True)


if __name__ == "__main__":
    main()
