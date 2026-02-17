#!/usr/bin/env python3
"""Parallel WHNC extremal scan for d_leaf <= 1 trees.

For each n, scans all trees (via geng partitions), filters d_leaf <= 1, and
computes the WHNC margin:

  margin = sum_{u in N(H)} (1/3 - P(u)) - sum_{h in H} (P(h) - 1/3),
  H = {v : P(v) > 1/3}.

Tracks minimum-margin trees and classifies whether minimizers are spiders,
including whether they match S(2^k,1)-type arm pattern.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from multiprocessing import Pool


def parse_graph6_bytes(line: bytes) -> tuple[int, list[list[int]]]:
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


def parse_graph6_str(g6: str) -> tuple[int, list[list[int]]]:
    return parse_graph6_bytes((g6.strip() + "\n").encode("ascii"))


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


def whnc_margin(probs: list[float], adj: list[list[int]], tol: float) -> float:
    one_third = 1.0 / 3.0
    heavy = {v for v, p in enumerate(probs) if p > one_third + tol}
    if not heavy:
        return 0.0

    nh = set()
    for h in heavy:
        nh.update(adj[h])
    nh -= heavy

    excess = sum(probs[h] - one_third for h in heavy)
    deficit = sum(one_third - probs[u] for u in nh)
    return deficit - excess


def push_top(
    top: list[tuple[float, str]], margin: float, g6: str, top_k: int
) -> None:
    if top_k <= 0:
        return
    if len(top) < top_k:
        top.append((margin, g6))
        top.sort(key=lambda x: x[0], reverse=True)
        return
    if margin < top[0][0]:
        top[0] = (margin, g6)
        top.sort(key=lambda x: x[0], reverse=True)


def spider_arms(n: int, adj: list[list[int]]) -> tuple[bool, list[int] | None]:
    deg = [len(adj[v]) for v in range(n)]
    hubs = [v for v in range(n) if deg[v] >= 3]
    if len(hubs) > 1:
        return False, None
    if len(hubs) == 0:
        return True, []  # path case

    hub = hubs[0]
    arms: list[int] = []
    for nb in adj[hub]:
        prev = hub
        cur = nb
        length = 1
        while True:
            d = len(adj[cur])
            if d == 1:
                arms.append(length)
                break
            if d != 2:
                return False, None
            a, b = adj[cur]
            nxt = a if b == prev else b
            prev, cur = cur, nxt
            length += 1
    arms.sort()
    return True, arms


def is_s2k1_pattern(arms: list[int] | None) -> bool:
    if arms is None or not arms:
        return False
    return arms.count(1) == 1 and all(a in (1, 2) for a in arms)


def worker(task: tuple[int, int, int, str, float, int]) -> dict:
    n, res, mod, geng, tol, top_k = task
    cmd = [geng, str(n), f"{n-1}:{n-1}", "-c", "-q", f"{res}/{mod}"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    total = 0
    dleaf = 0
    fail = 0
    min_margin = float("inf")
    min_g6 = ""
    top: list[tuple[float, str]] = []

    assert proc.stdout is not None
    for line in proc.stdout:
        g6 = line.decode("ascii").strip()
        total += 1
        nn, adj = parse_graph6_bytes(line)
        if not is_dleaf_le_1(nn, adj):
            continue
        dleaf += 1
        probs = hard_core_probs(nn, adj)
        margin = whnc_margin(probs, adj, tol)
        if margin < min_margin:
            min_margin = margin
            min_g6 = g6
        if margin < -tol:
            fail += 1
        push_top(top, margin, g6, top_k)

    proc.wait()
    return {
        "total": total,
        "dleaf": dleaf,
        "fail": fail,
        "min_margin": min_margin,
        "min_g6": min_g6,
        "top": top,
    }


def scan_n(
    n: int, workers: int, geng: str, tol: float, top_k: int
) -> dict:
    tasks = [(n, r, workers, geng, tol, top_k) for r in range(workers)]
    with Pool(workers) as pool:
        results = pool.map(worker, tasks)

    total = sum(r["total"] for r in results)
    dleaf = sum(r["dleaf"] for r in results)
    fail = sum(r["fail"] for r in results)

    min_margin = float("inf")
    min_g6 = ""
    top: list[tuple[float, str]] = []
    for r in results:
        if r["min_margin"] < min_margin:
            min_margin = r["min_margin"]
            min_g6 = r["min_g6"]
        for margin, g6 in r["top"]:
            push_top(top, margin, g6, top_k)

    top_sorted = sorted(top, key=lambda x: x[0])

    enriched = []
    for margin, g6 in top_sorted:
        nn, adj = parse_graph6_str(g6)
        spider, arms = spider_arms(nn, adj)
        enriched.append(
            {
                "margin": margin,
                "g6": g6,
                "is_spider": spider,
                "arms": arms,
                "is_s2k1_like": is_s2k1_pattern(arms),
            }
        )

    return {
        "n": n,
        "total_trees": total,
        "dleaf_trees": dleaf,
        "whnc_failures": fail,
        "min_margin": min_margin,
        "min_g6": min_g6,
        "top_minimizers": enriched,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parallel WHNC margin extremal scan."
    )
    parser.add_argument("--min-n", type=int, default=23)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument(
        "--out", default="results/whnc_extremal_scan.json"
    )
    args = parser.parse_args()

    all_results = {"config": vars(args), "per_n": {}}
    t_all = time.time()

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        out = scan_n(
            n=n,
            workers=args.workers,
            geng=args.geng,
            tol=args.tol,
            top_k=args.top_k,
        )
        all_results["per_n"][str(n)] = out

        top0 = out["top_minimizers"][0] if out["top_minimizers"] else None
        if top0 is None:
            top_str = "none"
        else:
            top_str = (
                f"margin={top0['margin']:.12f}, spider={top0['is_spider']}, "
                f"S(2^k,1)-like={top0['is_s2k1_like']}, arms={top0['arms']}"
            )
        print(
            f"n={n}: total={out['total_trees']:,} d_leaf<=1={out['dleaf_trees']:,} "
            f"fails={out['whnc_failures']} min={out['min_margin']:.12f} "
            f"({time.time() - t0:.1f}s)  {top_str}",
            flush=True,
        )

    all_results["wall_s"] = time.time() - t_all
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2)
    print(f"Saved {args.out} (wall {all_results['wall_s']:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
