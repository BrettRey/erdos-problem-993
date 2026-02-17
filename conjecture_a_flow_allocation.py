#!/usr/bin/env python3
"""Check edge-constrained deficit->excess transport feasibility for WHNC.

For each d_leaf <= 1 tree, define:
  H = {h : P(h) > 1/3}
  U = N(H) = neighbors of H outside H
  demand(h) = P(h) - 1/3
  supply(u) = 1/3 - P(u)

Question: does there exist nonnegative flow x_{h,u} on edges h-u (h in H, u in U)
such that:
  sum_{u~h} x_{h,u} >= demand(h)   for all h in H
  sum_{h~u} x_{h,u} <= supply(u)   for all u in U

This is a max-flow feasibility problem. If always feasible, it yields a local
allocation proof of WHNC stronger than global summation.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from collections import Counter
from multiprocessing import Pool


class Dinic:
    def __init__(self, n: int) -> None:
        self.n = n
        self.adj: list[list[int]] = [[] for _ in range(n)]
        self.to: list[int] = []
        self.cap: list[float] = []
        self.nxt: list[int] = []

    def add_edge(self, u: int, v: int, c: float) -> None:
        self.adj[u].append(len(self.to))
        self.to.append(v)
        self.cap.append(c)
        self.nxt.append(len(self.to))
        self.adj[v].append(len(self.to))
        self.to.append(u)
        self.cap.append(0.0)
        self.nxt.append(len(self.to) - 2)

    def max_flow(self, s: int, t: int, eps: float = 1e-12) -> float:
        flow = 0.0
        while True:
            level = [-1] * self.n
            q = [s]
            level[s] = 0
            qi = 0
            while qi < len(q):
                u = q[qi]
                qi += 1
                for ei in self.adj[u]:
                    if self.cap[ei] > eps:
                        v = self.to[ei]
                        if level[v] < 0:
                            level[v] = level[u] + 1
                            q.append(v)
            if level[t] < 0:
                break
            it = [0] * self.n

            def dfs(u: int, f: float) -> float:
                if u == t:
                    return f
                while it[u] < len(self.adj[u]):
                    ei = self.adj[u][it[u]]
                    it[u] += 1
                    v = self.to[ei]
                    if self.cap[ei] > eps and level[v] == level[u] + 1:
                        pushed = dfs(v, min(f, self.cap[ei]))
                        if pushed > eps:
                            self.cap[ei] -= pushed
                            self.cap[ei ^ 1] += pushed
                            return pushed
                return 0.0

            while True:
                pushed = dfs(s, 1e100)
                if pushed <= eps:
                    break
                flow += pushed
        return flow


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


def flow_feasible(
    adj: list[list[int]], probs: list[float], tol: float
) -> tuple[bool, float, float, int, int]:
    one_third = 1.0 / 3.0
    h_nodes = [v for v, pv in enumerate(probs) if pv > one_third + tol]
    if not h_nodes:
        return True, 0.0, 0.0, 0, 0

    h_set = set(h_nodes)
    u_set = set()
    for h in h_nodes:
        for u in adj[h]:
            if u not in h_set:
                u_set.add(u)
    u_nodes = sorted(u_set)

    demand = {h: probs[h] - one_third for h in h_nodes}
    supply = {u: one_third - probs[u] for u in u_nodes}

    total_demand = sum(demand.values())
    total_supply = sum(supply.values())

    # Flow network: source -> H -> U -> sink
    h_index = {h: i for i, h in enumerate(h_nodes)}
    u_index = {u: i for i, u in enumerate(u_nodes)}

    src = 0
    h_base = 1
    u_base = h_base + len(h_nodes)
    sink = u_base + len(u_nodes)
    g = Dinic(sink + 1)

    for h in h_nodes:
        g.add_edge(src, h_base + h_index[h], demand[h])
    for h in h_nodes:
        for u in adj[h]:
            if u in u_index:
                g.add_edge(h_base + h_index[h], u_base + u_index[u], 1e9)
    for u in u_nodes:
        cap = max(0.0, supply[u])
        g.add_edge(u_base + u_index[u], sink, cap)

    f = g.max_flow(src, sink, eps=tol)
    ok = (total_demand - f) <= max(tol, 1e-10)
    return ok, total_demand, f, len(h_nodes), len(u_nodes)


def worker(args: tuple[int, int, int, str, float]) -> dict:
    n, res, mod, geng, tol = args
    cmd = [geng, str(n), f"{n-1}:{n-1}", "-c", "-q", f"{res}/{mod}"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    total = 0
    dleaf = 0
    fail = 0
    max_gap = 0.0
    worst = None
    h_size_dist: Counter[int] = Counter()
    u_size_dist: Counter[int] = Counter()

    assert proc.stdout is not None
    for line in proc.stdout:
        n0, adj, g6 = parse_graph6(line)
        total += 1
        if not is_dleaf_le_1(n0, adj):
            continue
        dleaf += 1
        probs = hard_core_probs(n0, adj)

        ok, dem, flow, hs, us = flow_feasible(adj, probs, tol)
        h_size_dist[hs] += 1
        u_size_dist[us] += 1
        if not ok:
            fail += 1
            gap = dem - flow
            if gap > max_gap:
                max_gap = gap
                worst = {
                    "g6": g6,
                    "gap": gap,
                    "demand": dem,
                    "flow": flow,
                    "|H|": hs,
                    "|U|": us,
                }

    proc.wait()
    return {
        "total": total,
        "dleaf": dleaf,
        "fail": fail,
        "max_gap": max_gap,
        "worst": worst,
        "h_size_dist": dict(h_size_dist),
        "u_size_dist": dict(u_size_dist),
    }


def merge_counter_dicts(dicts: list[dict]) -> dict[int, int]:
    out: Counter[int] = Counter()
    for d in dicts:
        for k, v in d.items():
            out[int(k)] += int(v)
    return dict(out)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check WHNC transport feasibility via max-flow."
    )
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=23)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--geng", default="/opt/homebrew/bin/geng")
    parser.add_argument("--tol", type=float, default=1e-12)
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    total_trees = 0
    total_dleaf = 0
    total_fail = 0
    global_max_gap = 0.0
    global_worst = None
    all_h_dist: Counter[int] = Counter()
    all_u_dist: Counter[int] = Counter()
    per_n: dict[str, dict] = {}

    t_all = time.time()
    print("WHNC flow-allocation scan", flush=True)
    print("Feasibility of edge-constrained deficit->excess transport", flush=True)
    print("-" * 80, flush=True)

    for n in range(args.min_n, args.max_n + 1):
        t0 = time.time()
        tasks = [
            (n, r, args.workers, args.geng, args.tol) for r in range(args.workers)
        ]
        with Pool(args.workers) as pool:
            chunks = pool.map(worker, tasks)

        n_total = sum(c["total"] for c in chunks)
        n_dleaf = sum(c["dleaf"] for c in chunks)
        n_fail = sum(c["fail"] for c in chunks)

        total_trees += n_total
        total_dleaf += n_dleaf
        total_fail += n_fail

        n_max_gap = max(c["max_gap"] for c in chunks) if chunks else 0.0
        if n_max_gap > global_max_gap:
            global_max_gap = n_max_gap
            for c in chunks:
                if c["max_gap"] == n_max_gap and c["worst"] is not None:
                    global_worst = {"n": n, **c["worst"]}
                    break

        h_dist = merge_counter_dicts([c["h_size_dist"] for c in chunks])
        u_dist = merge_counter_dicts([c["u_size_dist"] for c in chunks])
        all_h_dist.update(h_dist)
        all_u_dist.update(u_dist)
        per_n[str(n)] = {
            "total": n_total,
            "dleaf": n_dleaf,
            "flow_fail": n_fail,
            "max_gap": n_max_gap,
            "h_size_dist": h_dist,
            "u_size_dist": u_dist,
            "elapsed_s": time.time() - t0,
        }

        print(
            f"n={n:2d}: total={n_total:9d} d_leaf<=1={n_dleaf:8d} "
            f"flow_fail={n_fail:6d} max_gap={n_max_gap:.6e} "
            f"({time.time() - t0:.1f}s)",
            flush=True,
        )

    print("-" * 80, flush=True)
    print(
        f"TOTAL trees={total_trees:,} d_leaf<=1={total_dleaf:,} "
        f"flow_fail={total_fail} wall={time.time() - t_all:.1f}s",
        flush=True,
    )
    print(f"Global max unmet demand gap: {global_max_gap:.6e}", flush=True)
    if global_worst is not None:
        print("Worst counterexample candidate:", global_worst, flush=True)

    print("Distribution of |H| among d_leaf<=1 trees:", flush=True)
    for k in sorted(all_h_dist):
        print(f"  |H|={k:2d}: {all_h_dist[k]:10d}", flush=True)
    print("Distribution of |U|=|N(H)| among d_leaf<=1 trees:", flush=True)
    for k in sorted(all_u_dist):
        print(f"  |U|={k:2d}: {all_u_dist[k]:10d}", flush=True)

    if args.out:
        payload = {
            "config": {
                "min_n": args.min_n,
                "max_n": args.max_n,
                "workers": args.workers,
                "geng": args.geng,
                "tol": args.tol,
            },
            "summary": {
                "total_trees": total_trees,
                "total_dleaf": total_dleaf,
                "flow_fail": total_fail,
                "global_max_gap": global_max_gap,
                "global_worst": global_worst,
                "wall_s": time.time() - t_all,
            },
            "dist": {
                "h_size": dict(all_h_dist),
                "u_size": dict(all_u_dist),
            },
            "per_n": per_n,
        }
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Saved {args.out}", flush=True)


if __name__ == "__main__":
    main()
