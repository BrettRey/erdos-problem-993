#!/usr/bin/env python3
"""
Check subdivision tail-ratio conditions on finite-core families.

For each tree in the finite-core class (b0, leaf_cap), check every edge uv:
  - compute d = first descent index of I(T),
  - compute p=P_uP_v, q=Q_uQ_v, r=P_uQ_v+Q_uP_v, j=q+r,
  - test:
      (i)  p_{d+2} <= p_d,
      (ii) j_{d+2} <= j_{d+1},
      (iii) tail log-concavity of p for t>=d+1,
      (iv) tail log-concavity of j for t>=d+2.

Reports the first counterexample for each condition.
"""

from __future__ import annotations

import argparse
import itertools
import json
from collections import deque
from typing import Iterable, List, Tuple

from indpoly import independence_poly
from trees import trees


def poly_add(a: List[int], b: List[int]) -> List[int]:
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


def poly_mul(a: List[int], b: List[int]) -> List[int]:
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return []
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def first_descent(seq: List[int]) -> int:
    for k in range(len(seq) - 1):
        if seq[k + 1] < seq[k]:
            return k
    return len(seq) - 1


def rooted_pair(adj: List[List[int]], root: int) -> Tuple[List[int], List[int]]:
    n = len(adj)
    parent = [-1] * n
    children: List[List[int]] = [[] for _ in range(n)]
    q = deque([root])
    parent[root] = root
    while q:
        v = q.popleft()
        for u in adj[v]:
            if parent[u] == -1:
                parent[u] = v
                children[v].append(u)
                q.append(u)
    # postorder
    order: List[int] = []
    stack = [(root, False)]
    while stack:
        v, processed = stack.pop()
        if processed:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    dp0: List[List[int]] = [[] for _ in range(n)]
    dp1: List[List[int]] = [[] for _ in range(n)]
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod0 = [1]
            for c in children[v]:
                prod0 = poly_mul(prod0, poly_add(dp0[c], dp1[c]))
            dp0[v] = prod0
            prod1 = [1]
            for c in children[v]:
                prod1 = poly_mul(prod1, dp0[c])
            dp1[v] = [0] + prod1
    return dp0[root], dp1[root]


def split_edge(adj: List[List[int]], u: int, v: int) -> Tuple[List[List[int]], int, List[List[int]], int]:
    n = len(adj)
    comp_u = []
    map_u = [-1] * n
    q = deque([u])
    map_u[u] = 0
    while q:
        x = q.popleft()
        comp_u.append(x)
        for y in adj[x]:
            if (x == u and y == v) or (x == v and y == u):
                continue
            if map_u[y] == -1:
                map_u[y] = len(comp_u)
                q.append(y)

    comp_v = []
    map_v = [-1] * n
    for i in range(n):
        if map_u[i] == -1:
            map_v[i] = len(comp_v)
            comp_v.append(i)

    def build(comp: List[int], mapping: List[int]) -> List[List[int]]:
        m = len(comp)
        out = [[] for _ in range(m)]
        for x in comp:
            ix = mapping[x]
            for y in adj[x]:
                if mapping[y] != -1:
                    out[ix].append(mapping[y])
        return out

    adj_u = build(comp_u, map_u)
    adj_v = build(comp_v, map_v)
    return adj_u, map_u[u], adj_v, map_v[v]


def tail_log_concave(seq: List[int], start: int) -> bool:
    # check seq[t]^2 >= seq[t-1]*seq[t+1] for all t>=start
    for t in range(start, len(seq) - 1):
        if seq[t - 1] == 0 or seq[t] == 0 or seq[t + 1] == 0:
            continue
        if seq[t] * seq[t] < seq[t - 1] * seq[t + 1]:
            return False
    return True


def encode_graph6_small(adj: List[List[int]]) -> str:
    n = len(adj)
    if n >= 63:
        raise ValueError("n too large for small graph6 encoder")
    aset = [set(nei) for nei in adj]
    bits: List[int] = []
    for j in range(1, n):
        sj = aset[j]
        for i in range(j):
            bits.append(1 if i in sj else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k : k + 6]:
            v = (v << 1) | b
        out.append(chr(v + 63))
    return "".join(out)


def build_tree_from_core(core_adj: List[List[int]], loads: Tuple[int, ...]) -> List[List[int]]:
    b = len(core_adj)
    n = b + sum(loads)
    adj: List[List[int]] = [nbrs[:] for nbrs in core_adj]
    adj.extend([] for _ in range(n - b))
    nxt = b
    for u, l in enumerate(loads):
        for _ in range(l):
            adj[u].append(nxt)
            adj[nxt].append(u)
            nxt += 1
    return adj


def candidate_adjs(b0: int, leaf_cap: int, core_backend: str) -> Iterable[List[List[int]]]:
    for b in range(1, b0 + 1):
        for _, core in trees(b, backend=core_backend):
            ranges: List[range] = []
            ok = True
            for u in range(b):
                lb = max(0, 3 - len(core[u]))
                if lb > leaf_cap:
                    ok = False
                    break
                ranges.append(range(lb, leaf_cap + 1))
            if not ok:
                continue
            for loads in itertools.product(*ranges):
                yield build_tree_from_core(core, loads)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--b0", type=int, required=True)
    ap.add_argument("--leaf-cap", type=int, required=True)
    ap.add_argument("--core-backend", default="networkx", choices=["networkx", "geng", "auto"])
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    stats = {
        "b0": args.b0,
        "leaf_cap": args.leaf_cap,
        "core_backend": args.core_backend,
        "checked_edges": 0,
        "b1_fail": None,
        "b2_fail": None,
        "p_ratio_fail": None,
        "j_ratio_fail": None,
    }

    for adj in candidate_adjs(args.b0, args.leaf_cap, args.core_backend):
        n = len(adj)
        g6 = encode_graph6_small(adj)
        I = independence_poly(n, adj)
        d = first_descent(I)

        for u in range(n):
            for v in adj[u]:
                if v < u:
                    continue
                stats["checked_edges"] += 1
                adj_u, ru, adj_v, rv = split_edge(adj, u, v)
                P_u, Q_u = rooted_pair(adj_u, ru)
                P_v, Q_v = rooted_pair(adj_v, rv)

                p = poly_mul(P_u, P_v)
                q = poly_mul(Q_u, Q_v)
                r = poly_add(poly_mul(P_u, Q_v), poly_mul(Q_u, P_v))
                j = poly_add(q, r)

                if d + 2 < len(p) and p[d + 2] > p[d] and stats["b1_fail"] is None:
                    stats["b1_fail"] = {"n": n, "edge": [u, v], "g6": g6, "d": d}
                if d + 2 < len(j) and j[d + 2] > j[d + 1] and stats["b2_fail"] is None:
                    stats["b2_fail"] = {"n": n, "edge": [u, v], "g6": g6, "d": d}

                if not tail_log_concave(p, max(1, d + 1)) and stats["p_ratio_fail"] is None:
                    stats["p_ratio_fail"] = {"n": n, "edge": [u, v], "g6": g6, "d": d}
                if not tail_log_concave(j, max(1, d + 2)) and stats["j_ratio_fail"] is None:
                    stats["j_ratio_fail"] = {"n": n, "edge": [u, v], "g6": g6, "d": d}

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2, sort_keys=True)
    else:
        print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
