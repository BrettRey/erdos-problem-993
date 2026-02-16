#!/usr/bin/env python3
"""
Check tail ratio-monotonicity conditions for subdivision decomposition.

For each edge uv in a tree T:
  - split into components A (u) and B (v),
  - compute rooted pairs (P_u,Q_u), (P_v,Q_v),
  - form p=P_uP_v, q=Q_uQ_v, r=P_uQ_v+Q_uP_v, j=q+r,
  - compute d = first descent index of I(T),
  - test:
      (i) p_{d+2} <= p_d,
      (ii) j_{d+2} <= j_{d+1},
      (iii) tail ratio monotonicity of p from t>=d+1,
      (iv) tail ratio monotonicity of j from t>=d+2.

Records the first counterexample for each condition.
"""

from __future__ import annotations

import argparse
import json
from collections import deque
from typing import List, Tuple

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
    # BFS from u without crossing v
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

    # remaining vertices are comp_v
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


def ratio_tail_nonincreasing(seq: List[int], start: int) -> bool:
    # Check ratios nonincreasing for t>=start (1-based ratio index).
    last = None
    for t in range(start, len(seq)):
        if seq[t - 1] == 0 or seq[t] == 0:
            continue
        r = seq[t] / seq[t - 1]
        if last is not None and r > last + 1e-12:
            return False
        last = r
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=14)
    ap.add_argument("--backend", default="networkx", choices=["networkx", "geng", "auto"])
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    stats = {
        "max_n": args.max_n,
        "backend": args.backend,
        "checked_edges": 0,
        "b1_fail": None,
        "b2_fail": None,
        "p_ratio_fail": None,
        "j_ratio_fail": None,
    }

    for n in range(2, args.max_n + 1):
        for g6, adj in trees(n, backend=args.backend):
            # compute I(T)
            p_root, q_root = rooted_pair(adj, 0)
            I = poly_add(p_root, q_root)
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

                    # boundary checks
                    if d + 2 < len(p) and p[d + 2] > p[d] and stats["b1_fail"] is None:
                        stats["b1_fail"] = {"n": n, "edge": [u, v], "g6": g6, "d": d}
                    if d + 2 < len(j) and j[d + 2] > j[d + 1] and stats["b2_fail"] is None:
                        stats["b2_fail"] = {"n": n, "edge": [u, v], "g6": g6, "d": d}

                    # ratio monotonicity tails
                    if not ratio_tail_nonincreasing(p, max(1, d + 1)) and stats["p_ratio_fail"] is None:
                        stats["p_ratio_fail"] = {"n": n, "edge": [u, v], "g6": g6, "d": d}
                    if not ratio_tail_nonincreasing(j, max(1, d + 2)) and stats["j_ratio_fail"] is None:
                        stats["j_ratio_fail"] = {"n": n, "edge": [u, v], "g6": g6, "d": d}

            if (
                stats["b1_fail"]
                and stats["b2_fail"]
                and stats["p_ratio_fail"]
                and stats["j_ratio_fail"]
            ):
                break

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2, sort_keys=True)
    else:
        print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
