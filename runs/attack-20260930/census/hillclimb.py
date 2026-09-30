#!/usr/bin/env python3
"""Adversarial kill-test beyond the exhaustive range: local search over trees on n
vertices minimising one of
  delta   : n * min_{k in W5} delta_k            (exact ints; star value is the conjectured floor)
  newton  : min_{k in W5} delta_k (k+1)(alpha-k+1)/(alpha+1)
  Vd      : min_{k in W5} V_k delta_k           (FLOAT)
  Q4      : min over k in W5 with V_k >= 4 of Q_k = V^{3/2}(2p_k - p_{k-1} - p_{k+1})  (FLOAT)
Moves: re-hang a random leaf at a random vertex; subtree prune-and-regraft.
Simulated-annealing acceptance, restarts from random Pruefer trees and from
near-stars. Every candidate best is re-evaluated exactly at the end.

Usage: python3 hillclimb.py OBJ N SECONDS SEED
"""
import json
import math
import os
import random
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from treepoly import ipoly_from_adj, margins, q_of, tilt, tilted_Q  # noqa: E402


def prufer_tree(n, rng):
    seq = [rng.randrange(n) for _ in range(n - 2)]
    deg = [1] * n
    for x in seq:
        deg[x] += 1
    edges = []
    import heapq
    leaves = [i for i in range(n) if deg[i] == 1]
    heapq.heapify(leaves)
    for x in seq:
        leaf = heapq.heappop(leaves)
        edges.append((leaf, x))
        deg[x] -= 1
        if deg[x] == 1:
            heapq.heappush(leaves, x)
    a = heapq.heappop(leaves); b = heapq.heappop(leaves)
    edges.append((a, b))
    return edges


def adj_of(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    return adj


def objective(obj, n, edges):
    p = ipoly_from_adj(adj_of(n, edges))
    ms = margins(p, n, "n5")
    if obj == "delta":
        return float(n * min(ms.values())), p
    alpha = len(p) - 1
    if obj == "newton":
        return float(min(d * Fraction((k + 1) * (alpha - k + 1), alpha + 1) for k, d in ms.items())), p
    best = 9e9
    for k, d in ms.items():
        lam, V = tilt(p, k, iters=70)
        if obj == "Vd":
            best = min(best, V * float(d))
        elif obj == "Q4" and V >= 4:
            best = min(best, tilted_Q(p, k, lam, V))
    return best, p


def mutate(n, edges, rng):
    adj = adj_of(n, edges)
    if rng.random() < 0.6:
        leaves = [v for v in range(n) if len(adj[v]) == 1]
        v = rng.choice(leaves)
        u = adj[v][0]
        w = rng.randrange(n)
        if w == v or w == u:
            return None
        new = [e for e in edges if not (v in e and u in e)]
        new.append((v, w))
        return new
    # subtree prune and regraft: cut an edge, reattach the part not containing r
    a, b = rng.choice(edges)
    rest = [e for e in edges if e != (a, b)]
    radj = adj_of(n, rest)
    side = {b}
    st = [b]
    while st:
        x = st.pop()
        for y in radj[x]:
            if y not in side:
                side.add(y); st.append(y)
    other = [v for v in range(n) if v not in side]
    x = rng.choice(list(side)); y = rng.choice(other)
    rest.append((x, y))
    return rest


def star_edges(n):
    return [(0, i) for i in range(1, n)]


def main():
    obj, n, secs, seed = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed)
    t0 = time.time()
    best_val, best_edges = 9e9, None
    evals = 0
    restarts = 0
    while time.time() - t0 < secs:
        restarts += 1
        r = rng.random()
        if r < 0.5:
            cur = prufer_tree(n, rng)
        else:
            cur = star_edges(n)
            for _ in range(rng.randrange(1, 6)):
                m = mutate(n, cur, rng)
                if m:
                    cur = m
        cv, _ = objective(obj, n, cur); evals += 1
        if cv < best_val:
            best_val, best_edges = cv, list(cur)
        T = 0.3 * abs(cv) + 1e-3
        steps = 0
        while steps < 400 and time.time() - t0 < secs:
            steps += 1
            m = mutate(n, cur, rng)
            if m is None:
                continue
            v, _ = objective(obj, n, m); evals += 1
            if v < cv or rng.random() < math.exp(-(v - cv) / T):
                cur, cv = m, v
                if cv < best_val:
                    best_val, best_edges = cv, list(cur)
            T *= 0.985
    p = ipoly_from_adj(adj_of(n, best_edges))
    alpha = len(p) - 1
    ms = margins(p, n, "n5")
    k = min(ms, key=lambda kk: ms[kk])
    star_val = Fraction(4, n + 2) if n % 2 == 0 else Fraction(4 * n, (n + 1) ** 2)
    degs = sorted((sum(1 for e in best_edges if v in e) for v in range(n)), reverse=True)
    out = dict(objective=obj, n=n, seconds=secs, seed=seed, evals=evals, restarts=restarts,
               best_value=best_val, alpha=alpha, q=q_of(alpha), exact_min_delta=str(ms[k]), argmin_k=k,
               n_min_delta=float(n * ms[k]), star_n_min_delta=float(n * star_val),
               beats_star=bool(ms[k] < star_val), degree_sequence_top=degs[:8], edges=best_edges)
    os.makedirs(os.path.join(HERE, "data", "hillclimb"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "data", "hillclimb", f"{obj}_n{n}_s{seed}.json"), "w"), indent=0)
    print(f"{obj} n={n} seed={seed} evals={evals} best={best_val:.5f} alpha={alpha} n*minD={out['n_min_delta']:.4f} "
          f"star={out['star_n_min_delta']:.4f} beats_star={out['beats_star']} degs={degs[:6]}", flush=True)


if __name__ == "__main__":
    main()
