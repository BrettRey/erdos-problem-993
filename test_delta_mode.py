#!/usr/bin/env python3
"""Test: mode(I(T)) - mode(IS(T-{v})) <= 1 for all trees T and vertices v.

This is a fundamental property that would imply mode(xRR) >= d(I(T)),
which is the key condition for the subdivision lemma proof.

The claim: removing any single vertex from a tree drops the IS polynomial
mode by at most 1. This is δ(T,v) <= 1 in the notation of the analysis.

Also test: mode(I(T)) - mode(IS(T-N[v])) for comparison (removing closed nbhd).
"""
import subprocess
import sys
import time
from collections import Counter, deque

from indpoly import _polyadd, _polymul, independence_poly

GENG = "/opt/homebrew/bin/geng"


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode("ascii")]
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


def first_descent(seq):
    for k in range(len(seq) - 1):
        if seq[k] > seq[k + 1]:
            return k
    return len(seq) - 1


def forest_is_poly(adj, vertices):
    """IS polynomial of a subgraph induced by vertices (may be disconnected)."""
    if not vertices:
        return [1]
    vset = set(vertices)
    visited = set()
    result = [1]

    for start in vset:
        if start in visited:
            continue
        # BFS to find component
        comp = []
        queue = deque([start])
        visited.add(start)
        while queue:
            x = queue.popleft()
            comp.append(x)
            for y in adj[x]:
                if y in vset and y not in visited:
                    visited.add(y)
                    queue.append(y)
        # Compute IS poly of this component (tree)
        if len(comp) == 1:
            comp_poly = [1, 1]
        else:
            # Root at first vertex and do DP
            root = comp[0]
            cset = set(comp)
            parent = {root: -1}
            order = [root]
            q = deque([root])
            while q:
                x = q.popleft()
                for y in adj[x]:
                    if y in cset and y not in parent:
                        parent[y] = x
                        q.append(y)
                        order.append(y)
            dp_in = {}
            dp_out = {}
            for v in reversed(order):
                children = [y for y in adj[v] if y in cset and parent.get(y) == v]
                if not children:
                    dp_in[v] = [0, 1]
                    dp_out[v] = [1]
                else:
                    prod_out = [1]
                    for c in children:
                        prod_out = _polymul(prod_out, dp_out[c])
                    dp_in[v] = [0] + prod_out
                    prod_both = [1]
                    for c in children:
                        prod_both = _polymul(prod_both, _polyadd(dp_in[c], dp_out[c]))
                    dp_out[v] = prod_both
            comp_poly = _polyadd(dp_in[root], dp_out[root])
        result = _polymul(result, comp_poly)
    return result


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    print(f"Testing δ(T,v) = mode(I(T)) - mode(IS(T-v)) <= 1, n up to {max_n}", flush=True)
    print("=" * 90, flush=True)

    t0 = time.time()
    total_verts = 0
    delta_dist = Counter()  # mode(I(T)) - mode(IS(T-{v}))
    delta_Nv_dist = Counter()  # mode(I(T)) - mode(IS(T-N[v]))
    max_delta = 0
    delta_ge2 = 0

    for n in range(3, max_n + 1):
        tn = time.time()
        n_verts = 0
        n_delta_ge2 = 0

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)
            d = first_descent(I_T)

            for v in range(nn):
                n_verts += 1
                total_verts += 1

                # IS(T - {v})
                remaining = set(range(nn)) - {v}
                R_v = forest_is_poly(adj, remaining)
                m_Rv = first_descent(R_v)

                delta = d - m_Rv
                delta_dist[delta] += 1
                if delta > max_delta:
                    max_delta = delta
                if delta >= 2:
                    n_delta_ge2 += 1
                    delta_ge2 += 1
                    if delta_ge2 <= 5:
                        print(f"  DELTA>=2 FOUND: n={nn}, v={v}, d={d}, m_Rv={m_Rv}, "
                              f"I_T={I_T[:6]}..., R_v={R_v[:6]}...", flush=True)

                # IS(T - N[v])
                Nv = {v} | set(adj[v])
                remaining_Nv = set(range(nn)) - Nv
                R_Nv = forest_is_poly(adj, remaining_Nv)
                m_RNv = first_descent(R_Nv)
                delta_Nv = d - m_RNv
                delta_Nv_dist[delta_Nv] += 1

        proc.wait()
        elapsed = time.time() - tn
        print(f"n={n:2d}: verts={n_verts:>10,}  δ>=2={n_delta_ge2}  ({elapsed:.1f}s)", flush=True)

    total_time = time.time() - t0
    print("=" * 90, flush=True)
    print(f"Total: {total_verts:,} (tree, vertex) pairs in {total_time:.1f}s", flush=True)
    print(flush=True)

    print(f"δ(T,v) = mode(I(T)) - mode(IS(T-v)):", flush=True)
    print(f"  Max δ: {max_delta}", flush=True)
    print(f"  δ >= 2: {'NONE' if delta_ge2 == 0 else f'{delta_ge2} ({100*delta_ge2/total_verts:.4f}%)'}",
          flush=True)
    print(flush=True)

    print("δ(T,v) distribution:", flush=True)
    for g in sorted(delta_dist.keys()):
        count = delta_dist[g]
        pct = 100 * count / total_verts
        print(f"  {g:+d}: {count:>12,} ({pct:.2f}%)", flush=True)

    print(flush=True)
    print("δ_N(T,v) = mode(I(T)) - mode(IS(T-N[v])) distribution:", flush=True)
    for g in sorted(delta_Nv_dist.keys()):
        count = delta_Nv_dist[g]
        pct = 100 * count / total_verts
        print(f"  {g:+d}: {count:>12,} ({pct:.2f}%)", flush=True)


if __name__ == "__main__":
    main()
