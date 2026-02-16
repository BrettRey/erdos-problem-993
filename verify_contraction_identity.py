#!/usr/bin/env python3
"""Verify the Subdivision-Contraction Identity:

    I(T') = I(T) + x · I(T/e)

where T' is T with edge e subdivided, and T/e is T with edge e contracted.

This means A = I(T') - I(T) = x · I(T/e), so the "extra" polynomial from
subdivision is just a shift of the IS polynomial of the contracted tree.

This reduces the delta bound to: mode(I(T/e)) >= mode(I(T)) - 1,
i.e., edge contraction drops the IS mode by at most 1.
"""
import subprocess
import sys
import time
from collections import deque

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


def contract_edge(n, adj, u, v):
    """Contract edge (u,v): merge v into u, return new adj list on n-1 vertices.
    Vertex v is removed; vertices > v are renumbered down by 1."""
    # Map old vertices to new
    new_id = {}
    for i in range(n):
        if i == v:
            new_id[i] = new_id[u]  # v maps to u
        elif i < v:
            new_id[i] = i
        else:
            new_id[i] = i - 1
    # Rebuild adjacency (n-1 vertices)
    new_n = n - 1
    new_adj = [[] for _ in range(new_n)]
    seen_edges = set()
    for i in range(n):
        for j in adj[i]:
            ni, nj = new_id[i], new_id[j]
            if ni != nj and (ni, nj) not in seen_edges:
                seen_edges.add((ni, nj))
                seen_edges.add((nj, ni))
                new_adj[ni].append(nj)
                new_adj[nj].append(ni)
    return new_n, new_adj


def split_at_edge(n, adj, u, v):
    A = set()
    queue = deque([u])
    A.add(u)
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y not in A and y != v:
                A.add(y)
                queue.append(y)
    return A, set(range(n)) - A


def rooted_is_poly(adj, vertices, root):
    vset = set(vertices)
    parent = {root: -1}
    order = [root]
    queue = deque([root])
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y in vset and y not in parent:
                parent[y] = x
                queue.append(y)
                order.append(y)
    dp_in = {}
    dp_out = {}
    for v in reversed(order):
        children = [y for y in adj[v] if y in vset and parent.get(y) == v]
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
    return dp_in[root], dp_out[root]


def subdivide_edge(n, adj, u, v):
    """Subdivide edge (u,v): insert new vertex w between u and v."""
    w = n
    new_adj = [list(nbrs) for nbrs in adj]
    new_adj.append([])
    # Remove edge u-v
    new_adj[u] = [x for x in new_adj[u] if x != v]
    new_adj[v] = [x for x in new_adj[v] if x != u]
    # Add edges u-w and w-v
    new_adj[u].append(w)
    new_adj[w].append(u)
    new_adj[w].append(v)
    new_adj[v].append(w)
    return n + 1, new_adj


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    print(f"Verifying Subdivision-Contraction Identity, n up to {max_n}", flush=True)
    print("I(T') = I(T) + x·I(T/e)", flush=True)
    print("=" * 80, flush=True)

    t0 = time.time()
    total = 0
    identity_fail = 0

    for n in range(3, max_n + 1):
        tn = time.time()
        n_edges = 0
        n_fail = 0

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                n_edges += 1
                total += 1

                # Compute I(T') by direct subdivision
                nn2, adj2 = subdivide_edge(nn, adj, u, v)
                I_Tp = independence_poly(nn2, adj2)

                # Compute I(T/e) by edge contraction
                nn3, adj3 = contract_edge(nn, adj, u, v)
                I_Te = independence_poly(nn3, adj3)

                # Check: I(T') = I(T) + x·I(T/e)
                x_I_Te = [0] + list(I_Te)  # multiply by x
                rhs = _polyadd(list(I_T), x_I_Te)

                # Compare
                if list(I_Tp) != list(rhs):
                    n_fail += 1
                    identity_fail += 1
                    if identity_fail <= 5:
                        print(f"  FAIL: n={nn}, edge=({u},{v})", flush=True)
                        print(f"    I(T') = {I_Tp[:8]}", flush=True)
                        print(f"    I(T)+xI(T/e) = {rhs[:8]}", flush=True)

        proc.wait()
        elapsed = time.time() - tn
        print(f"n={n:2d}: edges={n_edges:>10,}  identity_fail={n_fail}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 80, flush=True)
    print(f"Total: {total:,} edges in {total_time:.1f}s", flush=True)
    print(f"Identity I(T') = I(T) + x·I(T/e): "
          f"{'VERIFIED' if identity_fail == 0 else f'FAILED ({identity_fail})'}",
          flush=True)


if __name__ == "__main__":
    main()
