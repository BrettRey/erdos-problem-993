#!/usr/bin/env python3
"""Verify the formula A(x) = I(T') - I(T) against the code's computation.

The code computes A = Q_u*Q_v + x*P_u*P_v where Q_u = x*R_u.
This gives A_code = x^2*R_u*R_v + x*P_u*P_v.

But from first principles, I(T') = I(T) + P_u*P_v + x*R_u*R_v.
So A_theory = P_u*P_v + x*R_u*R_v.

Check: does A_code = A_theory?
"""
import subprocess
from collections import deque
from indpoly import independence_poly, _polymul, _polyadd

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


def subdivide_tree(n, adj, u, v):
    """Build the subdivided tree T' (add vertex w between u and v)."""
    w = n  # new vertex
    new_adj = [list(nbrs) for nbrs in adj] + [[]]
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
    print("Verifying subdivision formula")
    print("=" * 60)

    total = 0
    mismatches = 0

    for n in range(3, 13):
        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        n_edges = 0
        n_mismatch = 0

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                n_edges += 1
                total += 1

                # Direct: compute I(T') from actual subdivided tree
                n_new, adj_new = subdivide_tree(nn, adj, u, v)
                I_Tp_direct = independence_poly(n_new, adj_new)
                A_direct = [0] * max(len(I_Tp_direct), len(I_T))
                for k in range(len(I_Tp_direct)):
                    A_direct[k] += I_Tp_direct[k]
                for k in range(len(I_T)):
                    A_direct[k] -= I_T[k]
                while len(A_direct) > 1 and A_direct[-1] == 0:
                    A_direct.pop()

                # Code formula: A = x^2 R_u R_v + x P_u P_v
                side_u, side_v = split_at_edge(nn, adj, u, v)
                P_u, R_u = rooted_is_poly(adj, side_u, u)
                P_v, R_v = rooted_is_poly(adj, side_v, v)
                Q_u = [0] + R_u
                Q_v = [0] + R_v
                A_code = _polyadd(_polymul(Q_u, Q_v), [0] + _polymul(P_u, P_v))

                # Theory formula: A = P_u P_v + x R_u R_v
                A_theory = _polyadd(_polymul(P_u, P_v), [0] + _polymul(R_u, R_v))

                if A_direct != A_theory:
                    n_mismatch += 1
                    mismatches += 1
                    if n_mismatch <= 3:
                        print(f"  THEORY MISMATCH: n={nn} edge=({u},{v})")
                        print(f"    A_direct = {A_direct}")
                        print(f"    A_theory = {A_theory}")

                if A_direct == A_code:
                    # This would mean my analysis above is wrong
                    pass
                else:
                    if n_mismatch == 0 and n <= 6:
                        pass  # Don't print for every edge

        proc.wait()
        match_code = "CODE_OK" if n_mismatch == 0 else f"CODE_FAIL({n_mismatch})"
        print(f"n={n:2d}: edges={n_edges:5d}  theory_mismatches={n_mismatch}", flush=True)

    print()
    print(f"Total: {total} edges, {mismatches} theory mismatches")
    if mismatches == 0:
        print("A_theory = P_u*P_v + x*R_u*R_v is CORRECT")
    else:
        print("A_theory has mismatches!")

    # Verify A_code vs A_direct for a small example
    print()
    print("=== Detailed comparison for P4 edge 0-1 ===")
    adj_p4 = [[1], [0, 2], [1, 3], [2]]
    I_p4 = independence_poly(4, adj_p4)
    print(f"I(P4) = {I_p4}")

    n5, adj5 = subdivide_tree(4, adj_p4, 0, 1)
    I_p5 = independence_poly(n5, adj5)
    print(f"I(P5) = {I_p5}")

    side_u, side_v = split_at_edge(4, adj_p4, 0, 1)
    P_u, R_u = rooted_is_poly(adj_p4, side_u, 0)
    P_v, R_v = rooted_is_poly(adj_p4, side_v, 1)
    print(f"P_u={P_u}, R_u={R_u}, P_v={P_v}, R_v={R_v}")

    A_code = _polyadd(_polymul([0]+R_u, [0]+R_v), [0]+_polymul(P_u, P_v))
    A_theory = _polyadd(_polymul(P_u, P_v), [0]+_polymul(R_u, R_v))
    A_direct = [I_p5[k] - (I_p4[k] if k < len(I_p4) else 0) for k in range(len(I_p5))]

    print(f"A_code   = {A_code}")
    print(f"A_theory = {A_theory}")
    print(f"A_direct = {A_direct}")
    print(f"I + A_code   = {_polyadd(I_p4, A_code)}")
    print(f"I + A_theory = {_polyadd(I_p4, A_theory)}")


if __name__ == "__main__":
    main()
