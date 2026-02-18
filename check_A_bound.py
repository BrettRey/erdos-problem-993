#!/usr/bin/env python3
"""Check the bound A <= (1+x)I with the CORRECT A formula.

A = P_u*P_v + x*R_u*R_v
(1+x)I = I + xI where I = P_u*R_v + R_u*P_v + R_u*R_v

B = (1+x)I - A >= 0 coefficientwise?

Algebraically:
B = (P_u*R_v + R_u*P_v + R_u*R_v) + x(P_u*R_v + R_u*P_v + R_u*R_v)
    - (P_u*P_v + x*R_u*R_v)
  = P_u*R_v + R_u*P_v + R_u*R_v - P_u*P_v
    + x(P_u*R_v + R_u*P_v + R_u*R_v - R_u*R_v)
  = P_u*R_v + R_u*P_v + R_u*R_v - P_u*P_v
    + x(P_u*R_v + R_u*P_v)
  = (P_u + R_u)(P_v + R_v) - P_u*P_v - P_u*P_v + R_u*R_v - R_u*R_v
    ... let me just expand properly.

I = P_u*R_v + R_u*P_v + R_u*R_v

(1+x)I = P_u*R_v + R_u*P_v + R_u*R_v
         + x*P_u*R_v + x*R_u*P_v + x*R_u*R_v

A = P_u*P_v + x*R_u*R_v

B = (1+x)I - A
  = P_u*R_v + R_u*P_v + R_u*R_v - P_u*P_v
    + x*(P_u*R_v + R_u*P_v + R_u*R_v - R_u*R_v)
  = P_u*R_v + R_u*P_v + R_u*R_v - P_u*P_v
    + x*(P_u*R_v + R_u*P_v)
  = P_u*(R_v - P_v) + R_u*P_v + R_u*R_v + x*P_u*R_v + x*R_u*P_v

Check: is R_v - P_v >= 0 coefficientwise? For a single vertex: R = [1], P = [0,1].
R - P = [1, -1]. NEGATIVE! So R_v - P_v is NOT nonneg.

So the bound A <= (1+x)I is NOT obvious algebraically. Let's check computationally.
"""
import subprocess
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


def main():
    print("Checking A <= (1+x)I with correct formula", flush=True)
    print("=" * 60, flush=True)

    total = 0
    violations = 0

    for n in range(3, 17):
        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        n_edges = 0
        n_viol = 0

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                n_edges += 1
                total += 1

                side_u, side_v = split_at_edge(nn, adj, u, v)
                P_u, R_u = rooted_is_poly(adj, side_u, u)
                P_v, R_v = rooted_is_poly(adj, side_v, v)

                A = _polyadd(_polymul(P_u, P_v), [0] + _polymul(R_u, R_v))

                # (1+x)I
                opxI = _polyadd(I_T, [0] + I_T)

                # Check A <= (1+x)I coefficientwise
                bad = False
                for k in range(len(A)):
                    ak = A[k]
                    bk = opxI[k] if k < len(opxI) else 0
                    if ak > bk:
                        bad = True
                        if n_viol < 3:
                            print(f"  VIOLATION n={nn} edge=({u},{v}) k={k}: "
                                  f"A[k]={ak} > (1+x)I[k]={bk}", flush=True)
                            print(f"    I = {I_T}", flush=True)
                            print(f"    A = {A}", flush=True)
                            print(f"    (1+x)I = {opxI}", flush=True)
                        break

                if bad:
                    n_viol += 1
                    violations += 1

        proc.wait()
        print(f"n={n:2d}: edges={n_edges:7d}  violations={n_viol}", flush=True)

    print()
    print(f"Total: {total} edges, {violations} violations")
    if violations == 0:
        print("A <= (1+x)I HOLDS for all edges (correct formula)")
    else:
        print(f"A <= (1+x)I FAILS in {violations} cases!")


if __name__ == "__main__":
    main()
