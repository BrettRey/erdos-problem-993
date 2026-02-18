#!/usr/bin/env python3
"""Check mode decomposition conditions for subdivision lemma.

For each tree T and edge uv, we test whether:
  (a) mode(P_u) + mode(P_v) <= d(I(T))        [P_u P_v mode bound]
  (b) mode(R_u) + mode(R_v) <= d(I(T)) - 1    [R_u R_v mode bound]
  (1') mode(A) >= d(I(T))                      [A ascends long enough]
  (2)  A nonincreasing for k >= d+1             [A descends in tail]

If (a), (b), (1') all hold, subdivision preserves unimodality.
"""
import subprocess, sys
from indpoly import independence_poly, _polymul, _polyadd, is_unimodal


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode('ascii')]
    n = data[0]; idx = 1
    adj = [[] for _ in range(n)]
    bit = 5; word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5; idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj


def first_descent(seq):
    """Index of first strict decrease: smallest k with seq[k] > seq[k+1]."""
    for k in range(len(seq) - 1):
        if seq[k] > seq[k + 1]:
            return k
    return len(seq) - 1  # monotone non-decreasing


def mode_index(seq):
    """Index of the (first) maximum."""
    return max(range(len(seq)), key=lambda k: seq[k])


def subforest_poly(n, adj, exclude_set):
    """IS polynomial of the graph induced on V \ exclude_set."""
    keep = [v for v in range(n) if v not in exclude_set]
    if not keep:
        return [1]
    # Build subgraph
    remap = {v: i for i, v in enumerate(keep)}
    m = len(keep)
    sub_adj = [[] for _ in range(m)]
    for v in keep:
        for u in adj[v]:
            if u in remap:
                sub_adj[remap[v]].append(remap[u])
    return independence_poly(m, sub_adj)


def split_at_edge(n, adj, u, v):
    """Split tree T into components A (containing u) and B (containing v)
    after removing edge uv. Returns (A_vertices, B_vertices)."""
    # BFS from u avoiding v
    A = set()
    queue = [u]
    A.add(u)
    while queue:
        x = queue.pop()
        for y in adj[x]:
            if y not in A and y != v:
                A.add(y)
                queue.append(y)
    B = set(range(n)) - A
    return A, B


def main():
    geng = "/opt/homebrew/bin/geng"
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 14

    print("Subdivision mode decomposition check")
    print("=" * 70)
    print(f"Conditions: (a) mode(Pu)+mode(Pv) <= d")
    print(f"            (b) mode(Ru)+mode(Rv) <= d-1")
    print(f"            (1') mode(A) >= d")
    print(f"            (2) A nonincreasing for k >= d+1")
    print()

    total_edges = 0
    fail_a = 0
    fail_b = 0
    fail_1p = 0
    fail_2 = 0
    fail_uni = 0  # I(T') not unimodal (should never happen)
    tightest_a = None  # closest to failing (a)
    tightest_b = None

    for n in range(4, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        n_edges_this = 0
        n_fail_a = 0
        n_fail_b = 0
        n_fail_1p = 0
        n_fail_2 = 0

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj_list = parse_graph6(g6)

            I_T = independence_poly(nn, adj_list)
            d = first_descent(I_T)

            # Collect edges (avoid double counting)
            edges = []
            for u in range(nn):
                for v in adj_list[u]:
                    if u < v:
                        edges.append((u, v))

            for u, v in edges:
                n_edges_this += 1
                total_edges += 1

                A_verts, B_verts = split_at_edge(nn, adj_list, u, v)

                # P_u = I(A - {u}), P_v = I(B - {v})
                P_u = subforest_poly(nn, adj_list, {u})
                # But we need ONLY the A-{u} component, not all of V\{u}
                # Actually, I(A-{u}) is the IS poly of the forest A\{u}
                # with edges restricted to A. Similarly for B.
                # Let me do this properly with restricted adjacency.

                # Build A-restricted adjacency
                A_list = sorted(A_verts)
                B_list = sorted(B_verts)

                # P_u = I(A - {u}) with A's own edges
                Pu_verts = [x for x in A_list if x != u]
                Pv_verts = [x for x in B_list if x != v]

                # N[u] in A
                Nu_in_A = {u} | {x for x in adj_list[u] if x in A_verts}
                Nv_in_B = {v} | {x for x in adj_list[v] if x in B_verts}

                Ru_verts = [x for x in A_list if x not in Nu_in_A]
                Rv_verts = [x for x in B_list if x not in Nv_in_B]

                P_u = subforest_poly(nn, adj_list, set(range(nn)) - set(Pu_verts))
                P_v = subforest_poly(nn, adj_list, set(range(nn)) - set(Pv_verts))
                R_u = subforest_poly(nn, adj_list, set(range(nn)) - set(Ru_verts))
                R_v = subforest_poly(nn, adj_list, set(range(nn)) - set(Rv_verts))

                mode_Pu = mode_index(P_u)
                mode_Pv = mode_index(P_v)
                mode_Ru = mode_index(R_u)
                mode_Rv = mode_index(R_v)

                # Compute A(x) = Q_u Q_v + x P_u P_v
                # Q_u = x * R_u, Q_v = x * R_v
                Q_u = [0] + R_u  # multiply by x
                Q_v = [0] + R_v

                QQ = _polymul(Q_u, Q_v)
                PP = _polymul(P_u, P_v)
                xPP = [0] + PP  # x * P_u * P_v

                A_poly = _polyadd(QQ, xPP)
                mode_A = mode_index(A_poly)
                d_A = first_descent(A_poly)

                # Check conditions
                cond_a = (mode_Pu + mode_Pv <= d)
                cond_b = (mode_Ru + mode_Rv <= d - 1) if d >= 1 else True
                cond_1p = (d_A >= d)

                # Condition (2): A nonincreasing for k >= d+1
                cond_2 = True
                for k in range(d + 1, len(A_poly) - 1):
                    if A_poly[k + 1] > A_poly[k]:
                        cond_2 = False
                        break

                if not cond_a:
                    n_fail_a += 1
                    fail_a += 1
                if not cond_b:
                    n_fail_b += 1
                    fail_b += 1
                if not cond_1p:
                    n_fail_1p += 1
                    fail_1p += 1
                if not cond_2:
                    n_fail_2 += 1
                    fail_2 += 1

                # Check I(T') unimodality
                I_Tp = _polyadd(I_T, A_poly)
                if not is_unimodal(I_Tp):
                    fail_uni += 1

                # Track tightest (a) case
                gap_a = d - (mode_Pu + mode_Pv)
                if tightest_a is None or gap_a < tightest_a[0]:
                    tightest_a = (gap_a, g6, u, v, mode_Pu, mode_Pv, d, n)

                gap_b = (d - 1) - (mode_Ru + mode_Rv) if d >= 1 else 999
                if tightest_b is None or gap_b < tightest_b[0]:
                    tightest_b = (gap_b, g6, u, v, mode_Ru, mode_Rv, d, n)

        proc.wait()
        print(f"n={n:2d}: edges={n_edges_this:7d}  "
              f"fail(a)={n_fail_a:4d}  fail(b)={n_fail_b:4d}  "
              f"fail(1')={n_fail_1p:4d}  fail(2)={n_fail_2:4d}", flush=True)

    print()
    print(f"TOTAL: {total_edges} edge subdivisions checked")
    print(f"  Condition (a) failures: {fail_a}")
    print(f"  Condition (b) failures: {fail_b}")
    print(f"  Condition (1') failures: {fail_1p}")
    print(f"  Condition (2) failures: {fail_2}")
    print(f"  Unimodality failures: {fail_uni}")
    print()

    if tightest_a:
        gap, g6, u, v, mPu, mPv, d, nn = tightest_a
        print(f"Tightest (a): gap={gap}, n={nn}, edge=({u},{v}), "
              f"mode(Pu)={mPu}, mode(Pv)={mPv}, d={d}, tree={g6}")
    if tightest_b:
        gap, g6, u, v, mRu, mRv, d, nn = tightest_b
        print(f"Tightest (b): gap={gap}, n={nn}, edge=({u},{v}), "
              f"mode(Ru)={mRu}, mode(Rv)={mRv}, d-1={d-1}, tree={g6}")


if __name__ == "__main__":
    main()
