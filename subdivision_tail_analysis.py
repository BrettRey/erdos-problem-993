#!/usr/bin/env python3
"""Analyze subdivision tail behavior in detail.

Focus on:
1. Is d(I(T')) >= d(I(T)) always? (first descent doesn't move earlier)
2. When condition (2) fails (A bumps in tail), how big is the bump vs I's descent?
3. What structural features cause tight cases?
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
    for k in range(len(seq) - 1):
        if seq[k] > seq[k + 1]:
            return k
    return len(seq) - 1


def mode_index(seq):
    return max(range(len(seq)), key=lambda k: seq[k])


def subforest_poly(n, adj, exclude_set):
    keep = [v for v in range(n) if v not in exclude_set]
    if not keep:
        return [1]
    remap = {v: i for i, v in enumerate(keep)}
    m = len(keep)
    sub_adj = [[] for _ in range(m)]
    for v in keep:
        for u in adj[v]:
            if u in remap:
                sub_adj[remap[v]].append(remap[u])
    return independence_poly(m, sub_adj)


def split_at_edge(n, adj, u, v):
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
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 16

    print("Subdivision tail analysis")
    print("=" * 70)

    total_edges = 0
    descent_earlier = 0   # d(I') < d(I) ?
    descent_same = 0      # d(I') = d(I) ?
    descent_later = 0     # d(I') > d(I) ?
    tail_bump_cases = []  # cases where A has a bump in tail

    # Track combined tail: ΔI_k + ΔA_k for k >= d+1
    combined_fail = 0     # (I+A)_{k+1} > (I+A)_k for some k >= d+1

    # Track the tightest combined tail ratio
    worst_tail_ratio = 0.0
    worst_tail_info = None

    for n in range(4, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        n_edges = 0
        n_earlier = 0
        n_same = 0
        n_later = 0
        n_bump = 0
        n_combined_fail = 0

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj_list = parse_graph6(g6)

            I_T = independence_poly(nn, adj_list)
            d = first_descent(I_T)

            edges = []
            for u in range(nn):
                for v in adj_list[u]:
                    if u < v:
                        edges.append((u, v))

            for u, v in edges:
                n_edges += 1
                total_edges += 1

                A_verts, B_verts = split_at_edge(nn, adj_list, u, v)
                A_list = sorted(A_verts)
                B_list = sorted(B_verts)

                Pu_verts = [x for x in A_list if x != u]
                Pv_verts = [x for x in B_list if x != v]
                Nu_in_A = {u} | {x for x in adj_list[u] if x in A_verts}
                Nv_in_B = {v} | {x for x in adj_list[v] if x in B_verts}
                Ru_verts = [x for x in A_list if x not in Nu_in_A]
                Rv_verts = [x for x in B_list if x not in Nv_in_B]

                P_u = subforest_poly(nn, adj_list, set(range(nn)) - set(Pu_verts))
                P_v = subforest_poly(nn, adj_list, set(range(nn)) - set(Pv_verts))
                R_u = subforest_poly(nn, adj_list, set(range(nn)) - set(Ru_verts))
                R_v = subforest_poly(nn, adj_list, set(range(nn)) - set(Rv_verts))

                Q_u = [0] + R_u
                Q_v = [0] + R_v
                QQ = _polymul(Q_u, Q_v)
                PP = _polymul(P_u, P_v)
                xPP = [0] + PP
                A_poly = _polyadd(QQ, xPP)

                I_Tp = _polyadd(I_T, A_poly)
                d_Tp = first_descent(I_Tp)

                if d_Tp < d:
                    n_earlier += 1
                    descent_earlier += 1
                elif d_Tp == d:
                    n_same += 1
                    descent_same += 1
                else:
                    n_later += 1
                    descent_later += 1

                # Check A bumps in tail
                has_bump = False
                for k in range(d + 1, len(A_poly) - 1):
                    if A_poly[k + 1] > A_poly[k]:
                        has_bump = True
                        break
                if has_bump:
                    n_bump += 1

                # Check combined tail: I(T')_{k+1} <= I(T')_k for k >= d+1
                c_fail = False
                for k in range(d + 1, len(I_Tp) - 1):
                    if I_Tp[k + 1] > I_Tp[k]:
                        c_fail = True
                        n_combined_fail += 1
                        combined_fail += 1
                        break

                # Track tightest tail ratio in I(T')
                for k in range(d + 1, len(I_Tp) - 1):
                    if I_Tp[k] > 0:
                        ratio = I_Tp[k + 1] / I_Tp[k]
                        if ratio > worst_tail_ratio:
                            worst_tail_ratio = ratio
                            worst_tail_info = (g6, nn, u, v, k, ratio, d, d_Tp)

                # Detailed info for bump cases
                if has_bump and nn <= 14:
                    # Find the bump details
                    for k in range(d + 1, len(A_poly) - 1):
                        if A_poly[k + 1] > A_poly[k]:
                            bump_size = A_poly[k + 1] - A_poly[k]
                            i_drop = I_T[k] - I_T[k + 1]
                            tail_bump_cases.append({
                                'g6': g6, 'n': nn, 'u': u, 'v': v,
                                'k': k, 'd': d,
                                'bump': bump_size, 'i_drop': i_drop,
                                'A_k': A_poly[k], 'A_k1': A_poly[k + 1],
                                'I_k': I_T[k], 'I_k1': I_T[k + 1],
                                'deg_u': len(adj_list[u]),
                                'deg_v': len(adj_list[v]),
                                'deg_seq': sorted([len(adj_list[x]) for x in range(nn)], reverse=True),
                            })
                            break

        proc.wait()
        print(f"n={n:2d}: edges={n_edges:7d}  "
              f"d'<d={n_earlier:3d}  d'=d={n_same:5d}  d'>d={n_later:5d}  "
              f"A_bump={n_bump:3d}  combined_fail={n_combined_fail:3d}", flush=True)

    print()
    print(f"TOTAL: {total_edges} edges")
    print(f"  d(I') < d(I): {descent_earlier}")
    print(f"  d(I') = d(I): {descent_same}")
    print(f"  d(I') > d(I): {descent_later}")
    print(f"  A bump in tail: {len(tail_bump_cases)} (through n=14)")
    print(f"  Combined tail fail: {combined_fail}")
    print()

    if worst_tail_info:
        g6, nn, u, v, k, ratio, d, d_Tp = worst_tail_info
        print(f"Worst tail ratio in I(T'): {ratio:.6f} at k={k}")
        print(f"  tree={g6}, n={nn}, edge=({u},{v}), d(I)={d}, d(I')={d_Tp}")
    print()

    # Show bump case details
    if tail_bump_cases:
        print(f"A-bump cases (first {min(10, len(tail_bump_cases))}):")
        for case in tail_bump_cases[:10]:
            print(f"  n={case['n']} tree={case['g6']} edge=({case['u']},{case['v']}) "
                  f"d={case['d']} k={case['k']}")
            print(f"    A[k]={case['A_k']} A[k+1]={case['A_k1']} bump={case['bump']}")
            print(f"    I[k]={case['I_k']} I[k+1]={case['I_k1']} I_drop={case['i_drop']}")
            print(f"    deg(u)={case['deg_u']} deg(v)={case['deg_v']} "
                  f"deg_seq={case['deg_seq']}")


if __name__ == "__main__":
    main()
