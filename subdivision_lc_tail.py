#!/usr/bin/env python3
"""Test whether LC of I(T) + the bound A <= (1+x)I implies the combined tail.

Key question: if I is log-concave and unimodal with first descent d,
and A is nonneg with A <= (1+x)I coefficientwise and mode(A) >= d,
does I + A have to be unimodal?

We construct ADVERSARIAL examples:
- Take a real LC unimodal I
- Construct A that maximizes bumps in the tail while staying <= (1+x)I
- Check if I + A is unimodal

Also test: for each subdivision, how does a_k compare to i_k in the tail?
If a_k / i_k is decreasing in the tail, that combined with I's descent
would automatically give the combined tail condition.
"""
import subprocess, sys
from indpoly import independence_poly, _polymul, _polyadd, is_unimodal, is_log_concave


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


def main():
    geng = "/opt/homebrew/bin/geng"
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 16

    print("Tail ratio analysis: a_k / i_k behavior in descent")
    print("=" * 70)

    total_edges = 0
    ratio_increasing = 0   # Cases where a_k/i_k increases in tail
    worst_ratio_increase = 0.0
    worst_info = None

    # Track: is a_k <= i_k always in the tail?
    a_exceeds_i = 0

    # Track: the key inequality
    # In descent region, is a_{k+1}/i_{k+1} <= a_k/i_k?
    # This would give: a_{k+1}/a_k <= i_{k+1}/i_k,
    # hence a_k decreases at least as fast as i_k, giving combined tail.
    ratio_mono_fail = 0

    # NEW: check if (1+x)I - A has LC-like properties that help
    # Specifically: is B = (1+x)I - A log-concave?
    b_lc_fail = 0
    # Is B unimodal?
    b_uni_fail = 0

    for n in range(5, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        n_edges = 0
        n_ratio_fail = 0
        n_a_exceeds = 0
        n_b_lc = 0
        n_b_uni = 0

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

                # (1+x)I
                I_shifted = [0] + I_T  # x * I(T)
                one_plus_x_I = _polyadd(I_T, I_shifted)

                # B = (1+x)I - A
                B_poly = [0] * max(len(one_plus_x_I), len(A_poly))
                for k in range(len(one_plus_x_I)):
                    B_poly[k] += one_plus_x_I[k]
                for k in range(len(A_poly)):
                    B_poly[k] -= A_poly[k]
                # Trim trailing zeros
                while len(B_poly) > 1 and B_poly[-1] == 0:
                    B_poly.pop()

                # Check B properties
                if not is_log_concave(B_poly):
                    n_b_lc += 1
                    b_lc_fail += 1
                if not is_unimodal(B_poly):
                    n_b_uni += 1
                    b_uni_fail += 1

                # Check a_k / i_k in tail
                for k in range(d + 1, min(len(A_poly), len(I_T)) - 1):
                    if I_T[k] > 0 and I_T[k + 1] > 0:
                        r_k = A_poly[k] / I_T[k] if k < len(A_poly) else 0
                        r_k1 = A_poly[k + 1] / I_T[k + 1] if k + 1 < len(A_poly) else 0
                        if r_k1 > r_k + 1e-12:
                            n_ratio_fail += 1
                            ratio_mono_fail += 1
                            increase = r_k1 - r_k
                            if increase > worst_ratio_increase:
                                worst_ratio_increase = increase
                                worst_info = (g6, nn, u, v, k, r_k, r_k1, d)
                            break  # one failure per edge

                # Check a_k > i_k in tail
                for k in range(d + 1, len(A_poly)):
                    ak = A_poly[k] if k < len(A_poly) else 0
                    ik = I_T[k] if k < len(I_T) else 0
                    if ak > ik:
                        n_a_exceeds += 1
                        a_exceeds_i += 1
                        break

        proc.wait()
        print(f"n={n:2d}: edges={n_edges:7d}  "
              f"ratio_fail={n_ratio_fail:4d}  a>i={n_a_exceeds:4d}  "
              f"B_LC_fail={n_b_lc:5d}  B_uni_fail={n_b_uni:4d}", flush=True)

    print()
    print(f"TOTAL: {total_edges} edges")
    print(f"  a_k/i_k ratio increasing in tail: {ratio_mono_fail}")
    print(f"  a_k > i_k somewhere in tail: {a_exceeds_i}")
    print(f"  B = (1+x)I - A not LC: {b_lc_fail}")
    print(f"  B = (1+x)I - A not unimodal: {b_uni_fail}")
    print()

    if worst_info:
        g6, nn, u, v, k, r_k, r_k1, d = worst_info
        print(f"Worst ratio increase: {worst_ratio_increase:.6f}")
        print(f"  tree={g6}, n={nn}, edge=({u},{v}), k={k}")
        print(f"  a[k]/i[k]={r_k:.6f}, a[k+1]/i[k+1]={r_k1:.6f}, d={d}")


if __name__ == "__main__":
    main()
