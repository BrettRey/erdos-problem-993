#!/usr/bin/env python3
"""
Verify STP2(I, E) diagonal: E(m+1)*I(m-1) <= E(m)*I(m) for all m.

Using I(k) = E(k) + J(k-1), this rewrites as:
    c_m(E) + correction >= 0
where
    c_m(E) = E[m]^2 - E[m-1]*E[m+1]        (LC gap of E)
    correction = E[m]*J[m-1] - E[m+1]*J[m-2]

We check every tree n=3..18, every rooting, every valid m.

All arithmetic is exact (Python integers).
"""

import subprocess
import sys
from collections import deque
from fractions import Fraction

GENG = "/opt/homebrew/bin/geng"
MAX_N = 18


# ---------- graph6 decode ----------

def graph6_to_adj(g6: str) -> list:
    """Decode a graph6 string into an adjacency list."""
    s = g6.strip()
    data = [c - 63 for c in s.encode("ascii")]
    idx = 0
    if data[0] < 63:
        n = data[0]
        idx = 1
    elif data[1] < 63:
        n = (data[1] << 12) | (data[2] << 6) | data[3]
        idx = 4
    else:
        n = (data[2] << 30) | (data[3] << 24) | (data[4] << 18) | (data[5] << 12) | (data[6] << 6) | data[7]
        idx = 8

    bits = []
    for val in data[idx:]:
        for shift in range(5, -1, -1):
            bits.append((val >> shift) & 1)

    adj = [[] for _ in range(n)]
    bit_idx = 0
    for j in range(1, n):
        for i in range(j):
            if bit_idx < len(bits) and bits[bit_idx]:
                adj[i].append(j)
                adj[j].append(i)
            bit_idx += 1
    return adj


# ---------- polynomial operations (exact integer lists) ----------

def coeff(poly, k):
    """Return coefficient at index k, 0 if out of range."""
    if 0 <= k < len(poly):
        return poly[k]
    return 0


def polymul(a, b):
    """Multiply two polynomials (convolution). Exact integers."""
    if not a or not b:
        return [0]
    la, lb = len(a), len(b)
    result = [0] * (la + lb - 1)
    for i in range(la):
        if a[i] == 0:
            continue
        for j in range(lb):
            result[i + j] += a[i] * b[j]
    return result


# ---------- tree DP ----------

def tree_dp(adj, root):
    """
    Compute (E, J) for tree rooted at `root`.
    E = dp[root][0]: exclude-root polynomial.
    J = dp[root][1] / x: include-root polynomial (shifted down by 1).
    I = E + x*J is the full independence polynomial.

    For a leaf:  E = [1], J = [1]  =>  I = 1 + x
    For internal vertex with children c1..cd:
        E(v) = product of I(ci)
        J(v) = product of E(ci)
    """
    n = len(adj)
    parent = [-1] * n
    order = []
    visited = [False] * n
    queue = deque([root])
    visited[root] = True
    while queue:
        v = queue.popleft()
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                queue.append(u)

    children = [[] for _ in range(n)]
    for v in order[1:]:
        children[parent[v]].append(v)

    E = [None] * n
    J = [None] * n

    for v in reversed(order):
        if not children[v]:
            E[v] = [1]
            J[v] = [1]
        else:
            e_prod = [1]
            j_prod = [1]
            for c in children[v]:
                # I(c) = E(c) + x*J(c)
                shifted_j = [0] + list(J[c])
                max_len = max(len(E[c]), len(shifted_j))
                i_c = [0] * max_len
                for k in range(len(E[c])):
                    i_c[k] += E[c][k]
                for k in range(len(shifted_j)):
                    i_c[k] += shifted_j[k]

                e_prod = polymul(e_prod, i_c)
                j_prod = polymul(j_prod, E[c])

            E[v] = e_prod
            J[v] = j_prod

    return E[root], J[root]


def is_support_vertex(adj, v):
    """Check if v is adjacent to at least one leaf (degree-1 vertex)."""
    for u in adj[v]:
        if len(adj[u]) == 1:
            return True
    return False


# ---------- main scan ----------

def main():
    total_checks = 0
    total_failures = 0

    correction_neg_support = 0
    correction_neg_nonsupport = 0

    min_rescue_overall = float('inf')
    min_rescue_support = float('inf')
    min_rescue_nonsupport = float('inf')

    worst_case_overall = None
    worst_case_support = None
    worst_case_nonsupport = None

    total_trees = 0
    total_rootings = 0

    for n in range(3, MAX_N + 1):
        proc = subprocess.run(
            [GENG, "-q", "-c", str(n), f"{n-1}:{n-1}"],
            capture_output=True, text=True
        )
        trees_g6 = proc.stdout.strip().split("\n")
        if trees_g6 == [""]:
            trees_g6 = []

        n_trees = len(trees_g6)
        total_trees += n_trees
        n_checks = 0
        n_fails = 0
        n_corr_neg = 0

        for g6 in trees_g6:
            adj = graph6_to_adj(g6)
            nv = len(adj)

            for root in range(nv):
                total_rootings += 1
                E_poly, J_poly = tree_dp(adj, root)
                support = is_support_vertex(adj, root)

                deg_E = len(E_poly) - 1

                for m in range(1, deg_E + 1):
                    em = coeff(E_poly, m)
                    em1 = coeff(E_poly, m - 1)
                    emp1 = coeff(E_poly, m + 1)

                    c_m = em * em - em1 * emp1

                    jm1 = coeff(J_poly, m - 1)
                    jm2 = coeff(J_poly, m - 2)

                    correction = em * jm1 - emp1 * jm2

                    total_val = c_m + correction

                    total_checks += 1
                    n_checks += 1

                    if total_val < 0:
                        total_failures += 1
                        n_fails += 1
                        print(f"  FAILURE: n={n}, g6={g6}, root={root}, m={m}, "
                              f"c_m={c_m}, corr={correction}, total={total_val}, "
                              f"support={support}")

                    if correction < 0:
                        if support:
                            correction_neg_support += 1
                        else:
                            correction_neg_nonsupport += 1
                        n_corr_neg += 1

                        # rescue ratio = c_m / |correction|
                        abs_corr = -correction
                        if abs_corr > 0:
                            # Use Fraction for exact comparison
                            rescue = Fraction(c_m, abs_corr)
                            rescue_f = float(rescue)
                        else:
                            rescue_f = float('inf')

                        if rescue_f < min_rescue_overall:
                            min_rescue_overall = rescue_f
                            worst_case_overall = (n, g6, root, m, c_m, correction, total_val, support)

                        if support and rescue_f < min_rescue_support:
                            min_rescue_support = rescue_f
                            worst_case_support = (n, g6, root, m, c_m, correction, total_val, support)

                        if not support and rescue_f < min_rescue_nonsupport:
                            min_rescue_nonsupport = rescue_f
                            worst_case_nonsupport = (n, g6, root, m, c_m, correction, total_val, support)

        print(f"n={n:2d}: {n_trees:>8d} trees, {n_checks:>10d} checks, "
              f"{n_fails} failures, {n_corr_neg} correction<0")

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Trees scanned:    {total_trees}")
    print(f"Rootings checked: {total_rootings}")
    print(f"Total (m, root) checks: {total_checks}")
    print(f"STP2(I,E) failures:     {total_failures}")
    print()
    print(f"Correction < 0 at SUPPORT vertices:     {correction_neg_support}")
    print(f"Correction < 0 at NON-SUPPORT vertices: {correction_neg_nonsupport}")
    print(f"Correction < 0 total:                   {correction_neg_support + correction_neg_nonsupport}")
    print()

    if min_rescue_overall < float('inf'):
        print(f"Min rescue ratio overall:     {min_rescue_overall:.6f}")
        if worst_case_overall:
            n, g6, root, m, c_m, corr, tot, sup = worst_case_overall
            print(f"  at n={n}, g6={g6}, root={root}, m={m}, c_m={c_m}, corr={corr}, "
                  f"total={tot}, support={sup}")
    else:
        print("Correction never negative (rescue ratio N/A)")

    print()
    if min_rescue_support < float('inf'):
        print(f"Min rescue ratio (support):   {min_rescue_support:.6f}")
        if worst_case_support:
            n, g6, root, m, c_m, corr, tot, sup = worst_case_support
            print(f"  at n={n}, g6={g6}, root={root}, m={m}, c_m={c_m}, corr={corr}, "
                  f"total={tot}")
    else:
        print("Correction never negative at support vertices")

    print()
    if min_rescue_nonsupport < float('inf'):
        print(f"Min rescue ratio (non-supp):  {min_rescue_nonsupport:.6f}")
        if worst_case_nonsupport:
            n, g6, root, m, c_m, corr, tot, sup = worst_case_nonsupport
            print(f"  at n={n}, g6={g6}, root={root}, m={m}, c_m={c_m}, corr={corr}, "
                  f"total={tot}")
    else:
        print("Correction never negative at non-support vertices")

    print()
    if total_failures == 0:
        print(f"RESULT: STP2(I,E) holds at ALL rootings of ALL trees n=3..{MAX_N}.")
    else:
        print(f"RESULT: STP2(I,E) FAILS in {total_failures} cases!")


if __name__ == "__main__":
    main()
