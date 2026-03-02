#!/usr/bin/env python3
"""Verify the Id1 identity and 3-term integer form (IF) from GPT 5.2 Pro.

Id1: U_k V_k - U_{k-1} V_{k+1} = (V_{k+1}/V_k) d_{k-1}(U,V) + (U_k/V_k) c_k(V)
where d_{k-1}(U,V) = U_k V_{k-1} - U_{k-1} V_k
      c_k(V)       = V_k^2 - V_{k-1} V_{k+1}

Integer form (IF): J_k * Delta_k(A,J) + J_{k+1} * d_{k-1}(B,J) + B_k * c_k(J) >= 0

Also verify: g_k = c_k(E) + d_k(I,E) = d_k(I+xE, E) (leaf-augmentation)
"""
import sys
import subprocess
from collections import defaultdict
from indpoly import _polymul, _polyadd


def coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0


def delta_k(F, G, k):
    """LR minor: F_{k+1} G_k - F_k G_{k+1}"""
    return coeff(F, k+1) * coeff(G, k) - coeff(F, k) * coeff(G, k+1)


def d_km1(U, V, k):
    """Standard LR minor at k-1: U_k V_{k-1} - U_{k-1} V_k"""
    return coeff(U, k) * coeff(V, k-1) - coeff(U, k-1) * coeff(V, k)


def c_k(V, k):
    """LC gap: V_k^2 - V_{k-1} V_{k+1}"""
    return coeff(V, k)**2 - coeff(V, k-1) * coeff(V, k+1)


def xshift(p):
    return [0] + list(p)


def parse_g6(g6):
    s = g6.strip()
    n = ord(s[0]) - 63
    adj = [[] for _ in range(n)]
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        for sh in range(5, -1, -1):
            bits.append((v >> sh) & 1)
    k = 0
    for j in range(n):
        for i in range(j):
            if k < len(bits) and bits[k]:
                adj[i].append(j)
                adj[j].append(i)
            k += 1
    return n, adj


def dp_rooted(n, adj, root):
    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[root] = True
    queue = [root]
    head = 0
    while head < len(queue):
        v = queue[head]; head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                children[v].append(u)
                queue.append(u)
    order = []
    stack = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))
    dp0, dp1s = {}, {}
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            pS, pE = [1], [1]
            for c in children[v]:
                sc = _polyadd(dp0[c], xshift(dp1s[c]))
                pS = _polymul(pS, sc)
                pE = _polymul(pE, dp0[c])
            dp0[v] = pS
            dp1s[v] = pE
    return children[root], dp0, dp1s


def main():
    max_n = 18

    # Counters
    id1_checks = 0
    id1_fails = 0
    IF_checks = 0
    IF_fails = 0  # J_k * Delta_k(A,J) + J_{k+1} * d_{k-1}(B,J) + B_k * c_k(J) < 0
    W_checks = 0
    W_fails = 0  # weaker: J_k * Delta_k(A,J) + J_{k+1} * d_{k-1}(B,J) < 0
    d_neg_count = 0  # d_{k-1}(B,J) < 0
    total_stage_checks = 0

    # Profile: when d_{k-1}(B,J) < 0, how big is c_k(J) rescue?
    min_curvature_ratio = float('inf')  # min c_k(J) * B_k / |J_{k+1} * d_{k-1}(B,J)|
    min_cr_info = None

    for nn in range(3, max_n + 1):
        cmd = ["/opt/homebrew/bin/geng", "-q", str(nn),
               "%d:%d" % (nn-1, nn-1), "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0
        n_id1_fail = 0
        n_IF_fail = 0
        n_W_fail = 0

        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            n, adj = parse_g6(line)
            n_trees += 1
            degree = [len(adj[v]) for v in range(n)]

            for r in range(n):
                leaf_ch = [u for u in adj[r] if degree[u] == 1]
                if not leaf_ch:
                    continue
                children_r, dp0, dp1s = dp_rooted(n, adj, r)

                rl = []
                rnl = []
                for c in children_r:
                    nbrs = [u for u in adj[c] if u != r]
                    if len(nbrs) == 0:
                        rl.append(c)
                    else:
                        rnl.append(c)

                ell = len(rl)
                E_old = [1]
                for _ in range(ell):
                    E_old = _polymul(E_old, [1, 1])
                J_old = [1]

                for c in rnl:
                    Ec, Jc = dp0[c], dp1s[c]
                    A = _polymul(E_old, Ec)
                    B = _polymul(E_old, Jc)  # = f * r in GPT notation
                    J = _polymul(J_old, Ec)   # = g * q = J^(t)

                    md = max(len(A), len(B), len(J)) + 1

                    for k in range(md):
                        total_stage_checks += 1

                        # === Id1 verification ===
                        # LHS: B_k * J_k - B_{k-1} * J_{k+1}
                        lhs = coeff(B, k) * coeff(J, k) - coeff(B, k-1) * coeff(J, k+1)

                        # RHS (integer form, multiply through by J_k):
                        # J_{k+1} * d_{k-1}(B,J) + B_k * c_k(J)
                        # = J_{k+1} * (B_k J_{k-1} - B_{k-1} J_k) + B_k * (J_k^2 - J_{k-1} J_{k+1})
                        dkm1_BJ = d_km1(B, J, k)
                        ck_J = c_k(J, k)
                        rhs_times_Jk = coeff(J, k+1) * dkm1_BJ + coeff(B, k) * ck_J
                        lhs_times_Jk = lhs * coeff(J, k)

                        id1_checks += 1
                        if lhs_times_Jk != rhs_times_Jk:
                            id1_fails += 1
                            n_id1_fail += 1

                        # === Integer Form (IF) verification ===
                        # J_k * Delta_k(A,J) + J_{k+1} * d_{k-1}(B,J) + B_k * c_k(J) >= 0
                        Dk_AJ = delta_k(A, J, k)
                        term1 = coeff(J, k) * Dk_AJ
                        term2 = coeff(J, k+1) * dkm1_BJ
                        term3 = coeff(B, k) * ck_J
                        IF_val = term1 + term2 + term3

                        IF_checks += 1
                        if IF_val < 0:
                            IF_fails += 1
                            n_IF_fail += 1

                        # === Weaker form (W) ===
                        W_val = term1 + term2
                        W_checks += 1
                        if W_val < 0:
                            W_fails += 1
                            n_W_fail += 1

                        # Track d_{k-1}(B,J) sign
                        if dkm1_BJ < 0:
                            d_neg_count += 1
                            # curvature rescue ratio
                            if term2 != 0:
                                cr = term3 / abs(term2)
                                if cr < min_curvature_ratio:
                                    min_curvature_ratio = cr
                                    min_cr_info = (n, r, k, term1, term2, term3, cr)

                    E_old = _polyadd(A, xshift(B))
                    J_old = J

        proc.wait()
        print("n=%2d: %8d trees | id1_fail=%d IF_fail=%d W_fail=%d" % (
            nn, n_trees, n_id1_fail, n_IF_fail, n_W_fail), flush=True)

    print("\n" + "="*70)
    print("SUMMARY n=3..%d" % max_n)
    print("="*70)
    print("Id1 identity checks: %d, fails: %d" % (id1_checks, id1_fails))
    print("IF (3-term) checks:  %d, fails: %d" % (IF_checks, IF_fails))
    print("W  (2-term) checks:  %d, fails: %d" % (W_checks, W_fails))
    print("d_{k-1}(B,J) < 0:   %d / %d (%.2f%%)" % (
        d_neg_count, total_stage_checks,
        100 * d_neg_count / max(1, total_stage_checks)))
    print()
    if min_cr_info:
        info = min_cr_info
        print("Min curvature rescue ratio (term3/|term2| when d<0):")
        print("  ratio=%.6f at n=%d, r=%d, k=%d" % (info[6], info[0], info[1], info[2]))
        print("  term1=%d (Karlin), term2=%d (bad LR), term3=%d (curvature)" % (
            info[3], info[4], info[5]))
        print("  IF_total = %d" % (info[3] + info[4] + info[5]))


if __name__ == '__main__':
    main()
