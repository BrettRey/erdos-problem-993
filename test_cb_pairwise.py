#!/usr/bin/env python3
"""
Test the ladder-minor CB expansion of trees at support vertices.

For each tree (n=3..15), at each support vertex, we incrementally build
E_acc and J_acc by folding in non-leaf children one at a time.  At each
step t >= 2 we run four batteries of checks on the CB (Cauchy-Binet)
pairwise expansion terms.

Metrics
-------
(a) X_k >= 0         cross-term nonnegativity
(b) S(i,j,k) >= 0    pairwise symmetric sum
(c) F(i,j) >= 0      symmetrised bracket (pure-E part)
(d) Delta_{k-i,k-j}(P,Q) >= 0  for i>j where Delta_{i,j}(A,B)<0  (STP2 check)

All arithmetic is exact (Python ints).

Important indexing note
-----------------------
The CB cross-term index set includes boundary index -1 due to the shifted
term in Delta_{i,j}(A,B).  This script uses the full boundary-inclusive range.
"""

import subprocess


# ---------- polynomial helpers (exact integer arithmetic) ----------

def poly_mul(a, b):
    """Convolve two coefficient lists (exact ints)."""
    if not a or not b:
        return [0]
    la, lb = len(a), len(b)
    c = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            c[i + j] += ai * bj
    return c


def coeff(p, i):
    """Safe coefficient access: 0 outside valid range."""
    if 0 <= i < len(p):
        return p[i]
    return 0


# ---------- tree DP --------------------------------------------------

def tree_dp(adj, root):
    """
    Root the tree at *root* and compute (I, E, J) polynomials for every
    vertex via DFS.  Returns (I, E, J, children, parent).

    Convention (coefficient lists, index = power of x):
      E(v) = prod_c I(c)
      J(v) = prod_c E(c)
      I(v) = E(v) + x * J(v)
    """
    n = len(adj)
    parent = [-1] * n
    children = [[] for _ in range(n)]
    order = []

    # BFS to set parent / children / order
    visited = [False] * n
    stack = [root]
    visited[root] = True
    while stack:
        v = stack.pop()
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                stack.append(u)

    I = [None] * n
    E = [None] * n
    J = [None] * n

    for v in reversed(order):
        if not children[v]:          # leaf
            E[v] = [1]               # 1
            J[v] = [1]               # 1
            I[v] = [1, 1]            # 1 + x
        else:
            ev = [1]
            jv = [1]
            for c in children[v]:
                ev = poly_mul(ev, I[c])
                jv = poly_mul(jv, E[c])
            E[v] = ev
            J[v] = jv
            # I = E + x*J  (shift J right by 1)
            deg = max(len(ev), len(jv) + 1)
            iv = [0] * deg
            for i in range(len(ev)):
                iv[i] += ev[i]
            for i in range(len(jv)):
                iv[i + 1] += jv[i]
            I[v] = iv

    return I, E, J, children, parent


# ---------- parse graph6 --------------------------------------------

def graph6_to_adj(s):
    """Minimal graph6 decoder (ASCII, n < 63)."""
    s = s.strip()
    idx = 0
    n = ord(s[idx]) - 63
    idx += 1
    bits = []
    for ch in s[idx:]:
        val = ord(ch) - 63
        for b in range(5, -1, -1):
            bits.append((val >> b) & 1)
    adj = [[] for _ in range(n)]
    k = 0
    for j in range(1, n):
        for i in range(j):
            if k < len(bits) and bits[k]:
                adj[i].append(j)
                adj[j].append(i)
            k += 1
    return adj


# ---------- main scan ------------------------------------------------

def main():
    # Counters
    total_Xk = 0
    fail_Xk = 0
    total_Sijk = 0
    fail_Sijk = 0
    total_Fij = 0
    fail_Fij = 0
    total_stp2 = 0
    fail_stp2 = 0

    first_fail_Xk = None
    first_fail_Sijk = None
    first_fail_Fij = None
    first_fail_stp2 = None

    trees_checked = 0
    support_rootings = 0
    steps_checked = 0

    for n in range(3, 16):
        # Generate trees: connected graphs on n vertices with exactly n-1 edges
        cmd = ["/opt/homebrew/bin/geng", "-cq", str(n), "%d:%d" % (n-1, n-1)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for g6 in lines:
            adj = graph6_to_adj(g6)
            trees_checked += 1

            # Try each vertex as root; skip non-support vertices
            for root in range(n):
                # Support vertex: adjacent to at least one leaf
                is_support = False
                for nb in adj[root]:
                    if len(adj[nb]) == 1:   # nb is a leaf
                        is_support = True
                        break

                if not is_support:
                    continue

                support_rootings += 1

                # Compute full tree DP
                I, E, J, children, parent = tree_dp(adj, root)

                # Identify leaf vs non-leaf children of root
                leaf_kids = []
                nonleaf_kids = []
                for c in children[root]:
                    if not children[c]:  # c is a leaf (no children in rooted tree)
                        leaf_kids.append(c)
                    else:
                        nonleaf_kids.append(c)

                ell = len(leaf_kids)  # number of leaf children

                # Need >= 2 non-leaf children to reach step t >= 2
                if len(nonleaf_kids) < 2:
                    continue

                # E_acc starts as (1+x)^ell
                e_acc = [1]
                one_plus_x = [1, 1]
                for _ in range(ell):
                    e_acc = poly_mul(e_acc, one_plus_x)
                j_acc = [1]

                # Fold non-leaf children one by one
                for t_idx, c in enumerate(nonleaf_kids):
                    # The factor polynomials for this child
                    P_c = I[c]   # I(c)
                    Q_c = E[c]   # E(c)

                    # Update accumulators
                    e_new = poly_mul(e_acc, P_c)
                    j_new = poly_mul(j_acc, Q_c)

                    step_number = t_idx + 1  # 1-indexed

                    if step_number >= 2:
                        steps_checked += 1

                        # A = E_acc (old), B = J_acc (old)
                        A = e_acc
                        B = j_acc

                        max_new = max(len(e_new), len(j_new))

                        # Compute Lambda_k^{old}(A,B) = A[k]*B[k] - A[k-1]*B[k+1]
                        max_ab = max(len(A), len(B))
                        Lambda_old = []
                        for k in range(max_ab + 1):
                            val = coeff(A, k) * coeff(B, k) - coeff(A, k - 1) * coeff(B, k + 1)
                            Lambda_old.append(val)

                        # Compute Lambda_k^{new} = E_new[k]*J_new[k] - E_new[k-1]*J_new[k+1]
                        Lambda_new = []
                        for k in range(max_new + 1):
                            val = coeff(e_new, k) * coeff(j_new, k) - coeff(e_new, k - 1) * coeff(j_new, k + 1)
                            Lambda_new.append(val)

                        # Compute D_k = sum_i Lambda_i^{old} * P[k-i] * Q[k-i]
                        D = [0] * (max_new + 1)
                        for k in range(max_new + 1):
                            s = 0
                            for i in range(len(Lambda_old)):
                                s += Lambda_old[i] * coeff(P_c, k - i) * coeff(Q_c, k - i)
                            D[k] = s

                        # Check (a): X_k = Lambda_k^{new} - D_k >= 0
                        for k in range(max_new + 1):
                            xk = coeff(Lambda_new, k) - D[k]
                            total_Xk += 1
                            if xk < 0:
                                fail_Xk += 1
                                if first_fail_Xk is None:
                                    first_fail_Xk = (n, root, step_number, k, xk)

                        # Pairwise checks over (i, j) with i < j.
                        # Include boundary index -1 (from the shifted Delta term).
                        min_idx = -1
                        max_idx = max(len(A), len(B))
                        for i in range(min_idx, max_idx + 1):
                            for j in range(i + 1, max_idx + 1):
                                # Delta_{i,j}(A,B) = A[i]*B[j] - A[i-1]*B[j+1]
                                delta_ij = coeff(A, i) * coeff(B, j) - coeff(A, i - 1) * coeff(B, j + 1)
                                # Delta_{j,i}(A,B) = A[j]*B[i] - A[j-1]*B[i+1]
                                delta_ji = coeff(A, j) * coeff(B, i) - coeff(A, j - 1) * coeff(B, i + 1)

                                # (c) F(i,j) = Delta_{i,j} + Delta_{j,i}
                                fij = delta_ij + delta_ji
                                total_Fij += 1
                                if fij < 0:
                                    fail_Fij += 1
                                    if first_fail_Fij is None:
                                        first_fail_Fij = (n, root, step_number, i, j, fij)

                                # (b) S(i,j,k) for each valid k
                                for k in range(max_new + 1):
                                    s_val = (delta_ij * coeff(P_c, k - i) * coeff(Q_c, k - j)
                                             + delta_ji * coeff(P_c, k - j) * coeff(Q_c, k - i))
                                    total_Sijk += 1
                                    if s_val < 0:
                                        fail_Sijk += 1
                                        if first_fail_Sijk is None:
                                            first_fail_Sijk = (n, root, step_number, i, j, k, s_val)

                                # (d) STP2 check on (P,Q) side when (A,B) delta is negative
                                # If delta_ij < 0 (i < j), check
                                #   Delta_{k-i,k-j}(P,Q) = P[k-i]*Q[k-j] - P[k-i-1]*Q[k-j+1]
                                #   Here k-i < k-j is FALSE (k-i > k-j), so this tests TP2 of (P,Q)
                                if delta_ij < 0:
                                    for k in range(max_new + 1):
                                        ki = k - i
                                        kj = k - j
                                        d_pq = (coeff(P_c, ki) * coeff(Q_c, kj)
                                                - coeff(P_c, ki - 1) * coeff(Q_c, kj + 1))
                                        total_stp2 += 1
                                        if d_pq < 0:
                                            fail_stp2 += 1
                                            if first_fail_stp2 is None:
                                                first_fail_stp2 = (n, root, step_number, i, j, k, d_pq)

                                # If delta_ji < 0 (j > i), check
                                #   Delta_{k-j,k-i}(P,Q) = P[k-j]*Q[k-i] - P[k-j-1]*Q[k-i+1]
                                if delta_ji < 0:
                                    for k in range(max_new + 1):
                                        kj2 = k - j
                                        ki2 = k - i
                                        d_pq = (coeff(P_c, kj2) * coeff(Q_c, ki2)
                                                - coeff(P_c, kj2 - 1) * coeff(Q_c, ki2 + 1))
                                        total_stp2 += 1
                                        if d_pq < 0:
                                            fail_stp2 += 1
                                            if first_fail_stp2 is None:
                                                first_fail_stp2 = (n, root, step_number, j, i, k, d_pq)

                    # Update accumulators for next step
                    e_acc = e_new
                    j_acc = j_new

        print("n=%2d done  (trees so far: %d, support rootings: %d, steps: %d)" %
              (n, trees_checked, support_rootings, steps_checked),
              flush=True)

    # ---------- Report ----------
    print()
    print("=" * 70)
    print("CB PAIRWISE EXPANSION SCAN  n = 3 .. 15")
    print("=" * 70)
    print("Trees checked:       %12s" % format(trees_checked, ","))
    print("Support rootings:    %12s" % format(support_rootings, ","))
    print("Steps (t >= 2):      %12s" % format(steps_checked, ","))
    print()
    print("(a) X_k >= 0 checks: %12s   failures: %d" % (format(total_Xk, ","), fail_Xk))
    if first_fail_Xk:
        print("    First failure: n=%d, root=%d, step=%d, k=%d, val=%d" % first_fail_Xk)
    print()
    print("(b) S(i,j,k) >= 0:  %12s   failures: %d" % (format(total_Sijk, ","), fail_Sijk))
    if first_fail_Sijk:
        print("    First failure: n=%d, root=%d, step=%d, i=%d, j=%d, k=%d, val=%d" % first_fail_Sijk)
    print()
    print("(c) F(i,j) >= 0:    %12s   failures: %d" % (format(total_Fij, ","), fail_Fij))
    if first_fail_Fij:
        print("    First failure: n=%d, root=%d, step=%d, i=%d, j=%d, val=%d" % first_fail_Fij)
    else:
        print("    ** F(i,j) >= 0 holds universally: the symmetrised bracket is always nonneg **")
    print()
    print("(d) STP2(P,Q) >= 0: %12s   failures: %d" % (format(total_stp2, ","), fail_stp2))
    if first_fail_stp2:
        print("    First failure: n=%d, root=%d, step=%d, i=%d, j=%d, k=%d, val=%d" % first_fail_stp2)
    print()


if __name__ == "__main__":
    main()
