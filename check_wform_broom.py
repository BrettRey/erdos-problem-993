#!/usr/bin/env python3
"""W-form check on T_{3,4} broom rooted at support vertex x_0_0.

Tree: T_{3,4} — root v with 3 children w_i, each w_i has 4 children x_{i,j},
each x_{i,j} has 1 leaf y_{i,j}. n = 1 + 3 + 12 + 12 = 28.

Root at x_0_0 (vertex 4). This is a support vertex with:
  - 1 leaf child: y_0_0
  - 1 non-leaf child: the subtree through w_0

Since s = 1 (one non-leaf child), the W-form was PROVED for this case.
This script verifies that computationally and prints all the intermediate
quantities.

W form at step 1 (incorporating child c):
  A = E_old * E_c
  B = E_old * J_c
  C = J_old * E_c
  E_new = A + x*B
  J_new = C

  W_k = J_k * Delta_k(A,C) + J_{k+1} * d_{k-1}(B,C) >= 0

Also checks:
  - E_new ≽ J_new (ratio dominance: e_{k+1}*j_k >= e_k*j_{k+1})
  - SCC: Delta_k((1+x)*I_new, E_new) >= 0
"""
import sys
sys.path.insert(0, "/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993")
from indpoly import _polymul, _polyadd


def xshift(p):
    """Multiply polynomial by x (prepend 0)."""
    return [0] + list(p)


def coeff(p, k):
    """Safely get coefficient at index k."""
    return p[k] if 0 <= k < len(p) else 0


def delta_k(F, G, k):
    """Compute F_{k+1}*G_k - F_k*G_{k+1}."""
    return coeff(F, k+1) * coeff(G, k) - coeff(F, k) * coeff(G, k+1)


def poly_str(p, name="P"):
    """Pretty-print a polynomial."""
    return f"{name} = {p}"


def build_broom(m, t):
    """Build T_{m,t} broom tree.

    Vertex numbering:
      0: root v
      1..m: w_0, ..., w_{m-1}
      m+1..m+m*t: x vertices (x_{i,j} = 1 + m + i*t + j)
      1+m+m*t..1+m+2*m*t-1: y vertices (y_{i,j} = 1 + m + m*t + i*t + j)
    """
    n = 1 + m + m*t + m*t
    adj = [[] for _ in range(n)]
    labels = {}

    v = 0
    labels[v] = 'v'

    next_id = 1
    w_ids = []
    for i in range(m):
        w = next_id; next_id += 1
        w_ids.append(w)
        labels[w] = f'w_{i}'
        adj[v].append(w)
        adj[w].append(v)

    x_ids = []
    for i in range(m):
        for j in range(t):
            x = next_id; next_id += 1
            x_ids.append(x)
            labels[x] = f'x_{i}_{j}'
            adj[w_ids[i]].append(x)
            adj[x].append(w_ids[i])

    y_ids = []
    for i in range(m):
        for j in range(t):
            y = next_id; next_id += 1
            y_ids.append(y)
            labels[y] = f'y_{i}_{j}'
            x_idx = i * t + j
            adj[x_ids[x_idx]].append(y)
            adj[y].append(x_ids[x_idx])

    assert next_id == n
    return n, adj, labels, x_ids, w_ids, y_ids


def dp_rooted_full(n, adj, root):
    """Rooted tree DP returning per-vertex (dp0, dp1_stripped).

    dp0[v] = E_v (IS poly of subtree(v) excluding v)
    dp1_stripped[v] = J_v (IS poly of subtree(v) including v, divided by x)

    So I(subtree(v)) = E_v + x * J_v.
    """
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

    dp0 = [None] * n
    dp1s = [None] * n

    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            prod_S = [1]
            prod_E = [1]
            for c in children[v]:
                s_c = _polyadd(dp0[c], xshift(dp1s[c]))
                prod_S = _polymul(prod_S, s_c)
                prod_E = _polymul(prod_E, dp0[c])
            dp0[v] = prod_S
            dp1s[v] = prod_E

    return dp0, dp1s, children


def get_mode(poly):
    """Return leftmost mode index."""
    max_val = max(poly)
    for k, v in enumerate(poly):
        if v == max_val:
            return k
    return 0


def main():
    m, t = 3, 4
    n, adj, labels, x_ids, w_ids, y_ids = build_broom(m, t)
    print(f"Tree T_{{{m},{t}}}: n = {n}")
    print(f"Vertex degrees: v(0)={len(adj[0])}, w_i={len(adj[1])}, "
          f"x_i_j={len(adj[4])}, y_i_j={len(adj[16])}")

    # Root at x_0_0 = vertex 4
    root = x_ids[0]  # = 4
    print(f"\nRooting at vertex {root} (label: {labels[root]})")
    print(f"  degree = {len(adj[root])}")
    print(f"  neighbors: {[(u, labels[u]) for u in adj[root]]}")

    dp0, dp1s, children = dp_rooted_full(n, adj, root)

    # Children of the root when rooted at x_0_0
    root_children = children[root]
    print(f"  children of root: {[(c, labels[c]) for c in root_children]}")

    leaf_children = [c for c in root_children if not children[c]]
    nonleaf_children = [c for c in root_children if children[c]]
    ell = len(leaf_children)
    s = len(nonleaf_children)
    print(f"  leaf children (ell={ell}): {[(c, labels[c]) for c in leaf_children]}")
    print(f"  non-leaf children (s={s}): {[(c, labels[c]) for c in nonleaf_children]}")

    # ---- Step 0: Initialize from leaf children ----
    # E_old = (1+x)^ell (product of S_c = (1+x) for each leaf child)
    # J_old = [1] (product of E_c = [1] for each leaf child)
    #
    # Actually: E for the root so far = product of (E_c + x*J_c) over processed children
    # J for the root so far = product of E_c over processed children
    # Leaf child: E_c = [1], J_c = [1], S_c = [1,1]

    E_old = [1]
    J_old = [1]
    for lc in leaf_children:
        Ec = dp0[lc]   # = [1]
        Jc = dp1s[lc]  # = [1]
        Sc = _polyadd(Ec, xshift(Jc))  # = [1, 1]
        E_old = _polymul(E_old, Sc)
        J_old = _polymul(J_old, Ec)

    print(f"\n{'='*70}")
    print(f"After processing {ell} leaf child(ren):")
    print(f"  E_old = {E_old}")
    print(f"  J_old = {J_old}")

    # ---- Step 1: Incorporate non-leaf child c ----
    assert s == 1, f"Expected s=1, got s={s}"
    c = nonleaf_children[0]
    Ec = dp0[c]
    Jc = dp1s[c]
    Sc = _polyadd(Ec, xshift(Jc))
    Ic = _polyadd(Ec, xshift(Jc))

    print(f"\n{'='*70}")
    print(f"Non-leaf child c = {c} ({labels[c]})")
    print(f"  E_c (len={len(Ec)}): {Ec}")
    print(f"  J_c (len={len(Jc)}): {Jc}")
    print(f"  I_c = E_c + x*J_c (len={len(Ic)}): {Ic}")

    mode_Ic = get_mode(Ic)
    print(f"  mode(I_c) = {mode_Ic}")

    # Product build
    A = _polymul(E_old, Ec)
    B = _polymul(E_old, Jc)
    C = _polymul(J_old, Ec)

    E_new = _polyadd(A, xshift(B))
    J_new = list(C)

    # Final root poly
    I_root = _polyadd(E_new, xshift(J_new))
    mode_I = get_mode(I_root)

    print(f"\n{'='*70}")
    print(f"Product build (step 1, child c = {labels[c]}):")
    print(f"  A = E_old * E_c (len={len(A)}): {A}")
    print(f"  B = E_old * J_c (len={len(B)}): {B}")
    print(f"  C = J_old * E_c (len={len(C)}): {C}")
    print(f"  E_new = A + x*B (len={len(E_new)}): {E_new}")
    print(f"  J_new = C       (len={len(J_new)}): {J_new}")
    print(f"  I(T) = E_new + x*J_new (len={len(I_root)}): {I_root}")
    print(f"  mode(I(T)) = {mode_I}")

    # Verify E_new and J_new match the full DP
    E_root = dp0[root]
    J_root = dp1s[root]
    print(f"\nVerification against full DP:")
    print(f"  E_new == dp0[root]: {E_new == E_root}")
    print(f"  J_new == dp1s[root]: {J_new == J_root}")

    # ---- W-form check ----
    # W_k = J_k * Delta_k(A,C) + J_{k+1} * d_{k-1}(B,C)
    # where J here refers to J_new = C, and Delta_k(F,G) = F_{k+1}*G_k - F_k*G_{k+1}
    #
    # Wait: the W-form is about the STEP. The "J" in the W-form is J_old (the
    # partial product BEFORE this step), not J_new.
    #
    # Actually let me be precise. The P2 condition on (E_new, J_new) is:
    #   e_{k+1}^{new} * j_k^{new} >= e_k^{new} * j_{k+1}^{new}  for k < mode
    #
    # With E_new = A + xB, J_new = C, this becomes:
    #   (A_{k+1} + B_k) * C_k >= (A_k + B_{k-1}) * C_{k+1}
    # = A_{k+1}*C_k - A_k*C_{k+1} + B_k*C_k - B_{k-1}*C_{k+1}
    # = Delta_k(A,C) + B_k*C_k - B_{k-1}*C_{k+1}
    #
    # But A = E_old * Ec, B = E_old * Jc, C = J_old * Ec.
    # For s=1 and ell=1: E_old = [1,1], J_old = [1].
    # So A = [1,1]*Ec, B = [1,1]*Jc, C = Ec.

    print(f"\n{'='*70}")
    print("W-FORM CHECK")
    print(f"  W_k = J_new_k * Delta_k(A, C) + J_new_{{k+1}} * (B_k*C_k - B_{{k-1}}*C_{{k+1}})")
    print(f"  (But really: checking P2 on (E_new, J_new) directly)")
    print()

    # Direct P2 check on (E_new, J_new)
    print("P2 check: e_{k+1}^new * j_k^new >= e_k^new * j_{k+1}^new")
    print(f"  (for k = 0, ..., {mode_I - 1})")
    any_p2_fail = False
    for k in range(mode_I):
        lhs = coeff(E_new, k+1) * coeff(J_new, k)
        rhs = coeff(E_new, k) * coeff(J_new, k+1)
        diff = lhs - rhs
        marker = " <-- NEGATIVE (P2 FAILS)!" if diff < 0 else ""
        if diff < 0:
            any_p2_fail = True
        print(f"  k={k:2d}: e_{k+1}*j_{k} - e_{k}*j_{{k+1}} = {diff}{marker}")

    print(f"\nP2 overall: {'FAIL' if any_p2_fail else 'PASS'}")

    # P3 check: e_k >= j_{k-1} for k >= mode
    print(f"\nP3 check: e_k^new >= j_{{k-1}}^new for k >= {mode_I}")
    max_k = max(len(E_new), len(J_new) + 1)
    any_p3_fail = False
    for k in range(mode_I, max_k):
        ek = coeff(E_new, k)
        jk1 = coeff(J_new, k - 1)
        diff = ek - jk1
        marker = " <-- NEGATIVE (P3 FAILS)!" if diff < 0 else ""
        if diff < 0:
            any_p3_fail = True
        print(f"  k={k:2d}: e_{k} - j_{{k-1}} = {ek} - {jk1} = {diff}{marker}")

    print(f"\nP3 overall: {'FAIL' if any_p3_fail else 'PASS'}")

    # ---- W-form expansion ----
    # The W-form for the s=1 step. When we go from (E_old, J_old) to (E_new, J_new):
    #   E_new = A + xB where A = E_old*Ec, B = E_old*Jc
    #   J_new = C = J_old*Ec
    #
    # P2 at index k requires:
    #   (A_{k+1} + B_k)*C_k >= (A_k + B_{k-1})*C_{k+1}
    #   = Delta_k(A,C) + (B_k*C_k - B_{k-1}*C_{k+1}) >= 0
    #
    # Let's call:
    #   Term1_k = Delta_k(A,C) = A_{k+1}*C_k - A_k*C_{k+1}
    #   Term2_k = B_k*C_k - B_{k-1}*C_{k+1}
    #   W_k = Term1_k + Term2_k >= 0

    print(f"\n{'='*70}")
    print("W-FORM DECOMPOSITION (P2 = Term1 + Term2 >= 0)")
    print(f"  Term1_k = Delta_k(A,C) = A_{{k+1}}*C_k - A_k*C_{{k+1}}")
    print(f"  Term2_k = B_k*C_k - B_{{k-1}}*C_{{k+1}}")
    print()

    any_wfail = False
    for k in range(mode_I):
        t1 = delta_k(A, C, k)
        t2 = coeff(B, k) * coeff(C, k) - coeff(B, k-1) * coeff(C, k+1)
        wk = t1 + t2
        marker = " <-- NEGATIVE!" if wk < 0 else ""
        if wk < 0:
            any_wfail = True
        print(f"  k={k:2d}: Term1={t1:>15,d}  Term2={t2:>15,d}  W_k={wk:>15,d}{marker}")

    print(f"\nW-form overall: {'FAIL' if any_wfail else 'PASS'}")

    # ---- E_new ≽ J_new (full ratio dominance) ----
    print(f"\n{'='*70}")
    print("E_new ≽ J_new (ratio dominance: e_{k+1}*j_k >= e_k*j_{k+1} for ALL k)")
    any_rd_fail = False
    max_k_rd = max(len(E_new), len(J_new))
    for k in range(max_k_rd):
        lhs = coeff(E_new, k+1) * coeff(J_new, k)
        rhs = coeff(E_new, k) * coeff(J_new, k+1)
        diff = lhs - rhs
        marker = " <-- NEGATIVE (ratio dominance FAILS)!" if diff < 0 else ""
        if diff < 0:
            any_rd_fail = True
        print(f"  k={k:2d}: {diff:>15,d}{marker}")

    print(f"\nRatio dominance overall: {'FAIL' if any_rd_fail else 'PASS'}")

    # ---- SCC: Delta_k((1+x)*I, E) >= 0 ----
    print(f"\n{'='*70}")
    print("SCC check: Delta_k((1+x)*I, E) >= 0 for all k")
    I_tilde = _polyadd(I_root, xshift(E_new))  # (1+x)*I... wait
    # Actually (1+x)*I != I + x*E. Let me think.
    # I = E + x*J
    # (1+x)*I = I + x*I = (E + xJ) + x(E + xJ) = E + xJ + xE + x^2*J
    #         = E + x(E+J) + x^2*J
    # That's not the same as I + xE = E + xJ + xE = E + x(E+J)
    #
    # The "I_tilde" from verify_gpt_counterexample.py is I + xE, which is
    # NOT (1+x)*I. Let me check what was meant.
    #
    # From the scan: the P* check is on (E, J) via P2 and P3.
    # SCC (strong coefficient condition) from the paper is:
    #   a_k^2 >= a_{k-1}*a_{k+1} for all k (log-concavity of I)
    # That's a different thing. Let me just check LC of the full I(T).

    print("\nLog-concavity check on I(T):")
    any_lc_fail = False
    for k in range(1, len(I_root) - 1):
        lhs = I_root[k] * I_root[k]
        rhs = I_root[k-1] * I_root[k+1]
        diff = lhs - rhs
        marker = " <-- LC FAILS!" if diff < 0 else ""
        if diff < 0:
            any_lc_fail = True
        print(f"  k={k:2d}: i_k^2 - i_{{k-1}}*i_{{k+1}} = {diff:>15,d}{marker}")

    print(f"\nLC overall: {'FAIL' if any_lc_fail else 'PASS'}")

    # ---- Child-level checks: P2, P3 on (Ec, Jc) ----
    print(f"\n{'='*70}")
    print(f"Child subtree check for c = {c} ({labels[c]})")
    mode_c = get_mode(Ic)
    print(f"  mode(I_c) = {mode_c}")

    print(f"\n  P2 on (E_c, J_c):")
    any_c_p2_fail = False
    for k in range(mode_c):
        lhs = coeff(Ec, k+1) * coeff(Jc, k)
        rhs = coeff(Ec, k) * coeff(Jc, k+1)
        diff = lhs - rhs
        marker = " <-- FAIL" if diff < 0 else ""
        if diff < 0:
            any_c_p2_fail = True
        print(f"    k={k:2d}: {diff:>12,d}{marker}")

    print(f"  P2 on child: {'FAIL' if any_c_p2_fail else 'PASS'}")

    print(f"\n  P3 on (E_c, J_c):")
    any_c_p3_fail = False
    max_k_c = max(len(Ec), len(Jc) + 1)
    for k in range(mode_c, max_k_c):
        ek = coeff(Ec, k)
        jk1 = coeff(Jc, k - 1)
        diff = ek - jk1
        marker = " <-- FAIL" if diff < 0 else ""
        if diff < 0:
            any_c_p3_fail = True
        print(f"    k={k:2d}: e_{k} - j_{{k-1}} = {diff:>12,d}{marker}")

    print(f"  P3 on child: {'FAIL' if any_c_p3_fail else 'PASS'}")

    # ---- LC on child I_c ----
    print(f"\n  LC on I_c:")
    any_c_lc_fail = False
    for k in range(1, len(Ic) - 1):
        lhs = Ic[k] * Ic[k]
        rhs = Ic[k-1] * Ic[k+1]
        diff = lhs - rhs
        marker = " <-- LC FAILS!" if diff < 0 else ""
        if diff < 0:
            any_c_lc_fail = True
            print(f"    k={k:2d}: {diff:>15,d}{marker}")
    if not any_c_lc_fail:
        print(f"    All pass.")
    print(f"  LC on child: {'FAIL' if any_c_lc_fail else 'PASS'}")

    # ---- Summary ----
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"  Tree: T_{{3,4}} broom, n={n}")
    print(f"  Root: vertex {root} ({labels[root]}), support vertex, s=1, ell=1")
    print(f"  mode(I(T)) = {mode_I}")
    print(f"  P2 (prefix ratio dominance up to mode): {'FAIL' if any_p2_fail else 'PASS'}")
    print(f"  P3 (tail domination from mode):         {'FAIL' if any_p3_fail else 'PASS'}")
    print(f"  W-form (Term1+Term2 >= 0):              {'FAIL' if any_wfail else 'PASS'}")
    print(f"  E_new >= J_new (full ratio dominance):   {'FAIL' if any_rd_fail else 'PASS'}")
    print(f"  I(T) log-concave:                       {'FAIL' if any_lc_fail else 'PASS'}")
    print(f"  Child P2:                               {'FAIL' if any_c_p2_fail else 'PASS'}")
    print(f"  Child P3:                               {'FAIL' if any_c_p3_fail else 'PASS'}")
    print(f"  Child LC:                               {'FAIL' if any_c_lc_fail else 'PASS'}")


if __name__ == '__main__':
    main()
