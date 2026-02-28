"""Extended verification: Condition C at ALL rootings for n up to 18.

Also: verify the two key lift properties:
  1. (A+xB, A) satisfies Cond C  =>  ((1+x)^ell*A + xB, (1+x)^ell*A) satisfies Cond C
  2. Product closure: (I1,E1) and (I2,E2) satisfy Cond C with J<=E
     => (I1*I2, E1*E2) satisfies Cond C

And: check the critical gap -- does the induction need to PRODUCE
Condition C at all rootings, or only consume it at specific rootings?
"""

import subprocess
import sys
import time
from multiprocessing import Pool
from indpoly import _polymul, _polyadd


def parse_g6(g6):
    s = g6.strip()
    n = ord(s[0]) - 63
    adj = [[] for _ in range(n)]
    bits = []
    for ch in s[1:]:
        val = ord(ch) - 63
        for shift in range(5, -1, -1):
            bits.append((val >> shift) & 1)
    k = 0
    for j in range(n):
        for i in range(j):
            if k < len(bits) and bits[k]:
                adj[i].append(j)
                adj[j].append(i)
            k += 1
    return n, adj


def coeff(poly, k):
    return poly[k] if 0 <= k < len(poly) else 0


def check_condC(I_poly, E_poly):
    """Check Strong Condition C. Returns list of (k, val) for failures."""
    d_seq = []
    L = max(len(I_poly), len(E_poly)) + 1
    for k in range(L):
        dk = coeff(I_poly, k+1)*coeff(E_poly, k) - coeff(I_poly, k)*coeff(E_poly, k+1)
        d_seq.append(dk)
    fails = []
    for k in range(1, L):
        bkm1 = coeff(E_poly, k-1)
        if bkm1 == 0:
            continue
        bk = coeff(E_poly, k)
        bkp1 = coeff(E_poly, k+1)
        akm1 = coeff(I_poly, k-1)
        ck = bk*bk - bkm1*bkp1
        val = bkm1*d_seq[k] + bk*d_seq[k-1] + akm1*ck
        if val < 0:
            fails.append((k, val))
    return fails


def do_dp(n, adj, root):
    """Full DP."""
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
                s_c = _polyadd(dp0[c], [0] + dp1s[c])
                prod_S = _polymul(prod_S, s_c)
                prod_E = _polymul(prod_E, dp0[c])
            dp0[v] = prod_S
            dp1s[v] = prod_E
    return dp0, dp1s, children


def check_tree_all_rootings(g6_str):
    """Check Condition C at ALL rootings of a tree.
    Returns (n, n_rootings, n_fails, first_fail_info).
    """
    n, adj = parse_g6(g6_str)
    if n <= 2:
        return (n, 0, 0, None)

    n_fails = 0
    first_fail = None
    for root in range(n):
        dp0, dp1s, children = do_dp(n, adj, root)
        I_T = _polyadd(dp0[root], [0] + dp1s[root])
        E_root = dp0[root]
        fails = check_condC(I_T, E_root)
        if fails:
            n_fails += 1
            if first_fail is None:
                first_fail = (root, fails[0])

    return (n, n, n_fails, first_fail)


def check_tree_gap(g6_str):
    """Check the critical gap: at each non-leaf child c of each SV r,
    does the induction need to PRODUCE Cond C at rooting c?

    Specifically: c is a child of r in T. In a LARGER tree T', c could
    be used as a factor. But what's the rooting of c?
    Answer: c is rooted at c (i.e., the subtree T_c is rooted at c).
    This is a SPECIFIC rooting -- the one induced by rooting T at r.

    So the induction only needs Cond C at "child rootings" -- rootings
    induced by being a child of some vertex in some rooting.

    But EVERY vertex can be a child of every neighbor, so every vertex
    gets every possible rooting as a child. Hence we need ALL rootings.

    UNLESS we can show: we only need Cond C at rootings where the vertex
    is a child of a SUPPORT vertex. A vertex v is a child of a SV r only
    in rootings where r is a support vertex. But the choice of SV is up
    to us -- we can pick any SV of T.

    The constraint is: the induction at T picks SV r. Then it needs
    Cond C for (I_c, E_c) at each non-leaf child c. The rooting of
    subtree(c) is at c. We need Cond C to hold at this rooting.

    Can we CHOOSE r so that all non-leaf children c have Cond C at c?
    If Cond C holds at ALL rootings, this is automatic. If not, we'd
    need a clever choice of r.

    This function checks: for each tree, for each SV r, are there any
    non-leaf children c where Cond C fails at rooting c in subtree(c)?
    """
    n, adj = parse_g6(g6_str)
    if n <= 2:
        return (n, 0, 0, 0, None)

    degree = [len(adj[v]) for v in range(n)]
    leaf_count = [0] * n
    for v in range(n):
        for u in adj[v]:
            if degree[u] == 1:
                leaf_count[v] += 1

    total_factors = 0
    factor_fails = 0
    # Also track: does every SV have at least one non-leaf child where Cond C fails?
    sv_blocked = 0  # SVs where some non-leaf child fails Cond C

    for r in range(n):
        if leaf_count[r] == 0:
            continue

        dp0, dp1s, children = do_dp(n, adj, r)
        nonleaf_children = [c for c in children[r] if children[c]]
        any_fail = False

        for c in nonleaf_children:
            I_c = _polyadd(dp0[c], [0] + dp1s[c])
            E_c = dp0[c]
            fails = check_condC(I_c, E_c)
            total_factors += 1
            if fails:
                factor_fails += 1
                any_fail = True

        if any_fail:
            sv_blocked += 1

    return (n, total_factors, factor_fails, sv_blocked, None)


# ============================================================
# Run all-rootings check for n up to 18
# ============================================================
geng = '/opt/homebrew/bin/geng'

print("=" * 80)
print("EXTENDED VERIFICATION: Condition C at ALL rootings")
print("=" * 80)
print()

t0 = time.time()
total_checks = 0
total_fails = 0

for nn in range(3, 19):
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)
    trees = [line.strip() for line in proc.stdout if line.strip()]
    proc.wait()

    nn_checks = 0
    nn_fails = 0
    nn_trees_with_fail = 0

    # Use pool for larger n
    if nn >= 12 and len(trees) > 100:
        with Pool(8) as pool:
            results = pool.map(check_tree_all_rootings, trees,
                             chunksize=max(1, len(trees) // 32))
    else:
        results = [check_tree_all_rootings(g6) for g6 in trees]

    for (n_v, n_rootings, n_f, first_fail) in results:
        nn_checks += n_rootings
        nn_fails += n_f
        if n_f > 0:
            nn_trees_with_fail += 1

    total_checks += nn_checks
    total_fails += nn_fails
    elapsed = time.time() - t0
    print(f"n={nn:2d}: {len(trees):>10,d} trees, "
          f"{nn_checks:>12,d} (tree,root) checks, "
          f"Cond C fails: {nn_fails:>6d} "
          f"({nn_trees_with_fail} trees) | {elapsed:.1f}s")

print()
print(f"TOTAL: {total_checks:,d} checks, {total_fails} failures")
print()

# ============================================================
# Run factor-level check (non-leaf children of SV)
# ============================================================
print("=" * 80)
print("FACTOR-LEVEL CHECK: Condition C at factor level (n up to 18)")
print("=" * 80)
print()

t0 = time.time()
total_factor_checks = 0
total_factor_fails = 0

for nn in range(3, 19):
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)
    trees = [line.strip() for line in proc.stdout if line.strip()]
    proc.wait()

    if nn >= 12 and len(trees) > 100:
        with Pool(8) as pool:
            results = pool.map(check_tree_gap, trees,
                             chunksize=max(1, len(trees) // 32))
    else:
        results = [check_tree_gap(g6) for g6 in trees]

    nn_factors = 0
    nn_fails = 0
    nn_sv_blocked = 0
    for (n_v, tf, ff, sb, _) in results:
        nn_factors += tf
        nn_fails += ff
        nn_sv_blocked += sb

    total_factor_checks += nn_factors
    total_factor_fails += nn_fails
    elapsed = time.time() - t0
    print(f"n={nn:2d}: {len(trees):>10,d} trees, "
          f"{nn_factors:>12,d} factor checks, "
          f"fails: {nn_fails:>6d}, SVs blocked: {nn_sv_blocked:>6d} | "
          f"{elapsed:.1f}s")

print()
print(f"TOTAL: {total_factor_checks:,d} factor checks, {total_factor_fails} failures")
print()

# ============================================================
# CRITICAL: Show why the induction closes despite the gap
# ============================================================
print("=" * 80)
print("CRITICAL ANALYSIS: Why the induction closes")
print("=" * 80)
print()
print("The subtle issue is that the induction has two jobs:")
print("  JOB 1: Prove I(T) is unimodal (the theorem)")
print("  JOB 2: Prove Condition C for (I(T), dp0[v]) at all rootings v")
print("         (the strengthened inductive hypothesis)")
print()
print("JOB 1 follows from JOB 2 via:")
print("  Pick SV r, apply Cond C to get (A,B) satisfying it,")
print("  use identity to get P2, combine with P3 for unimodality.")
print()
print("But JOB 2 is the hard part. The inductive step only directly")
print("proves Cond C at the (A,B) product level inside a specific SV.")
print("It does NOT directly prove Cond C at an arbitrary rooting v.")
print()
print("KEY REALIZATION: Condition C for (I(T), dp0[v]) at rooting v")
print("is a statement about the TREE T and the vertex v.")
print("dp0[v] = product over children-of-v of I(subtree_c).")
print("So (I(T), dp0[v]) = (dp0[v] + x*dp1s[v], dp0[v]).")
print()
print("The inductive step at SV r produces:")
print("  dp0[r] = (1+x)^ell * prod(I_c over non-leaf c)")
print("  dp1s[r] = prod(E_c over non-leaf c)")
print()
print("To get Cond C at an ARBITRARY rooting v (not just r), we'd need")
print("a re-rooting argument: show that Cond C is preserved when the")
print("root moves from r to a neighbor of r.")
print()
print("ALTERNATIVE: Prove Condition C directly for (I_c, E_c) where c")
print("is any vertex and E_c = dp0[c]. This IS what product closure does:")
print("  dp0[c] = prod over children-of-c of I(subtree_gc)")
print("  I(subtree_c) = dp0[c] + x * dp1s[c]")
print("  So (I(subtree_c), dp0[c]) has I = E + xJ with J = dp1s[c] = prod E_gc")
print("  Each factor (I_gc, E_gc) satisfies Cond C by IH")
print("  Product closure: (prod I_gc, prod E_gc) = (dp0[c], dp1s[c]) ...")
print()
print("  WAIT. Product closure gives Cond C for (prod I_gc, prod E_gc).")
print("  But we need Cond C for (I(subtree_c), dp0[c]) = (E_c + x*J_c, E_c)")
print("  where E_c = prod I_gc and J_c = prod E_gc.")
print()
print("  So I = E + xJ where E = prod I_gc, J = prod E_gc.")
print("  And Cond C for (I, E) = Cond C for (E + xJ, E).")
print("  But product closure gives Cond C for (prod I_gc, prod E_gc) = (E, J).")
print("  We need Cond C for (E + xJ, E), which is DIFFERENT from Cond C for (E, J)!")
print()
print("  Cond C for (E+xJ, E) says:")
print("    e_{k-1} * d'_k + e_k * d'_{k-1} + (e+xj)_{k-1} * c_k >= 0")
print("  where d'_k = (e+xj)_{k+1}*e_k - (e+xj)_k*e_{k+1}")
print("            = (e_{k+1}+j_k)*e_k - (e_k+j_{k-1})*e_{k+1}")
print("            = e_{k+1}*e_k + j_k*e_k - e_k*e_{k+1} - j_{k-1}*e_{k+1}")
print("            = j_k*e_k - j_{k-1}*e_{k+1}")
print("            = d_k(J,E) (the LR minor of J vs E!)")
print()
print("  And c_k = e_k^2 - e_{k-1}*e_{k+1} (LC gap of E, same as before)")
print("  And (e+xj)_{k-1} = e_{k-1} + j_{k-2}")
print()
print("  So Cond C for (I=E+xJ, E) becomes:")
print("    e_{k-1}*d_k(J,E) + e_k*d_{k-1}(J,E) + (e_{k-1}+j_{k-2})*c_k >= 0")
print()
print("  This is NOT the same as Cond C for (E, J) which would be:")
print("    j_{k-1}*d_k(E,J) + j_k*d_{k-1}(E,J) + e_{k-1}*c'_k >= 0")
print("  where c'_k = j_k^2 - j_{k-1}*j_{k+1}")
print()
print("  So product closure of Cond C for (E, J) is IRRELEVANT.")
print("  What we need is Cond C for (E+xJ, E), which is a different object.")
print()
print("  But this is exactly what the factor-level scan checks!")
print("  Each factor (I_c, E_c) = (E_c + x*J_c, E_c), and the scan")
print("  verifies Cond C for this pair.")
print()
print("  THE INDUCTION SHOULD BE:")
print("  IH: For every tree T' on < n vertices, for every rooting v,")
print("       Cond C holds for (I(T'), dp0[v]) = (dp0[v]+x*dp1s[v], dp0[v]).")
print()
print("  Step: At rooting v of T with children c_1,...,c_d:")
print("    dp0[v] = prod_c I_c  where I_c = dp0[c]+x*dp1s[c]")
print("    dp1s[v] = prod_c dp0[c] = prod_c E_c")
print("    I(T) at v = dp0[v] + x*dp1s[v] = (prod I_c) + x*(prod E_c)")
print()
print("    By IH, each (I_c, E_c) satisfies Cond C.")
print("    Need: (prod I_c + x*prod E_c, prod I_c) satisfies Cond C.")
print()
print("    This requires: E = prod I_c, J = prod E_c, then")
print("    Cond C for (E+xJ, E).")
print()
print("    Note: E is a product of IS polys, J is a product of exclude polys.")
print("    The factors (I_c, E_c) satisfy Cond C and J_c <= E_c.")
print("    Product closure for Cond C of (E+xJ, E) is the key claim.")
print()
print("  This is a DIFFERENT product closure from '(I1*I2, E1*E2) preserves")
print("  Cond C'. It's: given factors satisfying Cond C, the AGGREGATE")
print("  (prod I_c + x*prod E_c, prod I_c) satisfies Cond C.")
print()
print("  These are different because the aggregate introduces a NEW x shift")
print("  (from the dp1 = x*dp1s structure), while the factors have their OWN")
print("  x shifts baked in.")
