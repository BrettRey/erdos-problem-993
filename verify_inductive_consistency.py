"""Verify inductive consistency of the Condition C framework.

THE QUESTION:
The proposed induction proves unimodality by:
  1. Pick support vertex r (adjacent to >= 1 leaf)
  2. P3 proved at r (leaf-swap injection)
  3. For P2: E = (1+x)^ell * A, J = B where
       A = prod_c I_c over non-leaf children c
       B = prod_c E_c over non-leaf children c
  4. Each factor (I_c, E_c) must satisfy Condition C + the side constraints
  5. Product closure => (A,B) satisfies Condition C => P2 holds

But what IS (I_c, E_c) at the factor level?

At non-leaf child c of r:
  - I_c = dp0[c] + x*dp1s[c] = IS polynomial of subtree(c)
  - E_c = dp0[c] = "exclude c" polynomial

Meanwhile, if we apply the induction to subtree(c) at c's own support vertex r_c:
  - c might or might not be a support vertex of subtree(c)
  - If c IS a support vertex: E_c_root = (1+x)^ell_c * A_c, J_c_root = B_c
  - The induction at subtree(c) would prove P* for (E_c_root, J_c_root)
    which is about the full IS poly of subtree(c)

THE GAP: The induction at subtree(c) proves that I(subtree_c) is unimodal.
But at the FACTOR level, we need Condition C for the pair (I_c, E_c).
The induction doesn't directly prove Condition C for (I_c, E_c) --
it proves P* (prefix ratio dom + tail dom) for (E_c_root, J_c_root).

So the question is: does proving P* for a subtree imply that the factor
pair (I_c, E_c) satisfies Condition C?

More precisely:
  Level 0: Leaf factor: I = 1+x, E = 1. Condition C trivially holds.
  Level 1: Non-leaf child c of r.
    I_c = dp0[c] + x*dp1s[c]
    E_c = dp0[c]
    J_c = dp1s[c]  (= dp1[c]/x)

    Now dp0[c] = prod over children gc of c: I_gc = dp0[gc] + x*dp1s[gc]
    And dp1s[c] = prod over children gc of c: E_gc = dp0[gc]

    So I_c = prod(I_gc) + x*prod(E_gc) = A_c + x*B_c
    where A_c = prod(I_gc), B_c = prod(E_gc)

    And E_c = A_c (not (1+x)^ell * A_c!)

    If c is a support vertex of subtree(c), with ell_c leaves among c's children:
      E_c_inductive = (1+x)^{ell_c} * A_c'  where A_c' = prod over non-leaf children of c
      But A_c = (1+x)^{ell_c} * A_c'  (since leaf children contribute (1+x) each)
      So E_c = A_c = (1+x)^{ell_c} * A_c'

    And B_c = prod(E_gc) = (1)^{ell_c} * prod_{gc non-leaf} E_gc = B_c'
      (since E_gc = dp0[gc] = [1] for leaf gc)

    So at the factor level: (I_c, E_c) = (A_c + x*B_c, A_c)
    And at the inductive level inside subtree(c): (E_c, J_c) = (A_c, B_c)

    THE CONDITION C CHECK at factor level is on (I_c, E_c) = (A_c + x*B_c, A_c)
    = (E_c + x*J_c, E_c) where J_c = B_c.

    The Condition C at the inductive level (inside subtree c, at support vertex c)
    would be on (A_c, B_c) -- but wait, the induction proves P* for the FULL IS poly
    I(subtree_c) = E_c + x*J_c at rooting c. P* is (P2 + P3) for the (E_c, J_c) pair.

    But the FACTOR-LEVEL Condition C is a STRONGER statement about (I_c, E_c):
    it asks whether certain 3-term sums involving I_c and E_c are nonneg.

    So the key question: does P* for (E_c, J_c) at support vertex c
    IMPLY Condition C for (I_c, E_c) = (E_c + x*J_c, E_c)?

This script checks this empirically.
"""

import subprocess
import sys
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


def get_mode(poly):
    if not poly:
        return 0
    mx = max(poly)
    for k, v in enumerate(poly):
        if v == mx:
            return k
    return 0


def check_condC(I_poly, E_poly):
    """Check Strong Condition C for pair (I, E).
    b_{k-1}*d_k + b_k*d_{k-1} + a_{k-1}*c_k >= 0 for k >= 1 with b_{k-1} > 0.
    Here a = I coefficients, b = E coefficients.
    d_k = a_{k+1}*b_k - a_k*b_{k+1} (LR minor)
    c_k = b_k^2 - b_{k-1}*b_{k+1} (LC gap of E)
    """
    d_seq = []
    L = max(len(I_poly), len(E_poly)) + 1
    for k in range(L):
        dk = coeff(I_poly, k+1)*coeff(E_poly, k) - coeff(I_poly, k)*coeff(E_poly, k+1)
        d_seq.append(dk)

    fail_indices = []
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
            fail_indices.append((k, val))
    return fail_indices  # empty = pass


def check_p2(E, J, m):
    """Check P2: e_{k+1}*j_k >= e_k*j_{k+1} for k=0,...,m-1."""
    fails = []
    for k in range(m):
        lhs = coeff(E, k+1) * coeff(J, k)
        rhs = coeff(E, k) * coeff(J, k+1)
        if lhs < rhs:
            fails.append(k)
    return fails


def check_p3(E, J, m):
    """Check P3: e_k >= j_{k-1} for k >= m."""
    fails = []
    max_k = max(len(E), len(J) + 1)
    for k in range(m, max_k):
        ek = coeff(E, k)
        jkm1 = coeff(J, k-1)
        if ek < jkm1:
            fails.append(k)
    return fails


def is_lc(seq):
    """Check log-concavity."""
    for k in range(1, len(seq) - 1):
        if seq[k]*seq[k] < seq[k-1]*seq[k+1]:
            return False
    return True


def do_dp(n, adj, root):
    """Full DP returning (dp0, dp1s, children) for all vertices."""
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


def extract_subtree(n, adj, root, subtree_root):
    """Extract the subtree rooted at subtree_root (when tree is rooted at root).
    Returns (sub_n, sub_adj) with relabeled vertices.
    """
    # First root at root to get parent/children
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

    # BFS from subtree_root using only tree-children
    sub_nodes = []
    q = [subtree_root]
    head = 0
    while head < len(q):
        v = q[head]; head += 1
        sub_nodes.append(v)
        for c in children[v]:
            q.append(c)

    # Relabel
    mapping = {old: new for new, old in enumerate(sub_nodes)}
    sub_n = len(sub_nodes)
    sub_adj = [[] for _ in range(sub_n)]
    for old_v in sub_nodes:
        for old_u in adj[old_v]:
            if old_u in mapping:
                new_v = mapping[old_v]
                new_u = mapping[old_u]
                if new_u not in sub_adj[new_v]:
                    sub_adj[new_v].append(new_u)
    return sub_n, sub_adj, mapping


def polymul_1px_pow(poly, ell):
    """Multiply polynomial by (1+x)^ell."""
    result = list(poly)
    for _ in range(ell):
        result = _polyadd(result, [0] + result)
    return result


# ============================================================
# MAIN ANALYSIS
# ============================================================

geng = '/opt/homebrew/bin/geng'
MAX_N = 16

# Counters
total_factors = 0
total_sv_factors = 0  # factors where c is a support vertex of subtree(c)
total_nonsv_factors = 0  # factors where c is NOT a support vertex

# Category 1: c is a support vertex of its subtree
cat1_condC_factor_fails = 0  # Condition C fails at factor level (I_c, E_c)
cat1_condC_inner_fails = 0   # Condition C fails for (A_c, B_c) inside subtree
cat1_pstar_fails = 0         # P* fails at c in subtree(c)
cat1_agree = 0               # factor condC and inner condC both pass
cat1_disagree = 0            # they disagree

# Category 2: c is NOT a support vertex (no leaf children)
cat2_condC_factor_fails = 0
cat2_has_sv_in_subtree = 0   # subtree(c) has some other support vertex

# Key relationship checks
factor_condC_needs_inner_condC = True   # can factor condC hold without inner condC?
inner_condC_needs_factor_condC = True   # can inner condC hold without factor condC?

# Track J <= E at factor level
jle_e_fails = 0

# Detailed examples
examples = []

print("=" * 80)
print("INDUCTIVE CONSISTENCY VERIFICATION")
print("Checking whether Condition C at the factor level matches the inductive claim")
print("=" * 80)
print()

for nn in range(3, MAX_N + 1):
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)

    nn_factors = 0
    nn_sv = 0
    nn_nonsv = 0
    nn_factor_condC_fail = 0
    nn_inner_condC_fail = 0

    for line in proc.stdout:
        g6 = line.strip()
        n, adj = parse_g6(g6)
        if n <= 2:
            continue

        # Find support vertices (vertices adjacent to at least one leaf)
        leaf_count = [0] * n
        degree = [len(adj[v]) for v in range(n)]
        for v in range(n):
            for u in adj[v]:
                if degree[u] == 1:
                    leaf_count[v] += 1

        for r in range(n):
            if leaf_count[r] == 0:
                continue

            dp0, dp1s, children = do_dp(n, adj, r)
            ell_r = sum(1 for c in children[r] if not children[c])
            nonleaf_children = [c for c in children[r] if children[c]]

            for c in nonleaf_children:
                total_factors += 1
                nn_factors += 1

                # Factor pair at the parent level
                I_c = _polyadd(dp0[c], [0] + dp1s[c])
                E_c = dp0[c]
                J_c = dp1s[c]

                # Check J_c <= E_c coefficientwise
                jle = all(coeff(J_c, k) <= coeff(E_c, k) for k in range(max(len(J_c), len(E_c))))
                if not jle:
                    jle_e_fails += 1

                # Check Condition C at factor level: (I_c, E_c)
                factor_condC = check_condC(I_c, E_c)
                factor_ok = len(factor_condC) == 0
                if not factor_ok:
                    nn_factor_condC_fail += 1

                # Now look at c's own structure within subtree(c)
                # c's children in the rooted tree
                c_children = children[c]
                c_leaf_children = [gc for gc in c_children if not children[gc]]
                c_nonleaf_children = [gc for gc in c_children if children[gc]]
                ell_c = len(c_leaf_children)

                if ell_c > 0:
                    # c IS a support vertex of subtree(c)
                    total_sv_factors += 1
                    nn_sv += 1

                    # The inductive decomposition at c:
                    # A_c = prod over gc non-leaf children: I_gc
                    # B_c = prod over gc non-leaf children: E_gc
                    # Then E_c = dp0[c] = (1+x)^ell_c * A_c
                    # And J_c = dp1s[c] = B_c (leaf E_gc = [1])

                    # Verify: E_c should equal (1+x)^ell_c * A_c
                    if c_nonleaf_children:
                        A_c = [1]
                        B_c = [1]
                        for gc in c_nonleaf_children:
                            I_gc = _polyadd(dp0[gc], [0] + dp1s[gc])
                            E_gc = dp0[gc]
                            A_c = _polymul(A_c, I_gc)
                            B_c = _polymul(B_c, E_gc)
                    else:
                        A_c = [1]
                        B_c = [1]

                    # Verify E_c = (1+x)^ell_c * A_c
                    expected_E = polymul_1px_pow(A_c, ell_c)
                    E_c_match = (E_c == expected_E)
                    if not E_c_match and len(examples) < 3:
                        examples.append(('E_MISMATCH', g6, n, r, c, ell_c,
                                         E_c, expected_E, A_c))

                    # Verify J_c = B_c
                    J_c_match = (J_c == B_c)
                    if not J_c_match and len(examples) < 3:
                        examples.append(('J_MISMATCH', g6, n, r, c, ell_c,
                                         J_c, B_c))

                    # Check Condition C for (A_c, B_c) -- the "inner" level
                    # This is what the next level of induction would need:
                    # We need Condition C for (A_c + x*B_c, A_c) = (I_sub, E_sub)
                    # where I_sub = A_c + x*B_c and E_sub = A_c
                    # Wait... that's exactly (I_c_sub, E_c_sub) if we rewrite
                    # I_c = (1+x)^ell_c * A_c + x*B_c
                    #      = E_c + x*J_c = E_c + x*B_c
                    # So (I_c, E_c) at the factor level includes the (1+x)^ell_c
                    # absorbed into A_c.

                    # The INDUCTIVE level would produce Condition C for (A_c + x*B_c, A_c)
                    # where A_c is the product of NON-LEAF grandchild I_gc's.
                    # But (I_c, E_c) = (E_c + x*J_c, E_c) = ((1+x)^ell_c * A_c + x*B_c, (1+x)^ell_c * A_c)
                    # This is NOT the same as (A_c + x*B_c, A_c).

                    # Check Condition C for (A_c + x*B_c, A_c) -- the "naked" inner pair
                    I_inner = _polyadd(A_c, [0] + B_c)
                    E_inner = A_c
                    inner_condC = check_condC(I_inner, E_inner)
                    inner_ok = len(inner_condC) == 0

                    if not inner_ok:
                        cat1_condC_inner_fails += 1
                        nn_inner_condC_fail += 1

                    if not factor_ok:
                        cat1_condC_factor_fails += 1

                    if factor_ok and inner_ok:
                        cat1_agree += 1
                    elif not factor_ok and not inner_ok:
                        cat1_agree += 1
                    else:
                        cat1_disagree += 1
                        if len(examples) < 5:
                            examples.append(('DISAGREE_SV', g6, n, r, c, ell_c,
                                             factor_ok, inner_ok,
                                             I_c, E_c, A_c, B_c,
                                             factor_condC, inner_condC))

                    # Also check: P* at c in subtree(c)
                    # E_at_c = dp0[c] = E_c, J_at_c = dp1s[c] = J_c
                    poly_c = _polyadd(E_c, [0] + J_c)
                    m_c = get_mode(poly_c)
                    p2_fails = check_p2(E_c, J_c, m_c)
                    p3_fails = check_p3(E_c, J_c, m_c)
                    if p2_fails or p3_fails:
                        cat1_pstar_fails += 1

                else:
                    # c is NOT a support vertex (has no leaf children)
                    total_nonsv_factors += 1
                    nn_nonsv += 1
                    if not factor_ok:
                        cat2_condC_factor_fails += 1

                    # Does subtree(c) have any support vertex?
                    # Check if any vertex in subtree has a leaf child
                    has_sv = False
                    q = [c]
                    head = 0
                    while head < len(q):
                        v = q[head]; head += 1
                        for gc in children[v]:
                            if children[gc]:
                                # gc has children
                                gc_has_leaf = any(not children[ggc] for ggc in children[gc])
                                if gc_has_leaf:
                                    has_sv = True
                            q.append(gc)
                    if has_sv:
                        cat2_has_sv_in_subtree += 1

    proc.wait()

    print(f"n={nn:2d}: factors={nn_factors:>8,d} "
          f"(sv={nn_sv:>7,d} nonsv={nn_nonsv:>7,d}) "
          f"factor_condC_fail={nn_factor_condC_fail} "
          f"inner_condC_fail={nn_inner_condC_fail}")

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print(f"Total factor pairs examined: {total_factors:,d}")
print(f"  Category 1 (c is support vertex of subtree): {total_sv_factors:,d}")
print(f"  Category 2 (c is NOT support vertex):        {total_nonsv_factors:,d}")
print()
print("--- Category 1: c is a support vertex ---")
print(f"  Factor-level Cond C failures:  {cat1_condC_factor_fails}")
print(f"  Inner (A_c,B_c) Cond C failures: {cat1_condC_inner_fails}")
print(f"  Factor & inner agree:  {cat1_agree}")
print(f"  Factor & inner DISAGREE: {cat1_disagree}")
print(f"  P* fails at c in subtree(c): {cat1_pstar_fails}")
print()
print("--- Category 2: c is NOT a support vertex ---")
print(f"  Factor-level Cond C failures: {cat2_condC_factor_fails}")
print(f"  Subtree has some other SV:    {cat2_has_sv_in_subtree}")
print()
print(f"J <= E failures at factor level: {jle_e_fails}")
print()

if examples:
    print("DETAILED EXAMPLES:")
    for ex in examples:
        print(f"  Type: {ex[0]}")
        if ex[0] == 'E_MISMATCH':
            _, g6, n, r, c, ell_c, E_c, expected, A_c = ex
            print(f"    g6={g6} n={n} root={r} child={c} ell_c={ell_c}")
            print(f"    E_c (actual)   = {E_c}")
            print(f"    (1+x)^ell * Ac = {expected}")
            print(f"    A_c = {A_c}")
        elif ex[0] == 'J_MISMATCH':
            _, g6, n, r, c, ell_c, J_c, B_c = ex
            print(f"    g6={g6} n={n} root={r} child={c} ell_c={ell_c}")
            print(f"    J_c = {J_c}")
            print(f"    B_c = {B_c}")
        elif ex[0] == 'DISAGREE_SV':
            _, g6, n, r, c, ell_c, fok, iok, I_c, E_c, A_c, B_c, ff, iif = ex
            print(f"    g6={g6} n={n} root={r} child={c} ell_c={ell_c}")
            print(f"    Factor (I_c,E_c) condC: {'PASS' if fok else 'FAIL at '+str(ff)}")
            print(f"    Inner  (Ac+xBc,Ac) condC: {'PASS' if iok else 'FAIL at '+str(iif)}")
            print(f"    I_c = {I_c}")
            print(f"    E_c = {E_c}")
            print(f"    A_c = {A_c}")
            print(f"    B_c = {B_c}")
        print()

# ============================================================
# KEY ANALYSIS: The inductive structure
# ============================================================
print("=" * 80)
print("ANALYSIS: Does the induction close?")
print("=" * 80)
print()
print("The factor pair (I_c, E_c) at a non-leaf child c of support vertex r is:")
print("  I_c = dp0[c] + x * dp1s[c]  = IS poly of subtree(c)")
print("  E_c = dp0[c]                 = exclude-c poly")
print()
print("If c is a support vertex of subtree(c) with ell_c leaf children:")
print("  E_c = (1+x)^ell_c * A_c     where A_c = prod(I_gc, gc non-leaf)")
print("  J_c = B_c                    where B_c = prod(E_gc, gc non-leaf)")
print("  I_c = E_c + x*J_c = (1+x)^ell_c * A_c + x * B_c")
print()
print("The INDUCTION at subtree(c) proves P* for (E_c, J_c).")
print("The FACTOR-LEVEL requirement is Condition C for (I_c, E_c).")
print()
print("These are RELATED but DIFFERENT statements:")
print("  P* for (E_c, J_c) is about ratio j_k/e_k being nonincreasing")
print("  Condition C for (I_c, E_c) is about b_{k-1}*d_k + b_k*d_{k-1} + a_{k-1}*c_k >= 0")
print("  where a = I_c coeffs, b = E_c coeffs")
print()
print("The KEY QUESTION: does the induction need to prove Condition C")
print("at EVERY subtree, or just P* (which implies unimodality)?")
print()
print("ANSWER: The induction needs Condition C, which is STRONGER than P*.")
print("P* => unimodality, but product closure requires Condition C on factors.")
print("So the inductive hypothesis must be Condition C (not just unimodality).")
print()

# ============================================================
# PART 2: Check whether the (1+x)^ell wrapping helps or hurts
# ============================================================
print("=" * 80)
print("PART 2: Effect of (1+x)^ell wrapping on Condition C")
print("=" * 80)
print()
print("At factor level: (I_c, E_c) has E_c = (1+x)^ell_c * A_c")
print("'Naked' inner level: (A_c + x*B_c, A_c)")
print()
print("Does (1+x)^ell wrapping make Condition C easier or harder?")
print()

# Systematic check: for each factor with ell_c >= 1,
# check Condition C for both (I_c, E_c) and (A_c + x*B_c, A_c)
wrap_helps = 0
wrap_hurts = 0
wrap_neutral = 0
naked_fail_wrap_ok = 0
naked_ok_wrap_fail = 0

for nn in range(3, min(MAX_N + 1, 13)):  # smaller range for this analysis
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)

    for line in proc.stdout:
        g6 = line.strip()
        n, adj = parse_g6(g6)
        if n <= 2:
            continue

        degree = [len(adj[v]) for v in range(n)]
        leaf_count = [0] * n
        for v in range(n):
            for u in adj[v]:
                if degree[u] == 1:
                    leaf_count[v] += 1

        for r in range(n):
            if leaf_count[r] == 0:
                continue

            dp0, dp1s, children = do_dp(n, adj, r)
            nonleaf_children = [c for c in children[r] if children[c]]

            for c in nonleaf_children:
                c_children = children[c]
                c_nonleaf = [gc for gc in c_children if children[gc]]
                ell_c = sum(1 for gc in c_children if not children[gc])

                if ell_c == 0:
                    continue  # no wrapping, same pair

                # Factor level
                I_c = _polyadd(dp0[c], [0] + dp1s[c])
                E_c = dp0[c]
                J_c = dp1s[c]

                # Inner (naked) level
                if c_nonleaf:
                    A_c = [1]
                    B_c = [1]
                    for gc in c_nonleaf:
                        I_gc = _polyadd(dp0[gc], [0] + dp1s[gc])
                        E_gc = dp0[gc]
                        A_c = _polymul(A_c, I_gc)
                        B_c = _polymul(B_c, E_gc)
                else:
                    A_c = [1]
                    B_c = [1]

                I_naked = _polyadd(A_c, [0] + B_c)

                factor_fails = check_condC(I_c, E_c)
                naked_fails = check_condC(I_naked, A_c)

                factor_ok = len(factor_fails) == 0
                naked_ok = len(naked_fails) == 0

                if factor_ok and naked_ok:
                    wrap_neutral += 1
                elif factor_ok and not naked_ok:
                    wrap_helps += 1
                    naked_fail_wrap_ok += 1
                elif not factor_ok and naked_ok:
                    wrap_hurts += 1
                    naked_ok_wrap_fail += 1
                else:
                    wrap_neutral += 1  # both fail

    proc.wait()

print(f"For factors where c is SV (ell_c >= 1), n <= {min(MAX_N, 12)}:")
print(f"  Both pass:              {wrap_neutral:,d}")
print(f"  Wrapping HELPS:         {wrap_helps} (naked fails, wrapped passes)")
print(f"  Wrapping HURTS:         {wrap_hurts} (naked passes, wrapped fails)")
print(f"  naked_fail, wrap_ok:    {naked_fail_wrap_ok}")
print(f"  naked_ok, wrap_fail:    {naked_ok_wrap_fail}")
print()

# ============================================================
# PART 3: What does the induction ACTUALLY need?
# ============================================================
print("=" * 80)
print("PART 3: Recursive structure of the inductive claim")
print("=" * 80)
print()
print("The induction proves: for every tree T, at every support vertex r,")
print("the factor pair (I_c, E_c) for each non-leaf child c satisfies")
print("Condition C.")
print()
print("But (I_c, E_c) = (E_c + x*J_c, E_c) where:")
print("  E_c = dp0[c] = prod over gc of I_gc")
print("  J_c = dp1s[c] = prod over gc of E_gc")
print()
print("And each I_gc = E_gc + x*J_gc is the IS poly of subtree(gc).")
print("Each E_gc = dp0[gc] = the exclude-gc poly.")
print()
print("So the inductive hypothesis at child c IS Condition C for (I_c, E_c),")
print("which is the SAME form as what we need at the root level.")
print()
print("THE INDUCTIVE CLAIM IS: For every tree T rooted at support vertex r,")
print("Condition C holds for (I(T), dp0[r]) where I(T) is the IS poly.")
print()
print("At the root level, we need Condition C for each factor (I_c, E_c).")
print("Each factor has the same form: (I(subtree_c), dp0[c]).")
print()
print("So the induction IS self-consistent IF Condition C for (I_c, E_c)")
print("can be proved from Condition C for the grandchild factors (I_gc, E_gc).")
print("This is exactly the product closure question.")
print()

# ============================================================
# PART 4: But c might not be a support vertex!
# ============================================================
print("=" * 80)
print("PART 4: Non-support-vertex children")
print("=" * 80)
print()
print("If c is NOT a support vertex of subtree(c), then the induction can't")
print("directly apply at c. We need Condition C for (I_c, E_c) but there's")
print("no support vertex at c to anchor the proof.")
print()
print("However, subtree(c) still has support vertices (every tree with >= 2")
print("vertices has at least one support vertex). The induction can be applied")
print("at some OTHER vertex r_c in subtree(c).")
print()
print("But that proves P* for I(subtree_c) at rooting r_c, not Condition C")
print("for (I_c, E_c) at rooting c.")
print()
print("KEY INSIGHT: The inductive claim should be:")
print("'For every tree T, for every vertex v, Condition C holds for")
print(" (I(T), dp0[v]) where dp0[v] is computed when T is rooted at v.'")
print()
print("This is stronger than 'Condition C at support vertices only'.")
print("Let me check if it holds universally...")
print()

# Check: does Condition C for (I_c, E_c) hold for ALL rootings, not just support vertices?
condC_all_rootings_fails = 0
condC_total_all = 0

for nn in range(3, min(MAX_N + 1, 14)):
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)
    nn_fails = 0
    nn_checks = 0

    for line in proc.stdout:
        g6 = line.strip()
        n, adj = parse_g6(g6)
        if n <= 2:
            continue

        for root in range(n):
            dp0, dp1s, children = do_dp(n, adj, root)
            I_T = _polyadd(dp0[root], [0] + dp1s[root])
            E_root = dp0[root]

            condC_result = check_condC(I_T, E_root)
            condC_total_all += 1
            nn_checks += 1
            if len(condC_result) > 0:
                condC_all_rootings_fails += 1
                nn_fails += 1

    proc.wait()
    print(f"n={nn:2d}: {nn_checks:>10,d} (tree,rooting) pairs, "
          f"Cond C fails: {nn_fails}")

print()
print(f"Total (tree,rooting) checks: {condC_total_all:,d}")
print(f"Condition C failures at ALL rootings: {condC_all_rootings_fails}")
print()
if condC_all_rootings_fails == 0:
    print("RESULT: Condition C holds at ALL rootings, not just support vertices!")
    print("This means the induction can use ANY rooting, including non-SV rootings.")
    print("The inductive claim 'Condition C for (I(T), dp0[r]) at every rooting r'")
    print("is self-consistent because:")
    print("  1. At each non-leaf child c, we need Condition C for (I_c, E_c)")
    print("  2. I_c = I(subtree_c) and E_c = dp0[c] at rooting c")
    print("  3. By induction on subtree_c at rooting c, this holds")
    print("  4. Product closure gives Condition C for the aggregate (A, B)")
    print("  5. The identity gives P2, which with P3 gives unimodality")
else:
    print(f"WARNING: Condition C FAILS at some rootings ({condC_all_rootings_fails} failures)")
    print("The induction must be restricted to specific rootings.")
    print("Need to check: do the failure rootings include any that the")
    print("induction would actually need?")

print()
print("=" * 80)
print("PART 5: Factor-level Condition C vs root-level Condition C")
print("=" * 80)
print()
print("Clarifying the levels:")
print()
print("ROOT LEVEL: At support vertex r of T:")
print("  E_r = (1+x)^ell * A, J_r = B")
print("  Condition C for (A, B) => P2 for (E_r, J_r) => unimodality")
print()
print("  Where A = prod(I_c), B = prod(E_c) over non-leaf children c.")
print()
print("FACTOR LEVEL: Each factor (I_c, E_c) is:")
print("  I_c = I(subtree_c) rooted at c")
print("  E_c = dp0[c] at rooting c")
print()
print("INDUCTIVE STRUCTURE:")
print("  Need: Condition C for (A, B) = (prod I_c, prod E_c)")
print("  Have: Condition C for each (I_c, E_c)")
print("  Bridge: PRODUCT CLOSURE of Condition C")
print()
print("  But wait -- Condition C for (prod I_c, prod E_c) is NOT the same as")
print("  Condition C at the root level, which is for (A, B) where E_r = (1+x)^ell * A.")
print()
print("  Actually, the identity decomposes P2 in terms of (A, B), and P2 follows")
print("  if the RHS (involving d_k, c_k, a_k of the (A,B) pair) is nonneg.")
print("  THAT is what 'Condition C for (A,B)' means.")
print()
print("  But (A, B) is NOT of the form (I(tree), dp0[root]) for any tree!")
print("  A = prod(I_c) is the product of IS polynomials, not itself an IS poly.")
print("  B = prod(E_c) is the product of exclude polys.")
print()
print("  However, (A, B) still has the structure that I_c = E_c + x*J_c for each")
print("  factor, so A = prod(E_c + x*J_c) and B = prod(E_c).")
print()
print("  The product closure question is: starting from factors satisfying")
print("  Condition C, does the product also satisfy Condition C?")
print("  This is purely an algebraic question about the product of pairs.")
print()

# ============================================================
# PART 6: Verify that (I(T), dp0[v]) has the same Condition C
# as the factor pair (I_c, E_c) when looking at T rooted at v
# ============================================================
print("=" * 80)
print("PART 6: Self-consistency check -- is the factor pair exactly")
print("what the induction produces?")
print("=" * 80)
print()

# For each support vertex r of a tree T:
#   For each non-leaf child c:
#     The factor pair (I_c, E_c) = (I(subtree_c), dp0_c_in_subtree)
#     This is exactly (I(T'), dp0[root(T')]) where T' = subtree(c), root = c
#     So Condition C for this pair IS the inductive claim applied to T' at c.
#
#     The induction at T' picks a support vertex r' of T'.
#     If r' = c, great. If r' != c, we need Condition C at rooting c,
#     but the induction only gives it at r'.
#
#     UNLESS Condition C holds at ALL rootings (checked in Part 4).

# Final verification: at a non-leaf child c that IS a SV,
# is Condition C for (I_c, E_c) = Condition C for ((1+x)^ell * A_c + x*B_c, (1+x)^ell * A_c)?
# And the induction inside subtree(c) would give Condition C for each (I_gc, E_gc).
# Product closure gives Condition C for (A_c, B_c).
# Then we need: Condition C for (A_c, B_c) => Condition C for
#   ((1+x)^ell * A_c + x*B_c, (1+x)^ell * A_c)

# Let's check: if Condition C holds for (A, B), does it hold for
#   ((1+x)^ell * A + x*B, (1+x)^ell * A)?
# This is the "(1+x)^ell lifting" step.

print("Check: does Condition C for (A, B) imply Condition C for")
print("  ((1+x)^ell * A + x*B, (1+x)^ell * A)?")
print()

lift_total = 0
lift_condC_base_ok_lift_fail = 0
lift_condC_base_fail = 0

for nn in range(3, min(MAX_N + 1, 13)):
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)

    for line in proc.stdout:
        g6 = line.strip()
        n, adj = parse_g6(g6)
        if n <= 2:
            continue

        degree = [len(adj[v]) for v in range(n)]
        leaf_count = [0] * n
        for v in range(n):
            for u in adj[v]:
                if degree[u] == 1:
                    leaf_count[v] += 1

        for r in range(n):
            if leaf_count[r] == 0:
                continue

            dp0, dp1s, children = do_dp(n, adj, r)
            nonleaf_children = [c for c in children[r] if children[c]]

            for c in nonleaf_children:
                c_children = children[c]
                c_nonleaf = [gc for gc in c_children if children[gc]]
                ell_c = sum(1 for gc in c_children if not children[gc])

                if ell_c == 0:
                    continue

                if c_nonleaf:
                    A_c = [1]
                    B_c = [1]
                    for gc in c_nonleaf:
                        I_gc = _polyadd(dp0[gc], [0] + dp1s[gc])
                        E_gc = dp0[gc]
                        A_c = _polymul(A_c, I_gc)
                        B_c = _polymul(B_c, E_gc)
                else:
                    A_c = [1]
                    B_c = [1]

                # Check Condition C for (A_c + x*B_c, A_c)
                I_base = _polyadd(A_c, [0] + B_c)
                base_fails = check_condC(I_base, A_c)

                # Check Condition C for ((1+x)^ell*A_c + x*B_c, (1+x)^ell*A_c)
                lifted_A = polymul_1px_pow(A_c, ell_c)
                I_lifted = _polyadd(lifted_A, [0] + B_c)
                lifted_fails = check_condC(I_lifted, lifted_A)

                lift_total += 1
                if len(base_fails) > 0:
                    lift_condC_base_fail += 1
                elif len(lifted_fails) > 0:
                    lift_condC_base_ok_lift_fail += 1

    proc.wait()

print(f"Total (A_c, B_c) pairs with ell >= 1, n <= {min(MAX_N, 12)}: {lift_total:,d}")
print(f"Base (A+xB, A) fails:  {lift_condC_base_fail}")
print(f"Base OK but lift FAILS: {lift_condC_base_ok_lift_fail}")
print()
if lift_condC_base_ok_lift_fail == 0:
    print("GOOD: (1+x)^ell lifting never breaks Condition C.")
    print("If Condition C holds for (A+xB, A), it also holds for")
    print("((1+x)^ell*A + xB, (1+x)^ell*A).")
else:
    print(f"WARNING: {lift_condC_base_ok_lift_fail} cases where lifting breaks Cond C!")

print()
print("=" * 80)
print("FINAL CONCLUSIONS")
print("=" * 80)
print()
print("The inductive proof structure is:")
print()
print("  INDUCTIVE HYPOTHESIS (IH):")
print("    For every tree T' on < n vertices, for every vertex v of T',")
print("    Condition C holds for (I(T'), dp0[v]) where dp0[v] is computed")
print("    when T' is rooted at v.")
print()
print("  INDUCTIVE STEP:")
print("    Let T be a tree on n vertices. Pick support vertex r.")
print("    Non-leaf children c_1,...,c_s have subtrees T_1,...,T_s.")
print("    By IH, Condition C holds for (I(T_i), dp0[c_i]) at rooting c_i.")
print("    By product closure, Condition C holds for (prod I(T_i), prod dp0[c_i]).")
print("    The identity + (1+x)^ell smoothing gives P2.")
print("    P3 is proved. So I(T) is unimodal.")
print()
print("  BUT WAIT: IH says Condition C at ALL rootings.")
print("  The inductive step only proves P* (unimodality) at one rooting (support vertex).")
print("  It does NOT directly prove Condition C at all rootings of T!")
print()
print("  This is a GAP unless either:")
print("  (a) Condition C holds at all rootings (verified up to n=13), or")
print("  (b) The induction only needs Condition C at rootings used as children")
print("      of support vertices -- which means at ALL vertices, since any vertex")
print("      can be a child of a support vertex in a larger tree.")
print()
if condC_all_rootings_fails == 0:
    print("  Since Condition C holds at ALL rootings (verified), option (a) works.")
    print("  The induction IS self-consistent.")
    print()
    print("  REMAINING ALGEBRAIC CHALLENGES:")
    print("  1. Prove product closure of Condition C (with J <= E on factors)")
    print("  2. Prove that the (1+x)^ell lifting preserves Condition C")
    print("  3. Prove Condition C at ALL rootings (not just support vertices)")
    print("     -- or restructure to only need specific rootings")
else:
    print("  Condition C FAILS at some rootings. Need to investigate further.")
