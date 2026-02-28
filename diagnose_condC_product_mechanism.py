#!/usr/bin/env python3
"""Diagnose the product closure mechanism of Strong Condition C.

At a support vertex r (adjacent to at least one leaf), the DP gives:
    E = (1+x)^ell * prod_c E_c,   J = prod_c E_c
where ell = number of leaf children, and the products run over non-leaf
children c of r.  For each such child c:
    I_c = E_c + x * J_c     (IS poly of subtree rooted at c)
    E_c = dp[c][0]           (c excluded)
    J_c = dp[c][1] / x       (c included, stripped of leading x)

Strong Condition C for a pair (I, E) with a = coeff(I), b = coeff(E):
    d_k = a_{k+1} b_k - a_k b_{k+1}       (LR minor / ratio-dominance)
    c_k = b_k^2 - b_{k-1} b_{k+1}         (LC curvature of E)
    CondC_k = b_{k-1} d_k + b_k d_{k-1} + a_{k-1} c_k   for k >= 1

This script:
1. Extracts factor pairs (I_c, E_c) from tree DP at support vertices
2. Computes all pairwise products (I1*I2, E1*E2)
3. Decomposes Condition C into its three terms when d_k < 0
4. Decomposes d_k^{prod} using Cauchy-Binet identity
5. Finds the most adversarial product pairs
"""

from __future__ import annotations

import subprocess
import sys
import time
from collections import defaultdict

from graph6 import parse_graph6
from indpoly import _polymul, _polyadd, is_log_concave


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def coeff(poly: list[int], k: int) -> int:
    """Safe coefficient access."""
    return poly[k] if 0 <= k < len(poly) else 0


def mode_index(poly: list[int]) -> int:
    """Leftmost mode index."""
    if not poly:
        return 0
    mx = max(poly)
    for k, v in enumerate(poly):
        if v == mx:
            return k
    return 0


def is_dleaf_le_1(n: int, adj: list[list[int]]) -> bool:
    """Check d_leaf <= 1: no vertex has more than one leaf neighbour."""
    leaves = {v for v in range(n) if len(adj[v]) == 1}
    for v in range(n):
        if sum(1 for u in adj[v] if u in leaves) > 1:
            return False
    return True


def dp_tree(n: int, adj: list[list[int]], root: int):
    """Full DP returning dp0[v], dp1_stripped[v] for every vertex.

    dp0[v] = E_v = IS poly of subtree(v) with v excluded
    dp1_stripped[v] = J_v = dp[v][1]/x = product of dp0[c] over children c

    I(subtree(v)) = dp0[v] + x * dp1_stripped[v]
    """
    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[root] = True
    queue = [root]
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
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
    dp1s = [None] * n  # dp1_stripped

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

    return dp0, dp1s, children, parent


# ---------------------------------------------------------------------------
# Condition C computation
# ---------------------------------------------------------------------------

def compute_d(I_poly: list[int], E_poly: list[int], k: int) -> int:
    """d_k = a_{k+1} * b_k - a_k * b_{k+1} where a=I, b=E."""
    return coeff(I_poly, k + 1) * coeff(E_poly, k) - coeff(I_poly, k) * coeff(E_poly, k + 1)


def compute_c(E_poly: list[int], k: int) -> int:
    """c_k = b_k^2 - b_{k-1} * b_{k+1} (LC curvature)."""
    bk = coeff(E_poly, k)
    return bk * bk - coeff(E_poly, k - 1) * coeff(E_poly, k + 1)


def condition_C(I_poly: list[int], E_poly: list[int], k: int) -> tuple[int, int, int, int]:
    """Compute the three terms of Condition C at index k.

    Returns (T1, T2, T3, total) where:
        T1 = b_{k-1} * d_k       (current LR minor, weighted by prev coeff)
        T2 = b_k * d_{k-1}       (memory term: previous LR minor)
        T3 = a_{k-1} * c_k       (curvature bonus: LC of E weighted by I)
        total = T1 + T2 + T3
    """
    dk = compute_d(I_poly, E_poly, k)
    dk_prev = compute_d(I_poly, E_poly, k - 1)
    ck = compute_c(E_poly, k)
    bk = coeff(E_poly, k)
    bk_prev = coeff(E_poly, k - 1)
    ak_prev = coeff(I_poly, k - 1)

    T1 = bk_prev * dk
    T2 = bk * dk_prev
    T3 = ak_prev * ck
    return T1, T2, T3, T1 + T2 + T3


# ---------------------------------------------------------------------------
# d_k product decomposition via Cauchy-Binet
# ---------------------------------------------------------------------------

def decompose_dk_product(I1: list[int], E1: list[int],
                         I2: list[int], E2: list[int], k: int) -> dict:
    """Decompose d_k^{prod} using the Cauchy-Binet identity.

    For products A = I1*I2, B = E1*E2:
    d_k^{prod} = A_{k+1}*B_k - A_k*B_{k+1}

    By Cauchy-Binet (2x2 minors of products of matrices):
    d_k^{prod} = sum_{(i,j): i+j=k} sum_{(i',j'): i'+j'=k}
                 det[[a1_{i+1}*a2_{j+1}, a1_{i'+1}*a2_{j'+1}],
                     [b1_i*b2_j,         b1_{i'}*b2_{j'}]]  ... (wrong form)

    Actually, let's use the direct expansion. Define:
      r_f(i) = a_f(i+1) * b_f(i)   ("up-cross")
      s_f(i) = a_f(i) * b_f(i+1)   ("down-cross")
    so d_f(i) = r_f(i) - s_f(i).

    Then d_k^{prod} = sum_{i+j=k} sum_{i'+j'=k}
        [a1_{i+1}*a2_{j+1} * b1_{i'}*b2_{j'} - a1_{i'}*a2_{j'} * b1_{i+1}*b2_{j+1}]
        ... but this double sum is the Cauchy-Binet form.

    The Cauchy-Binet identity for the determinant of a product gives:
    d_k^{prod} = sum_{0<=p<q, p+q=2k+1???} ...

    Let me use a simpler approach: just compute the "diagonal" and
    "off-diagonal" contributions.

    Diagonal (i=i', j=j'): sum_{i+j=k} [r1(i)*r2(j) - s1(i)*s2(j)]
      = sum_{i+j=k} [d1(i)*d2(j) + d1(i)*s2(j) + s1(i)*d2(j)]

    Off-diagonal (i!=i'): sum_{i<i', j=k-i, j'=k-i'}
      = sum_{i<i'} [a1_{i+1}*b1_{i'} - a1_{i'}*b1_{i+1}] * [a2_{j+1}*b2_{j'} - a2_{j'}*b2_{j+1}]
    where j = k-i, j' = k-i', so j > j'.
    The first bracket is det[[a1_{i+1}, a1_{i'}], [b1_{i+1}, b1_{i'}]]
    ... actually it's a1_{i+1}*b1_{i'} - a1_{i'*b1_{i+1}} which is a
    "shifted" minor of factor 1.

    Let's just compute diagonal vs off-diagonal numerically.
    """
    I_prod = _polymul(I1, I2)
    E_prod = _polymul(E1, E2)
    dk_direct = compute_d(I_prod, E_prod, k)

    # Diagonal contribution: sum_{i+j=k} [r1_i*r2_j - s1_i*s2_j]
    # where r_f(i) = a_f(i+1)*b_f(i), s_f(i) = a_f(i)*b_f(i+1)
    diag = 0
    diag_parts_pure = 0    # sum d1(i)*d2(j)
    diag_parts_cross1 = 0  # sum d1(i)*s2(j)
    diag_parts_cross2 = 0  # sum s1(i)*d2(j)

    for i in range(k + 1):
        j = k - i
        r1i = coeff(I1, i + 1) * coeff(E1, i)
        s1i = coeff(I1, i) * coeff(E1, i + 1)
        r2j = coeff(I2, j + 1) * coeff(E2, j)
        s2j = coeff(I2, j) * coeff(E2, j + 1)
        d1i = r1i - s1i
        d2j = r2j - s2j
        diag += r1i * r2j - s1i * s2j
        diag_parts_pure += d1i * d2j
        diag_parts_cross1 += d1i * s2j
        diag_parts_cross2 += s1i * d2j

    # Off-diagonal: sum over i < i' with i+j=k, i'+j'=k (so j > j')
    # Each off-diagonal pair contributes:
    #   [a1_{i+1}*a2_{j+1}*b1_{i'}*b2_{j'} + a1_{i'+1}*a2_{j'+1}*b1_i*b2_j]
    # - [a1_{i'+1}*a2_{j'+1}*b1_{i'}*b2_{j'} ... ] wait, let me be more careful.
    #
    # d_k^prod = A_{k+1}*B_k - A_k*B_{k+1}
    # A_{k+1} = sum_{i+j=k+1} a1_i*a2_j
    # B_k = sum_{p+q=k} b1_p*b2_q
    # A_k = sum_{i+j=k} a1_i*a2_j
    # B_{k+1} = sum_{p+q=k+1} b1_p*b2_q
    #
    # d_k^prod = sum_{i+j=k+1} sum_{p+q=k} a1_i*a2_j*b1_p*b2_q
    #          - sum_{i+j=k} sum_{p+q=k+1} a1_i*a2_j*b1_p*b2_q
    #
    # Reindex first sum: let i'=i-1 so i'+j=k, and the sum becomes
    # sum_{i'+j=k} sum_{p+q=k} a1_{i'+1}*a2_j*b1_p*b2_q
    #
    # Hmm, this is messy. Let's just compute the off-diagonal as total - diag.
    offdiag = dk_direct - diag

    return {
        "dk_direct": dk_direct,
        "diag": diag,
        "offdiag": offdiag,
        "diag_pure": diag_parts_pure,
        "diag_cross1": diag_parts_cross1,
        "diag_cross2": diag_parts_cross2,
        "diag_check": diag_parts_pure + diag_parts_cross1 + diag_parts_cross2 == diag,
    }


# ---------------------------------------------------------------------------
# Factor pair extraction
# ---------------------------------------------------------------------------

def extract_support_factor_pairs(n: int, adj: list[list[int]]) -> list[dict]:
    """Extract (I_c, E_c) factor pairs from support vertices.

    For each support vertex r (has at least one leaf child), root the
    tree at r and collect the (I_c, E_c) pair from each non-leaf child c.

    Returns list of dicts with keys: I, E, subtree_size, tree_n.
    """
    deg = [len(adj[v]) for v in range(n)]
    leaves = {v for v in range(n) if deg[v] == 1}

    pairs = []
    seen_polys = set()

    for r in range(n):
        leaf_children_of_r = [v for v in adj[r] if v in leaves]
        if not leaf_children_of_r:
            continue

        dp0, dp1s, children, parent = dp_tree(n, adj, r)

        non_leaf_children = [c for c in children[r] if c not in leaves]
        for c in non_leaf_children:
            E_c = dp0[c]
            I_c = _polyadd(dp0[c], [0] + dp1s[c])
            key = (tuple(I_c), tuple(E_c))
            if key not in seen_polys:
                seen_polys.add(key)
                # Count subtree size
                sub_size = 0
                stk = [c]
                while stk:
                    v = stk.pop()
                    sub_size += 1
                    for ch in children[v]:
                        stk.append(ch)
                pairs.append({
                    "I": I_c,
                    "E": E_c,
                    "subtree_size": sub_size,
                    "tree_n": n,
                })

    return pairs


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def check_condC_all_k(I_poly, E_poly) -> list[dict]:
    """Check Condition C for a pair at all k, return per-k data."""
    max_k = max(len(I_poly), len(E_poly))
    results = []
    for k in range(1, max_k):
        T1, T2, T3, total = condition_C(I_poly, E_poly, k)
        dk = compute_d(I_poly, E_poly, k)
        ck = compute_c(E_poly, k)
        if T1 == 0 and T2 == 0 and T3 == 0:
            continue
        results.append({
            "k": k, "dk": dk, "ck": ck,
            "T1": T1, "T2": T2, "T3": T3, "total": total,
        })
    return results


def main():
    GENG = "/opt/homebrew/bin/geng"
    MAX_N_FACTORS = 13   # collect factors from trees up to this size
    MAX_PAIRS = 2000     # cap on unique factor pairs to keep runtime sane

    print("=" * 90)
    print("DIAGNOSTIC: Strong Condition C Product Closure Mechanism")
    print(f"Collecting factor pairs from support vertices of all trees n <= {MAX_N_FACTORS}")
    print(f"Cap on factor pairs: {MAX_PAIRS}")
    print("=" * 90)

    # Phase 1: collect unique factor pairs
    t0 = time.time()
    all_pairs: list[dict] = []
    seen_keys = set()
    trees_seen = 0

    for n in range(3, MAX_N_FACTORS + 1):
        proc = subprocess.Popen(
            [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"],
            stdout=subprocess.PIPE,
        )
        count_n = 0
        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            trees_seen += 1
            pairs = extract_support_factor_pairs(nn, adj)
            for p in pairs:
                key = (tuple(p["I"]), tuple(p["E"]))
                if key not in seen_keys:
                    seen_keys.add(key)
                    all_pairs.append(p)
            count_n += 1
        proc.wait()
        print(f"  n={n:2d}: {count_n:>6d} trees, {len(all_pairs):>6d} unique factor pairs so far",
              flush=True)

    elapsed_phase1 = time.time() - t0
    print(f"\nPhase 1 complete: {trees_seen:,} trees scanned, "
          f"{len(all_pairs):,} unique factor pairs ({elapsed_phase1:.1f}s)")

    # Sort by subtree size descending (larger subtrees = more interesting)
    # and keep only MAX_PAIRS
    all_pairs.sort(key=lambda p: -p["subtree_size"])
    if len(all_pairs) > MAX_PAIRS:
        print(f"Trimming to {MAX_PAIRS} largest-subtree factor pairs "
              f"(dropped {len(all_pairs) - MAX_PAIRS})")
        all_pairs = all_pairs[:MAX_PAIRS]

    # Phase 1.5: verify Condition C on individual factors
    print(f"\n{'-'*90}")
    print("Phase 1.5: Verify Condition C on individual factor pairs")
    factor_failures = 0
    factor_dk_neg_total = 0
    for idx, p in enumerate(all_pairs):
        per_k = check_condC_all_k(p["I"], p["E"])
        for entry in per_k:
            if entry["dk"] < 0:
                factor_dk_neg_total += 1
            if entry["total"] < 0:
                factor_failures += 1
                print(f"  FACTOR FAIL pair {idx}, k={entry['k']}: "
                      f"total={entry['total']}")
                break

    n_E_lc = sum(1 for p in all_pairs if is_log_concave(p["E"]))
    n_I_lc = sum(1 for p in all_pairs if is_log_concave(p["I"]))
    print(f"  Factor CondC failures: {factor_failures} / {len(all_pairs)}")
    print(f"  Factor (pair,k) with d_k < 0: {factor_dk_neg_total}")
    print(f"  E is LC: {n_E_lc}/{len(all_pairs)}, I is LC: {n_I_lc}/{len(all_pairs)}")

    # Phase 2: pairwise products
    n_pairs = len(all_pairs)
    n_product_pairs = n_pairs * (n_pairs + 1) // 2
    print(f"\n{'-'*90}")
    print(f"Phase 2: Pairwise products ({n_pairs} pairs -> "
          f"{n_product_pairs:,} products)")
    print("-" * 90)

    t_phase2 = time.time()

    # Accumulators
    total_condC_checks = 0
    n_dk_negative = 0
    n_condC_fail_pairs = 0

    # Compensation breakdown when d_k < 0
    compensation = defaultdict(int)

    # Tightest margins
    tightest_margins = []  # (margin_ratio, i, j, k, details)

    # d_k decomposition stats
    n_decomp_diag_neg = 0   # diagonal contribution to d_k is negative
    n_decomp_offdiag_neg = 0

    # Curvature vs memory ratio accumulator
    curv_fracs = []
    mem_fracs = []

    progress_interval = max(1, n_pairs // 20)

    for i in range(n_pairs):
        if i % progress_interval == 0:
            elapsed = time.time() - t_phase2
            pct = 100 * i / n_pairs
            print(f"  Progress: {pct:5.1f}% (i={i}/{n_pairs}, elapsed={elapsed:.1f}s)",
                  flush=True)

        I1, E1 = all_pairs[i]["I"], all_pairs[i]["E"]

        for j in range(i, n_pairs):
            I2, E2 = all_pairs[j]["I"], all_pairs[j]["E"]
            I_prod = _polymul(I1, I2)
            E_prod = _polymul(E1, E2)

            max_k = max(len(I_prod), len(E_prod))
            pair_fail = False
            pair_min_margin = None
            pair_worst_k = -1
            pair_dk_neg = 0

            for k in range(1, max_k):
                T1, T2, T3, total = condition_C(I_prod, E_prod, k)
                dk = compute_d(I_prod, E_prod, k)
                if T1 == 0 and T2 == 0 and T3 == 0:
                    continue
                total_condC_checks += 1

                if dk < 0:
                    n_dk_negative += 1
                    pair_dk_neg += 1

                    # Compensation
                    if total >= 0:
                        t2_pos = T2 > 0
                        t3_pos = T3 > 0
                        if t2_pos and t3_pos:
                            compensation["T2+T3"] += 1
                        elif t2_pos:
                            compensation["T2_only"] += 1
                        elif t3_pos:
                            compensation["T3_only"] += 1
                        else:
                            compensation["T1_nonneg"] += 1

                        # Curvature vs memory ratio
                        pos_total = max(T1, 0) + max(T2, 0) + max(T3, 0)
                        if pos_total > 0:
                            curv_fracs.append(max(T3, 0) / pos_total)
                            mem_fracs.append(max(T2, 0) / pos_total)
                    else:
                        compensation["FAIL"] += 1

                    # d_k decomposition
                    decomp = decompose_dk_product(I1, E1, I2, E2, k)
                    if decomp["diag"] < 0:
                        n_decomp_diag_neg += 1
                    if decomp["offdiag"] < 0:
                        n_decomp_offdiag_neg += 1

                # Margin tracking
                norm = abs(T1) + abs(T2) + abs(T3)
                if norm > 0:
                    margin = total / norm
                    if pair_min_margin is None or margin < pair_min_margin:
                        pair_min_margin = margin
                        pair_worst_k = k

                if total < 0:
                    pair_fail = True

            if pair_fail:
                n_condC_fail_pairs += 1

            if pair_dk_neg > 0 and pair_min_margin is not None:
                tightest_margins.append((
                    pair_min_margin, i, j, pair_worst_k,
                    {
                        "dk_neg_count": pair_dk_neg,
                        "sub1": all_pairs[i]["subtree_size"],
                        "sub2": all_pairs[j]["subtree_size"],
                    }
                ))

    tightest_margins.sort(key=lambda x: x[0])
    elapsed_phase2 = time.time() - t_phase2

    print(f"\nPhase 2 complete: {n_product_pairs:,} products, "
          f"{total_condC_checks:,} checks ({elapsed_phase2:.1f}s)\n")

    # ======================================================================
    # REPORT
    # ======================================================================
    print("=" * 90)
    print("RESULTS")
    print("=" * 90)

    print(f"\n1. PRODUCT-LEVEL CONDITION C")
    print(f"   Total product pairs:           {n_product_pairs:>12,}")
    print(f"   CondC FAIL pairs:              {n_condC_fail_pairs:>12,}")
    print(f"   Total (pair, k) with d_k < 0:  {n_dk_negative:>12,}")

    print(f"\n2. COMPENSATION WHEN d_k < 0")
    for label in ["T2_only", "T3_only", "T2+T3", "T1_nonneg", "FAIL"]:
        cnt = compensation.get(label, 0)
        pct = 100 * cnt / max(n_dk_negative, 1)
        desc = {
            "T2_only":   "Memory (b_k * d_{k-1}) alone",
            "T3_only":   "Curvature (a_{k-1} * c_k) alone",
            "T2+T3":     "Both memory + curvature",
            "T1_nonneg": "T1 non-negative (b_{k-1}*d_k >= 0 despite d_k < 0)",
            "FAIL":      "NOTHING (Condition C fails)",
        }[label]
        print(f"   {label:12s}: {cnt:>10,} ({pct:5.1f}%)  {desc}")

    t2_participates = compensation.get("T2_only", 0) + compensation.get("T2+T3", 0)
    t3_participates = compensation.get("T3_only", 0) + compensation.get("T2+T3", 0)
    print(f"\n   Memory participates:     {t2_participates:>10,} / {n_dk_negative}")
    print(f"   Curvature participates:  {t3_participates:>10,} / {n_dk_negative}")

    print(f"\n3. CURVATURE vs MEMORY DISTRIBUTION")
    if curv_fracs:
        curv_fracs.sort()
        mem_fracs.sort()
        nc = len(curv_fracs)
        print(f"   {nc} rescued (pair,k) events:")
        print(f"   Curvature fraction of positive load:")
        print(f"     min={curv_fracs[0]:.4f}  p10={curv_fracs[nc//10]:.4f}  "
              f"p25={curv_fracs[nc//4]:.4f}  median={curv_fracs[nc//2]:.4f}  "
              f"p75={curv_fracs[3*nc//4]:.4f}  max={curv_fracs[-1]:.4f}")
        print(f"   Memory fraction of positive load:")
        print(f"     min={mem_fracs[0]:.4f}  p10={mem_fracs[nc//10]:.4f}  "
              f"p25={mem_fracs[nc//4]:.4f}  median={mem_fracs[nc//2]:.4f}  "
              f"p75={mem_fracs[3*nc//4]:.4f}  max={mem_fracs[-1]:.4f}")
        curv_dom = sum(1 for c, m in zip(curv_fracs, mem_fracs) if c > m)
        mem_dom = sum(1 for c, m in zip(curv_fracs, mem_fracs) if m > c)
        print(f"   Curvature dominant: {curv_dom} ({100*curv_dom/nc:.1f}%)")
        print(f"   Memory dominant:    {mem_dom} ({100*mem_dom/nc:.1f}%)")
        print(f"   Equal:              {nc - curv_dom - mem_dom}")
    else:
        print("   No rescued events to analyze.")

    print(f"\n4. d_k PRODUCT DECOMPOSITION (diagonal vs off-diagonal)")
    print(f"   When d_k^prod < 0:")
    print(f"     diagonal < 0:    {n_decomp_diag_neg:>10,} / {n_dk_negative}")
    print(f"     off-diagonal < 0: {n_decomp_offdiag_neg:>10,} / {n_dk_negative}")

    print(f"\n5. TIGHTEST MARGINS")
    n_show = min(20, len(tightest_margins))
    if tightest_margins:
        print(f"   {'Rank':>4s}  {'Margin':>10s}  {'i':>4s}  {'j':>4s}  {'k':>3s}  "
              f"{'dk_neg':>6s}  {'Sub1':>4s}  {'Sub2':>4s}")
        for rank, (margin, ii, jj, kk, details) in enumerate(tightest_margins[:n_show]):
            print(f"   {rank+1:4d}  {margin:10.6f}  {ii:4d}  {jj:4d}  {kk:3d}  "
                  f"{details['dk_neg_count']:6d}  "
                  f"{details['sub1']:4d}  {details['sub2']:4d}")

        # Deep dive on the tightest
        print(f"\n   === Deep dive: tightest pair ===")
        margin, ii, jj, kk, _ = tightest_margins[0]
        I1, E1 = all_pairs[ii]["I"], all_pairs[ii]["E"]
        I2, E2 = all_pairs[jj]["I"], all_pairs[jj]["E"]
        I_prod = _polymul(I1, I2)
        E_prod = _polymul(E1, E2)

        print(f"   Factor 1 (subtree {all_pairs[ii]['subtree_size']}):")
        print(f"     I = {I1}")
        print(f"     E = {E1}")
        print(f"     LC(I)={is_log_concave(I1)}, LC(E)={is_log_concave(E1)}")
        print(f"   Factor 2 (subtree {all_pairs[jj]['subtree_size']}):")
        print(f"     I = {I2}")
        print(f"     E = {E2}")
        print(f"     LC(I)={is_log_concave(I2)}, LC(E)={is_log_concave(E2)}")
        print(f"   Product:")
        print(f"     I_prod = {I_prod}")
        print(f"     E_prod = {E_prod}")
        print(f"     LC(I_prod)={is_log_concave(I_prod)}, LC(E_prod)={is_log_concave(E_prod)}")
        print(f"     mode(I_prod) = {mode_index(I_prod)}")

        print(f"\n   Per-k breakdown for product (tightest pair):")
        print(f"   {'k':>3s}  {'dk':>14s}  {'ck':>14s}  {'T1':>14s}  "
              f"{'T2':>14s}  {'T3':>14s}  {'Total':>14s}  {'margin':>10s}")
        max_k = max(len(I_prod), len(E_prod))
        for k in range(1, max_k):
            T1, T2, T3, total = condition_C(I_prod, E_prod, k)
            dk = compute_d(I_prod, E_prod, k)
            ck = compute_c(E_prod, k)
            if T1 == 0 and T2 == 0 and T3 == 0:
                continue
            norm = abs(T1) + abs(T2) + abs(T3)
            mr = total / norm if norm > 0 else float("inf")
            flag = " <-- dk<0" if dk < 0 else ""
            print(f"   {k:3d}  {dk:14d}  {ck:14d}  {T1:14d}  "
                  f"{T2:14d}  {T3:14d}  {total:14d}  {mr:10.6f}{flag}")

        # d_k decomposition at tightest k
        print(f"\n   d_k decomposition at k={kk}:")
        decomp = decompose_dk_product(I1, E1, I2, E2, kk)
        print(f"     d_k^prod (direct) = {decomp['dk_direct']}")
        print(f"     diagonal          = {decomp['diag']}")
        print(f"       pure (Σ d1*d2)  = {decomp['diag_pure']}")
        print(f"       cross1 (Σ d1*s2)= {decomp['diag_cross1']}")
        print(f"       cross2 (Σ s1*d2)= {decomp['diag_cross2']}")
        print(f"       check diag sum  = {decomp['diag_check']}")
        print(f"     off-diagonal      = {decomp['offdiag']}")

        # Show factor-level d values
        print(f"\n   Factor-level d_k values:")
        for label, Ix, Ex in [("Factor 1", I1, E1), ("Factor 2", I2, E2)]:
            vals = []
            for ki in range(max(len(Ix), len(Ex)) + 1):
                dv = compute_d(Ix, Ex, ki)
                if dv != 0:
                    vals.append(f"d[{ki}]={dv}")
            print(f"     {label}: {', '.join(vals)}")

        # Deep dive on #2 and #3 if available
        for extra_rank in range(1, min(3, len(tightest_margins))):
            margin2, ii2, jj2, kk2, det2 = tightest_margins[extra_rank]
            I1b, E1b = all_pairs[ii2]["I"], all_pairs[ii2]["E"]
            I2b, E2b = all_pairs[jj2]["I"], all_pairs[jj2]["E"]
            I_prod2 = _polymul(I1b, I2b)
            E_prod2 = _polymul(E1b, E2b)
            T1r, T2r, T3r, totalr = condition_C(I_prod2, E_prod2, kk2)
            dkr = compute_d(I_prod2, E_prod2, kk2)
            print(f"\n   === Rank {extra_rank+1} (margin={margin2:.6f}, k={kk2}) ===")
            print(f"   Factor 1 (sub={all_pairs[ii2]['subtree_size']}): "
                  f"I={I1b[:6]}{'...' if len(I1b)>6 else ''}, "
                  f"E={E1b[:6]}{'...' if len(E1b)>6 else ''}")
            print(f"   Factor 2 (sub={all_pairs[jj2]['subtree_size']}): "
                  f"I={I2b[:6]}{'...' if len(I2b)>6 else ''}, "
                  f"E={E2b[:6]}{'...' if len(E2b)>6 else ''}")
            print(f"   At k={kk2}: dk={dkr}, T1={T1r}, T2={T2r}, T3={T3r}, total={totalr}")
            decomp2 = decompose_dk_product(I1b, E1b, I2b, E2b, kk2)
            print(f"   Decomp: diag={decomp2['diag']} (pure={decomp2['diag_pure']}, "
                  f"c1={decomp2['diag_cross1']}, c2={decomp2['diag_cross2']}), "
                  f"offdiag={decomp2['offdiag']}")

    else:
        print("   No d_k < 0 events in any product pair.")

    # Phase 3: factor-pair patterns for tightest product pairs
    print(f"\n6. FACTOR-LEVEL CONDITION C (full range for each factor)")
    for rank in range(min(5, len(tightest_margins))):
        margin, ii, jj, kk, _ = tightest_margins[rank]
        I1, E1 = all_pairs[ii]["I"], all_pairs[ii]["E"]
        I2, E2 = all_pairs[jj]["I"], all_pairs[jj]["E"]
        print(f"\n   Rank {rank+1} (margin={margin:.6f}, product k={kk}):")
        for label, Ix, Ex in [("Factor 1", I1, E1), ("Factor 2", I2, E2)]:
            max_fk = max(len(Ix), len(Ex))
            any_printed = False
            for k in range(1, max_fk):
                T1, T2, T3, total = condition_C(Ix, Ex, k)
                dk = compute_d(Ix, Ex, k)
                if T1 == 0 and T2 == 0 and T3 == 0:
                    continue
                any_printed = True
                flag = " *dk<0*" if dk < 0 else ""
                norm = abs(T1) + abs(T2) + abs(T3)
                mr = total / norm if norm > 0 else float("inf")
                print(f"     {label} k={k}: dk={dk:>10d}  T1={T1:>12d}  T2={T2:>12d}  "
                      f"T3={T3:>12d}  total={total:>12d}  m={mr:.4f}{flag}")
            if not any_printed:
                print(f"     {label}: (trivial, all terms zero)")

    # Phase 4: analyse WHERE d_k < 0 relative to modes
    print(f"\n7. WHERE d_k < 0 OCCURS IN THE PRODUCT")
    # For each product with d_k < 0, record k relative to mode of I_prod
    dk_neg_relative_mode = defaultdict(int)
    dk_neg_abs_position = defaultdict(int)
    sample_count = 0
    for i in range(min(500, n_pairs)):
        for j in range(i, min(500, n_pairs)):
            I1, E1 = all_pairs[i]["I"], all_pairs[i]["E"]
            I2, E2 = all_pairs[j]["I"], all_pairs[j]["E"]
            I_prod = _polymul(I1, I2)
            E_prod = _polymul(E1, E2)
            m = mode_index(I_prod)
            max_k = max(len(I_prod), len(E_prod))
            sample_count += 1
            for k in range(1, max_k):
                dk = compute_d(I_prod, E_prod, k)
                if dk < 0:
                    dk_neg_relative_mode[k - m] += 1
                    dk_neg_abs_position[k] += 1
    if dk_neg_relative_mode:
        print(f"   Distribution of (k - mode) when d_k < 0 ({sample_count} products sampled):")
        for offset in sorted(dk_neg_relative_mode.keys()):
            print(f"     k - mode = {offset:+3d}: {dk_neg_relative_mode[offset]:>8,}")
    if dk_neg_abs_position:
        print(f"   Distribution of absolute k when d_k < 0:")
        for k in sorted(dk_neg_abs_position.keys()):
            print(f"     k = {k:3d}: {dk_neg_abs_position[k]:>8,}")

    total_elapsed = time.time() - t0
    print(f"\n{'='*90}")
    print(f"TOTAL TIME: {total_elapsed:.1f}s")
    if n_condC_fail_pairs == 0:
        print("VERDICT: PASS -- Condition C closed under all tested products")
    else:
        print(f"VERDICT: FAIL -- {n_condC_fail_pairs} product pairs violate Condition C")
    print("=" * 90)


if __name__ == "__main__":
    main()
