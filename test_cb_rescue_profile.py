#!/usr/bin/env python3
"""
Profile the relationship between D_k and X_k in the CB decomposition,
focusing on cases where X_k < 0.

For each tree n=3..17 (exhaustive via geng), at each support vertex,
at each step t >= 2, at each k:
  1. Compute D_k and X_k (as in test_cb_pairwise.py)
  2. When X_k < 0, record D_k, X_k, and the ratio D_k/|X_k|
  3. Group the cross terms by "gap" g = |j-i|:
       gap_sum(g) = sum_{j-i=g} [Delta_{i,j}(A,B)*P(k-i)*Q(k-j)
                                + Delta_{j,i}(A,B)*P(k-j)*Q(k-i)]
     where i,j range over the full CB support indices, including boundary
     index -1 from the shifted term in Delta.
     Check if gap_sum(g) >= 0 for each gap g.

Reports:
  - Distribution of D_k/|X_k| when X_k < 0
  - Number of X_k < 0 cases per n
  - Gap-sum analysis: for which gaps g does gap_sum(g) < 0 occur?
  - Is there a gap threshold G such that sum_{g=1..G} gap_sum(g) >= 0 always?

All arithmetic is exact (Python ints); ratios use Fraction for exactness.
"""

import subprocess
from fractions import Fraction
from collections import defaultdict


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
    MAX_N = 17

    # Per-n counters for X_k < 0
    xk_neg_count_by_n = defaultdict(int)
    total_xk_by_n = defaultdict(int)

    # Collection of (D_k, X_k, ratio) when X_k < 0
    rescue_ratios = []  # list of Fraction(D_k, |X_k|)
    rescue_Dk = []
    rescue_Xk = []

    # Gap-sum analysis
    # gap_neg_count[g] = number of times gap_sum(g) < 0
    # gap_total_count[g] = total checks for gap g
    gap_neg_count = defaultdict(int)
    gap_total_count = defaultdict(int)

    # For each (step, k) where X_k < 0, record max G such that
    # partial_sum(1..G) >= 0 for all G up to max_gap
    # Actually: record smallest G where cumulative becomes >= 0
    # and whether cumulative is always >= 0 (up to max gap)
    # We'll track: for each case, the list of cumulative gap sums,
    # and find threshold behavior.
    cumul_always_nonneg = 0  # count of X_k<0 cases where cumul gap sum >= 0 for all G
    cumul_nonneg_by_G = defaultdict(int)  # cumul_nonneg_by_G[G] = count where sum_{g=1..G} >= 0

    # Sanity: does sum_g gap_sum(g) reconstruct X_k exactly?
    gap_reconstruct_failures = 0
    first_gap_reconstruct_failure = None

    trees_checked = 0
    support_rootings = 0
    steps_checked = 0

    for n in range(3, MAX_N + 1):
        # Generate trees: connected graphs on n vertices with exactly n-1 edges
        cmd = ["/opt/homebrew/bin/geng", "-cq", str(n), "%d:%d" % (n - 1, n - 1)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for g6 in lines:
            adj = graph6_to_adj(g6)
            trees_checked += 1

            # Try each vertex as root; skip non-support vertices
            for root in range(n):
                # Support vertex: adjacent to at least one leaf
                is_support = any(len(adj[nb]) == 1 for nb in adj[root])
                if not is_support:
                    continue

                support_rootings += 1

                # Compute full tree DP
                I, E, J, children, parent = tree_dp(adj, root)

                # Identify leaf vs non-leaf children of root
                leaf_kids = []
                nonleaf_kids = []
                for c in children[root]:
                    if not children[c]:  # c is a leaf
                        leaf_kids.append(c)
                    else:
                        nonleaf_kids.append(c)

                ell = len(leaf_kids)

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
                    P_c = I[c]
                    Q_c = E[c]

                    e_new = poly_mul(e_acc, P_c)
                    j_new = poly_mul(j_acc, Q_c)

                    step_number = t_idx + 1

                    if step_number >= 2:
                        steps_checked += 1

                        A = e_acc
                        B = j_acc

                        max_new = max(len(e_new), len(j_new))

                        # Compute Lambda_k^{old}(A,B)
                        max_ab = max(len(A), len(B))
                        Lambda_old = []
                        for k in range(max_ab + 1):
                            val = coeff(A, k) * coeff(B, k) - coeff(A, k - 1) * coeff(B, k + 1)
                            Lambda_old.append(val)

                        # Compute Lambda_k^{new}
                        Lambda_new = []
                        for k in range(max_new + 1):
                            val = coeff(e_new, k) * coeff(j_new, k) - coeff(e_new, k - 1) * coeff(j_new, k + 1)
                            Lambda_new.append(val)

                        # Compute D_k
                        D = [0] * (max_new + 1)
                        for k in range(max_new + 1):
                            s = 0
                            for i in range(len(Lambda_old)):
                                s += Lambda_old[i] * coeff(P_c, k - i) * coeff(Q_c, k - i)
                            D[k] = s

                        # Index window for Delta_{i,j}(A,B):
                        # include boundary index -1 (from shifted CB term), and
                        # go up to max degree index where Delta can be nonzero.
                        min_idx = -1
                        max_idx = max(len(A), len(B))

                        # For each k, check X_k and compute gap sums
                        for k in range(max_new + 1):
                            xk = coeff(Lambda_new, k) - D[k]
                            total_xk_by_n[n] += 1

                            # Compute gap sums for this k
                            # gap_sums_k[g] = sum over all (i,j) with j-i=g of
                            #   Delta_{i,j}(A,B)*P(k-i)*Q(k-j) + Delta_{j,i}(A,B)*P(k-j)*Q(k-i)
                            gap_sums_k = defaultdict(int)
                            max_gap_k = 0

                            for i in range(min_idx, max_idx + 1):
                                for j in range(i + 1, max_idx + 1):
                                    g = j - i
                                    delta_ij = (coeff(A, i) * coeff(B, j)
                                                - coeff(A, i - 1) * coeff(B, j + 1))
                                    delta_ji = (coeff(A, j) * coeff(B, i)
                                                - coeff(A, j - 1) * coeff(B, i + 1))
                                    term = (delta_ij * coeff(P_c, k - i) * coeff(Q_c, k - j)
                                            + delta_ji * coeff(P_c, k - j) * coeff(Q_c, k - i))
                                    gap_sums_k[g] += term
                                    if g > max_gap_k:
                                        max_gap_k = g

                            gap_total = sum(gap_sums_k.values())
                            if gap_total != xk:
                                gap_reconstruct_failures += 1
                                if first_gap_reconstruct_failure is None:
                                    first_gap_reconstruct_failure = (
                                        n, root, step_number, k, xk, gap_total
                                    )

                            # Record gap-sum sign data (for all k, not just X_k<0)
                            for g in range(1, max_gap_k + 1):
                                gs = gap_sums_k.get(g, 0)
                                gap_total_count[g] += 1
                                if gs < 0:
                                    gap_neg_count[g] += 1

                            if xk < 0:
                                xk_neg_count_by_n[n] += 1
                                dk = D[k]
                                rescue_Dk.append(dk)
                                rescue_Xk.append(xk)
                                if xk != 0:
                                    rescue_ratios.append(Fraction(dk, abs(xk)))

                                # Cumulative gap sum analysis for this X_k<0 case
                                cumul = 0
                                all_nonneg = True
                                for g in range(1, max_gap_k + 1):
                                    cumul += gap_sums_k.get(g, 0)
                                    if cumul >= 0:
                                        cumul_nonneg_by_G[g] += 1
                                    else:
                                        all_nonneg = False
                                if all_nonneg and max_gap_k >= 1:
                                    cumul_always_nonneg += 1

                    # Update accumulators
                    e_acc = e_new
                    j_acc = j_new

        print("n=%2d done  (trees: %d, support rootings: %d, steps: %d, X_k<0 this n: %d)"
              % (n, trees_checked, support_rootings, steps_checked,
                 xk_neg_count_by_n.get(n, 0)),
              flush=True)

    # ---------- Report ----------
    total_xk_neg = sum(xk_neg_count_by_n.values())
    total_xk_all = sum(total_xk_by_n.values())

    print()
    print("=" * 72)
    print("CB RESCUE PROFILE  n = 3 .. %d" % MAX_N)
    print("=" * 72)
    print("Trees checked:       %12s" % format(trees_checked, ","))
    print("Support rootings:    %12s" % format(support_rootings, ","))
    print("Steps (t >= 2):      %12s" % format(steps_checked, ","))
    print("Total X_k checks:    %12s" % format(total_xk_all, ","))
    print("Total X_k < 0:       %12s" % format(total_xk_neg, ","))
    print()

    # --- X_k < 0 counts per n ---
    print("-" * 72)
    print("X_k < 0 counts per n:")
    print("-" * 72)
    print("  n   total_Xk    Xk<0    rate")
    for n in range(3, MAX_N + 1):
        tot = total_xk_by_n.get(n, 0)
        neg = xk_neg_count_by_n.get(n, 0)
        rate = ("%.6f" % (neg / tot)) if tot > 0 else "n/a"
        print("  %2d  %9s  %7s  %s" % (n, format(tot, ","), format(neg, ","), rate))
    print()

    # --- D_k / |X_k| distribution ---
    print("-" * 72)
    print("D_k / |X_k| distribution when X_k < 0  (%d cases)" % len(rescue_ratios))
    print("-" * 72)
    if rescue_ratios:
        sorted_ratios = sorted(rescue_ratios)
        N = len(sorted_ratios)

        def percentile(arr, p):
            idx = int(p * (len(arr) - 1))
            return arr[idx]

        stats = {
            "min": sorted_ratios[0],
            "p5": percentile(sorted_ratios, 0.05),
            "p25": percentile(sorted_ratios, 0.25),
            "median": percentile(sorted_ratios, 0.50),
            "p75": percentile(sorted_ratios, 0.75),
            "p95": percentile(sorted_ratios, 0.95),
            "max": sorted_ratios[-1],
        }
        for label, val in stats.items():
            print("  %-8s %s  (%.6f)" % (label, val, float(val)))

        # How many have ratio >= 1 (D_k >= |X_k|, i.e., D_k dominates the deficit)?
        ge1 = sum(1 for r in rescue_ratios if r >= 1)
        print()
        print("  D_k >= |X_k| (ratio >= 1): %d / %d = %.4f"
              % (ge1, N, ge1 / N))
        # How many have D_k >= 0?
        dk_nonneg = sum(1 for d in rescue_Dk if d >= 0)
        print("  D_k >= 0 when X_k < 0:     %d / %d = %.4f"
              % (dk_nonneg, total_xk_neg, dk_nonneg / total_xk_neg if total_xk_neg else 0))

        # |X_k| distribution
        sorted_abs_xk = sorted(abs(x) for x in rescue_Xk)
        print()
        print("  |X_k| distribution:")
        print("    min: %d" % sorted_abs_xk[0])
        print("    median: %d" % percentile(sorted_abs_xk, 0.50))
        print("    max: %d" % sorted_abs_xk[-1])
    else:
        print("  No X_k < 0 cases found.")
    print()

    # --- Reconstruction sanity ---
    print("-" * 72)
    print("Reconstruction sanity:  X_k == sum_g gap_sum(g)")
    print("-" * 72)
    print("  total checks: %12s" % format(total_xk_all, ","))
    print("  failures:     %12s" % format(gap_reconstruct_failures, ","))
    if first_gap_reconstruct_failure:
        n, root, step_number, k, xk, gap_total = first_gap_reconstruct_failure
        print("  first failure: n=%d, root=%d, step=%d, k=%d, X_k=%d, gap_sum=%d"
              % (n, root, step_number, k, xk, gap_total))
    print()

    # --- Gap-sum analysis ---
    print("-" * 72)
    print("Gap-sum analysis: gap_sum(g) sign across ALL (step, k) checks")
    print("-" * 72)
    max_gap_seen = max(gap_total_count.keys()) if gap_total_count else 0
    print("  gap g   total_checks   neg_count   neg_rate")
    for g in range(1, max_gap_seen + 1):
        tot = gap_total_count.get(g, 0)
        neg = gap_neg_count.get(g, 0)
        if tot > 0:
            print("  %5d   %12s   %9s   %.6f" % (g, format(tot, ","), format(neg, ","), neg / tot))
    print()

    # --- Gap threshold analysis for X_k < 0 cases ---
    print("-" * 72)
    print("Gap threshold analysis (X_k < 0 cases only): cumulative gap sum >= 0")
    print("-" * 72)
    print("Total X_k < 0 cases: %d" % total_xk_neg)
    print("Cases where cumul gap sum >= 0 for ALL g: %d" % cumul_always_nonneg)
    if total_xk_neg > 0:
        print()
        print("  Cumul sum_{g=1..G} gap_sum(g) >= 0:")
        print("  G     count_nonneg / total     rate")
        for G in range(1, max_gap_seen + 1):
            cnt = cumul_nonneg_by_G.get(G, 0)
            print("  %3d   %8d / %8d    %.6f" % (G, cnt, total_xk_neg, cnt / total_xk_neg))
    print()

    # --- Verification: X_k = sum of all gap_sums ---
    # (Sanity: X_k should equal sum over g of gap_sum(g) for that k.)
    # We already know this from the algebra, but let's flag it.
    print("-" * 72)
    print("SUMMARY")
    print("-" * 72)
    if total_xk_neg == 0:
        print("X_k >= 0 everywhere -- no rescue needed.")
    else:
        print("X_k < 0 in %d / %d checks (%.4f%%)."
              % (total_xk_neg, total_xk_all, 100.0 * total_xk_neg / total_xk_all))
        if rescue_ratios:
            print("D_k / |X_k| ratio: min = %.6f, median = %.6f"
                  % (float(stats["min"]), float(stats["median"])))
        if cumul_always_nonneg == total_xk_neg:
            print("Cumulative gap sum >= 0 for ALL gaps in EVERY X_k<0 case.")
            print("==> A gap-cumulative rescue is viable.")
        else:
            deficit = total_xk_neg - cumul_always_nonneg
            print("Cumulative gap sum fails to be nonneg in %d / %d X_k<0 cases."
                  % (deficit, total_xk_neg))


if __name__ == "__main__":
    main()
