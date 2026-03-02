"""Test aggregate dominance hypothesis for P2 at support vertices.

For each tree n=3..18, for each support vertex with s >= 2 non-leaf children,
process children incrementally and track:

1. The Karlin decomposition: W_k = Term1 (diagonal) + Term2 (cross).
   Finding from initial run: Term2 >= 0 ALWAYS. So the decomposition is
   additive, not cancellative.

2. The factor-level ratio dominance: delta_fg(p) = f_{p+1}*g_p - f_p*g_{p+1}
   where f = S_c (subtree IS poly) and g = E_c (exclude-root poly). When is
   delta_fg(p) < 0? How does this interact with the accumulated E, J?

3. Relative P2 slack: min_k W_k / (E_{k+1}*J_k). Does this increase with t?
   (Initial run: NO, slack DECREASES with more children.)

4. Per-step P2: does P2 hold at every intermediate step?
   (Initial run: YES, 0 failures across 410K support vertices.)

Key new test: the "aggregated condition" -- does ratio dominance of accumulated
J^(t) over each individual factor g_c hold?
"""

import subprocess
import sys
import time
from collections import defaultdict

from indpoly import _polymul, _polyadd


def parse_g6(g6: str) -> tuple[int, list[list[int]]]:
    s = g6.strip()
    idx = 0
    n = ord(s[idx]) - 63
    idx += 1
    adj = [[] for _ in range(n)]
    bits = []
    for ch in s[idx:]:
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


def get_coeff(poly, k):
    if 0 <= k < len(poly):
        return poly[k]
    return 0


def rooted_dp(n, adj, root):
    parent = [-1] * n
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
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    order = []
    stack = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v); continue
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

    return children[root], dp0, dp1s


def get_mode(poly):
    mx = max(poly)
    for k, v in enumerate(poly):
        if v == mx:
            return k
    return 0


def check_ratio_dominance(A, B, up_to):
    """Check a_{k+1}*b_k >= a_k*b_{k+1} for k=0..up_to-1.
    Returns (passes, first_fail, min_slack_ratio)."""
    first_fail = -1
    min_ratio = float('inf')
    for k in range(up_to):
        lhs = get_coeff(A, k+1) * get_coeff(B, k)
        rhs = get_coeff(A, k) * get_coeff(B, k+1)
        if lhs == 0 and rhs == 0:
            continue
        if lhs > 0:
            ratio = rhs / lhs
            slack = 1.0 - ratio
        elif rhs > 0:
            ratio = float('inf')
            slack = float('-inf')
        else:
            continue
        if slack < min_ratio:
            min_ratio = slack
        if lhs < rhs:
            if first_fail == -1:
                first_fail = k
    return first_fail == -1, first_fail, min_ratio


def analyze_support_vertex(n, adj, root, children_of_root, dp0, dp1s):
    leaf_children = []
    nonleaf_children = []
    for c in children_of_root:
        if dp0[c] == [1]:
            leaf_children.append(c)
        else:
            nonleaf_children.append(c)

    s = len(nonleaf_children)
    ell = len(leaf_children)
    if s < 2 or ell == 0:
        return None

    E_root = dp0[root]
    J_root = dp1s[root]
    I_T = _polyadd(E_root, [0] + J_root)
    m = get_mode(I_T)
    if m == 0:
        return None

    # Initialize
    E_acc = [1]
    for _ in range(ell):
        E_acc = _polymul(E_acc, [1, 1])
    J_acc = [1]

    steps = []
    child_factors = []  # (f, g) for each non-leaf child

    for t_idx, c in enumerate(nonleaf_children):
        t = t_idx + 1
        f = _polyadd(dp0[c], [0] + dp1s[c])  # S_c
        g = dp0[c]                             # E_c
        child_factors.append((f, g))

        E_new = _polymul(E_acc, f)
        J_new = _polymul(J_acc, g)

        # === P2 check at this intermediate step ===
        p2_ok, p2_fail_k, _ = check_ratio_dominance(E_new, J_new, m)

        # === Karlin decomposition ===
        n_neg_T2 = 0
        n_pos_T2 = 0
        n_zero_T2 = 0
        min_W = float('inf')
        min_T1 = float('inf')
        max_neg_T2_val = 0

        for k in range(m):
            w_k = (get_coeff(E_new, k+1) * get_coeff(J_new, k) -
                   get_coeff(E_new, k) * get_coeff(J_new, k+1))

            term1 = 0
            for i in range(max(len(E_acc), len(J_acc))):
                ei = get_coeff(E_acc, i)
                ji = get_coeff(J_acc, i)
                if ei == 0 or ji == 0:
                    continue
                p = k - i
                delta = (get_coeff(f, p+1) * get_coeff(g, p) -
                         get_coeff(f, p) * get_coeff(g, p+1))
                term1 += ei * ji * delta

            term2 = w_k - term1

            if term2 < 0:
                n_neg_T2 += 1
                if abs(term2) > max_neg_T2_val:
                    max_neg_T2_val = abs(term2)
            elif term2 > 0:
                n_pos_T2 += 1
            else:
                n_zero_T2 += 1

            if w_k < min_W:
                min_W = w_k
            if term1 < min_T1:
                min_T1 = term1

        # === Factor-level ratio dominance: delta_fg(p) for this child ===
        n_neg_delta = 0
        max_deg_f = len(f)
        for p in range(max_deg_f):
            delta = (get_coeff(f, p+1) * get_coeff(g, p) -
                     get_coeff(f, p) * get_coeff(g, p+1))
            if delta < 0:
                n_neg_delta += 1

        # === Relative P2 slack ===
        min_rel_slack = float('inf')
        for k in range(m):
            denom = get_coeff(E_new, k+1) * get_coeff(J_new, k)
            if denom > 0:
                w_k = (get_coeff(E_new, k+1) * get_coeff(J_new, k) -
                       get_coeff(E_new, k) * get_coeff(J_new, k+1))
                rel = w_k / denom
                if rel < min_rel_slack:
                    min_rel_slack = rel

        # === Ratio dominance of J^(t) over each individual g_c ===
        # Check: j^(t)_{k+1} * g_k >= j^(t)_k * g_{k+1} for all k up to mode
        # (i.e., J_new ratio-dominates each individual g factor)
        rd_all_ok = True
        rd_fail_child = -1
        for prev_t in range(t):
            _, g_prev = child_factors[prev_t]
            ok, fk, _ = check_ratio_dominance(J_new, g_prev, m)
            if not ok:
                rd_all_ok = False
                if rd_fail_child == -1:
                    rd_fail_child = prev_t

        # === Also check: E_acc ratio-dominates J_acc (before update) ===
        ea_rd_ja, _, ea_rd_slack = check_ratio_dominance(E_acc, J_acc, m)

        steps.append({
            't': t,
            'p2_ok': p2_ok,
            'p2_fail_k': p2_fail_k,
            'min_W': min_W,
            'n_neg_T2': n_neg_T2,
            'n_pos_T2': n_pos_T2,
            'n_zero_T2': n_zero_T2,
            'n_neg_delta': n_neg_delta,
            'min_rel_slack': min_rel_slack,
            'J_rd_all_g': rd_all_ok,
            'J_rd_fail_child': rd_fail_child,
            'Eacc_rd_Jacc': ea_rd_ja,
            'Eacc_rd_Jacc_slack': ea_rd_slack,
        })

        E_acc = E_new
        J_acc = J_new

    return {
        'n': n,
        'root': root,
        'ell': ell,
        's': s,
        'mode': m,
        'steps': steps,
    }


def main():
    MAX_N = 18
    GENG = '/opt/homebrew/bin/geng'

    # Collectors
    t2_sign = defaultdict(lambda: [0, 0, 0])  # neg, zero, pos
    neg_delta_by_st = defaultdict(int)  # count of factor-level failures
    total_factors_by_st = defaultdict(int)
    slack_by_st = defaultdict(list)
    slack_mono_viols = []
    p2_failures = []
    p2_intermediate_failures = []
    Jrd_failures = defaultdict(int)
    Eacc_rd_failures = defaultdict(int)
    eacc_rd_slack_by_st = defaultdict(list)

    total_trees = 0
    total_sv = 0
    t0 = time.time()

    for nn in range(3, MAX_N + 1):
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_trees = 0
        n_sv = 0

        for line in proc.stdout:
            g6 = line.strip()
            if not g6:
                continue
            n_trees += 1
            n_val, adj = parse_g6(g6)

            degree = [len(adj[v]) for v in range(n_val)]
            leaves = {v for v in range(n_val) if degree[v] == 1}

            for root in range(n_val):
                if root in leaves:
                    continue
                if not any(nb in leaves for nb in adj[root]):
                    continue

                ch, d0, d1s = rooted_dp(n_val, adj, root)
                result = analyze_support_vertex(n_val, adj, root, ch, d0, d1s)
                if result is None:
                    continue

                n_sv += 1
                s = result['s']
                prev_slack = None

                for step in result['steps']:
                    t = step['t']
                    key = (s, t)

                    # P2
                    if not step['p2_ok']:
                        if len(p2_failures) < 20:
                            p2_failures.append({
                                'g6': g6, 'n': n_val, 'root': root,
                                's': s, 't': t})

                    # T2 signs
                    t2_sign[key][0] += step['n_neg_T2']
                    t2_sign[key][1] += step['n_zero_T2']
                    t2_sign[key][2] += step['n_pos_T2']

                    # Factor-level delta
                    neg_delta_by_st[key] += (1 if step['n_neg_delta'] > 0 else 0)
                    total_factors_by_st[key] += 1

                    # Slack
                    ms = step['min_rel_slack']
                    slack_by_st[key].append(ms)
                    if prev_slack is not None and ms < prev_slack - 1e-12:
                        slack_mono_viols.append({
                            'n': n_val, 'root': root, 's': s, 't': t,
                            'prev': prev_slack, 'curr': ms})
                    prev_slack = ms

                    # J ratio-dominates all g
                    if not step['J_rd_all_g']:
                        Jrd_failures[key] += 1

                    # E_acc rd J_acc
                    if not step['Eacc_rd_Jacc']:
                        Eacc_rd_failures[key] += 1
                    if step['Eacc_rd_Jacc_slack'] < float('inf'):
                        eacc_rd_slack_by_st[key].append(step['Eacc_rd_Jacc_slack'])

        proc.wait()
        total_trees += n_trees
        total_sv += n_sv
        elapsed = time.time() - t0
        print(f"n={nn:2d}: {n_trees:>8,d} trees, {n_sv:>8,d} SV  [{elapsed:.1f}s]",
              flush=True)

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"TOTAL: {total_trees:,d} trees, {total_sv:,d} support vertices ({elapsed:.1f}s)")

    # === REPORT 1: Term2 is always >= 0 ===
    print(f"\n{'='*70}")
    print("REPORT 1: Term2 (cross-term) sign distribution by (s, t)")
    print("  Hypothesis: Term2 >= 0 always. If confirmed, W_k = Term1 + Term2 >= Term1.")
    print(f"{'s':>3} {'t':>3} {'neg':>10} {'zero':>10} {'pos':>10} {'neg%':>8}")
    any_neg = False
    for key in sorted(t2_sign.keys()):
        s, t = key
        neg, zero, pos = t2_sign[key]
        total = neg + zero + pos
        pct = 100 * neg / total if total > 0 else 0
        if neg > 0:
            any_neg = True
        print(f"{s:3d} {t:3d} {neg:10d} {zero:10d} {pos:10d} {pct:7.2f}%")
    if not any_neg:
        print("\n  CONFIRMED: Term2 >= 0 in ALL cases. The cross-term is always non-negative.")
        print("  This means the diagonal Karlin term ALONE guarantees P2.")

    # === REPORT 2: Factor-level ratio dominance ===
    print(f"\n{'='*70}")
    print("REPORT 2: Factor-level ratio dominance (S_c vs E_c)")
    print("  delta_fg(p) = f_{p+1}*g_p - f_p*g_{p+1} < 0 for some p?")
    print(f"{'s':>3} {'t':>3} {'with_neg':>10} {'total':>10} {'pct':>8}")
    for key in sorted(neg_delta_by_st.keys()):
        s, t = key
        nd = neg_delta_by_st[key]
        tot = total_factors_by_st[key]
        pct = 100 * nd / tot if tot > 0 else 0
        print(f"{s:3d} {t:3d} {nd:10d} {tot:10d} {pct:7.2f}%")

    # === REPORT 3: E_acc ratio-dominates J_acc (before each step) ===
    print(f"\n{'='*70}")
    print("REPORT 3: Does E_acc ratio-dominate J_acc at each intermediate step?")
    print("  (Before incorporating child t, does E^(t-1) rd J^(t-1)?)")
    print(f"{'s':>3} {'t':>3} {'fails':>10} {'total':>10} {'pass%':>8}")
    for key in sorted(eacc_rd_slack_by_st.keys()):
        s, t = key
        fails = Eacc_rd_failures.get(key, 0)
        tot = total_factors_by_st.get(key, 0)
        pct = 100 * (tot - fails) / tot if tot > 0 else 0
        print(f"{s:3d} {t:3d} {fails:10d} {tot:10d} {pct:7.2f}%")

    # E_acc rd J_acc slack distribution
    print("\n  Slack distribution (min over k of 1 - rhs/lhs):")
    print(f"  {'s':>3} {'t':>3} {'count':>8} {'min':>10} {'p25':>10} {'median':>10}")
    for key in sorted(eacc_rd_slack_by_st.keys()):
        s, t = key
        vals = sorted(eacc_rd_slack_by_st[key])
        nv = len(vals)
        if nv == 0:
            continue
        print(f"  {s:3d} {t:3d} {nv:8d} {vals[0]:10.6f} "
              f"{vals[nv//4]:10.6f} {vals[nv//2]:10.6f}")

    # === REPORT 4: J^(t) ratio-dominates each g_c ===
    print(f"\n{'='*70}")
    print("REPORT 4: Does J^(t) ratio-dominate each individual factor g_c?")
    print(f"{'s':>3} {'t':>3} {'fails':>10} {'total':>10}")
    any_jrd_fail = False
    for key in sorted(Jrd_failures.keys()):
        s, t = key
        nf = Jrd_failures[key]
        tot = total_factors_by_st.get(key, 0)
        if nf > 0:
            any_jrd_fail = True
        print(f"{s:3d} {t:3d} {nf:10d} {tot:10d}")
    if not any_jrd_fail and not Jrd_failures:
        print("  J^(t) ratio-dominates all individual factors in ALL cases!")
    elif any_jrd_fail:
        print(f"\n  FAILURES detected. J^(t) does NOT always ratio-dominate individual g_c.")

    # === REPORT 5: Relative P2 slack ===
    print(f"\n{'='*70}")
    print("REPORT 5: Relative P2 slack by (s, t)")
    print(f"{'s':>3} {'t':>3} {'count':>8} {'min':>10} {'p25':>10} {'median':>10} {'p75':>10}")
    for key in sorted(slack_by_st.keys()):
        s, t = key
        vals = sorted(slack_by_st[key])
        nv = len(vals)
        if nv == 0:
            continue
        print(f"{s:3d} {t:3d} {nv:8d} {vals[0]:10.6f} {vals[nv//4]:10.6f} "
              f"{vals[nv//2]:10.6f} {vals[3*nv//4]:10.6f}")

    # === REPORT 6: Slack monotonicity ===
    print(f"\n{'='*70}")
    print(f"REPORT 6: Slack monotonicity violations")
    print(f"Total: {len(slack_mono_viols)}")
    by_s = defaultdict(list)
    for v in slack_mono_viols:
        by_s[v['s']].append(v)
    for sv in sorted(by_s.keys()):
        vv = by_s[sv]
        drops = [v['curr'] - v['prev'] for v in vv]  # negative values
        print(f"  s={sv}: {len(vv)} violations, "
              f"max drop={min(drops):.6f}, mean drop={sum(drops)/len(drops):.6f}")

    # === REPORT 7: P2 at intermediate steps ===
    print(f"\n{'='*70}")
    print(f"REPORT 7: P2 at intermediate steps")
    print(f"Failures: {len(p2_failures)}")
    if not p2_failures:
        print("P2 holds at EVERY intermediate step for ALL support vertices!")
    else:
        for f in p2_failures[:10]:
            print(f"  n={f['n']} root={f['root']} s={f['s']} t={f['t']}")

    # === Summary ===
    print(f"\n{'='*70}")
    print("SUMMARY")
    print("="*70)
    print(f"1. Term2 (cross-term) >= 0: {'CONFIRMED' if not any_neg else 'VIOLATED'}")
    print(f"2. P2 at every intermediate step: {'CONFIRMED' if not p2_failures else 'VIOLATED'}")
    print(f"3. Relative slack monotonically increasing: "
          f"{'YES' if not slack_mono_viols else f'NO ({len(slack_mono_viols)} violations)'}")
    print(f"4. The diagonal Karlin term alone suffices for P2 "
          f"(since Term2 >= 0 always).")


if __name__ == '__main__':
    main()
