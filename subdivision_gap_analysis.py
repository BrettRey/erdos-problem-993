#!/usr/bin/env python3
"""Analyze the 0.1% gap cases where mode(A) = d+2.

For these cases, the sign pattern of Δ(I+A) is: +...+, ?₁, ?₂, -...-
A valley at d+1 requires ?₁ < 0 AND ?₂ > 0. We check:
1. What are the actual signs of ?₁ and ?₂?
2. How tight is the margin?
3. What structural features produce mode(A) = d+2?
4. Can we bound |ΔA_{d+1}| vs |ΔI_{d+1}|?

Also: extend to n=19 to get more gap cases and test tightness.
"""
import subprocess
import sys
import time
from collections import Counter, deque

from indpoly import (
    _polyadd,
    _polymul,
    independence_poly,
    is_log_concave,
)

GENG = "/opt/homebrew/bin/geng"


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode("ascii")]
    n = data[0]
    idx = 1
    adj = [[] for _ in range(n)]
    bit = 5
    word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5
                idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj


def first_descent(seq):
    for k in range(len(seq) - 1):
        if seq[k] > seq[k + 1]:
            return k
    return len(seq) - 1


def mode_index(seq):
    return max(range(len(seq)), key=lambda k: seq[k])


def split_at_edge(n, adj, u, v):
    A = set()
    queue = deque([u])
    A.add(u)
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y not in A and y != v:
                A.add(y)
                queue.append(y)
    return A, set(range(n)) - A


def rooted_is_poly(adj, vertices, root):
    """Compute (P, R) for forest on `vertices` rooted at `root`.
    P = IS poly with root included, R = IS poly with root excluded.
    """
    vset = set(vertices)
    parent = {root: -1}
    order = [root]
    queue = deque([root])
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y in vset and y not in parent:
                parent[y] = x
                queue.append(y)
                order.append(y)

    dp_in = {}
    dp_out = {}
    for v in reversed(order):
        children = [y for y in adj[v] if y in vset and parent.get(y) == v]
        if not children:
            dp_in[v] = [0, 1]
            dp_out[v] = [1]
        else:
            prod_out = [1]
            for c in children:
                prod_out = _polymul(prod_out, dp_out[c])
            dp_in[v] = [0] + prod_out

            prod_both = [1]
            for c in children:
                prod_both = _polymul(prod_both, _polyadd(dp_in[c], dp_out[c]))
            dp_out[v] = prod_both

    return dp_in[root], dp_out[root]


def compute_A_detailed(n, adj, u, v):
    """Compute A(x) and return (A, P_u, P_v, R_u, R_v, side sizes)."""
    side_u, side_v = split_at_edge(n, adj, u, v)
    P_u, R_u = rooted_is_poly(adj, side_u, u)
    P_v, R_v = rooted_is_poly(adj, side_v, v)

    Q_u = [0] + R_u
    Q_v = [0] + R_v
    term1 = _polymul(Q_u, Q_v)
    term2 = [0] + _polymul(P_u, P_v)
    A = _polyadd(term1, term2)
    return A, P_u, P_v, R_u, R_v, len(side_u), len(side_v)


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    print(f"Gap analysis: mode(A) = d+2 cases, n up to {max_n}", flush=True)
    print("=" * 78, flush=True)

    t0 = time.time()
    gap_cases = []
    total_edges = 0
    mode_dist = Counter()  # mode(A) - d distribution

    for n in range(4, max_n + 1):
        tn = time.time()
        n_edges = 0
        n_gap = 0

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)
            d = first_descent(I_T)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                n_edges += 1
                A, P_u, P_v, R_u, R_v, su, sv = compute_A_detailed(nn, adj, u, v)
                d_A = first_descent(A)
                diff = d_A - d
                mode_dist[diff] += 1

                if diff >= 2:
                    n_gap += 1
                    # Compute the key quantities at d and d+1
                    delta_I_d = I_T[d] - I_T[d + 1] if d + 1 < len(I_T) else I_T[d]
                    delta_A_d = (A[d + 1] if d + 1 < len(A) else 0) - (A[d] if d < len(A) else 0)

                    I_Tp = _polyadd(I_T, A)
                    delta_sum_d = I_Tp[d + 1] - I_Tp[d] if d + 1 < len(I_Tp) else -I_Tp[d]

                    # At d+1
                    delta_I_d1 = (I_T[d + 1] - I_T[d + 2]) if d + 2 < len(I_T) else I_T[d + 1]
                    delta_A_d1 = ((A[d + 2] if d + 2 < len(A) else 0) -
                                  (A[d + 1] if d + 1 < len(A) else 0))
                    delta_sum_d1 = (I_Tp[d + 2] - I_Tp[d + 1]) if d + 2 < len(I_Tp) else -I_Tp[d + 1]

                    # Ratio: how much of I's descent does A's ascent consume?
                    margin_d = delta_I_d + delta_A_d  # positive means A doesn't overcome I
                    margin_d1 = delta_I_d1 + delta_A_d1  # at d+1

                    # Actually we want: at position d, Δ(I+A)_d = ΔI_d + ΔA_d
                    # ΔI_d = I[d+1] - I[d] < 0 (I descending)
                    # ΔA_d = A[d+1] - A[d] (A still ascending since mode(A) >= d+2)
                    # For unimodality we need: the "?" entries don't create a valley

                    deg_u = len(adj[u])
                    deg_v = len(adj[v])
                    deg_seq = sorted([len(adj[x]) for x in range(nn)], reverse=True)

                    a_lc = is_log_concave(A)

                    # Compute LC ratio of A at d and d+1
                    a_d = A[d] if d < len(A) else 0
                    a_d1 = A[d + 1] if d + 1 < len(A) else 0
                    a_d2 = A[d + 2] if d + 2 < len(A) else 0
                    a_dm1 = A[d - 1] if d - 1 >= 0 and d - 1 < len(A) else 0

                    # LC at d: a_d^2 >= a_{d-1} * a_{d+1}
                    lc_ratio_d = (a_d ** 2) / (a_dm1 * a_d1) if a_dm1 > 0 and a_d1 > 0 else float('inf')
                    # LC at d+1: a_{d+1}^2 >= a_d * a_{d+2}
                    lc_ratio_d1 = (a_d1 ** 2) / (a_d * a_d2) if a_d > 0 and a_d2 > 0 else float('inf')

                    gap_cases.append({
                        'g6': g6, 'n': nn, 'u': u, 'v': v,
                        'd': d, 'd_A': d_A, 'diff': diff,
                        'deg_u': deg_u, 'deg_v': deg_v,
                        'su': su, 'sv': sv,
                        'deg_seq': deg_seq,
                        'delta_I_d': -delta_I_d,  # store as positive drop
                        'delta_A_d': delta_A_d,    # positive = A ascending
                        'delta_sum_d': delta_sum_d,
                        'delta_I_d1': -delta_I_d1,
                        'delta_A_d1': delta_A_d1,
                        'delta_sum_d1': delta_sum_d1,
                        'margin_d': -delta_sum_d,  # positive = I+A still descending
                        'margin_d1': -delta_sum_d1,
                        'A_lc': a_lc,
                        'lc_ratio_d': lc_ratio_d,
                        'lc_ratio_d1': lc_ratio_d1,
                        'I_d': I_T[d], 'I_d1': I_T[d + 1] if d + 1 < len(I_T) else 0,
                        'I_d2': I_T[d + 2] if d + 2 < len(I_T) else 0,
                        'A_d': a_d, 'A_d1': a_d1, 'A_d2': a_d2,
                        'A': A[:min(len(A), d + 5)],
                        'I': I_T[:min(len(I_T), d + 5)],
                    })

        proc.wait()
        elapsed = time.time() - tn
        total_edges += n_edges
        print(f"n={n:2d}: edges={n_edges:>10,}  gap_cases={n_gap:4d}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 78, flush=True)
    print(f"Total: {total_edges:,} edges, {len(gap_cases)} gap cases, {total_time:.1f}s", flush=True)
    print()

    # Mode distribution
    print("mode(A) - d(I) distribution:", flush=True)
    for diff in sorted(mode_dist.keys()):
        count = mode_dist[diff]
        pct = 100.0 * count / total_edges
        print(f"  +{diff}: {count:>10,} ({pct:.3f}%)", flush=True)
    print()

    if not gap_cases:
        print("No gap cases found.", flush=True)
        return

    # Analyze gap cases
    print(f"=== GAP CASE ANALYSIS ({len(gap_cases)} cases) ===", flush=True)
    print()

    # Key question: is ?₁ = Δ(I+A)_d ever negative when ?₂ = Δ(I+A)_{d+1} is positive?
    # That would create a valley.
    valley_risk = 0
    q1_neg = 0
    q2_pos = 0
    for c in gap_cases:
        if c['delta_sum_d'] < 0:
            q1_neg += 1
        if c['delta_sum_d1'] > 0:
            q2_pos += 1
        if c['delta_sum_d'] < 0 and c['delta_sum_d1'] > 0:
            valley_risk += 1

    print(f"?₁ = Δ(I+A)_d   negative: {q1_neg}/{len(gap_cases)}", flush=True)
    print(f"?₂ = Δ(I+A)_{{d+1}} positive: {q2_pos}/{len(gap_cases)}", flush=True)
    print(f"VALLEY RISK (?₁<0 and ?₂>0): {valley_risk}/{len(gap_cases)}", flush=True)
    print()

    # Sign patterns
    sign_patterns = Counter()
    for c in gap_cases:
        s1 = '+' if c['delta_sum_d'] > 0 else ('-' if c['delta_sum_d'] < 0 else '0')
        s2 = '+' if c['delta_sum_d1'] > 0 else ('-' if c['delta_sum_d1'] < 0 else '0')
        sign_patterns[(s1, s2)] += 1

    print("Sign patterns (?₁, ?₂):", flush=True)
    for (s1, s2), count in sign_patterns.most_common():
        print(f"  ({s1}, {s2}): {count}", flush=True)
    print()

    # Margin analysis: how tight is the tightest case?
    # margin_d = I's drop at d minus A's rise at d (positive = safe)
    # Actually margin_d = -(delta_sum_d) = I_drop - A_rise at position d
    min_margin_d = min(c['margin_d'] for c in gap_cases)
    min_margin_d1 = min(c['margin_d1'] for c in gap_cases)
    print(f"Min margin at d:   {min_margin_d} (positive = I+A nonincreasing)", flush=True)
    print(f"Min margin at d+1: {min_margin_d1}", flush=True)
    print()

    # Tightest cases
    tightest = sorted(gap_cases, key=lambda c: c['margin_d1'])[:10]
    print("10 tightest cases (by margin at d+1):", flush=True)
    for c in tightest:
        print(f"  n={c['n']} edge=({c['u']},{c['v']}) d={c['d']} d_A={c['d_A']} "
              f"deg=({c['deg_u']},{c['deg_v']}) sides=({c['su']},{c['sv']})", flush=True)
        print(f"    At d:   I_drop={c['delta_I_d']:6d}  A_rise={c['delta_A_d']:6d}  "
              f"margin={c['margin_d']:6d}", flush=True)
        print(f"    At d+1: I_drop={c['delta_I_d1']:6d}  A_rise={c['delta_A_d1']:6d}  "
              f"margin={c['margin_d1']:6d}", flush=True)
        print(f"    I[d..d+2] = {c['I_d']}, {c['I_d1']}, {c['I_d2']}", flush=True)
        print(f"    A[d..d+2] = {c['A_d']}, {c['A_d1']}, {c['A_d2']}", flush=True)
        print(f"    LC ratio at d: {c['lc_ratio_d']:.4f}  at d+1: {c['lc_ratio_d1']:.4f}",
              flush=True)
        print(f"    deg_seq: {c['deg_seq']}", flush=True)
        print()

    # Structural patterns
    print("Structural analysis:", flush=True)
    # What degrees do the subdivided edge endpoints have?
    deg_pairs = Counter()
    for c in gap_cases:
        du, dv = sorted([c['deg_u'], c['deg_v']])
        deg_pairs[(du, dv)] += 1
    print("  Edge endpoint degrees (sorted):", flush=True)
    for (du, dv), count in deg_pairs.most_common(10):
        print(f"    ({du},{dv}): {count}", flush=True)
    print()

    # Side size distribution
    print("  Side size balance (min, max):", flush=True)
    balance = Counter()
    for c in gap_cases:
        s_min, s_max = sorted([c['su'], c['sv']])
        balance[(s_min, s_max)] += 1
    for (sm, sx), count in balance.most_common(10):
        print(f"    ({sm},{sx}): {count}", flush=True)
    print()

    # Key ratios
    print("  A_rise / I_drop ratios at d+1:", flush=True)
    ratios = []
    for c in gap_cases:
        if c['delta_I_d1'] > 0:
            r = c['delta_A_d1'] / c['delta_I_d1']
            ratios.append(r)
    if ratios:
        ratios.sort()
        print(f"    max: {max(ratios):.6f}", flush=True)
        print(f"    p90: {ratios[int(0.9 * len(ratios))]:.6f}", flush=True)
        print(f"    p99: {ratios[int(0.99 * len(ratios))]:.6f}", flush=True)
        print(f"    median: {ratios[len(ratios) // 2]:.6f}", flush=True)
        print(f"    min: {min(ratios):.6f}", flush=True)

    # Check bound: A ≤ (1+x)I gives a_{d+1} ≤ i_{d+1} + i_d
    # and a_{d+2} ≤ i_{d+2} + i_{d+1}
    # Does A's LC (a_{d+1}^2 ≥ a_d * a_{d+2}) + these bounds close the gap?
    print()
    print("  Bound check: using A ≤ (1+x)I + LC of A:", flush=True)
    for c in tightest[:5]:
        # From A ≤ (1+x)I: a_{d+1} ≤ i_{d+1} + i_d
        bound_ad1 = c['I_d1'] + c['I_d']
        # Actual a_{d+1}
        actual_ad1 = c['A_d1']
        # From LC of A: a_{d+2} ≤ a_{d+1}^2 / a_d
        lc_bound_ad2 = c['A_d1'] ** 2 // c['A_d'] if c['A_d'] > 0 else 0
        # Actual a_{d+2}
        actual_ad2 = c['A_d2']
        # ΔA_{d+1} = a_{d+2} - a_{d+1}
        # Using LC bound: ΔA_{d+1} ≤ a_{d+1}(a_{d+1}/a_d - 1) = a_{d+1}(a_{d+1} - a_d)/a_d
        delta_A_bound = actual_ad1 * (actual_ad1 - c['A_d']) // c['A_d'] if c['A_d'] > 0 else 0
        # I's drop at d+1: ΔI_{d+1} = i_{d+1} - i_{d+2} (positive)
        i_drop = c['delta_I_d1']

        print(f"    n={c['n']} edge=({c['u']},{c['v']}):", flush=True)
        print(f"      a[d+1]={actual_ad1} ≤ {bound_ad1} = i[d+1]+i[d]", flush=True)
        print(f"      a[d+2]={actual_ad2} ≤ {lc_bound_ad2} (LC bound)", flush=True)
        print(f"      ΔA[d+1] actual={c['delta_A_d1']}, LC bound={delta_A_bound}", flush=True)
        print(f"      I drop at d+1={i_drop}", flush=True)
        print(f"      Safe? LC_bound < I_drop: {delta_A_bound < i_drop}", flush=True)


if __name__ == "__main__":
    main()
