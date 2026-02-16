#!/usr/bin/env python3
"""Comprehensive subdivision analysis with CORRECT formula.

CORRECT: A(x) = P_u * P_v + x * R_u * R_v
where P_u = IS poly of u-side with u INCLUDED,
      R_u = IS poly of u-side with u EXCLUDED.

Previous code had A = x^2 R_u R_v + x P_u P_v (wrong: extra factor of x).

Checks:
  1. A(x) log-concave?
  2. A(x) unimodal?
  3. I(T') = I(T) + A unimodal?
  4. I(T') log-concave?
  5. first_descent(A) vs first_descent(I) (mode gap)
  6. Combined tail: I(T') nonincreasing from d(I)+1?
  7. B = (1+x)I - A unimodal?
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
    is_unimodal,
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


def compute_A(n, adj, u, v):
    """Compute CORRECT A(x) = P_u*P_v + x*R_u*R_v."""
    side_u, side_v = split_at_edge(n, adj, u, v)
    P_u, R_u = rooted_is_poly(adj, side_u, u)
    P_v, R_v = rooted_is_poly(adj, side_v, v)
    pp = _polymul(P_u, P_v)
    rr = _polymul(R_u, R_v)
    A = _polyadd(pp, [0] + rr)  # P_u*P_v + x*R_u*R_v
    return A


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    print(f"CORRECTED subdivision analysis, n up to {max_n}", flush=True)
    print("Formula: A = P_u*P_v + x*R_u*R_v", flush=True)
    print("=" * 80, flush=True)

    t0 = time.time()
    grand_totals = {
        'edges': 0,
        'A_lc_fail': 0, 'A_uni_fail': 0,
        'Ip_lc_fail': 0, 'Ip_uni_fail': 0,
        'B_uni_fail': 0, 'B_lc_fail': 0,
        'combined_tail_fail': 0,
        'descent_earlier': 0, 'descent_same': 0, 'descent_later': 0,
    }
    mode_dist = Counter()

    # Track tightest tail cases
    worst_tail_cases = []

    for n in range(3, max_n + 1):
        tn = time.time()
        counts = {k: 0 for k in grand_totals}
        counts['edges'] = 0
        n_mode_dist = Counter()

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)
            d = first_descent(I_T)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                counts['edges'] += 1

                A = compute_A(nn, adj, u, v)
                I_Tp = _polyadd(I_T, A)
                d_A = first_descent(A)
                d_Tp = first_descent(I_Tp)

                diff = d_A - d
                n_mode_dist[diff] += 1

                # Check A properties
                if not is_log_concave(A):
                    counts['A_lc_fail'] += 1
                if not is_unimodal(A):
                    counts['A_uni_fail'] += 1

                # Check I(T') properties
                if not is_log_concave(I_Tp):
                    counts['Ip_lc_fail'] += 1
                if not is_unimodal(I_Tp):
                    counts['Ip_uni_fail'] += 1

                # Check descent monotonicity
                if d_Tp < d:
                    counts['descent_earlier'] += 1
                elif d_Tp == d:
                    counts['descent_same'] += 1
                else:
                    counts['descent_later'] += 1

                # Check combined tail: I(T') nonincreasing from d+1
                for k in range(d + 1, len(I_Tp) - 1):
                    if I_Tp[k + 1] > I_Tp[k]:
                        counts['combined_tail_fail'] += 1
                        break

                # Check B = (1+x)I - A
                I_shifted = [0] + I_T
                opxI = _polyadd(I_T, I_shifted)
                B = [0] * max(len(opxI), len(A))
                for k in range(len(opxI)):
                    B[k] += opxI[k]
                for k in range(len(A)):
                    B[k] -= A[k]
                while len(B) > 1 and B[-1] == 0:
                    B.pop()

                if any(b < 0 for b in B):
                    pass  # A > (1+x)I would be a problem
                else:
                    if not is_unimodal(B):
                        counts['B_uni_fail'] += 1
                    if not is_log_concave(B):
                        counts['B_lc_fail'] += 1

                # Track tightest combined tail
                if diff >= 1 and d + 1 < len(I_Tp) - 1:
                    for k in range(d, min(d + diff + 1, len(I_Tp) - 1)):
                        delta = I_Tp[k + 1] - I_Tp[k]
                        if delta < 0 and (not worst_tail_cases or delta > worst_tail_cases[-1][0]):
                            pass  # only track risky ones
                    # Track sign pattern in uncertain zone
                    signs = []
                    for k in range(d, min(d + diff, len(I_Tp) - 1)):
                        delta = I_Tp[k + 1] - I_Tp[k]
                        signs.append('+' if delta > 0 else ('-' if delta < 0 else '0'))

        proc.wait()
        elapsed = time.time() - tn
        mode_dist.update(n_mode_dist)
        for k in grand_totals:
            grand_totals[k] += counts[k]

        print(f"n={n:2d} | edges={counts['edges']:>10,} "
              f"| A_lc={counts['A_lc_fail']} A_uni={counts['A_uni_fail']} "
              f"| I'_lc={counts['Ip_lc_fail']} I'_uni={counts['Ip_uni_fail']} "
              f"| tail={counts['combined_tail_fail']} "
              f"| d'<d={counts['descent_earlier']} "
              f"| B_uni={counts['B_uni_fail']} B_lc={counts['B_lc_fail']} "
              f"| {elapsed:.1f}s",
              flush=True)

    total_time = time.time() - t0
    print("=" * 80, flush=True)
    e = grand_totals['edges']
    print(f"Total: {e:,} edge subdivisions in {total_time:.1f}s", flush=True)
    print()
    print("Results:", flush=True)
    print(f"  A(x) LC failures:      {grand_totals['A_lc_fail']}", flush=True)
    print(f"  A(x) unimodal failures: {grand_totals['A_uni_fail']}", flush=True)
    print(f"  I(T') LC failures:     {grand_totals['Ip_lc_fail']}", flush=True)
    print(f"  I(T') unimodal failures: {grand_totals['Ip_uni_fail']}", flush=True)
    print(f"  Combined tail failures: {grand_totals['combined_tail_fail']}", flush=True)
    print(f"  d(I') < d(I):          {grand_totals['descent_earlier']}", flush=True)
    print(f"  d(I') = d(I):          {grand_totals['descent_same']}", flush=True)
    print(f"  d(I') > d(I):          {grand_totals['descent_later']}", flush=True)
    print(f"  B=(1+x)I-A uni fail:   {grand_totals['B_uni_fail']}", flush=True)
    print(f"  B=(1+x)I-A LC fail:    {grand_totals['B_lc_fail']}", flush=True)
    print()
    print("mode(A) - d(I) distribution:", flush=True)
    for diff in sorted(mode_dist.keys()):
        count = mode_dist[diff]
        pct = 100.0 * count / e if e > 0 else 0
        label = f"+{diff}" if diff >= 0 else str(diff)
        print(f"  {label}: {count:>10,} ({pct:.3f}%)", flush=True)


if __name__ == "__main__":
    main()
