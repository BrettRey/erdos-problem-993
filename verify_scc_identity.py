"""Verify the SCC identity and its implications.

The identity (from Condition C framework):
  b_{k-1} * Delta_k = b_{k-1} * d_k + b_k * d_{k-1} + a_{k-1} * c_k

where:
  a_k = I_k (IS poly coefficients)
  b_k = E_k (exclude-root poly)
  e_k = ((1+x)I)_k = I_k + I_{k-1}
  d_k = a_{k+1}*b_k - a_k*b_{k+1} = I_{k+1}*E_k - I_k*E_{k+1}  (LR minor of I vs E)
  c_k = b_k^2 - b_{k-1}*b_{k+1}  (LC gap of E)
  Delta_k = e_{k+1}*b_k - e_k*b_{k+1}  (SCC)

Note: d_k = E_k*J_k - E_{k+1}*J_{k-1}  (since I_k = E_k + J_{k-1}).
d_k CAN be negative (~14% of checks).

Key consequence: if c_k >= 0 (E is LC) AND d_k, d_{k-1} >= 0, then
ALL THREE terms are nonneg, so Delta_k >= 0 (SCC holds).
When d_k < 0, the curvature term a_{k-1}*c_k must compensate.
"""

import subprocess
import sys
import time

sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul, _polyadd

GENG = '/opt/homebrew/bin/geng'


def parse_g6(g6: str):
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
    dp0 = [None] * n
    dp1s = [None] * n
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            prod_I = [1]
            prod_E = [1]
            for c in children[v]:
                Ic = _polyadd(dp0[c], [0] + dp1s[c])
                prod_I = _polymul(prod_I, Ic)
                prod_E = _polymul(prod_E, dp0[c])
            dp0[v] = prod_I
            dp1s[v] = prod_E
    return dp0, dp1s, children


def coeff(poly, k):
    if k < 0 or poly is None or k >= len(poly):
        return 0
    return poly[k]


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    args = parser.parse_args()

    t0 = time.time()
    identity_failures = 0
    dk_neg = 0
    ck_neg = 0
    scc_neg = 0
    total_checks = 0
    total_trees = 0
    total_support_checks = 0

    term_stats = {
        'all_three_pos': 0,
        'dk_neg_only': 0,
        'dkm1_neg_only': 0,
        'both_d_neg': 0,
        'ck_neg': 0,
    }
    min_margin = None
    min_margin_info = None

    for nn in range(args.min_n, args.max_n + 1):
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0

        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            n, adj = parse_g6(line)
            n_trees += 1

            leaf_count = [0] * n
            for v in range(n):
                for u in adj[v]:
                    if len(adj[u]) == 1:
                        leaf_count[v] += 1

            for root in range(n):
                if leaf_count[root] == 0:
                    continue

                total_support_checks += 1
                dp0, dp1s, children = dp_rooted(n, adj, root)

                E = dp0[root]
                J = dp1s[root]
                I_poly = _polyadd(E, [0] + J)
                e = _polyadd(I_poly, [0] + I_poly)  # (1+x)I

                max_k = max(len(e), len(E))

                for k in range(max_k):
                    total_checks += 1
                    a_k = coeff(I_poly, k)
                    a_km1 = coeff(I_poly, k-1)
                    a_kp1 = coeff(I_poly, k+1)
                    b_k = coeff(E, k)
                    b_km1 = coeff(E, k-1)
                    b_kp1 = coeff(E, k+1)
                    e_k = coeff(e, k)
                    e_kp1 = coeff(e, k+1)

                    # SCC
                    delta_k = e_kp1 * b_k - e_k * b_kp1

                    # d_k = I_{k+1}*E_k - I_k*E_{k+1} (LR minor of I vs E)
                    d_k = a_kp1 * b_k - a_k * b_kp1
                    d_km1 = a_k * b_km1 - a_km1 * b_k

                    # c_k = b_k^2 - b_{k-1}*b_{k+1} (LC gap of E)
                    c_k = b_k * b_k - b_km1 * b_kp1

                    if c_k < 0:
                        ck_neg += 1

                    if d_k < 0:
                        dk_neg += 1

                    # Identity: b_{k-1} * delta_k = b_{k-1} * d_k + b_k * d_{k-1} + a_{k-1} * c_k
                    lhs = b_km1 * delta_k
                    rhs = b_km1 * d_k + b_k * d_km1 + a_km1 * c_k

                    if lhs != rhs:
                        identity_failures += 1
                        if identity_failures <= 3:
                            print(f"IDENTITY FAILURE at n={n}, root={root}, k={k}")
                            print(f"  LHS = {lhs}, RHS = {rhs}, diff = {lhs - rhs}")

                    if delta_k < 0:
                        scc_neg += 1

                    # Track term analysis
                    if b_km1 > 0:
                        if d_k < 0 and d_km1 < 0:
                            term_stats['both_d_neg'] += 1
                        elif d_k < 0:
                            term_stats['dk_neg_only'] += 1
                        elif d_km1 < 0:
                            term_stats['dkm1_neg_only'] += 1
                        else:
                            term_stats['all_three_pos'] += 1

                        if c_k < 0:
                            term_stats['ck_neg'] += 1

                        if min_margin is None or delta_k < min_margin:
                            min_margin = delta_k
                            min_margin_info = (n, root, k, delta_k, d_k, d_km1, c_k,
                                               a_km1, b_km1, b_k, line)

        proc.wait()
        total_trees += n_trees
        elapsed = time.time() - t0
        print(f"n={nn}: {n_trees} trees, {total_checks} checks, "
              f"id_fails={identity_failures}, dk<0={dk_neg}, ck<0={ck_neg}, "
              f"scc<0={scc_neg}, {elapsed:.1f}s", flush=True)

    print(f"\n{'='*60}")
    print(f"SUMMARY (n={args.min_n}..{args.max_n})")
    print(f"{'='*60}")
    print(f"Trees:            {total_trees}")
    print(f"Support vertices: {total_support_checks}")
    print(f"Total checks:     {total_checks}")
    print(f"Identity fails:   {identity_failures}")
    print(f"d_k < 0:          {dk_neg} ({100*dk_neg/max(1,total_checks):.2f}%)")
    print(f"c_k < 0 (LC gap): {ck_neg}")
    print(f"SCC < 0:          {scc_neg}")
    print(f"Min SCC margin:   {min_margin}")
    if min_margin_info:
        n, root, k, delta, dk, dkm, ck, akm, bkm, bk, g6 = min_margin_info
        print(f"  at n={n}, root={root}, k={k}")
        print(f"  d_k={dk}, d_{{k-1}}={dkm}, c_k={ck}")
        print(f"  a_{{k-1}}={akm}, b_{{k-1}}={bkm}, b_k={bk}")
    print(f"\nTerm analysis (b_km1 > 0):")
    for key, val in term_stats.items():
        print(f"  {key}: {val}")
    print(f"\nElapsed: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
