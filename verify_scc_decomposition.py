"""Verify SCC decomposition via ★-identity adapted for (1+x)I.

Key observation: at incremental product stages,
    b_new = b_old * P                          (multiplicative)
    e_new = e_old * Q + x*(1+x)*b_old * R      (NOT multiplicative, but clean)

where e = (1+x)*I, b = E, P = I_c, Q = E_c, R = J_c, P = Q + xR.

This gives SCC Δ_k^new = e_new_{k+1}*b_new_k - e_new_k*b_new_{k+1} as a sum of:

    Δ_k(e*Q, b*Q)           : main TP2 term (≥ 0 by old SCC + Q PF2)
    + (e*Q, b*R) correction  : from P = Q + xR in b
    + (b*R, b*Q) correction  : from x*(1+x)*b*R in e
    + LC(b*R)                : from cross terms

Check which terms are nonneg and what margin structure emerges.

Usage:
    python3 verify_scc_decomposition.py --max-n 18 --workers 8
"""

import argparse
import subprocess
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul, _polyadd

GENG = '/opt/homebrew/bin/geng'


def parse_g6(g6: str):
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
    if poly is None or k < 0 or k >= len(poly):
        return 0
    return poly[k]


def poly_shift(p, s):
    """Return x^s * p as a list (prepend s zeros)."""
    return [0]*s + list(p)


def process_tree(g6):
    n, adj = parse_g6(g6)

    result = {
        'n': n,
        'total_checks': 0,
        'identity_fails': 0,
        'scc_fails': 0,
        # Signs of the 5 terms
        'term1_neg': 0,  # Δ_k(e*Q, b*Q) < 0
        'term2_neg': 0,  # α_{k+1}*γ_{k-1} - α_k*γ_k < 0
        'term3_neg': 0,  # γ_k*β_k - γ_{k-1}*β_{k+1} < 0
        'term4_neg': 0,  # γ_{k-1}*β_k - γ_{k-2}*β_{k+1} < 0
        'term5_neg': 0,  # γ_{k-1}^2 - γ_{k-2}*γ_k < 0 (LC of b*R)
        # Combined terms
        'terms23_neg': 0,
        'terms234_neg': 0,
        'terms2345_neg': 0,
        # Tightest
        'min_scc': float('inf'),
        'min_term1': float('inf'),
        'min_terms2345': float('inf'),
        # Ratio tracking when correction < 0
        'min_ratio_t1_over_corr': float('inf'),  # min term1/|corr| when corr < 0
        'tightest_ratios': [],  # (ratio, k, stage, nfac, term1, corr, g6)
    }

    leaf_count = [0] * n
    for v in range(n):
        for u in adj[v]:
            if len(adj[u]) == 1:
                leaf_count[v] += 1

    for root in range(n):
        if leaf_count[root] == 0:
            continue

        dp0, dp1s, children = dp_rooted(n, adj, root)
        non_leaf_children = [c for c in children[root] if len(adj[c]) > 1]
        s = len(non_leaf_children)
        if s == 0:
            continue

        factors = []
        for c in non_leaf_children:
            Ic = _polyadd(dp0[c], [0] + dp1s[c])
            Ec = dp0[c]
            Rc = dp1s[c]
            factors.append((Ic, Ec, Rc))

        # Accumulated pair
        E_acc = [1]  # b
        J_acc = [1]  # j

        for stage, (P, Q, R) in enumerate(factors, 1):
            # Old I = E + xJ, e = (1+x)*I, b = E
            I_old = _polyadd(E_acc, [0] + J_acc)
            e_old = _polyadd(I_old, [0] + I_old)  # (1+x)*I
            b_old = E_acc

            # New accumulated
            E_new = _polymul(E_acc, P)   # b_new
            J_new = _polymul(J_acc, Q)   # j_new
            I_new = _polyadd(E_new, [0] + J_new)  # a_new
            e_new = _polyadd(I_new, [0] + I_new)  # (1+x)*I_new

            # Convolutions for decomposition
            eQ = _polymul(e_old, Q)       # α = e*Q
            bQ = _polymul(b_old, Q)       # β = b*Q
            bR = _polymul(b_old, R)       # γ = b*R

            max_k = max(len(e_new), len(E_new)) - 1

            for k in range(max_k):
                # Direct SCC
                delta_direct = coeff(e_new, k+1) * coeff(E_new, k) - \
                               coeff(e_new, k) * coeff(E_new, k+1)

                # 5-term decomposition
                ak1 = coeff(eQ, k+1)   # α_{k+1}
                ak0 = coeff(eQ, k)     # α_k
                bk0 = coeff(bQ, k)     # β_k
                bk1 = coeff(bQ, k+1)   # β_{k+1}
                gk0 = coeff(bR, k)     # γ_k
                gkm1 = coeff(bR, k-1)  # γ_{k-1}
                gkm2 = coeff(bR, k-2)  # γ_{k-2}

                term1 = ak1*bk0 - ak0*bk1                  # Δ_k(e*Q, b*Q)
                term2 = ak1*gkm1 - ak0*gk0                 # e*Q vs b*R cross
                term3 = gk0*bk0 - gkm1*bk1                 # b*R vs b*Q cross
                term4 = gkm1*bk0 - gkm2*bk1                # shifted b*R vs b*Q
                term5 = gkm1*gkm1 - gkm2*gk0               # LC(b*R)

                decomp = term1 + term2 + term3 + term4 + term5

                if delta_direct == 0 and decomp == 0:
                    continue

                result['total_checks'] += 1

                # Identity check
                if delta_direct != decomp:
                    result['identity_fails'] += 1

                if delta_direct < 0:
                    result['scc_fails'] += 1

                if delta_direct < result['min_scc']:
                    result['min_scc'] = delta_direct

                # Term signs
                if term1 < 0:
                    result['term1_neg'] += 1
                if term2 < 0:
                    result['term2_neg'] += 1
                if term3 < 0:
                    result['term3_neg'] += 1
                if term4 < 0:
                    result['term4_neg'] += 1
                if term5 < 0:
                    result['term5_neg'] += 1

                if term2 + term3 < 0:
                    result['terms23_neg'] += 1
                if term2 + term3 + term4 < 0:
                    result['terms234_neg'] += 1

                correction = term2 + term3 + term4 + term5
                if correction < 0 and term1 > 0:
                    ratio = term1 / abs(correction)
                    if ratio < result['min_ratio_t1_over_corr']:
                        result['min_ratio_t1_over_corr'] = ratio
                    result['tightest_ratios'].append(
                        (ratio, k, stage, s, int(term1), int(correction), g6.strip())
                    )
                    result['tightest_ratios'].sort(key=lambda x: x[0])
                    result['tightest_ratios'] = result['tightest_ratios'][:5]
                if correction < 0:
                    result['terms2345_neg'] += 1
                if correction < result['min_terms2345']:
                    result['min_terms2345'] = correction

                if term1 < result['min_term1']:
                    result['min_term1'] = term1

            # Update accumulated
            E_acc = E_new
            J_acc = J_new

    return result


def process_batch(g6_lines):
    keys = ['trees', 'total_checks', 'identity_fails', 'scc_fails',
            'term1_neg', 'term2_neg', 'term3_neg', 'term4_neg', 'term5_neg',
            'terms23_neg', 'terms234_neg', 'terms2345_neg']
    batch = {k: 0 for k in keys}
    batch['min_scc'] = float('inf')
    batch['min_term1'] = float('inf')
    batch['min_terms2345'] = float('inf')
    batch['min_ratio'] = float('inf')
    batch['tightest_ratios'] = []

    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        for k in keys:
            if k != 'trees':
                batch[k] += r.get(k, 0)
        if r['min_scc'] < batch['min_scc']:
            batch['min_scc'] = r['min_scc']
        if r['min_term1'] < batch['min_term1']:
            batch['min_term1'] = r['min_term1']
        if r['min_terms2345'] < batch['min_terms2345']:
            batch['min_terms2345'] = r['min_terms2345']
        if r['min_ratio_t1_over_corr'] < batch['min_ratio']:
            batch['min_ratio'] = r['min_ratio_t1_over_corr']
        batch['tightest_ratios'].extend(r['tightest_ratios'])

    batch['tightest_ratios'].sort(key=lambda x: x[0])
    batch['tightest_ratios'] = batch['tightest_ratios'][:20]
    return batch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=500)
    args = parser.parse_args()

    totals = {
        'trees': 0, 'total_checks': 0, 'identity_fails': 0, 'scc_fails': 0,
        'term1_neg': 0, 'term2_neg': 0, 'term3_neg': 0, 'term4_neg': 0, 'term5_neg': 0,
        'terms23_neg': 0, 'terms234_neg': 0, 'terms2345_neg': 0,
    }
    g_min_scc = float('inf')
    g_min_term1 = float('inf')
    g_min_2345 = float('inf')
    g_min_ratio = float('inf')
    g_tightest = []

    t0 = time.time()

    for nn in range(args.min_n, args.max_n + 1):
        tn = time.time()
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_totals = {k: 0 for k in totals}
        n_min_scc = float('inf')

        if args.workers <= 1:
            for line in proc.stdout:
                line = line.strip()
                if not line:
                    continue
                r = process_tree(line)
                n_totals['trees'] += 1
                for k in totals:
                    if k != 'trees':
                        n_totals[k] += r.get(k, 0)
                if r['min_scc'] < n_min_scc:
                    n_min_scc = r['min_scc']
                if r['min_ratio_t1_over_corr'] < g_min_ratio:
                    g_min_ratio = r['min_ratio_t1_over_corr']
                g_tightest.extend(r['tightest_ratios'])
        else:
            batch = []
            with Pool(args.workers) as pool:
                futures = []
                for line in proc.stdout:
                    line = line.strip()
                    if not line:
                        continue
                    batch.append(line)
                    if len(batch) >= args.batch_size:
                        futures.append(pool.apply_async(process_batch, (batch,)))
                        batch = []
                if batch:
                    futures.append(pool.apply_async(process_batch, (batch,)))
                for fut in futures:
                    br = fut.get()
                    n_totals['trees'] += br['trees']
                    for k in totals:
                        if k != 'trees':
                            n_totals[k] += br.get(k, 0)
                    if br['min_scc'] < n_min_scc:
                        n_min_scc = br['min_scc']
                    if br['min_term1'] < g_min_term1:
                        g_min_term1 = br['min_term1']
                    if br['min_terms2345'] < g_min_2345:
                        g_min_2345 = br['min_terms2345']
                    if br.get('min_ratio', float('inf')) < g_min_ratio:
                        g_min_ratio = br['min_ratio']
                    g_tightest.extend(br.get('tightest_ratios', []))

        proc.wait()
        elapsed_n = time.time() - tn
        for k in totals:
            totals[k] += n_totals[k]
        if n_min_scc < g_min_scc:
            g_min_scc = n_min_scc

        print(f"n={nn:3d}: {n_totals['trees']:>10,} trees | "
              f"id:{n_totals['identity_fails']}  scc:{n_totals['scc_fails']}  "
              f"t1<0:{n_totals['term1_neg']:>7,}  t2<0:{n_totals['term2_neg']:>7,}  "
              f"t3<0:{n_totals['term3_neg']:>7,}  t4<0:{n_totals['term4_neg']:>7,}  "
              f"t5<0:{n_totals['term5_neg']:>7,}  "
              f"corr<0:{n_totals['terms2345_neg']:>7,}  [{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*90}")
    print(f"SCC 5-TERM DECOMPOSITION")
    print(f"{'='*90}")
    print(f"Range:                    n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:              {totals['trees']:>15,}")
    print(f"Total checks:             {totals['total_checks']:>15,}")
    print(f"Identity failures:        {totals['identity_fails']:>15,}")
    print(f"SCC failures:             {totals['scc_fails']:>15,}")
    print()
    print(f"Term sign counts (negative):")
    print(f"  Term 1 (Δ(e*Q,b*Q)):   {totals['term1_neg']:>15,}")
    print(f"  Term 2 (e*Q vs b*R):    {totals['term2_neg']:>15,}")
    print(f"  Term 3 (b*R vs b*Q):    {totals['term3_neg']:>15,}")
    print(f"  Term 4 (shifted):       {totals['term4_neg']:>15,}")
    print(f"  Term 5 (LC(b*R)):       {totals['term5_neg']:>15,}")
    print()
    print(f"  Terms 2+3 < 0:          {totals['terms23_neg']:>15,}")
    print(f"  Terms 2+3+4 < 0:        {totals['terms234_neg']:>15,}")
    print(f"  Terms 2+3+4+5 < 0:      {totals['terms2345_neg']:>15,}")
    print()
    print(f"Minimum values:")
    print(f"  Min SCC (Δ_k):          {g_min_scc}")
    print(f"  Min Term 1:             {g_min_term1}")
    print(f"  Min correction (2345):  {g_min_2345}")
    print(f"  Min Term1/|corr| ratio: {g_min_ratio:.6f}")

    g_tightest.sort(key=lambda x: x[0])
    g_tightest = g_tightest[:30]
    if g_tightest:
        print(f"\n{'='*90}")
        print(f"TIGHTEST Term1/|correction| RATIOS (when correction < 0)")
        print(f"{'='*90}")
        print(f"{'Rank':>4s}  {'Ratio':>10s}  {'k':>3s}  {'Stage':>5s}  {'#fac':>4s}  "
              f"{'Term1':>14s}  {'Corr':>14s}  g6")
        for rank, (ratio, k, stage, nfac, t1, corr, g6) in enumerate(g_tightest):
            print(f"{rank+1:4d}  {ratio:10.6f}  {k:3d}  {stage:5d}  {nfac:4d}  "
                  f"{t1:14d}  {corr:14d}  {g6}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
