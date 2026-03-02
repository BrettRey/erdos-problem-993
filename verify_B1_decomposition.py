"""Decompose B1 = Δ(e_new, E_old·E_c) into two TP2-closure sub-terms.

B1 = Δ_k(E_old·(1+x)I_c, E_old·E_c) + Δ_k(x(1+x)·J_old·E_c, E_old·E_c)
   = [TP2 of SCC_factor]           + [TP2 of C2_prev]

Term1 ≥ 0 by Karlin's theorem (E_old PF2, SCC of factor ≥ 0).
Term2 sign: same direction as C2 at previous stage (E_c PF2 preserves sign).

This script verifies these claims and profiles the ratio structure.

Usage:
    python3 verify_B1_decomposition.py --max-n 18 --workers 8
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


def process_tree(g6):
    n, adj = parse_g6(g6)

    result = {
        'n': n,
        'total_checks': 0,
        'identity_fails': 0,
        'B1_fails': 0,
        'T1_neg': 0,     # sub-term 1 (TP2 of SCC_factor) < 0
        'T2_neg': 0,     # sub-term 2 (TP2 of C2_prev) < 0
        'C2_prev_neg': 0,  # C2 at previous stage < 0 (for sign comparison)
        'sign_match': 0,   # T2 and C2_prev have same sign
        'sign_mismatch': 0,
        'min_T1': float('inf'),
        'min_ratio_T1_T2': float('inf'),  # min T1/|T2| when T2 < 0
        'tightest': [],
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

        E_acc = [1]
        J_acc = [1]

        for stage, (P, Q, R) in enumerate(factors, 1):
            # Previous stage C2: Δ_k(x(1+x)·J_acc, E_acc)
            I_old = _polyadd(E_acc, [0] + J_acc)
            e_old = _polyadd(I_old, [0] + I_old)

            # New accumulated
            E_new = _polymul(E_acc, P)
            J_new = _polymul(J_acc, Q)
            I_new = _polyadd(E_new, [0] + J_new)
            e_new = _polyadd(I_new, [0] + I_new)

            # B1 = Δ_k(e_new, E_acc·Q)
            EoldQ = _polymul(E_acc, Q)

            # Sub-term 1: Δ_k(E_old·(1+x)I_c, E_old·Q)
            # = Δ_k(E_old·(1+x)·P, E_old·Q) since I_c = P
            oneplus_P = _polyadd(P, [0] + P)  # (1+x)·P
            Eold_1xP = _polymul(E_acc, oneplus_P)  # E_old·(1+x)P

            # Sub-term 2: Δ_k(x(1+x)·J_old·E_c, E_old·E_c)
            # = Δ_k(x(1+x)·J_acc·Q, E_acc·Q)
            x1xJ = _polyadd([0] + J_acc, [0, 0] + J_acc)  # x(1+x)·J_acc
            x1xJQ = _polymul(x1xJ, Q)

            # C2 at previous stage: Δ_k(x(1+x)·J_acc, E_acc)
            x1xJ_prev = _polyadd([0] + J_acc, [0, 0] + J_acc)

            max_k = max(len(e_new), len(EoldQ)) - 1

            for k in range(max_k):
                # Direct B1
                B1_direct = coeff(e_new, k+1) * coeff(EoldQ, k) - \
                            coeff(e_new, k) * coeff(EoldQ, k+1)

                # Sub-term 1
                T1 = coeff(Eold_1xP, k+1) * coeff(EoldQ, k) - \
                     coeff(Eold_1xP, k) * coeff(EoldQ, k+1)

                # Sub-term 2
                T2 = coeff(x1xJQ, k+1) * coeff(EoldQ, k) - \
                     coeff(x1xJQ, k) * coeff(EoldQ, k+1)

                # C2 at previous stage
                C2_prev = coeff(x1xJ_prev, k+1) * coeff(E_acc, k) - \
                          coeff(x1xJ_prev, k) * coeff(E_acc, k+1)

                if B1_direct == 0 and T1 == 0 and T2 == 0:
                    continue

                result['total_checks'] += 1

                if T1 + T2 != B1_direct:
                    result['identity_fails'] += 1

                if B1_direct < 0:
                    result['B1_fails'] += 1

                if T1 < 0:
                    result['T1_neg'] += 1
                if T1 < result['min_T1']:
                    result['min_T1'] = T1

                if T2 < 0:
                    result['T2_neg'] += 1

                if C2_prev < 0:
                    result['C2_prev_neg'] += 1

                # Sign comparison (nonzero entries only)
                if T2 != 0 and C2_prev != 0:
                    if (T2 > 0) == (C2_prev > 0):
                        result['sign_match'] += 1
                    else:
                        result['sign_mismatch'] += 1

                # Ratio when T2 < 0
                if T2 < 0 and T1 > 0:
                    r = T1 / abs(T2)
                    if r < result['min_ratio_T1_T2']:
                        result['min_ratio_T1_T2'] = r
                    result['tightest'].append(
                        (r, k, stage, s, int(T1), int(T2), int(C2_prev), g6.strip()))
                    result['tightest'].sort(key=lambda x: x[0])
                    result['tightest'] = result['tightest'][:3]

            E_acc = E_new
            J_acc = J_new

    return result


def process_batch(g6_lines):
    keys = ['trees', 'total_checks', 'identity_fails', 'B1_fails',
            'T1_neg', 'T2_neg', 'C2_prev_neg', 'sign_match', 'sign_mismatch']
    batch = {k: 0 for k in keys}
    batch['min_T1'] = float('inf')
    batch['min_ratio'] = float('inf')
    batch['tightest'] = []

    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        for k in keys:
            if k != 'trees':
                batch[k] += r.get(k, 0)
        if r['min_T1'] < batch['min_T1']:
            batch['min_T1'] = r['min_T1']
        if r['min_ratio_T1_T2'] < batch['min_ratio']:
            batch['min_ratio'] = r['min_ratio_T1_T2']
        batch['tightest'].extend(r['tightest'])
        batch['tightest'].sort(key=lambda x: x[0])
        batch['tightest'] = batch['tightest'][:10]

    return batch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=500)
    args = parser.parse_args()

    totals = {
        'trees': 0, 'total_checks': 0, 'identity_fails': 0, 'B1_fails': 0,
        'T1_neg': 0, 'T2_neg': 0, 'C2_prev_neg': 0,
        'sign_match': 0, 'sign_mismatch': 0,
    }
    g_min_T1 = float('inf')
    g_min_ratio = float('inf')
    g_tight = []

    t0 = time.time()

    for nn in range(args.min_n, args.max_n + 1):
        tn = time.time()
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_totals = {k: 0 for k in totals}

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
                if r['min_T1'] < g_min_T1:
                    g_min_T1 = r['min_T1']
                if r['min_ratio_T1_T2'] < g_min_ratio:
                    g_min_ratio = r['min_ratio_T1_T2']
                g_tight.extend(r['tightest'])
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
                    if br['min_T1'] < g_min_T1:
                        g_min_T1 = br['min_T1']
                    if br.get('min_ratio', float('inf')) < g_min_ratio:
                        g_min_ratio = br['min_ratio']
                    g_tight.extend(br.get('tightest', []))

        proc.wait()
        elapsed_n = time.time() - tn
        for k in totals:
            totals[k] += n_totals[k]

        print(f"n={nn:3d}: {n_totals['trees']:>10,} trees | "
              f"id:{n_totals['identity_fails']}  B1<0:{n_totals['B1_fails']}  "
              f"T1<0:{n_totals['T1_neg']:>7,}  T2<0:{n_totals['T2_neg']:>7,}  "
              f"C2p<0:{n_totals['C2_prev_neg']:>7,}  "
              f"match:{n_totals['sign_match']:>7,}  mismatch:{n_totals['sign_mismatch']:>5,}  "
              f"[{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*90}")
    print(f"B1 SUB-TERM DECOMPOSITION")
    print(f"{'='*90}")
    print(f"Range:                    n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:              {totals['trees']:>15,}")
    print(f"Total checks:             {totals['total_checks']:>15,}")
    print(f"Identity failures:        {totals['identity_fails']:>15,}")
    print(f"B1 < 0:                   {totals['B1_fails']:>15,}")
    print()
    print(f"Sub-term 1 (TP2 of SCC_factor):")
    print(f"  T1 < 0:                {totals['T1_neg']:>15,}")
    print(f"  Min T1:                 {g_min_T1}")
    print()
    print(f"Sub-term 2 (TP2 of C2_prev):")
    print(f"  T2 < 0:                {totals['T2_neg']:>15,}")
    print(f"  C2_prev < 0:            {totals['C2_prev_neg']:>15,}")
    print(f"  Sign match (T2, C2p):   {totals['sign_match']:>15,}")
    print(f"  Sign mismatch:          {totals['sign_mismatch']:>15,}")
    if totals['sign_match'] + totals['sign_mismatch'] > 0:
        pct = 100 * totals['sign_match'] / (totals['sign_match'] + totals['sign_mismatch'])
        print(f"  Match rate:             {pct:.2f}%")
    print()
    print(f"Min T1/|T2| when T2 < 0: {g_min_ratio:.6f}")

    g_tight.sort(key=lambda x: x[0])
    g_tight = g_tight[:20]
    if g_tight:
        print(f"\n{'='*90}")
        print(f"TIGHTEST T1/|T2| ratios (when T2 < 0)")
        print(f"{'='*90}")
        print(f"{'Rank':>4s}  {'Ratio':>10s}  {'k':>3s}  {'Stage':>5s}  {'#fac':>4s}  "
              f"{'T1':>14s}  {'T2':>14s}  {'C2p':>10s}  g6")
        for rank, (ratio, k, stage, nfac, t1, t2, c2p, g6) in enumerate(g_tight):
            print(f"{rank+1:4d}  {ratio:10.6f}  {k:3d}  {stage:5d}  {nfac:4d}  "
                  f"{t1:14d}  {t2:14d}  {c2p:10d}  {g6}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
