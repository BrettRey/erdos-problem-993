"""2-term SCC decomposition via linearity in first/second argument.

Two decompositions of SCC = Δ_k(e_new, b_new):

(A) Linearity in 1st argument: e_new = α + (x+x²)γ
    SCC = Δ_k(α, b_new) + Δ_k((x+x²)γ, b_new)
    = (Term1+Term2) + (Term3+Term4+Term5)

(B) Linearity in 2nd argument: b_new = β + xγ
    SCC = Δ_k(e_new, β) + Δ_k(e_new, xγ)
    = (Term1+Term3+Term4) + (Term2+Term5)

Also checks the tree-level decomposition:
(C) SCC = c_k(b_new) + Δ_k(x(1+x)·J_acc, E_acc·P)
    where c_k(b) = LC gap = b_k² - b_{k-1}·b_{k+1}

Usage:
    python3 verify_scc_2term.py --max-n 22 --workers 8
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
        'scc_fails': 0,
        # Decomposition (A): by first argument
        'A1_neg': 0,  # Δ_k(α, b_new) < 0 (= Term1+Term2)
        'A2_neg': 0,  # Δ_k((x+x²)γ, b_new) < 0 (= Term3+Term4+Term5)
        # Decomposition (B): by second argument
        'B1_neg': 0,  # Δ_k(e_new, β) < 0 (= Term1+Term3+Term4)
        'B2_neg': 0,  # Δ_k(e_new, xγ) < 0 (= Term2+Term5)
        # Decomposition (C): LC + correction
        'C1_neg': 0,  # c_k(b_new) < 0
        'C2_neg': 0,  # Δ_k(x(1+x)J, b_new) < 0
        # Ratios when correction < 0
        'min_ratio_A': float('inf'),  # min A1/|A2| when A2 < 0
        'min_ratio_B': float('inf'),  # min B1/|B2| when B2 < 0
        'min_ratio_C': float('inf'),  # min C1/|C2| when C2 < 0
        'tightest_A': [],
        'tightest_B': [],
        'tightest_C': [],
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
            I_old = _polyadd(E_acc, [0] + J_acc)
            e_old = _polyadd(I_old, [0] + I_old)  # (1+x)*I
            b_old = E_acc

            E_new = _polymul(E_acc, P)
            J_new = _polymul(J_acc, Q)
            I_new = _polyadd(E_new, [0] + J_new)
            e_new = _polyadd(I_new, [0] + I_new)  # (1+x)*I_new

            # Convolutions
            alpha = _polymul(e_old, Q)    # e*Q
            beta = _polymul(b_old, Q)     # b*Q
            gamma = _polymul(b_old, R)    # b*R

            # For decomposition (C): x(1+x)*J_acc (not J_old·Q but full J_new)
            # Actually: e_new = (1+x)*b_new + x(1+x)*J_new
            # So SCC = c_k(b_new) + Δ_k(x(1+x)*J_new, b_new)
            xJ = [0] + list(J_new)       # x*J_new
            x1xJ = _polyadd(xJ, [0] + xJ)  # x(1+x)*J_new

            max_k = max(len(e_new), len(E_new)) - 1

            for k in range(max_k):
                # Direct SCC
                scc_direct = coeff(e_new, k+1) * coeff(E_new, k) - \
                             coeff(e_new, k) * coeff(E_new, k+1)

                # Decomposition (A): by first argument
                # A1 = Δ_k(α, b_new)
                A1 = coeff(alpha, k+1) * coeff(E_new, k) - \
                     coeff(alpha, k) * coeff(E_new, k+1)
                # A2 = Δ_k((x+x²)γ, b_new)
                # [(x+x²)γ]_k = γ_{k-1} + γ_{k-2}
                xr_k1 = coeff(gamma, k) + coeff(gamma, k-1)
                xr_k = coeff(gamma, k-1) + coeff(gamma, k-2)
                A2 = xr_k1 * coeff(E_new, k) - xr_k * coeff(E_new, k+1)

                # Decomposition (B): by second argument
                # B1 = Δ_k(e_new, β)
                B1 = coeff(e_new, k+1) * coeff(beta, k) - \
                     coeff(e_new, k) * coeff(beta, k+1)
                # B2 = Δ_k(e_new, xγ)
                B2 = coeff(e_new, k+1) * coeff(gamma, k-1) - \
                     coeff(e_new, k) * coeff(gamma, k)

                # Decomposition (C): LC + correction
                bk = coeff(E_new, k)
                bkp1 = coeff(E_new, k+1)
                bkm1 = coeff(E_new, k-1)
                C1 = bk * bk - bkm1 * bkp1  # LC gap
                C2 = coeff(x1xJ, k+1) * bk - coeff(x1xJ, k) * bkp1

                if scc_direct == 0 and A1 == 0 and A2 == 0:
                    continue

                result['total_checks'] += 1

                # Identity checks
                if A1 + A2 != scc_direct:
                    result['identity_fails'] += 1
                if B1 + B2 != scc_direct:
                    result['identity_fails'] += 1
                if C1 + C2 != scc_direct:
                    result['identity_fails'] += 1

                if scc_direct < 0:
                    result['scc_fails'] += 1

                # Signs
                if A1 < 0:
                    result['A1_neg'] += 1
                if A2 < 0:
                    result['A2_neg'] += 1
                if B1 < 0:
                    result['B1_neg'] += 1
                if B2 < 0:
                    result['B2_neg'] += 1
                if C1 < 0:
                    result['C1_neg'] += 1
                if C2 < 0:
                    result['C2_neg'] += 1

                # Ratios for (A)
                if A2 < 0 and A1 > 0:
                    r = A1 / abs(A2)
                    if r < result['min_ratio_A']:
                        result['min_ratio_A'] = r
                    result['tightest_A'].append(
                        (r, k, stage, s, int(A1), int(A2), g6.strip()))
                    result['tightest_A'].sort(key=lambda x: x[0])
                    result['tightest_A'] = result['tightest_A'][:3]
                elif A1 < 0 and A2 > 0:
                    r = A2 / abs(A1)
                    if r < result['min_ratio_A']:
                        result['min_ratio_A'] = r

                # Ratios for (B)
                if B2 < 0 and B1 > 0:
                    r = B1 / abs(B2)
                    if r < result['min_ratio_B']:
                        result['min_ratio_B'] = r
                    result['tightest_B'].append(
                        (r, k, stage, s, int(B1), int(B2), g6.strip()))
                    result['tightest_B'].sort(key=lambda x: x[0])
                    result['tightest_B'] = result['tightest_B'][:3]

                # Ratios for (C)
                if C2 < 0 and C1 > 0:
                    r = C1 / abs(C2)
                    if r < result['min_ratio_C']:
                        result['min_ratio_C'] = r
                    result['tightest_C'].append(
                        (r, k, stage, s, int(C1), int(C2), g6.strip()))
                    result['tightest_C'].sort(key=lambda x: x[0])
                    result['tightest_C'] = result['tightest_C'][:3]

            E_acc = E_new
            J_acc = J_new

    return result


def process_batch(g6_lines):
    keys = ['trees', 'total_checks', 'identity_fails', 'scc_fails',
            'A1_neg', 'A2_neg', 'B1_neg', 'B2_neg', 'C1_neg', 'C2_neg']
    batch = {k: 0 for k in keys}
    batch['min_ratio_A'] = float('inf')
    batch['min_ratio_B'] = float('inf')
    batch['min_ratio_C'] = float('inf')
    batch['tightest_A'] = []
    batch['tightest_B'] = []
    batch['tightest_C'] = []

    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        for k in keys:
            if k != 'trees':
                batch[k] += r.get(k, 0)
        for label in ['A', 'B', 'C']:
            rk = f'min_ratio_{label}'
            tk = f'tightest_{label}'
            if r[rk] < batch[rk]:
                batch[rk] = r[rk]
            batch[tk].extend(r[tk])
            batch[tk].sort(key=lambda x: x[0])
            batch[tk] = batch[tk][:10]

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
        'A1_neg': 0, 'A2_neg': 0, 'B1_neg': 0, 'B2_neg': 0, 'C1_neg': 0, 'C2_neg': 0,
    }
    g_min_A = float('inf')
    g_min_B = float('inf')
    g_min_C = float('inf')
    g_tight_A = []
    g_tight_B = []
    g_tight_C = []

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
                for label in ['A', 'B', 'C']:
                    rk = f'min_ratio_{label}'
                    tk = f'tightest_{label}'
                    gmin = {'A': g_min_A, 'B': g_min_B, 'C': g_min_C}
                    if r[rk] < gmin[label]:
                        if label == 'A': g_min_A = r[rk]
                        elif label == 'B': g_min_B = r[rk]
                        else: g_min_C = r[rk]
                    if label == 'A': g_tight_A.extend(r[tk])
                    elif label == 'B': g_tight_B.extend(r[tk])
                    else: g_tight_C.extend(r[tk])
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
                    for label in ['A', 'B', 'C']:
                        rk = f'min_ratio_{label}'
                        tk = f'tightest_{label}'
                        if br[rk] < {'A': g_min_A, 'B': g_min_B, 'C': g_min_C}[label]:
                            if label == 'A': g_min_A = br[rk]
                            elif label == 'B': g_min_B = br[rk]
                            else: g_min_C = br[rk]
                        if label == 'A': g_tight_A.extend(br[tk])
                        elif label == 'B': g_tight_B.extend(br[tk])
                        else: g_tight_C.extend(br[tk])

        proc.wait()
        elapsed_n = time.time() - tn
        for k in totals:
            totals[k] += n_totals[k]

        print(f"n={nn:3d}: {n_totals['trees']:>10,} trees | "
              f"id:{n_totals['identity_fails']}  scc:{n_totals['scc_fails']}  "
              f"A1<0:{n_totals['A1_neg']:>7,}  A2<0:{n_totals['A2_neg']:>7,}  "
              f"B1<0:{n_totals['B1_neg']:>7,}  B2<0:{n_totals['B2_neg']:>7,}  "
              f"C1<0:{n_totals['C1_neg']:>7,}  C2<0:{n_totals['C2_neg']:>7,}  "
              f"[{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*90}")
    print(f"SCC 2-TERM DECOMPOSITIONS")
    print(f"{'='*90}")
    print(f"Range:                    n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:              {totals['trees']:>15,}")
    print(f"Total checks:             {totals['total_checks']:>15,}")
    print(f"Identity failures:        {totals['identity_fails']:>15,}")
    print(f"SCC failures:             {totals['scc_fails']:>15,}")
    print()
    print(f"(A) Linearity in 1st arg: e_new = eQ + x(1+x)bR")
    print(f"    A1 = Δ(eQ, bP) < 0:  {totals['A1_neg']:>15,}")
    print(f"    A2 = Δ(x(1+x)bR, bP) < 0: {totals['A2_neg']:>12,}")
    print(f"    Min ratio (pos/neg):  {g_min_A:.6f}")
    print()
    print(f"(B) Linearity in 2nd arg: b_new = bQ + xbR")
    print(f"    B1 = Δ(e_new, bQ) < 0: {totals['B1_neg']:>13,}")
    print(f"    B2 = Δ(e_new, xbR) < 0: {totals['B2_neg']:>12,}")
    print(f"    Min ratio (pos/neg):  {g_min_B:.6f}")
    print()
    print(f"(C) LC + correction: SCC = c_k(E) + Δ(x(1+x)J, E)")
    print(f"    C1 = c_k(E) < 0:     {totals['C1_neg']:>15,}")
    print(f"    C2 = Δ(x(1+x)J, E) < 0: {totals['C2_neg']:>11,}")
    print(f"    Min ratio (pos/neg):  {g_min_C:.6f}")

    for label, tight, gmin in [('A', g_tight_A, g_min_A),
                                ('B', g_tight_B, g_min_B),
                                ('C', g_tight_C, g_min_C)]:
        tight.sort(key=lambda x: x[0])
        tight = tight[:15]
        if tight:
            print(f"\n{'='*90}")
            print(f"TIGHTEST ratios for decomposition ({label})")
            print(f"{'='*90}")
            print(f"{'Rank':>4s}  {'Ratio':>10s}  {'k':>3s}  {'Stage':>5s}  {'#fac':>4s}  "
                  f"{'Pos':>14s}  {'Neg':>14s}  g6")
            for rank, (ratio, k, stage, nfac, pos, neg, g6) in enumerate(tight):
                print(f"{rank+1:4d}  {ratio:10.6f}  {k:3d}  {stage:5d}  {nfac:4d}  "
                      f"{pos:14d}  {neg:14d}  {g6}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
