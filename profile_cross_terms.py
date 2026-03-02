"""Profile the 4-term bilinear decomposition of SCC.

SCC = Δ_k(e_old·Q + x(1+x)·E_old·R, E_old·Q + x·E_old·R)
    = Term_QQ + Term_QR + Term_RQ + Term_RR

where:
  Term_QQ = Δ_k(e_old·Q, E_old·Q)  ≥ 0 PROVED (Karlin + SCC_old)
  Term_QR = Δ_k(e_old·Q, x·E_old·R)  cross term
  Term_RQ = Δ_k(x(1+x)·E_old·R, E_old·Q)  cross term
  Term_RR = Δ_k(x(1+x)·E_old·R, x·E_old·R) = c_{k-1}(E_old·R)  ≥ 0 PROVED (LC)

Questions:
  1. Do the cross terms partially cancel?
  2. Is Term_QR + Term_RQ ≤ 0 always?
  3. Is |Term_QR + Term_RQ| ≤ Term_QQ + Term_RR?
  4. Is there a Cauchy-Schwarz-type bound?

Usage:
    python3 profile_cross_terms.py --max-n 20 --workers 8
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
        'QQ_neg': 0, 'QR_neg': 0, 'RQ_neg': 0, 'RR_neg': 0,
        'cross_neg': 0,      # Term_QR + Term_RQ < 0
        'cross_pos': 0,      # Term_QR + Term_RQ > 0
        'diag_dominates': 0, # Term_QQ + Term_RR ≥ |cross|
        'diag_fails': 0,     # Term_QQ + Term_RR < |cross|
        'cs_holds': 0,       # |Term_QR·Term_RQ| ≤ Term_QQ·Term_RR
        'cs_fails': 0,       # Cauchy-Schwarz fails
        'min_diag_ratio': float('inf'),  # min (QQ+RR)/|cross|
        'min_cs_ratio': float('inf'),    # min QQ·RR / (QR·RQ) when QR·RQ > 0
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
            I_old = _polyadd(E_acc, [0] + J_acc)
            e_old = _polyadd(I_old, [0] + I_old)  # (1+x)*I_old

            # Convolutions for the 4 terms
            eQ = _polymul(e_old, Q)           # e_old * Q
            EQ = _polymul(E_acc, Q)           # E_old * Q
            ER = _polymul(E_acc, R)           # E_old * R
            x1xER = _polyadd([0] + ER, [0, 0] + ER)  # x(1+x) * E_old * R
            xER = [0] + list(ER)              # x * E_old * R

            # New accumulated
            E_new = _polymul(E_acc, P)
            J_new = _polymul(J_acc, Q)
            I_new = _polyadd(E_new, [0] + J_new)
            e_new = _polyadd(I_new, [0] + I_new)

            max_k = max(len(e_new), len(E_new)) - 1

            for k in range(max_k):
                # Direct SCC
                scc = coeff(e_new, k+1) * coeff(E_new, k) - \
                      coeff(e_new, k) * coeff(E_new, k+1)

                # 4-term decomposition
                QQ = coeff(eQ, k+1) * coeff(EQ, k) - coeff(eQ, k) * coeff(EQ, k+1)
                QR = coeff(eQ, k+1) * coeff(xER, k) - coeff(eQ, k) * coeff(xER, k+1)
                RQ = coeff(x1xER, k+1) * coeff(EQ, k) - coeff(x1xER, k) * coeff(EQ, k+1)
                RR = coeff(x1xER, k+1) * coeff(xER, k) - coeff(x1xER, k) * coeff(xER, k+1)

                if scc == 0 and QQ == 0 and RR == 0:
                    continue

                result['total_checks'] += 1

                if QQ + QR + RQ + RR != scc:
                    result['identity_fails'] += 1

                if QQ < 0: result['QQ_neg'] += 1
                if QR < 0: result['QR_neg'] += 1
                if RQ < 0: result['RQ_neg'] += 1
                if RR < 0: result['RR_neg'] += 1

                cross = QR + RQ
                diag = QQ + RR

                if cross < 0:
                    result['cross_neg'] += 1
                elif cross > 0:
                    result['cross_pos'] += 1

                # Diagonal domination
                if cross != 0:
                    if diag >= abs(cross):
                        result['diag_dominates'] += 1
                        r = diag / abs(cross)
                        if r < result['min_diag_ratio']:
                            result['min_diag_ratio'] = r
                    else:
                        result['diag_fails'] += 1
                        r = diag / abs(cross)
                        result['tightest'].append(
                            (r, k, stage, s, int(QQ), int(QR), int(RQ), int(RR),
                             int(diag), int(cross), g6.strip()))
                        result['tightest'].sort(key=lambda x: x[0])
                        result['tightest'] = result['tightest'][:5]

                # Cauchy-Schwarz type: QQ·RR ≥ QR·RQ (when QR·RQ > 0)
                if QR * RQ > 0:  # same sign
                    if QQ * RR >= QR * RQ:
                        result['cs_holds'] += 1
                    else:
                        result['cs_fails'] += 1
                        if QQ * RR > 0:
                            r = (QQ * RR) / (QR * RQ)
                            if r < result['min_cs_ratio']:
                                result['min_cs_ratio'] = r

            E_acc = E_new
            J_acc = J_new

    return result


def process_batch(g6_lines):
    keys = ['trees', 'total_checks', 'identity_fails',
            'QQ_neg', 'QR_neg', 'RQ_neg', 'RR_neg',
            'cross_neg', 'cross_pos', 'diag_dominates', 'diag_fails',
            'cs_holds', 'cs_fails']
    batch = {k: 0 for k in keys}
    batch['min_diag_ratio'] = float('inf')
    batch['min_cs_ratio'] = float('inf')
    batch['tightest'] = []

    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        for k in keys:
            if k != 'trees':
                batch[k] += r.get(k, 0)
        if r['min_diag_ratio'] < batch['min_diag_ratio']:
            batch['min_diag_ratio'] = r['min_diag_ratio']
        if r['min_cs_ratio'] < batch['min_cs_ratio']:
            batch['min_cs_ratio'] = r['min_cs_ratio']
        batch['tightest'].extend(r['tightest'])
        batch['tightest'].sort(key=lambda x: x[0])
        batch['tightest'] = batch['tightest'][:20]

    return batch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=500)
    args = parser.parse_args()

    totals = {
        'trees': 0, 'total_checks': 0, 'identity_fails': 0,
        'QQ_neg': 0, 'QR_neg': 0, 'RQ_neg': 0, 'RR_neg': 0,
        'cross_neg': 0, 'cross_pos': 0, 'diag_dominates': 0, 'diag_fails': 0,
        'cs_holds': 0, 'cs_fails': 0,
    }
    g_min_diag = float('inf')
    g_min_cs = float('inf')
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
                if r['min_diag_ratio'] < g_min_diag:
                    g_min_diag = r['min_diag_ratio']
                if r['min_cs_ratio'] < g_min_cs:
                    g_min_cs = r['min_cs_ratio']
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
                    if br['min_diag_ratio'] < g_min_diag:
                        g_min_diag = br['min_diag_ratio']
                    if br['min_cs_ratio'] < g_min_cs:
                        g_min_cs = br['min_cs_ratio']
                    g_tight.extend(br.get('tightest', []))

        proc.wait()
        elapsed_n = time.time() - tn
        for k in totals:
            totals[k] += n_totals[k]

        print(f"n={nn:3d}: {n_totals['trees']:>10,} trees | "
              f"QQ<0:{n_totals['QQ_neg']}  QR<0:{n_totals['QR_neg']:>7,}  "
              f"RQ<0:{n_totals['RQ_neg']:>7,}  RR<0:{n_totals['RR_neg']}  "
              f"cross<0:{n_totals['cross_neg']:>7,}  "
              f"diag≥|cr|:{n_totals['diag_dominates']:>7,}  "
              f"diag<|cr|:{n_totals['diag_fails']:>5,}  "
              f"[{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*90}")
    print(f"4-TERM BILINEAR DECOMPOSITION")
    print(f"{'='*90}")
    print(f"Range:                    n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:              {totals['trees']:>15,}")
    print(f"Total checks:             {totals['total_checks']:>15,}")
    print(f"Identity failures:        {totals['identity_fails']:>15,}")
    print()
    print(f"Term signs:")
    print(f"  QQ < 0 (should be 0):   {totals['QQ_neg']:>15,}")
    print(f"  QR < 0:                 {totals['QR_neg']:>15,}")
    print(f"  RQ < 0:                 {totals['RQ_neg']:>15,}")
    print(f"  RR < 0 (should be 0):   {totals['RR_neg']:>15,}")
    print()
    print(f"Cross terms (QR + RQ):")
    print(f"  Negative:               {totals['cross_neg']:>15,}")
    print(f"  Positive:               {totals['cross_pos']:>15,}")
    cross_total = totals['cross_neg'] + totals['cross_pos']
    if cross_total > 0:
        print(f"  % negative:             {100*totals['cross_neg']/cross_total:.1f}%")
    print()
    print(f"Diagonal domination (QQ+RR ≥ |QR+RQ|):")
    print(f"  Holds:                  {totals['diag_dominates']:>15,}")
    print(f"  Fails:                  {totals['diag_fails']:>15,}")
    if g_min_diag < float('inf'):
        print(f"  Min (QQ+RR)/|cross|:    {g_min_diag:.6f}")
    print()
    print(f"Cauchy-Schwarz (QQ·RR ≥ QR·RQ when same sign):")
    print(f"  Holds:                  {totals['cs_holds']:>15,}")
    print(f"  Fails:                  {totals['cs_fails']:>15,}")
    if g_min_cs < float('inf'):
        print(f"  Min QQ·RR/(QR·RQ):      {g_min_cs:.6f}")

    g_tight.sort(key=lambda x: x[0])
    g_tight = g_tight[:20]
    if g_tight:
        print(f"\n{'='*90}")
        print(f"TIGHTEST diagonal domination failures (diag/|cross| < 1)")
        print(f"{'='*90}")
        print(f"{'Rank':>4s}  {'Ratio':>8s}  {'k':>3s}  {'Stg':>3s}  {'#f':>3s}  "
              f"{'QQ':>10s}  {'QR':>10s}  {'RQ':>10s}  {'RR':>10s}  "
              f"{'diag':>10s}  {'cross':>10s}")
        for rank, item in enumerate(g_tight[:20]):
            ratio, k, stage, nfac, qq, qr, rq, rr, diag, cross, g6 = item
            print(f"{rank+1:4d}  {ratio:8.4f}  {k:3d}  {stage:3d}  {nfac:3d}  "
                  f"{qq:10d}  {qr:10d}  {rq:10d}  {rr:10d}  "
                  f"{diag:10d}  {cross:10d}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
