"""When cross = QR + RQ < 0, profile (QQ + RR) / |cross| margin.

Since SCC = QQ + cross + RR ≥ 0 always, when cross < 0 we must have
diag = QQ + RR ≥ |cross|. Question: how tight is this?

Usage:
    python3 profile_cross_negative.py --max-n 22 --workers 8
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
        'total': 0,
        'cross_neg': 0,
        'min_diag_over_cross': float('inf'),  # min (QQ+RR)/|cross| when cross < 0
        'min_QQ_over_cross': float('inf'),     # min QQ/|cross| when cross < 0
        'min_RR_over_cross': float('inf'),     # min RR/|cross| when cross < 0
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
            e_old = _polyadd(I_old, [0] + I_old)

            eQ = _polymul(e_old, Q)
            EQ = _polymul(E_acc, Q)
            ER = _polymul(E_acc, R)
            x1xER = _polyadd([0] + ER, [0, 0] + ER)
            xER = [0] + list(ER)

            E_new = _polymul(E_acc, P)
            J_new = _polymul(J_acc, Q)
            I_new = _polyadd(E_new, [0] + J_new)
            e_new = _polyadd(I_new, [0] + I_new)

            max_k = max(len(e_new), len(E_new)) - 1

            for k in range(max_k):
                QQ = coeff(eQ, k+1) * coeff(EQ, k) - coeff(eQ, k) * coeff(EQ, k+1)
                QR = coeff(eQ, k+1) * coeff(xER, k) - coeff(eQ, k) * coeff(xER, k+1)
                RQ = coeff(x1xER, k+1) * coeff(EQ, k) - coeff(x1xER, k) * coeff(EQ, k+1)
                RR = coeff(x1xER, k+1) * coeff(xER, k) - coeff(x1xER, k) * coeff(xER, k+1)

                cross = QR + RQ
                diag = QQ + RR
                scc = diag + cross

                if scc == 0 and cross == 0:
                    continue

                result['total'] += 1

                if cross < 0:
                    result['cross_neg'] += 1
                    r = diag / abs(cross) if cross != 0 else float('inf')
                    rqq = QQ / abs(cross) if cross != 0 else float('inf')
                    rrr = RR / abs(cross) if cross != 0 else float('inf')

                    if r < result['min_diag_over_cross']:
                        result['min_diag_over_cross'] = r
                    if rqq < result['min_QQ_over_cross']:
                        result['min_QQ_over_cross'] = rqq
                    if rrr < result['min_RR_over_cross']:
                        result['min_RR_over_cross'] = rrr

                    if r < 5.0:
                        result['tightest'].append(
                            (r, k, stage, s, int(QQ), int(QR), int(RQ), int(RR),
                             int(scc), g6.strip()))
                        result['tightest'].sort(key=lambda x: x[0])
                        result['tightest'] = result['tightest'][:5]

            E_acc = E_new
            J_acc = J_new

    return result


def process_batch(g6_lines):
    batch = {
        'trees': 0, 'total': 0, 'cross_neg': 0,
        'min_diag': float('inf'), 'min_QQ': float('inf'), 'min_RR': float('inf'),
        'tightest': [],
    }

    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        batch['total'] += r['total']
        batch['cross_neg'] += r['cross_neg']
        if r['min_diag_over_cross'] < batch['min_diag']:
            batch['min_diag'] = r['min_diag_over_cross']
        if r['min_QQ_over_cross'] < batch['min_QQ']:
            batch['min_QQ'] = r['min_QQ_over_cross']
        if r['min_RR_over_cross'] < batch['min_RR']:
            batch['min_RR'] = r['min_RR_over_cross']
        batch['tightest'].extend(r['tightest'])
        batch['tightest'].sort(key=lambda x: x[0])
        batch['tightest'] = batch['tightest'][:20]

    return batch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=20)
    parser.add_argument('--min-n', type=int, default=3)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=500)
    args = parser.parse_args()

    total_trees = 0
    total_checks = 0
    total_cross_neg = 0
    g_min_diag = float('inf')
    g_min_QQ = float('inf')
    g_min_RR = float('inf')
    g_tight = []

    t0 = time.time()

    for nn in range(args.min_n, args.max_n + 1):
        tn = time.time()
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_trees = 0
        n_checks = 0
        n_cross_neg = 0
        n_min_diag = float('inf')

        if args.workers <= 1:
            for line in proc.stdout:
                line = line.strip()
                if not line:
                    continue
                r = process_tree(line)
                n_trees += 1
                n_checks += r['total']
                n_cross_neg += r['cross_neg']
                if r['min_diag_over_cross'] < n_min_diag:
                    n_min_diag = r['min_diag_over_cross']
                if r['min_diag_over_cross'] < g_min_diag:
                    g_min_diag = r['min_diag_over_cross']
                if r['min_QQ_over_cross'] < g_min_QQ:
                    g_min_QQ = r['min_QQ_over_cross']
                if r['min_RR_over_cross'] < g_min_RR:
                    g_min_RR = r['min_RR_over_cross']
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
                    n_trees += br['trees']
                    n_checks += br['total']
                    n_cross_neg += br['cross_neg']
                    if br['min_diag'] < n_min_diag:
                        n_min_diag = br['min_diag']
                    if br['min_diag'] < g_min_diag:
                        g_min_diag = br['min_diag']
                    if br['min_QQ'] < g_min_QQ:
                        g_min_QQ = br['min_QQ']
                    if br['min_RR'] < g_min_RR:
                        g_min_RR = br['min_RR']
                    g_tight.extend(br.get('tightest', []))

        proc.wait()
        elapsed_n = time.time() - tn
        total_trees += n_trees
        total_checks += n_checks
        total_cross_neg += n_cross_neg

        min_str = f"{n_min_diag:.4f}" if n_min_diag < float('inf') else "n/a"
        print(f"n={nn:3d}: {n_trees:>10,} trees | "
              f"cross<0: {n_cross_neg:>8,} / {n_checks:>10,}  "
              f"min diag/|cross|: {min_str}  [{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*80}")
    print(f"CROSS-NEGATIVE MARGIN ANALYSIS")
    print(f"{'='*80}")
    print(f"Range:           n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:     {total_trees:>12,}")
    print(f"Total checks:    {total_checks:>12,}")
    print(f"Cross < 0:       {total_cross_neg:>12,}  ({100*total_cross_neg/total_checks:.2f}%)")
    print()
    print(f"When cross < 0:")
    print(f"  Min (QQ+RR)/|cross|: {g_min_diag:.6f}")
    print(f"  Min QQ/|cross|:      {g_min_QQ:.6f}")
    print(f"  Min RR/|cross|:      {g_min_RR:.6f}")

    g_tight.sort(key=lambda x: x[0])
    g_tight = g_tight[:30]
    if g_tight:
        print(f"\n{'='*80}")
        print(f"TIGHTEST (QQ+RR)/|cross| when cross < 0")
        print(f"{'='*80}")
        print(f"{'Rank':>4s}  {'Ratio':>8s}  {'k':>3s}  {'Stg':>3s}  {'#f':>3s}  "
              f"{'QQ':>12s}  {'QR':>12s}  {'RQ':>12s}  {'RR':>12s}  {'SCC':>12s}")
        for rank, item in enumerate(g_tight):
            ratio, k, stage, nfac, qq, qr, rq, rr, scc, g6 = item
            print(f"{rank+1:4d}  {ratio:8.4f}  {k:3d}  {stage:3d}  {nfac:3d}  "
                  f"{qq:12d}  {qr:12d}  {rq:12d}  {rr:12d}  {scc:12d}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
