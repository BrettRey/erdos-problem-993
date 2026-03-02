"""Profile T1/T2/T3 rescue mechanism at incremental product stages.

For each tree at each support vertex, builds E^(k) incrementally and
at every intermediate stage, decomposes SCC = T1 + T2 + T3 where:
  T1 = b_{k-1} * d_k    (current LR minor, can be negative)
  T2 = b_k * d_{k-1}    (memory of previous LR minor)
  T3 = a_{k-1} * c_k    (LC curvature bonus)

Key questions:
  1. Does T3 >= |T1| always when T1 < 0?  (curvature alone suffices)
  2. At which k values and tree structures is T3/|T1| tightest?
  3. Is T2 ever needed as sole rescue (T3 < |T1| but T2 + T3 >= |T1|)?

Usage:
    python3 profile_T123_incremental.py --max-n 18 --workers 8
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
    if k < 0 or k >= len(poly):
        return 0
    return poly[k]


def process_tree(g6):
    n, adj = parse_g6(g6)

    result = {
        'n': n,
        'total_checks': 0,          # total (stage, k) with nonzero terms
        'dk_neg_total': 0,           # total with d_k < 0
        'dk_neg_k1': 0,              # d_k < 0 at k=1
        'dk_neg_k2plus': 0,          # d_k < 0 at k >= 2
        'T3_alone_rescues': 0,       # T3 >= |T1| (curvature alone suffices)
        'T3_alone_rescues_k2plus': 0,
        'T2_needed': 0,              # T3 < |T1| but T2+T3 >= |T1| (memory needed)
        'T2_needed_k2plus': 0,
        'SCC_fails': 0,              # T1+T2+T3 < 0
        # Tightest T3/|T1| ratios at k>=2
        'tightest_k2plus': [],       # (ratio, k, stage, nfactors, g6)
        # Tightest SCC margin overall
        'min_scc': float('inf'),
        'min_scc_k': -1,
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
            factors.append((Ic, Ec))

        E_acc = [1]
        J_acc = [1]

        for stage, (Ic, Ec) in enumerate(factors, 1):
            E_acc = _polymul(E_acc, Ic)
            J_acc = _polymul(J_acc, Ec)
            I_acc = _polyadd(E_acc, [0] + J_acc)

            max_k = max(len(I_acc), len(E_acc))

            for k in range(1, max_k):
                bkm1 = coeff(E_acc, k - 1)
                if bkm1 == 0:
                    continue

                ak = coeff(I_acc, k)
                akp1 = coeff(I_acc, k + 1)
                akm1 = coeff(I_acc, k - 1)
                bk = coeff(E_acc, k)
                bkp1 = coeff(E_acc, k + 1)

                dk = akp1 * bk - ak * bkp1
                dkm1 = ak * bkm1 - akm1 * bk
                ck = bk * bk - bkm1 * bkp1

                T1 = bkm1 * dk
                T2 = bk * dkm1
                T3 = akm1 * ck
                S = T1 + T2 + T3

                if T1 == 0 and T2 == 0 and T3 == 0:
                    continue

                result['total_checks'] += 1

                if S < result['min_scc']:
                    result['min_scc'] = S
                    result['min_scc_k'] = k

                if S < 0:
                    result['SCC_fails'] += 1

                if dk < 0:
                    result['dk_neg_total'] += 1
                    is_k1 = (k == 1)
                    if is_k1:
                        result['dk_neg_k1'] += 1
                    else:
                        result['dk_neg_k2plus'] += 1

                    # Check rescue mechanism
                    if T3 >= abs(T1):
                        result['T3_alone_rescues'] += 1
                        if not is_k1:
                            result['T3_alone_rescues_k2plus'] += 1
                    else:
                        # T3 < |T1|, need T2
                        if T2 + T3 >= abs(T1):
                            result['T2_needed'] += 1
                            if not is_k1:
                                result['T2_needed_k2plus'] += 1
                        # else: SCC fails (counted above)

                    # Record ratio at k>=2
                    if not is_k1 and abs(T1) > 0:
                        ratio = T3 / abs(T1)
                        result['tightest_k2plus'].append(
                            (ratio, k, stage, s, g6.strip(), T1, T2, T3, S)
                        )

    # Keep only the 5 tightest
    result['tightest_k2plus'].sort(key=lambda x: x[0])
    result['tightest_k2plus'] = result['tightest_k2plus'][:5]

    return result


def process_batch(g6_lines):
    batch = {
        'trees': 0,
        'total_checks': 0,
        'dk_neg_total': 0,
        'dk_neg_k1': 0,
        'dk_neg_k2plus': 0,
        'T3_alone_rescues': 0,
        'T3_alone_rescues_k2plus': 0,
        'T2_needed': 0,
        'T2_needed_k2plus': 0,
        'SCC_fails': 0,
        'tightest_k2plus': [],
        'min_scc': float('inf'),
    }
    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        batch['total_checks'] += r['total_checks']
        batch['dk_neg_total'] += r['dk_neg_total']
        batch['dk_neg_k1'] += r['dk_neg_k1']
        batch['dk_neg_k2plus'] += r['dk_neg_k2plus']
        batch['T3_alone_rescues'] += r['T3_alone_rescues']
        batch['T3_alone_rescues_k2plus'] += r['T3_alone_rescues_k2plus']
        batch['T2_needed'] += r['T2_needed']
        batch['T2_needed_k2plus'] += r['T2_needed_k2plus']
        batch['SCC_fails'] += r['SCC_fails']
        batch['tightest_k2plus'].extend(r['tightest_k2plus'])
        if r['min_scc'] < batch['min_scc']:
            batch['min_scc'] = r['min_scc']

    batch['tightest_k2plus'].sort(key=lambda x: x[0])
    batch['tightest_k2plus'] = batch['tightest_k2plus'][:20]
    return batch


def main():
    parser = argparse.ArgumentParser(description='Profile T1/T2/T3 at incremental product stages')
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=500)
    args = parser.parse_args()

    totals = {
        'trees': 0, 'total_checks': 0,
        'dk_neg_total': 0, 'dk_neg_k1': 0, 'dk_neg_k2plus': 0,
        'T3_alone_rescues': 0, 'T3_alone_rescues_k2plus': 0,
        'T2_needed': 0, 'T2_needed_k2plus': 0,
        'SCC_fails': 0,
    }
    global_tightest = []  # tightest at k>=2

    t0 = time.time()

    for nn in range(args.min_n, args.max_n + 1):
        tn = time.time()
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_totals = {k: 0 for k in totals}
        n_tightest = []

        if args.workers <= 1:
            for line in proc.stdout:
                line = line.strip()
                if not line:
                    continue
                r = process_tree(line)
                n_totals['trees'] += 1
                n_totals['total_checks'] += r['total_checks']
                n_totals['dk_neg_total'] += r['dk_neg_total']
                n_totals['dk_neg_k1'] += r['dk_neg_k1']
                n_totals['dk_neg_k2plus'] += r['dk_neg_k2plus']
                n_totals['T3_alone_rescues'] += r['T3_alone_rescues']
                n_totals['T3_alone_rescues_k2plus'] += r['T3_alone_rescues_k2plus']
                n_totals['T2_needed'] += r['T2_needed']
                n_totals['T2_needed_k2plus'] += r['T2_needed_k2plus']
                n_totals['SCC_fails'] += r['SCC_fails']
                n_tightest.extend(r['tightest_k2plus'])
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
                    for k in n_totals:
                        n_totals[k] += br[k]
                    n_tightest.extend(br['tightest_k2plus'])

        proc.wait()
        elapsed_n = time.time() - tn

        for k in totals:
            totals[k] += n_totals[k]

        n_tightest.sort(key=lambda x: x[0])
        n_tightest = n_tightest[:10]
        global_tightest.extend(n_tightest)
        global_tightest.sort(key=lambda x: x[0])
        global_tightest = global_tightest[:50]

        # Per-n summary
        dk2 = n_totals['dk_neg_k2plus']
        t3r = n_totals['T3_alone_rescues_k2plus']
        t2n = n_totals['T2_needed_k2plus']
        min_r = n_tightest[0][0] if n_tightest else float('inf')
        print(f"n={nn:3d}: {n_totals['trees']:>10,} trees | "
              f"dk<0 k=1:{n_totals['dk_neg_k1']:>8,}  k≥2:{dk2:>8,} | "
              f"T3≥|T1|(k≥2):{t3r:>8,}  T2needed:{t2n:>6,}  "
              f"min_ratio={min_r:.4f}  [{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*90}")
    print(f"T1/T2/T3 RESCUE MECHANISM PROFILE")
    print(f"{'='*90}")
    print(f"Range:                    n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:              {totals['trees']:>15,}")
    print(f"Total (stage,k) checks:   {totals['total_checks']:>15,}")
    print(f"SCC failures:             {totals['SCC_fails']:>15,}")
    print()
    print(f"d_k < 0 events:")
    print(f"  Total:                  {totals['dk_neg_total']:>15,}")
    print(f"  At k=1:                 {totals['dk_neg_k1']:>15,}")
    print(f"  At k>=2:                {totals['dk_neg_k2plus']:>15,}")
    print()
    print(f"Rescue mechanism (all k):")
    pct_t3 = 100 * totals['T3_alone_rescues'] / max(1, totals['dk_neg_total'])
    pct_t2 = 100 * totals['T2_needed'] / max(1, totals['dk_neg_total'])
    print(f"  T3 >= |T1| (curvature alone):  {totals['T3_alone_rescues']:>12,}  ({pct_t3:.1f}%)")
    print(f"  T2 needed (memory required):   {totals['T2_needed']:>12,}  ({pct_t2:.1f}%)")
    print()
    print(f"Rescue mechanism (k >= 2 only):")
    pct_t3_k2 = 100 * totals['T3_alone_rescues_k2plus'] / max(1, totals['dk_neg_k2plus'])
    pct_t2_k2 = 100 * totals['T2_needed_k2plus'] / max(1, totals['dk_neg_k2plus'])
    print(f"  T3 >= |T1| (curvature alone):  {totals['T3_alone_rescues_k2plus']:>12,}  ({pct_t3_k2:.1f}%)")
    print(f"  T2 needed (memory required):   {totals['T2_needed_k2plus']:>12,}  ({pct_t2_k2:.1f}%)")

    print(f"\n{'='*90}")
    print(f"TIGHTEST T3/|T1| AT k >= 2")
    print(f"{'='*90}")
    print(f"{'Rank':>4s}  {'Ratio':>10s}  {'k':>3s}  {'Stage':>5s}  {'#fac':>4s}  "
          f"{'T1':>14s}  {'T2':>14s}  {'T3':>14s}  {'SCC':>14s}  g6")
    for rank, (ratio, k, stage, nfac, g6, T1, T2, T3, S) in enumerate(global_tightest[:30]):
        print(f"{rank+1:4d}  {ratio:10.6f}  {k:3d}  {stage:5d}  {nfac:4d}  "
              f"{T1:14d}  {T2:14d}  {T3:14d}  {S:14d}  {g6}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
