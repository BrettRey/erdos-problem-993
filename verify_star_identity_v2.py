"""Refined ★ identity check: separate P2/T_k failures by (nfac, stage).

Key question: do P2 failures occur at genuinely INTERMEDIATE stages
(stage < nfac) of multi-factor products, or only at final stages?

Usage:
    python3 verify_star_identity_v2.py --max-n 18 --workers 8
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
        # P2 failures by category
        'p2_fail_final': 0,         # stage == nfac (final stage, includes nfac=1)
        'p2_fail_intermediate': 0,  # stage < nfac (genuinely intermediate)
        # T_k < 0 by category
        'Tk_neg_final': 0,
        'Tk_neg_intermediate': 0,
        # D_k < 0 by category
        'Dk_neg_final': 0,
        'Dk_neg_intermediate': 0,
        # SCC check (b_{k-1} * Δ_k + ... form)
        'scc_fails': 0,
        # Minimum values at intermediate stages only
        'min_p2_intermediate': float('inf'),
        'min_Tk_intermediate': float('inf'),
        # Examples of intermediate failures
        'intermediate_examples': [],
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
            alpha_new = _polymul(E_acc, P)
            beta_new = _polymul(J_acc, Q)

            alphaQ = _polymul(E_acc, Q)
            betaQ = _polymul(J_acc, Q)
            alphaR = _polymul(E_acc, R)

            # Also compute I_acc and SCC
            I_acc = _polyadd(alpha_new, [0] + beta_new)

            is_final = (stage == s)
            is_intermediate = (stage < s)

            max_k = max(len(alpha_new), len(beta_new)) - 1

            for k in range(max_k):
                # P2: Δ_k = α_new_{k+1} · β_new_k − α_new_k · β_new_{k+1}
                delta_k = coeff(alpha_new, k+1) * coeff(beta_new, k) - \
                          coeff(alpha_new, k) * coeff(beta_new, k+1)

                # Main term: D_k
                D_k = coeff(alphaQ, k+1) * coeff(betaQ, k) - \
                      coeff(alphaQ, k) * coeff(betaQ, k+1)

                # Correction: T_k
                T_k = coeff(alphaR, k) * coeff(betaQ, k) - \
                      coeff(alphaR, k-1) * coeff(betaQ, k+1)

                if delta_k == 0 and D_k == 0 and T_k == 0:
                    continue

                result['total_checks'] += 1

                if delta_k < 0:
                    if is_intermediate:
                        result['p2_fail_intermediate'] += 1
                        if delta_k < result['min_p2_intermediate']:
                            result['min_p2_intermediate'] = delta_k
                        if len(result['intermediate_examples']) < 3:
                            result['intermediate_examples'].append({
                                'type': 'p2_inter', 'k': k, 'stage': stage,
                                'nfac': s, 'delta': delta_k, 'Dk': D_k, 'Tk': T_k,
                                'g6': g6.strip(),
                            })
                    else:
                        result['p2_fail_final'] += 1

                if T_k < 0:
                    if is_intermediate:
                        result['Tk_neg_intermediate'] += 1
                        if T_k < result['min_Tk_intermediate']:
                            result['min_Tk_intermediate'] = T_k
                    else:
                        result['Tk_neg_final'] += 1

                if D_k < 0:
                    if is_intermediate:
                        result['Dk_neg_intermediate'] += 1
                    else:
                        result['Dk_neg_final'] += 1

                # Also check SCC at this stage
                bkm1 = coeff(alpha_new, k - 1)  # Wait, SCC uses (I, E)
                # Let me compute it properly
                # a = I_acc, b = E_acc = alpha_new
                ak = coeff(I_acc, k)
                akp1 = coeff(I_acc, k + 1)
                akm1 = coeff(I_acc, k - 1)
                bk = coeff(alpha_new, k)
                bkp1 = coeff(alpha_new, k + 1)
                bkm1_scc = coeff(alpha_new, k - 1)

                if bkm1_scc > 0:
                    dk_scc = akp1 * bk - ak * bkp1
                    dkm1_scc = ak * bkm1_scc - akm1 * bk
                    ck_scc = bk * bk - bkm1_scc * bkp1
                    scc = bkm1_scc * dk_scc + bk * dkm1_scc + akm1 * ck_scc
                    if scc < 0:
                        result['scc_fails'] += 1

            E_acc = alpha_new
            J_acc = beta_new

    result['intermediate_examples'] = result['intermediate_examples'][:5]
    return result


def process_batch(g6_lines):
    keys = ['trees', 'total_checks', 'p2_fail_final', 'p2_fail_intermediate',
            'Tk_neg_final', 'Tk_neg_intermediate', 'Dk_neg_final', 'Dk_neg_intermediate',
            'scc_fails']
    batch = {k: 0 for k in keys}
    batch['min_p2_intermediate'] = float('inf')
    batch['min_Tk_intermediate'] = float('inf')
    batch['intermediate_examples'] = []

    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        for k in keys:
            if k != 'trees':
                batch[k] += r.get(k, 0)
        if r['min_p2_intermediate'] < batch['min_p2_intermediate']:
            batch['min_p2_intermediate'] = r['min_p2_intermediate']
        if r['min_Tk_intermediate'] < batch['min_Tk_intermediate']:
            batch['min_Tk_intermediate'] = r['min_Tk_intermediate']
        batch['intermediate_examples'].extend(r['intermediate_examples'])

    batch['intermediate_examples'] = batch['intermediate_examples'][:10]
    return batch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=500)
    args = parser.parse_args()

    totals = {
        'trees': 0, 'total_checks': 0,
        'p2_fail_final': 0, 'p2_fail_intermediate': 0,
        'Tk_neg_final': 0, 'Tk_neg_intermediate': 0,
        'Dk_neg_final': 0, 'Dk_neg_intermediate': 0,
        'scc_fails': 0,
    }
    global_min_p2_inter = float('inf')
    global_min_Tk_inter = float('inf')
    all_examples = []

    t0 = time.time()

    for nn in range(args.min_n, args.max_n + 1):
        tn = time.time()
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_totals = {k: 0 for k in totals}
        n_min_p2 = float('inf')
        n_min_Tk = float('inf')
        n_examples = []

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
                if r['min_p2_intermediate'] < n_min_p2:
                    n_min_p2 = r['min_p2_intermediate']
                if r['min_Tk_intermediate'] < n_min_Tk:
                    n_min_Tk = r['min_Tk_intermediate']
                n_examples.extend(r['intermediate_examples'])
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
                    if br['min_p2_intermediate'] < n_min_p2:
                        n_min_p2 = br['min_p2_intermediate']
                    if br['min_Tk_intermediate'] < n_min_Tk:
                        n_min_Tk = br['min_Tk_intermediate']
                    n_examples.extend(br['intermediate_examples'])

        proc.wait()
        elapsed_n = time.time() - tn

        for k in totals:
            totals[k] += n_totals[k]
        if n_min_p2 < global_min_p2_inter:
            global_min_p2_inter = n_min_p2
        if n_min_Tk < global_min_Tk_inter:
            global_min_Tk_inter = n_min_Tk
        all_examples.extend(n_examples)

        p2i = n_totals['p2_fail_intermediate']
        p2f = n_totals['p2_fail_final']
        tki = n_totals['Tk_neg_intermediate']
        tkf = n_totals['Tk_neg_final']

        print(f"n={nn:3d}: {n_totals['trees']:>10,} trees | "
              f"P2fail final:{p2f:>7,} inter:{p2i:>7,} | "
              f"Tk<0 final:{tkf:>7,} inter:{tki:>7,} | "
              f"SCC:{n_totals['scc_fails']}  [{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*90}")
    print(f"★ IDENTITY: FINAL vs INTERMEDIATE STAGE ANALYSIS")
    print(f"{'='*90}")
    print(f"Range:                     n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:               {totals['trees']:>15,}")
    print(f"Total checks:              {totals['total_checks']:>15,}")
    print()
    print(f"P2 failures (Δ_k < 0):")
    print(f"  At FINAL stage:          {totals['p2_fail_final']:>15,}")
    print(f"  At INTERMEDIATE stage:   {totals['p2_fail_intermediate']:>15,}")
    print()
    print(f"T_k < 0 events:")
    print(f"  At FINAL stage:          {totals['Tk_neg_final']:>15,}")
    print(f"  At INTERMEDIATE stage:   {totals['Tk_neg_intermediate']:>15,}")
    print()
    print(f"D_k < 0 events:")
    print(f"  At FINAL stage:          {totals['Dk_neg_final']:>15,}")
    print(f"  At INTERMEDIATE stage:   {totals['Dk_neg_intermediate']:>15,}")
    print()
    print(f"SCC failures:              {totals['scc_fails']:>15,}")
    print()
    print(f"Min P2 at intermediate:    {global_min_p2_inter}")
    print(f"Min T_k at intermediate:   {global_min_Tk_inter}")

    if all_examples:
        print(f"\nIntermediate P2 failure examples:")
        for i, ex in enumerate(all_examples[:15]):
            print(f"  [{i+1}] stage {ex['stage']}/{ex['nfac']} k={ex['k']}: "
                  f"Δ={ex['delta']}  D={ex['Dk']}  T={ex['Tk']}  {ex['g6']}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
