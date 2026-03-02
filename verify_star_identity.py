"""Verify Instance 1's ★ identity and Y_k monotonicity at incremental product stages.

At each stage the pair (E_acc, J_acc) evolves multiplicatively:
    E_acc_new = E_acc * I_c   (convolution with P = I_c)
    J_acc_new = J_acc * E_c   (convolution with Q = E_c)

With P = Q + xR where R = J_c (= dp1s[c]), the P2 condition
    Δ_k(E, J) = E_{k+1} · J_k − E_k · J_{k+1}
decomposes as:
    Δ_k^new = D_k + T_k

where:
    D_k = Δ_k(E*Q, J*Q) ≥ 0  (by TP2 closure: old P2 + Q is PF2)
    T_k = (E*R)_k · (J*Q)_k − (E*R)_{k-1} · (J*Q)_{k+1}

Key checks:
    1. Identity: Δ_k^new == D_k + T_k  (algebraic, should be exact)
    2. P2 at all intermediate stages: Δ_k^new ≥ 0
    3. D_k ≥ 0 (TP2 main term)
    4. T_k ≥ 0 (bridge lemma conjecture)
    5. Y_k = (E*R)_k / (J*Q)_{k+1} nondecreasing

Usage:
    python3 verify_star_identity.py --max-n 18 --workers 8
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
        'p2_fails': 0,
        'Dk_neg': 0,
        'Tk_neg': 0,
        'Yk_not_nondec': 0,
        'min_Dk': float('inf'),
        'min_Tk': float('inf'),
        'min_p2': float('inf'),
        # Tightest T_k failures (if any)
        'tightest_Tk': [],  # (T_k, k, stage, nfactors, g6)
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
            Ic = _polyadd(dp0[c], [0] + dp1s[c])   # P = I_c
            Ec = dp0[c]                              # Q = E_c
            Rc = dp1s[c]                             # R = J_c (so P = Q + xR)
            factors.append((Ic, Ec, Rc))

        E_acc = [1]  # α = E_acc
        J_acc = [1]  # β = J_acc

        for stage, (P, Q, R) in enumerate(factors, 1):
            # Old pair: (E_acc, J_acc) = (α, β)
            # New pair: α_new = α*P, β_new = β*Q

            alpha_new = _polymul(E_acc, P)  # E_acc_new
            beta_new = _polymul(J_acc, Q)   # J_acc_new

            # Compute auxiliary convolutions
            alphaQ = _polymul(E_acc, Q)     # E_acc * Q
            betaQ = _polymul(J_acc, Q)      # J_acc * Q
            alphaR = _polymul(E_acc, R)     # E_acc * R

            max_k = max(len(alpha_new), len(beta_new)) - 1

            # Track Y_k values for monotonicity check
            Y_values = []  # (k, Y_k) where Y_k = (E*R)_k / (J*Q)_{k+1}

            for k in range(max_k):
                # P2 at new stage: Δ_k = α_new_{k+1} · β_new_k − α_new_k · β_new_{k+1}
                delta_k = coeff(alpha_new, k+1) * coeff(beta_new, k) - \
                          coeff(alpha_new, k) * coeff(beta_new, k+1)

                # Main term: D_k = Δ_k(α*Q, β*Q)
                D_k = coeff(alphaQ, k+1) * coeff(betaQ, k) - \
                      coeff(alphaQ, k) * coeff(betaQ, k+1)

                # Correction: T_k = (α*R)_k · (β*Q)_k − (α*R)_{k-1} · (β*Q)_{k+1}
                T_k = coeff(alphaR, k) * coeff(betaQ, k) - \
                      coeff(alphaR, k-1) * coeff(betaQ, k+1)

                # Skip trivial zeros
                if delta_k == 0 and D_k == 0 and T_k == 0:
                    continue

                result['total_checks'] += 1

                # Check 1: Identity
                if delta_k != D_k + T_k:
                    result['identity_fails'] += 1

                # Check 2: P2
                if delta_k < 0:
                    result['p2_fails'] += 1
                if delta_k < result['min_p2']:
                    result['min_p2'] = delta_k

                # Check 3: D_k ≥ 0
                if D_k < 0:
                    result['Dk_neg'] += 1
                if D_k < result['min_Dk']:
                    result['min_Dk'] = D_k

                # Check 4: T_k ≥ 0
                if T_k < 0:
                    result['Tk_neg'] += 1
                    result['tightest_Tk'].append(
                        (T_k, k, stage, s, g6.strip(), D_k, delta_k)
                    )
                if T_k < result['min_Tk']:
                    result['min_Tk'] = T_k

                # Track Y_k for monotonicity
                denom = coeff(betaQ, k+1)
                if denom > 0:
                    Y_k = coeff(alphaR, k) / denom
                    Y_values.append((k, Y_k))

            # Check 5: Y_k nondecreasing
            for i in range(1, len(Y_values)):
                if Y_values[i][1] < Y_values[i-1][1] - 1e-12:
                    result['Yk_not_nondec'] += 1

            # Update accumulated pair
            E_acc = alpha_new
            J_acc = beta_new

    # Keep only 10 tightest T_k
    result['tightest_Tk'].sort(key=lambda x: x[0])
    result['tightest_Tk'] = result['tightest_Tk'][:10]

    return result


def process_batch(g6_lines):
    batch = {
        'trees': 0,
        'total_checks': 0,
        'identity_fails': 0,
        'p2_fails': 0,
        'Dk_neg': 0,
        'Tk_neg': 0,
        'Yk_not_nondec': 0,
        'min_Dk': float('inf'),
        'min_Tk': float('inf'),
        'min_p2': float('inf'),
        'tightest_Tk': [],
    }
    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        batch['total_checks'] += r['total_checks']
        batch['identity_fails'] += r['identity_fails']
        batch['p2_fails'] += r['p2_fails']
        batch['Dk_neg'] += r['Dk_neg']
        batch['Tk_neg'] += r['Tk_neg']
        batch['Yk_not_nondec'] += r['Yk_not_nondec']
        if r['min_Dk'] < batch['min_Dk']:
            batch['min_Dk'] = r['min_Dk']
        if r['min_Tk'] < batch['min_Tk']:
            batch['min_Tk'] = r['min_Tk']
        if r['min_p2'] < batch['min_p2']:
            batch['min_p2'] = r['min_p2']
        batch['tightest_Tk'].extend(r['tightest_Tk'])

    batch['tightest_Tk'].sort(key=lambda x: x[0])
    batch['tightest_Tk'] = batch['tightest_Tk'][:20]
    return batch


def main():
    parser = argparse.ArgumentParser(description='Verify ★ identity and Y_k monotonicity')
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=500)
    args = parser.parse_args()

    totals = {
        'trees': 0, 'total_checks': 0,
        'identity_fails': 0, 'p2_fails': 0,
        'Dk_neg': 0, 'Tk_neg': 0, 'Yk_not_nondec': 0,
    }
    global_min_Dk = float('inf')
    global_min_Tk = float('inf')
    global_min_p2 = float('inf')
    global_tightest = []

    t0 = time.time()

    for nn in range(args.min_n, args.max_n + 1):
        tn = time.time()
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_totals = {k: 0 for k in totals}
        n_min_Dk = float('inf')
        n_min_Tk = float('inf')
        n_min_p2 = float('inf')
        n_tightest = []

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
                if r['min_Dk'] < n_min_Dk:
                    n_min_Dk = r['min_Dk']
                if r['min_Tk'] < n_min_Tk:
                    n_min_Tk = r['min_Tk']
                if r['min_p2'] < n_min_p2:
                    n_min_p2 = r['min_p2']
                n_tightest.extend(r['tightest_Tk'])
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
                    if br['min_Dk'] < n_min_Dk:
                        n_min_Dk = br['min_Dk']
                    if br['min_Tk'] < n_min_Tk:
                        n_min_Tk = br['min_Tk']
                    if br['min_p2'] < n_min_p2:
                        n_min_p2 = br['min_p2']
                    n_tightest.extend(br['tightest_Tk'])

        proc.wait()
        elapsed_n = time.time() - tn

        for k in totals:
            totals[k] += n_totals[k]

        if n_min_Dk < global_min_Dk:
            global_min_Dk = n_min_Dk
        if n_min_Tk < global_min_Tk:
            global_min_Tk = n_min_Tk
        if n_min_p2 < global_min_p2:
            global_min_p2 = n_min_p2

        n_tightest.sort(key=lambda x: x[0])
        n_tightest = n_tightest[:20]
        global_tightest.extend(n_tightest)
        global_tightest.sort(key=lambda x: x[0])
        global_tightest = global_tightest[:50]

        tk_str = f"{n_min_Tk}" if n_min_Tk < float('inf') else "n/a"
        dk_str = f"{n_min_Dk}" if n_min_Dk < float('inf') else "n/a"

        print(f"n={nn:3d}: {n_totals['trees']:>10,} trees | "
              f"id:{n_totals['identity_fails']}  "
              f"P2:{n_totals['p2_fails']}  "
              f"Dk<0:{n_totals['Dk_neg']}  "
              f"Tk<0:{n_totals['Tk_neg']}  "
              f"Yk↓:{n_totals['Yk_not_nondec']}  "
              f"minTk={tk_str}  [{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*90}")
    print(f"★ IDENTITY AND BRIDGE LEMMA VERIFICATION")
    print(f"{'='*90}")
    print(f"Range:                    n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:              {totals['trees']:>15,}")
    print(f"Total (stage,k) checks:   {totals['total_checks']:>15,}")
    print()
    print(f"Identity failures (Δ ≠ D+T):  {totals['identity_fails']:>10,}")
    print(f"P2 failures (Δ_k < 0):        {totals['p2_fails']:>10,}")
    print(f"D_k < 0 events:               {totals['Dk_neg']:>10,}")
    print(f"T_k < 0 events:               {totals['Tk_neg']:>10,}")
    print(f"Y_k non-monotone events:       {totals['Yk_not_nondec']:>10,}")
    print()
    print(f"Min D_k: {global_min_Dk}")
    print(f"Min T_k: {global_min_Tk}")
    print(f"Min P2 (Δ_k): {global_min_p2}")

    if global_tightest:
        print(f"\n{'='*90}")
        print(f"TIGHTEST T_k (most negative or smallest positive)")
        print(f"{'='*90}")
        print(f"{'Rank':>4s}  {'T_k':>14s}  {'k':>3s}  {'Stage':>5s}  {'#fac':>4s}  "
              f"{'D_k':>14s}  {'Δ_k':>14s}  g6")
        for rank, (Tk, k, stage, nfac, g6, Dk, deltak) in enumerate(global_tightest[:30]):
            print(f"{rank+1:4d}  {Tk:14d}  {k:3d}  {stage:5d}  {nfac:4d}  "
                  f"{Dk:14d}  {deltak:14d}  {g6}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
