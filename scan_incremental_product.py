"""Incremental product scan for Strong Condition C.

For each tree at each support vertex, builds E^(k) and J^(k) incrementally
by adding one non-leaf child at a time. Checks at EVERY intermediate stage:
  1. Strong Condition C for (I^(k), E^(k))
  2. LC of E^(k)
  3. J^(k) <= E^(k) coefficientwise

Reports the tightest Condition C margin and which stage achieves it.

Usage:
    python3 scan_incremental_product.py --max-n 18 --workers 8
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
    """Compute rooted DP. Returns (dp0, dp1_shifted, children) for all vertices.

    dp0[v] = exclude-v polynomial
    dp1s[v] = include-v polynomial / x  (shifted, so dp1s[v][k] = coeff of x^k in dp1[v]/x)
    """
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
    dp1s = [None] * n  # dp1[v] / x

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


def check_strong_cond_c(a_poly, b_poly):
    """Check Strong Condition C for (I, E) where a=coeffs(I), b=coeffs(E).

    Returns (passes, min_margin, worst_k, min_nontrivial_margin, nontrivial_k, num_checks).
    min_margin: minimum over all k (including boundary zeros).
    min_nontrivial_margin: minimum over k where b_k > 0 (excluding trivial tail zeros).
    """
    max_k = max(len(a_poly), len(b_poly))
    min_margin = float('inf')
    worst_k = -1
    min_nt_margin = float('inf')  # nontrivial minimum
    nt_k = -1
    checks = 0

    for k in range(1, max_k):
        bkm1 = coeff(b_poly, k - 1)
        if bkm1 == 0:
            continue
        checks += 1

        ak = coeff(a_poly, k)
        akp1 = coeff(a_poly, k + 1)
        akm1 = coeff(a_poly, k - 1)
        bk = coeff(b_poly, k)
        bkp1 = coeff(b_poly, k + 1)

        dk = akp1 * bk - ak * bkp1
        dkm1 = ak * bkm1 - akm1 * bk
        ck = bk * bk - bkm1 * bkp1

        margin = bkm1 * dk + bk * dkm1 + akm1 * ck

        if margin < min_margin:
            min_margin = margin
            worst_k = k

        # Nontrivial: skip k where b_k = 0 (tail boundary, trivially 0)
        if bk > 0 and margin < min_nt_margin:
            min_nt_margin = margin
            nt_k = k

    if checks == 0:
        return True, float('inf'), -1, float('inf'), -1, 0

    return min_margin >= 0, min_margin, worst_k, min_nt_margin, nt_k, checks


def check_lc(poly):
    """Check log-concavity. Returns (is_lc, min_gap, worst_k)."""
    min_gap = float('inf')
    worst_k = -1
    for k in range(1, len(poly) - 1):
        gap = poly[k] * poly[k] - poly[k - 1] * poly[k + 1]
        if gap < min_gap:
            min_gap = gap
            worst_k = k
    if len(poly) <= 2:
        return True, float('inf'), -1
    return min_gap >= 0, min_gap, worst_k


def check_dominance(a_poly, b_poly):
    """Check a >= b coefficientwise. Returns (passes, worst_deficit)."""
    max_len = max(len(a_poly), len(b_poly))
    worst = float('inf')
    for k in range(max_len):
        diff = coeff(a_poly, k) - coeff(b_poly, k)
        if diff < worst:
            worst = diff
    return worst >= 0, worst


def process_tree(g6):
    n, adj = parse_g6(g6)
    result = {
        'n': n,
        'support_vertices': 0,
        'total_stages': 0,       # total intermediate stages checked
        'condC_fails': 0,
        'lc_fails': 0,
        'dom_fails': 0,          # J > E failures
        'min_condC_margin': float('inf'),       # includes boundary zeros
        'min_nt_margin': float('inf'),          # nontrivial (b_k > 0)
        'min_nt_stage': -1,
        'min_nt_nfactors': -1,
        'min_condC_stage': -1,
        'min_condC_nfactors': -1,
        'stage_dist': {},        # {(stage, nfactors): count}
        'examples': [],
    }

    # Find support vertices
    leaf_count = [0] * n
    for v in range(n):
        for u in adj[v]:
            if len(adj[u]) == 1:
                leaf_count[v] += 1

    for root in range(n):
        if leaf_count[root] == 0:
            continue
        result['support_vertices'] += 1

        dp0, dp1s, children = dp_rooted(n, adj, root)

        # Non-leaf children of root
        non_leaf_children = [c for c in children[root] if len(adj[c]) > 1]
        s = len(non_leaf_children)

        if s == 0:
            # All children are leaves: A=1, B=1, trivially Cond C
            result['total_stages'] += 1
            continue

        # Build factor pairs: (I_c, E_c) for each non-leaf child
        factors = []
        for c in non_leaf_children:
            Ic = _polyadd(dp0[c], [0] + dp1s[c])
            Ec = dp0[c]
            factors.append((Ic, Ec))

        # Incremental product
        E_acc = [1]  # E^(0) = 1
        J_acc = [1]  # J^(0) = 1

        for stage, (Ic, Ec) in enumerate(factors, 1):
            E_acc = _polymul(E_acc, Ic)
            J_acc = _polymul(J_acc, Ec)

            # I^(k) = E^(k) + x*J^(k)
            I_acc = _polyadd(E_acc, [0] + J_acc)

            result['total_stages'] += 1
            key = (stage, s)
            result['stage_dist'][key] = result['stage_dist'].get(key, 0) + 1

            # Check 1: Strong Condition C
            passes, margin, worst_k, nt_margin, nt_k, checks = check_strong_cond_c(I_acc, E_acc)
            if not passes:
                result['condC_fails'] += 1
                if len(result['examples']) < 5:
                    result['examples'].append({
                        'type': 'condC',
                        'root': root,
                        'stage': stage,
                        'nfactors': s,
                        'margin': int(margin),
                        'worst_k': worst_k,
                    })

            if margin < result['min_condC_margin']:
                result['min_condC_margin'] = margin
                result['min_condC_stage'] = stage
                result['min_condC_nfactors'] = s

            if nt_margin < result['min_nt_margin']:
                result['min_nt_margin'] = nt_margin
                result['min_nt_stage'] = stage
                result['min_nt_nfactors'] = s

            # Check 2: LC of E^(k)
            lc_ok, lc_gap, lc_k = check_lc(E_acc)
            if not lc_ok:
                result['lc_fails'] += 1
                if len(result['examples']) < 5:
                    result['examples'].append({
                        'type': 'lc',
                        'root': root,
                        'stage': stage,
                        'nfactors': s,
                        'lc_gap': int(lc_gap),
                        'lc_k': lc_k,
                    })

            # Check 3: J^(k) <= E^(k)
            dom_ok, dom_deficit = check_dominance(E_acc, J_acc)
            if not dom_ok:
                result['dom_fails'] += 1
                if len(result['examples']) < 5:
                    result['examples'].append({
                        'type': 'dom',
                        'root': root,
                        'stage': stage,
                        'nfactors': s,
                        'deficit': int(dom_deficit),
                    })

    return result


def process_batch(g6_lines):
    batch = {
        'trees': 0,
        'support_vertices': 0,
        'total_stages': 0,
        'condC_fails': 0,
        'lc_fails': 0,
        'dom_fails': 0,
        'min_condC_margin': float('inf'),
        'min_condC_info': None,
        'min_nt_margin': float('inf'),
        'min_nt_info': None,
        'examples': [],
    }
    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree(g6)
        batch['trees'] += 1
        batch['support_vertices'] += r['support_vertices']
        batch['total_stages'] += r['total_stages']
        batch['condC_fails'] += r['condC_fails']
        batch['lc_fails'] += r['lc_fails']
        batch['dom_fails'] += r['dom_fails']
        if r['min_condC_margin'] < batch['min_condC_margin']:
            batch['min_condC_margin'] = r['min_condC_margin']
            batch['min_condC_info'] = {
                'g6': g6.strip(),
                'n': r['n'],
                'stage': r['min_condC_stage'],
                'nfactors': r['min_condC_nfactors'],
            }
        if r['min_nt_margin'] < batch['min_nt_margin']:
            batch['min_nt_margin'] = r['min_nt_margin']
            batch['min_nt_info'] = {
                'g6': g6.strip(),
                'n': r['n'],
                'stage': r['min_nt_stage'],
                'nfactors': r['min_nt_nfactors'],
            }
        batch['examples'].extend(r['examples'][:2])
    return batch


def main():
    parser = argparse.ArgumentParser(description='Incremental product Condition C scan')
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=500)
    args = parser.parse_args()

    total_trees = 0
    total_sv = 0
    total_stages = 0
    total_condC_fails = 0
    total_lc_fails = 0
    total_dom_fails = 0
    global_min_margin = float('inf')
    global_min_info = None
    global_min_nt = float('inf')
    global_min_nt_info = None
    all_examples = []

    t0 = time.time()

    for nn in range(args.min_n, args.max_n + 1):
        tn = time.time()
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_trees = 0
        n_sv = 0
        n_stages = 0
        n_condC = 0
        n_lc = 0
        n_dom = 0
        n_min_margin = float('inf')
        n_min_info = None
        n_min_nt = float('inf')
        n_min_nt_info = None

        if args.workers <= 1:
            for line in proc.stdout:
                line = line.strip()
                if not line:
                    continue
                r = process_tree(line)
                n_trees += 1
                n_sv += r['support_vertices']
                n_stages += r['total_stages']
                n_condC += r['condC_fails']
                n_lc += r['lc_fails']
                n_dom += r['dom_fails']
                if r['min_nt_margin'] < n_min_nt:
                    n_min_nt = r['min_nt_margin']
                    n_min_nt_info = {
                        'g6': line,
                        'n': r['n'],
                        'stage': r['min_nt_stage'],
                        'nfactors': r['min_nt_nfactors'],
                    }
                if r['min_condC_margin'] < n_min_margin:
                    n_min_margin = r['min_condC_margin']
                    n_min_info = {
                        'g6': line,
                        'n': r['n'],
                        'stage': r['min_condC_stage'],
                        'nfactors': r['min_condC_nfactors'],
                    }
                all_examples.extend(r['examples'][:2])
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
                    n_sv += br['support_vertices']
                    n_stages += br['total_stages']
                    n_condC += br['condC_fails']
                    n_lc += br['lc_fails']
                    n_dom += br['dom_fails']
                    if br['min_condC_margin'] < n_min_margin:
                        n_min_margin = br['min_condC_margin']
                        n_min_info = br['min_condC_info']
                    if br['min_nt_margin'] < n_min_nt:
                        n_min_nt = br['min_nt_margin']
                        n_min_nt_info = br['min_nt_info']
                    all_examples.extend(br['examples'][:2])

        proc.wait()
        elapsed_n = time.time() - tn
        total_trees += n_trees
        total_sv += n_sv
        total_stages += n_stages
        total_condC_fails += n_condC
        total_lc_fails += n_lc
        total_dom_fails += n_dom
        if n_min_margin < global_min_margin:
            global_min_margin = n_min_margin
            global_min_info = n_min_info
        if n_min_nt < global_min_nt:
            global_min_nt = n_min_nt
            global_min_nt_info = n_min_nt_info

        condC_status = "PASS" if n_condC == 0 else f"FAIL ({n_condC})"
        lc_status = "PASS" if n_lc == 0 else f"FAIL ({n_lc})"
        dom_status = "PASS" if n_dom == 0 else f"FAIL ({n_dom})"

        nt_str = f"{n_min_nt}" if n_min_nt < float('inf') else "n/a"
        print(f"n={nn:3d}: {n_trees:>10,} trees, {n_stages:>12,} stages | "
              f"CondC: {condC_status:<12s} LC: {lc_status:<12s} J≤E: {dom_status:<10s} "
              f"min_nt_margin={nt_str}  [{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*80}")
    print(f"INCREMENTAL PRODUCT SCAN SUMMARY")
    print(f"{'='*80}")
    print(f"Vertex range:              n = {args.min_n} .. {args.max_n}")
    print(f"Total trees:               {total_trees:>15,}")
    print(f"Total support vertices:    {total_sv:>15,}")
    print(f"Total intermediate stages: {total_stages:>15,}")
    print(f"Condition C failures:      {total_condC_fails:>15,}")
    print(f"LC failures (E^(k)):       {total_lc_fails:>15,}")
    print(f"J^(k) > E^(k) failures:   {total_dom_fails:>15,}")
    print(f"Min nontrivial margin:     {global_min_nt}")
    if global_min_nt_info:
        print(f"  at n={global_min_nt_info['n']}, stage {global_min_nt_info['stage']}"
              f"/{global_min_nt_info['nfactors']}, g6={global_min_nt_info['g6']}")
    print(f"Elapsed time:              {elapsed:.1f}s")

    if all_examples:
        print(f"\nFirst {min(10, len(all_examples))} failure examples:")
        for i, ex in enumerate(all_examples[:10]):
            print(f"  [{i+1}] {ex}")


if __name__ == '__main__':
    main()
