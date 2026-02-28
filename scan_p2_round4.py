"""Round 4 P2 scans: SV1, canonical max-ell, and Delta identity profile.

Three scans in a single pass over all trees:

1. SV1 (sign variation): For each support vertex, compute d_k (LR minors of A,B)
   and count sign changes in the prefix up to mode. SV1 = at most 1 sign change.

2. Canonical max-ell: For each tree, pick the vertex with the most leaf neighbors.
   Check P* only at that vertex.

3. Delta identity profile: Decompose Delta_k = d_k + memory + curvature and report
   minimum Delta_k and worst-case contributions.

Usage:
    python3 scan_p2_round4.py --max-n 18 --workers 8
    python3 scan_p2_round4.py --max-n 22 --workers 8  # ~15 min
"""

import argparse
import json
import subprocess
import sys
import time
from multiprocessing import Pool

from indpoly import _polymul, _polyadd


def parse_g6(g6: str) -> tuple[int, list[list[int]]]:
    """Parse graph6 string to (n, adjacency list)."""
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
    """Rooted DP returning (E, J) where E=dp[root][0], J=dp[root][1]/x."""
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
            prod_S = [1]
            prod_E = [1]
            for c in children[v]:
                s_c = _polyadd(dp0[c], [0] + dp1s[c])
                prod_S = _polymul(prod_S, s_c)
                prod_E = _polymul(prod_E, dp0[c])
            dp0[v] = prod_S
            dp1s[v] = prod_E

    return dp0[root], dp1s[root], children


def dp_rooted_with_AB(n, adj, root):
    """Rooted DP returning (E, J, A, B, children, ell).

    At a support vertex root with leaf children:
      E = (1+x)^ell * A
      J = B
    where A = prod_{non-leaf children} I_c, B = prod_{non-leaf children} E_c.
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
    dp1s = [None] * n
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            prod_S = [1]
            prod_E = [1]
            for c in children[v]:
                s_c = _polyadd(dp0[c], [0] + dp1s[c])
                prod_S = _polymul(prod_S, s_c)
                prod_E = _polymul(prod_E, dp0[c])
            dp0[v] = prod_S
            dp1s[v] = prod_E

    E = dp0[root]
    J = dp1s[root]

    # Count leaf children and compute A, B (excluding leaf children)
    leaf_children = [c for c in children[root] if not children[c]]
    nonleaf_children = [c for c in children[root] if children[c]]
    ell = len(leaf_children)

    # A = product over non-leaf children of I_c = dp0[c] + x*dp1s[c]
    # B = product over non-leaf children of E_c = dp0[c]
    A = [1]
    B = [1]
    for c in nonleaf_children:
        I_c = _polyadd(dp0[c], [0] + dp1s[c])
        A = _polymul(A, I_c)
        B = _polymul(B, dp0[c])

    return E, J, A, B, children, ell


def get_mode(poly):
    max_val = max(poly)
    for k, v in enumerate(poly):
        if v == max_val:
            return k
    return 0


def coeff(poly, k):
    return poly[k] if 0 <= k < len(poly) else 0


def check_p2(E, J, m):
    for k in range(m):
        if coeff(E, k+1) * coeff(J, k) < coeff(E, k) * coeff(J, k+1):
            return False, k
    return True, -1


def check_p3(E, J, m):
    max_k = max(len(E), len(J) + 1)
    for k in range(m, max_k):
        if coeff(E, k) < coeff(J, k - 1):
            return False, k
    return True, -1


def count_sign_changes(seq):
    """Count sign changes in a sequence, ignoring zeros."""
    changes = 0
    last_sign = 0
    for v in seq:
        if v > 0:
            s = 1
        elif v < 0:
            s = -1
        else:
            continue
        if last_sign != 0 and s != last_sign:
            changes += 1
        last_sign = s
    return changes


def sign_pattern(seq):
    """Return sign pattern string for display."""
    return ''.join('+' if v > 0 else ('-' if v < 0 else '0') for v in seq)


def process_tree(args):
    """Process one tree for all three scans."""
    g6 = args
    n, adj = parse_g6(g6)
    if n <= 2:
        return None

    # Get I(T) and mode using root 0
    E0, J0 = dp_rooted(n, adj, 0)[:2]
    poly = _polyadd(E0, [0] + J0)
    m = get_mode(poly)

    result = {
        'g6': g6.strip(),
        'n': n,
        'm': m,
    }

    # --- Scan 1: SV1 at all support vertices ---
    sv1_worst_changes = 0
    sv1_fails = 0
    sv1_support_checked = 0
    sv1_d_negative_count = 0  # support vertices where any d_k < 0

    # --- Scan 2: Canonical max-ell ---
    max_ell = 0
    max_ell_vertices = []

    # --- Scan 3: Delta identity (check sign only: b_{k-1}*Delta_k >= 0?) ---
    delta_min = 0  # 0 means no failure found; -1 means failure
    delta_min_info = None

    # First pass: find leaf-neighbor counts for all vertices
    leaf_count = [0] * n
    for v in range(n):
        for u in adj[v]:
            if len(adj[u]) == 1:
                leaf_count[v] += 1

    max_ell = max(leaf_count)
    max_ell_vertices = [v for v in range(n) if leaf_count[v] == max_ell]

    # Check each support vertex
    for r in range(n):
        if leaf_count[r] == 0:
            continue  # not a support vertex

        E, J, A, B, children, ell = dp_rooted_with_AB(n, adj, r)
        sv1_support_checked += 1

        # Compute d_k = a_{k+1}*b_k - a_k*b_{k+1} for k = 0,...,m-1
        d_seq = []
        for k in range(m):
            dk = coeff(A, k+1) * coeff(B, k) - coeff(A, k) * coeff(B, k+1)
            d_seq.append(dk)

        changes = count_sign_changes(d_seq)
        if changes > sv1_worst_changes:
            sv1_worst_changes = changes

        has_negative = any(dk < 0 for dk in d_seq)
        if has_negative:
            sv1_d_negative_count += 1

        if changes > 1:
            sv1_fails += 1

        # Delta profile: Delta_k = d_k + memory + curvature
        # Check sign using integer arithmetic: b_{k-1}*Delta_k = b_{k-1}*d_k + b_k*d_{k-1} + a_{k-1}*c_k
        for k in range(1, m):
            bkm1 = coeff(B, k - 1)
            if bkm1 == 0:
                continue  # skip degenerate

            dk = d_seq[k] if k < len(d_seq) else 0
            dkm1 = d_seq[k - 1] if k - 1 < len(d_seq) else 0
            bk = coeff(B, k)
            bkp1 = coeff(B, k + 1)
            akm1 = coeff(A, k - 1)
            ck = bk * bk - bkm1 * bkp1  # LC gap of B

            delta_times_bkm1 = bkm1 * dk + bk * dkm1 + akm1 * ck

            if delta_times_bkm1 < 0:
                delta_min = -1
                delta_min_info = {
                    'root': r, 'k': k, 'ell': ell,
                    'dk': dk, 'dkm1': dkm1, 'ck': ck,
                    'bk': bk, 'bkm1': bkm1, 'akm1': akm1,
                }

    # Scan 2: check P* at the canonical max-ell vertex
    canonical_r = max_ell_vertices[0]
    E_c, J_c = dp_rooted(n, adj, canonical_r)[:2]
    p2_ok, _ = check_p2(E_c, J_c, m)
    p3_ok, _ = check_p3(E_c, J_c, m)

    result['sv1_worst_changes'] = sv1_worst_changes
    result['sv1_fails'] = sv1_fails
    result['sv1_support_checked'] = sv1_support_checked
    result['sv1_d_negative_count'] = sv1_d_negative_count
    result['max_ell'] = max_ell
    result['canonical_p2'] = p2_ok
    result['canonical_p3'] = p3_ok
    result['canonical_pstar'] = p2_ok and p3_ok
    result['delta_fail'] = (delta_min < 0)
    result['delta_min_info'] = delta_min_info

    return result


def run_scan(max_n, workers):
    geng = '/opt/homebrew/bin/geng'
    t0 = time.time()

    # Aggregates
    total_trees = 0
    sv1_total_fails = 0
    sv1_total_support = 0
    sv1_total_d_negative = 0
    sv1_max_changes = 0
    canonical_total_fail = 0
    ell_dist = {}
    delta_total_fail = 0
    delta_fail_examples = []

    sv1_fail_examples = []
    canonical_fail_examples = []

    for nn in range(3, max_n + 1):
        cmd = [geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        trees = [line.strip() for line in proc.stdout if line.strip()]
        proc.wait()

        n_trees = len(trees)
        if workers > 1 and n_trees > 100:
            with Pool(workers) as pool:
                results = pool.map(process_tree, trees,
                                   chunksize=max(1, n_trees // (workers * 4)))
        else:
            results = [process_tree(t) for t in trees]

        n_sv1_fail = 0
        n_canonical_fail = 0
        n_max_changes = 0

        for res in results:
            if res is None:
                continue
            total_trees += 1

            # SV1
            sv1_total_support += res['sv1_support_checked']
            sv1_total_d_negative += res['sv1_d_negative_count']
            if res['sv1_fails'] > 0:
                sv1_total_fails += 1
                n_sv1_fail += 1
                if len(sv1_fail_examples) < 10:
                    sv1_fail_examples.append(res)
            if res['sv1_worst_changes'] > n_max_changes:
                n_max_changes = res['sv1_worst_changes']
            if res['sv1_worst_changes'] > sv1_max_changes:
                sv1_max_changes = res['sv1_worst_changes']

            # Canonical max-ell
            ell = res['max_ell']
            ell_dist[ell] = ell_dist.get(ell, 0) + 1
            if not res['canonical_pstar']:
                canonical_total_fail += 1
                n_canonical_fail += 1
                if len(canonical_fail_examples) < 10:
                    canonical_fail_examples.append(res)

            # Delta identity
            if res['delta_fail']:
                delta_total_fail += 1
                if len(delta_fail_examples) < 10:
                    delta_fail_examples.append(res)

        n_delta_fail = sum(1 for r in results if r and r['delta_fail'])

        elapsed = time.time() - t0
        print(f"n={nn:2d}: {n_trees:>10,d} trees | "
              f"SV1 fail={n_sv1_fail} max_sc={n_max_changes} | "
              f"canon fail={n_canonical_fail} | "
              f"delta fail={n_delta_fail} | "
              f"{elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"TOTAL: {total_trees:,d} trees ({elapsed:.1f}s)")
    print(f"\n--- SV1 (sign variation of d_k) ---")
    print(f"  Support vertices checked: {sv1_total_support:,d}")
    print(f"  With any d_k < 0: {sv1_total_d_negative:,d}")
    print(f"  Max sign changes seen: {sv1_max_changes}")
    print(f"  Trees with SV1 violation (>1 change): {sv1_total_fails}")
    if sv1_fail_examples:
        print(f"  First failures:")
        for ex in sv1_fail_examples[:5]:
            print(f"    n={ex['n']} g6={ex['g6'][:30]} changes={ex['sv1_worst_changes']}")

    print(f"\n--- Canonical max-ell rooting ---")
    print(f"  P* failures at max-ell vertex: {canonical_total_fail}")
    print(f"  ell distribution: {dict(sorted(ell_dist.items()))}")
    if canonical_fail_examples:
        print(f"  First failures:")
        for ex in canonical_fail_examples[:5]:
            print(f"    n={ex['n']} g6={ex['g6'][:30]} ell={ex['max_ell']} "
                  f"P2={ex['canonical_p2']} P3={ex['canonical_p3']}")

    print(f"\n--- Delta identity (Condition C) ---")
    print(f"  b_{{k-1}}*Delta_k < 0 failures: {delta_total_fail}")
    if delta_fail_examples:
        print(f"  First failures:")
        for ex in delta_fail_examples[:5]:
            info = ex['delta_min_info']
            print(f"    n={ex['n']} g6={ex['g6'][:30]} root={info['root']} k={info['k']}")

    # Save results
    out = {
        'max_n': max_n,
        'total_trees': total_trees,
        'elapsed': elapsed,
        'sv1': {
            'support_checked': sv1_total_support,
            'd_negative_count': sv1_total_d_negative,
            'max_sign_changes': sv1_max_changes,
            'trees_with_sv1_violation': sv1_total_fails,
        },
        'canonical': {
            'pstar_failures': canonical_total_fail,
            'ell_distribution': dict(sorted(ell_dist.items())),
        },
        'delta_identity': {
            'condition_c_failures': delta_total_fail,
        },
    }
    fname = f"results/round4_scan_n{max_n}.json"
    with open(fname, 'w') as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\nResults saved to {fname}")


def main():
    parser = argparse.ArgumentParser(description='Round 4 P2 scans')
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    run_scan(args.max_n, args.workers)


if __name__ == '__main__':
    main()
