"""P⋆ invariant scanner for tree independence polynomials.

Tests the two-state (E, J) invariant proposed by GPT 5.2 Pro (Feb 2027):
For a tree T rooted at r, let E = dp[r][0] and J = dp[r][1]/x (i.e.,
J = product of dp[c][0] over children c of r). Then I(T) = E + x*J.

P2 (prefix TP2):  e_{k+1} * j_k >= e_k * j_{k+1}  for k = 0,...,m-1
    Equivalently: ratio j_k/e_k is nonincreasing up to mode.
P3 (tail domination):  e_k >= j_{k-1}  for k >= m
    Equivalently: past the mode, E termwise dominates the shifted J.

where m = mode(I(T;x)) (the actual leftmost mode index).

Usage:
    python3 scan_pstar_invariant.py --max-n 15 --all-rootings
    python3 scan_pstar_invariant.py --max-n 20 --workers 8
"""

import argparse
import json
import subprocess
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

from indpoly import _polymul, _polyadd


def parse_g6(g6: str) -> tuple[int, list[list[int]]]:
    """Parse graph6 string to (n, adjacency list)."""
    s = g6.strip()
    idx = 0
    n = ord(s[idx]) - 63
    idx += 1
    adj = [[] for _ in range(n)]
    bit_pos = 0
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


def dp_rooted(n: int, adj: list[list[int]], root: int) -> tuple[list[int], list[int]]:
    """Compute DP at a given root, returning (E, J).

    E = dp[root][0] = IS polynomial of T with root excluded.
    J = product of dp[c][0] over children c = dp[root][1] / x.

    So I(T;x) = E(x) + x * J(x), and:
      [x^k] I(T) = E[k] + J[k-1]
    """
    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[root] = True
    queue = [root]
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    # Post-order traversal
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
    dp1_stripped = [None] * n  # dp[v][1] / x

    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1_stripped[v] = [1]
        else:
            # dp0[v] = product of (dp0[c] + dp1[c]) = product of S_c
            # dp1_stripped[v] = product of dp0[c]
            prod_S = [1]
            prod_E = [1]
            for c in children[v]:
                # S_c = dp0[c] + x * dp1_stripped[c]
                s_c = _polyadd(dp0[c], [0] + dp1_stripped[c])
                prod_S = _polymul(prod_S, s_c)
                prod_E = _polymul(prod_E, dp0[c])
            dp0[v] = prod_S
            dp1_stripped[v] = prod_E

    return dp0[root], dp1_stripped[root]


def get_mode(poly: list[int]) -> int:
    """Return the leftmost mode index (smallest k with max coefficient)."""
    max_val = max(poly)
    for k, v in enumerate(poly):
        if v == max_val:
            return k
    return 0


def get_coeff(poly: list[int], k: int) -> int:
    """Safely get coefficient at index k, 0 if out of range."""
    if 0 <= k < len(poly):
        return poly[k]
    return 0


def check_p2(E: list[int], J: list[int], m: int) -> tuple[bool, int, float]:
    """Check P2: e_{k+1}*j_k >= e_k*j_{k+1} for k = 0,...,m-1.

    Returns (passes, first_fail_k, worst_ratio).
    worst_ratio = max over k of (e_k * j_{k+1}) / (e_{k+1} * j_k).
    Values > 1.0 indicate failure.
    """
    first_fail = -1
    worst_ratio = 0.0
    for k in range(m):
        lhs = get_coeff(E, k + 1) * get_coeff(J, k)
        rhs = get_coeff(E, k) * get_coeff(J, k + 1)
        if lhs == 0 and rhs == 0:
            continue
        if lhs > 0:
            ratio = rhs / lhs
        elif rhs > 0:
            ratio = float('inf')
        else:
            ratio = 0.0
        if ratio > worst_ratio:
            worst_ratio = ratio
        if lhs < rhs:
            if first_fail == -1:
                first_fail = k
    return first_fail == -1, first_fail, worst_ratio


def check_p3(E: list[int], J: list[int], m: int) -> tuple[bool, int, float]:
    """Check P3: e_k >= j_{k-1} for k >= m.

    Returns (passes, first_fail_k, worst_ratio).
    worst_ratio = max over k of j_{k-1} / e_k. Values > 1.0 indicate failure.
    """
    first_fail = -1
    worst_ratio = 0.0
    max_k = max(len(E), len(J) + 1)
    for k in range(m, max_k):
        e_k = get_coeff(E, k)
        j_km1 = get_coeff(J, k - 1)
        if e_k == 0 and j_km1 == 0:
            continue
        if e_k > 0:
            ratio = j_km1 / e_k
        elif j_km1 > 0:
            ratio = float('inf')
        else:
            ratio = 0.0
        if ratio > worst_ratio:
            worst_ratio = ratio
        if e_k < j_km1:
            if first_fail == -1:
                first_fail = k
    return first_fail == -1, first_fail, worst_ratio


def process_tree(args):
    """Process a single tree: check P⋆ for one or all rootings."""
    g6, all_rootings = args
    n, adj = parse_g6(g6)
    if n <= 1:
        return None

    # Compute I(T) and mode
    # Use any rooting to get I(T)
    E0, J0 = dp_rooted(n, adj, 0)
    poly = _polyadd(E0, [0] + J0)
    m = get_mode(poly)

    rootings_to_check = range(n) if all_rootings else [0]
    best_result = None

    for r in rootings_to_check:
        if r == 0:
            E, J = E0, J0
        else:
            E, J = dp_rooted(n, adj, r)

        p2_ok, p2_fail_k, p2_ratio = check_p2(E, J, m)
        p3_ok, p3_fail_k, p3_ratio = check_p3(E, J, m)

        result = {
            'root': r,
            'p2_ok': p2_ok,
            'p3_ok': p3_ok,
            'p2_fail_k': p2_fail_k,
            'p3_fail_k': p3_fail_k,
            'p2_worst': p2_ratio,
            'p3_worst': p3_ratio,
            'both_ok': p2_ok and p3_ok,
        }

        if result['both_ok']:
            return {
                'g6': g6.strip(),
                'n': n,
                'm': m,
                'status': 'pass',
                'best_root': r,
                'rootings_checked': len(list(rootings_to_check)) if all_rootings else 1,
            }

        # Track best (fewest failures) for reporting
        if best_result is None or (result['p2_ok'] and not best_result['p2_ok']):
            best_result = result

    # All rootings failed
    return {
        'g6': g6.strip(),
        'n': n,
        'm': m,
        'status': 'fail',
        'best': best_result,
        'rootings_checked': n if all_rootings else 1,
    }


def run_scan(max_n: int, all_rootings: bool, workers: int):
    """Run P⋆ scan for trees up to max_n."""
    geng = '/opt/homebrew/bin/geng'
    results_by_n = {}
    failures = []
    total_trees = 0
    total_pass = 0
    t0 = time.time()

    for nn in range(3, max_n + 1):
        cmd = [geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        trees = []
        for line in proc.stdout:
            line = line.strip()
            if line:
                trees.append((line, all_rootings))
        proc.wait()

        n_trees = len(trees)
        n_pass = 0
        n_fail = 0
        n_p2_fail = 0
        n_p3_fail = 0
        worst_p2 = 0.0
        worst_p3 = 0.0
        fail_examples = []

        if workers > 1 and n_trees > 100:
            with Pool(workers) as pool:
                results = pool.map(process_tree, trees, chunksize=max(1, n_trees // (workers * 4)))
        else:
            results = [process_tree(t) for t in trees]

        for res in results:
            if res is None:
                continue
            total_trees += 1
            if res['status'] == 'pass':
                n_pass += 1
                total_pass += 1
            else:
                n_fail += 1
                b = res['best']
                if not b['p2_ok']:
                    n_p2_fail += 1
                if not b['p3_ok']:
                    n_p3_fail += 1
                if b['p2_worst'] > worst_p2:
                    worst_p2 = b['p2_worst']
                if b['p3_worst'] > worst_p3:
                    worst_p3 = b['p3_worst']
                if len(fail_examples) < 5:
                    fail_examples.append(res)
                failures.append(res)

        elapsed = time.time() - t0
        tag = "ALL-ROOT" if all_rootings else "ROOT-0"
        print(f"n={nn:2d}: {n_trees:>10,d} trees | "
              f"pass={n_pass:>10,d} fail={n_fail:>6,d} "
              f"(P2:{n_p2_fail} P3:{n_p3_fail}) | "
              f"worst P2={worst_p2:.4f} P3={worst_p3:.4f} | "
              f"{elapsed:.1f}s [{tag}]")

        results_by_n[nn] = {
            'trees': n_trees,
            'pass': n_pass,
            'fail': n_fail,
            'p2_fail': n_p2_fail,
            'p3_fail': n_p3_fail,
            'worst_p2': worst_p2,
            'worst_p3': worst_p3,
        }
        if fail_examples:
            results_by_n[nn]['examples'] = fail_examples

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"TOTAL: {total_trees:,d} trees, {total_pass:,d} pass, "
          f"{total_trees - total_pass:,d} fail ({elapsed:.1f}s)")

    if failures:
        print(f"\nFirst 10 failures:")
        for f in failures[:10]:
            b = f['best']
            p2_str = 'OK' if b['p2_ok'] else 'FAIL@' + str(b['p2_fail_k'])
            p3_str = 'OK' if b['p3_ok'] else 'FAIL@' + str(b['p3_fail_k'])
            print(f"  n={f['n']} m={f['m']} g6={f['g6'][:30]} "
                  f"P2={p2_str} P3={p3_str} "
                  f"(checked {f['rootings_checked']} rootings)")
    else:
        print("\nNo failures found!")

    # Save results
    out = {
        'max_n': max_n,
        'all_rootings': all_rootings,
        'total_trees': total_trees,
        'total_pass': total_pass,
        'total_fail': total_trees - total_pass,
        'by_n': results_by_n,
        'fail_count': len(failures),
        'first_failures': [
            {k: v for k, v in f.items() if k != 'best' or not isinstance(v, dict)}
            for f in failures[:50]
        ],
    }
    fname = f"results/pstar_scan_n{max_n}_{'allroot' if all_rootings else 'root0'}.json"
    with open(fname, 'w') as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"Results saved to {fname}")


def main():
    parser = argparse.ArgumentParser(description='P⋆ invariant scanner')
    parser.add_argument('--max-n', type=int, default=15,
                        help='Maximum tree size (default 15)')
    parser.add_argument('--all-rootings', action='store_true',
                        help='Check all rootings (slower but definitive)')
    parser.add_argument('--workers', type=int, default=1,
                        help='Number of parallel workers')
    args = parser.parse_args()
    run_scan(args.max_n, args.all_rootings, args.workers)


if __name__ == '__main__':
    main()
