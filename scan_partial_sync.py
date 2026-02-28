"""Partial synchronicity scanner for tree IS polynomial factor pairs.

Tests the Hu-Wang-Zhao-Zhao partial synchronicity condition on (I_v, E_v)
for each tree T and each vertex v as root:

    I_v = dp[v][0] + dp[v][1]   (full IS polynomial, since v is root of T)
    E_v = dp[v][0]              (exclude-root polynomial)

Definition (Hu et al., arXiv:1507.08430):
Two sequences a, b are partially synchronized if for all m >= n >= 0:
    a_m * b_n + a_n * b_m >= a_{m+1} * b_{n-1} + a_{n-1} * b_{m+1}
where a_{-1} = b_{-1} = 0 and terms beyond degree are 0.

Usage:
    python3 scan_partial_sync.py --max-n 15
    python3 scan_partial_sync.py --max-n 15 --workers 4
"""

import argparse
import subprocess
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')

from indpoly import _polymul, _polyadd

GENG = '/opt/homebrew/bin/geng'


def parse_g6(g6: str) -> tuple[int, list[list[int]]]:
    """Parse graph6 string to (n, adjacency list)."""
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


def dp_rooted(n: int, adj: list[list[int]], root: int) -> tuple[list[int], list[int]]:
    """Compute DP at a given root, returning (dp0, dp1).

    dp0 = dp[root][0] (exclude root)
    dp1 = dp[root][1] (include root, with leading 0 for the x factor)

    I(T;x) = dp0 + dp1
    E(T;x) = dp0
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
    dp1 = [None] * n

    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod_S = [1]
            prod_E = [1]
            for c in children[v]:
                s_c = _polyadd(dp0[c], dp1[c])
                prod_S = _polymul(prod_S, s_c)
                prod_E = _polymul(prod_E, dp0[c])
            dp0[v] = prod_S
            dp1[v] = [0] + prod_E  # x * product of dp0[c]

    return dp0[root], dp1[root]


def coeff(poly: list[int], k: int) -> int:
    """Return coefficient at index k, with out-of-range convention."""
    if k < 0:
        return 0
    if k >= len(poly):
        return 0
    return poly[k]


def check_partial_sync(a: list[int], b: list[int]) -> list[tuple[int, int, int]]:
    """Check partial synchronicity of sequences a, b.

    Returns list of (m, n_idx, violation_amount) for each failure.
    Condition: for all m >= n >= 0,
        a_m * b_n + a_n * b_m >= a_{m+1} * b_{n-1} + a_{n-1} * b_{m+1}
    """
    max_idx = max(len(a), len(b)) - 1
    failures = []
    for m_idx in range(max_idx + 1):
        for n_idx in range(m_idx + 1):  # n <= m
            lhs = coeff(a, m_idx) * coeff(b, n_idx) + coeff(a, n_idx) * coeff(b, m_idx)
            rhs = coeff(a, m_idx + 1) * coeff(b, n_idx - 1) + coeff(a, n_idx - 1) * coeff(b, m_idx + 1)
            if lhs < rhs:
                failures.append((m_idx, n_idx, rhs - lhs))
    return failures


def edges_from_adj(adj: list[list[int]]) -> list[tuple[int, int]]:
    """Extract edge list from adjacency list."""
    edges = []
    for u, nbrs in enumerate(adj):
        for v in nbrs:
            if u < v:
                edges.append((u, v))
    return edges


def process_tree_g6(g6: str) -> dict:
    """Process a single tree: test partial sync at every rooting."""
    n, adj = parse_g6(g6)
    result = {
        'n': n,
        'g6': g6.strip(),
        'vertex_pairs': 0,
        'mn_pairs': 0,
        'failures': [],
    }

    for root in range(n):
        dp0, dp1 = dp_rooted(n, adj, root)
        I_v = _polyadd(dp0, dp1)  # full IS polynomial
        E_v = dp0                  # exclude-root polynomial

        fails = check_partial_sync(I_v, E_v)
        max_idx = max(len(I_v), len(E_v)) - 1
        # Number of (m, n) pairs: sum_{m=0}^{max_idx} (m+1) = (max_idx+1)(max_idx+2)/2
        result['mn_pairs'] += (max_idx + 1) * (max_idx + 2) // 2
        result['vertex_pairs'] += 1

        if fails:
            edges = edges_from_adj(adj)
            for (m_val, n_val, viol) in fails:
                result['failures'].append({
                    'edges': edges,
                    'root': root,
                    'm': m_val,
                    'n': n_val,
                    'violation': int(viol),
                })

    return result


def process_batch(g6_lines: list[str]) -> dict:
    """Process a batch of g6 strings."""
    batch_result = {
        'trees': 0,
        'vertex_pairs': 0,
        'mn_pairs': 0,
        'failures': [],
    }
    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = process_tree_g6(g6)
        batch_result['trees'] += 1
        batch_result['vertex_pairs'] += r['vertex_pairs']
        batch_result['mn_pairs'] += r['mn_pairs']
        batch_result['failures'].extend(r['failures'])
    return batch_result


def main():
    parser = argparse.ArgumentParser(description='Partial synchronicity scan for (I_v, E_v)')
    parser.add_argument('--max-n', type=int, default=15,
                        help='Maximum number of vertices (default: 15)')
    parser.add_argument('--min-n', type=int, default=3,
                        help='Minimum number of vertices (default: 3)')
    parser.add_argument('--workers', type=int, default=1,
                        help='Number of parallel workers (default: 1)')
    parser.add_argument('--batch-size', type=int, default=500,
                        help='Batch size for parallel processing (default: 500)')
    args = parser.parse_args()

    total_trees = 0
    total_vertex_pairs = 0
    total_mn_pairs = 0
    all_failures = []

    t0 = time.time()

    for n in range(args.min_n, args.max_n + 1):
        tn = time.time()
        # Generate trees: connected graphs with n-1 edges
        cmd = [GENG, '-q', str(n), f'{n-1}:{n-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_trees = 0
        n_vpairs = 0
        n_mnpairs = 0
        n_failures = 0

        if args.workers <= 1:
            # Sequential
            for line in proc.stdout:
                line = line.strip()
                if not line:
                    continue
                r = process_tree_g6(line)
                n_trees += 1
                n_vpairs += r['vertex_pairs']
                n_mnpairs += r['mn_pairs']
                if r['failures']:
                    n_failures += len(r['failures'])
                    all_failures.extend(r['failures'])
        else:
            # Parallel: collect in batches
            batch = []
            with Pool(args.workers) as pool:
                results_queue = []
                for line in proc.stdout:
                    line = line.strip()
                    if not line:
                        continue
                    batch.append(line)
                    if len(batch) >= args.batch_size:
                        results_queue.append(pool.apply_async(process_batch, (batch,)))
                        batch = []
                if batch:
                    results_queue.append(pool.apply_async(process_batch, (batch,)))

                for fut in results_queue:
                    br = fut.get()
                    n_trees += br['trees']
                    n_vpairs += br['vertex_pairs']
                    n_mnpairs += br['mn_pairs']
                    n_failures += len(br['failures'])
                    all_failures.extend(br['failures'])

        proc.wait()
        elapsed_n = time.time() - tn
        total_trees += n_trees
        total_vertex_pairs += n_vpairs
        total_mn_pairs += n_mnpairs

        status = "PASS" if n_failures == 0 else f"FAIL ({n_failures} violations)"
        print(f"n={n:3d}: {n_trees:>10,} trees, {n_vpairs:>12,} (tree,vertex) pairs, "
              f"{n_mnpairs:>15,} (m,n) pairs, {status}  [{elapsed_n:.1f}s]")

    elapsed = time.time() - t0

    print(f"\n{'='*72}")
    print(f"SUMMARY")
    print(f"{'='*72}")
    print(f"Vertex range:          n = {args.min_n} .. {args.max_n}")
    print(f"Total trees tested:    {total_trees:>15,}")
    print(f"Total (tree,v) pairs:  {total_vertex_pairs:>15,}")
    print(f"Total (m,n) checks:    {total_mn_pairs:>15,}")
    print(f"Partial sync failures: {len(all_failures):>15,}")
    print(f"Elapsed time:          {elapsed:>14.1f}s")

    if all_failures:
        print(f"\nFirst {min(5, len(all_failures))} failures:")
        for i, f in enumerate(all_failures[:5]):
            print(f"  [{i+1}] edges={f['edges']}, root={f['root']}, "
                  f"m={f['m']}, n={f['n']}, violation={f['violation']}")
    else:
        print("\nAll (I_v, E_v) pairs satisfy Hu-Wang-Zhao-Zhao partial synchronicity.")


if __name__ == '__main__':
    main()
