"""Check whether LC of E at support vertices follows from LC of subtree IS polys.

At a support vertex r with leaf children and non-leaf children c_1, ..., c_s:
  E = dp[r][0] = (1+x)^ell * prod_{j=1}^s I(T_{c_j})

Since (1+x) is LC and Cauchy products of nonneg LC polys are LC,
E is LC iff all I(T_{c_j}) are LC.

This script checks whether I(T_{c_j}) is LC for all subtrees at all support vertices
through n <= max_n. Since global LC of IS polys holds through n=25, this should always
hold for subtrees with < 26 vertices.

Key question: does LC of E at support vertices extend beyond n=25?
Could a non-LC tree (first at n=26) appear as a subtree T_{c_j}?
"""

import subprocess
import sys
import time

sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul, _polyadd

GENG = '/opt/homebrew/bin/geng'


def parse_g6(g6: str):
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
    if k < 0 or poly is None or k >= len(poly):
        return 0
    return poly[k]


def is_lc(poly):
    """Check if polynomial is log-concave (nonneg + c_k >= 0)."""
    for k in range(1, len(poly) - 1):
        if poly[k] * poly[k] < poly[k-1] * poly[k+1]:
            return False
    return True


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=22)
    parser.add_argument('--min-n', type=int, default=3)
    args = parser.parse_args()

    t0 = time.time()
    total_trees = 0
    total_support = 0
    total_factors = 0
    lc_E_failures = 0
    lc_factor_failures = 0
    max_subtree_size = 0

    for nn in range(args.min_n, args.max_n + 1):
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0

        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            n, adj = parse_g6(line)
            n_trees += 1

            leaf_count = [0] * n
            for v in range(n):
                for u in adj[v]:
                    if len(adj[u]) == 1:
                        leaf_count[v] += 1

            for root in range(n):
                if leaf_count[root] == 0:
                    continue
                total_support += 1

                dp0, dp1s, children = dp_rooted(n, adj, root)

                # E at root
                E = dp0[root]
                if not is_lc(E):
                    lc_E_failures += 1

                # Check each non-leaf child's subtree IS poly
                for c in children[root]:
                    if len(adj[c]) == 1:
                        continue  # leaf child, IS poly = 1+x, always LC
                    total_factors += 1
                    Ic = _polyadd(dp0[c], [0] + dp1s[c])

                    # Subtree size
                    # Count vertices in subtree of c
                    subtree_size = 0
                    stack = [c]
                    visited = set()
                    while stack:
                        v = stack.pop()
                        if v in visited:
                            continue
                        visited.add(v)
                        subtree_size += 1
                        for ch in children[v]:
                            stack.append(ch)
                    if subtree_size > max_subtree_size:
                        max_subtree_size = subtree_size

                    if not is_lc(Ic):
                        lc_factor_failures += 1
                        if lc_factor_failures <= 5:
                            print(f"FACTOR LC FAILURE: n={n}, root={root}, child={c}, "
                                  f"subtree_size={subtree_size}")
                            print(f"  I_c = {list(Ic)[:15]}...")

        proc.wait()
        total_trees += n_trees
        elapsed = time.time() - t0
        print(f"n={nn}: {n_trees} trees, E_lc_fails={lc_E_failures}, "
              f"factor_lc_fails={lc_factor_failures}, max_subtree={max_subtree_size}, "
              f"{elapsed:.1f}s", flush=True)

    print(f"\n{'='*60}")
    print(f"SUMMARY (n={args.min_n}..{args.max_n})")
    print(f"{'='*60}")
    print(f"Trees:              {total_trees}")
    print(f"Support vertices:   {total_support}")
    print(f"Non-leaf factors:   {total_factors}")
    print(f"Max subtree size:   {max_subtree_size}")
    print(f"LC(E) failures:     {lc_E_failures}")
    print(f"Factor LC failures: {lc_factor_failures}")
    print(f"\nConclusion: {'ALL factors LC => E is LC at all support vertices' if lc_factor_failures == 0 else 'FACTOR LC FAILURES FOUND'}")
    print(f"\nElapsed: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
