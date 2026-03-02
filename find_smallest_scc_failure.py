#!/usr/bin/env python3
"""Find the smallest tree where SCC fails at a support vertex.

SCC (Strong Condition C) at rooted tree (T, r):
  E = dp[r][0]           (exclude-root IS poly)
  J = dp[r][1] / x       (include-root IS poly, shifted down by 1)
  I = E + x*J = full IS poly  (equivalently I = E + dp1)
  Check: (1+x)*I ≽ E, i.e. all coefficients of (1+x)*I - E are >= 0.

A support vertex is one adjacent to at least one leaf (degree-1 vertex).

Uses multiprocessing with geng's res/mod partitioning for parallelism.
"""

import subprocess
import sys
import time
from multiprocessing import Pool

# ---------------------------------------------------------------------------
# Inline graph6 parser (avoid import issues in workers)
# ---------------------------------------------------------------------------

def parse_graph6(line: bytes):
    s = line.strip()
    if s.startswith(b">>graph6<<"):
        s = s[10:]
    pos = 0
    if s[pos] < 126:
        n = s[pos] - 63
        pos += 1
    elif s[pos] == 126 and s[pos + 1] < 126:
        n = 0
        for i in range(1, 4):
            n = (n << 6) | (s[pos + i] - 63)
        pos += 4
    else:
        n = 0
        for i in range(2, 8):
            n = (n << 6) | (s[pos + i] - 63)
        pos += 8
    bits = []
    for i in range(pos, len(s)):
        val = s[i] - 63
        for shift in (5, 4, 3, 2, 1, 0):
            bits.append((val >> shift) & 1)
    adj = [[] for _ in range(n)]
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                adj[i].append(j)
                adj[j].append(i)
            idx += 1
    return n, adj


# ---------------------------------------------------------------------------
# Polynomial arithmetic (inlined for worker pickling)
# ---------------------------------------------------------------------------

try:
    from numpy import convolve as _np_convolve
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

_INT64_SAFE = 2**62

def _polymul(a, b):
    if not a or not b:
        return []
    if _HAS_NUMPY:
        max_a = max(abs(x) for x in a)
        max_b = max(abs(x) for x in b)
        max_terms = min(len(a), len(b))
        if max_a > 0 and max_b > 0 and max_terms * max_a * max_b < _INT64_SAFE:
            return _np_convolve(a, b).astype(int).tolist()
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out

def _polyadd(a, b):
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def dp_rooted(n, adj, root):
    """Compute (dp0[root], dp1[root]) for tree rooted at root."""
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
            prod_sum = [1]
            prod_exc = [1]
            for c in children[v]:
                summand = _polyadd(dp0[c], dp1[c])
                prod_sum = _polymul(prod_sum, summand)
                prod_exc = _polymul(prod_exc, dp0[c])
            dp0[v] = prod_sum
            dp1[v] = [0] + prod_exc

    return dp0[root], dp1[root]


def check_scc(E, dp1):
    """Check (1+x)*I ≽ E where I = E + dp1.

    Returns (passes, failing_k, delta_k).
    """
    I = _polyadd(E, dp1)
    lenI = len(I)
    # (1+x)*I coefficients: prod[k] = I[k] + I[k-1]
    prod = [0] * (lenI + 1)
    for k in range(lenI):
        prod[k] += I[k]
        prod[k + 1] += I[k]

    lenE = len(E)
    max_len = max(len(prod), lenE)
    for k in range(max_len):
        p = prod[k] if k < len(prod) else 0
        e = E[k] if k < lenE else 0
        delta = p - e
        if delta < 0:
            return False, k, delta
    return True, -1, 0


def find_support_vertices(n, adj):
    """Return list of support vertices (adjacent to at least one degree-1 vertex)."""
    deg = [len(adj[v]) for v in range(n)]
    support = []
    for v in range(n):
        for u in adj[v]:
            if deg[u] == 1:
                support.append(v)
                break
    return support


def worker(args):
    """Worker: enumerate and check one geng partition.

    Returns (count, failure_or_None).
    failure = dict with tree info if SCC fails, else None.
    """
    n, res, mod = args
    geng = "/opt/homebrew/bin/geng"
    edges = n - 1
    cmd = [geng, "-cq", str(n), f"{edges}:{edges}", f"{res}/{mod}"]

    count = 0
    with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            tree_n, adj = parse_graph6(line)
            count += 1
            support = find_support_vertices(tree_n, adj)
            for r in support:
                E, dp1 = dp_rooted(tree_n, adj, r)
                passes, fail_k, delta_k = check_scc(E, dp1)
                if not passes:
                    J = dp1[1:]
                    I = _polyadd(E, dp1)
                    proc.terminate()
                    return count, {
                        "n": tree_n,
                        "adj": adj,
                        "root": r,
                        "E": E,
                        "J": J,
                        "I": I,
                        "fail_k": fail_k,
                        "delta_k": delta_k,
                        "degree_seq": sorted([len(adj[v]) for v in range(tree_n)], reverse=True),
                    }
    return count, None


def main():
    WORKERS = 10

    for n in range(23, 28):
        print(f"\n{'='*60}")
        print(f"Searching n = {n} with {WORKERS} workers")
        print(f"{'='*60}")
        sys.stdout.flush()

        t0 = time.time()
        tasks = [(n, res, WORKERS) for res in range(WORKERS)]

        total_count = 0
        failure = None

        with Pool(WORKERS) as pool:
            results = pool.map(worker, tasks)

        for cnt, result in results:
            total_count += cnt
            if result is not None and failure is None:
                failure = result

        elapsed = time.time() - t0

        if failure is not None:
            f = failure
            print(f"\n*** SCC FAILURE FOUND ***")
            print(f"  n = {f['n']}")
            print(f"  root (support vertex) = {f['root']}")
            print(f"  degree sequence: {f['degree_seq']}")
            print(f"  adjacency: {f['adj']}")
            print(f"  E = {f['E']}")
            print(f"  J = {f['J']}")
            print(f"  I = {f['I']}")
            print(f"  failing k = {f['fail_k']}, delta = {f['delta_k']}")
            print(f"  checked {total_count:,} trees total in {elapsed:.1f}s")
            print(f"\nFirst SCC failure found at n = {n}.")
            sys.stdout.flush()
            return
        else:
            rate = total_count / elapsed if elapsed > 0 else 0
            print(f"  Completed n = {n}: {total_count:,} trees, 0 SCC failures, "
                  f"{elapsed:.1f}s ({rate:.0f} trees/s)")
            sys.stdout.flush()

    print(f"\nNo SCC failures found for n = 23..27.")


if __name__ == "__main__":
    main()
