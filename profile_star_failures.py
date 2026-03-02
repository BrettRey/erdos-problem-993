"""
Profile the ⋆-failing SVs: what trees/structures trigger ⋆ failure?
Focus on the tightest w_k cases (smallest Term1 + Term2 + Term3).
"""

import subprocess
from itertools import permutations
from collections import defaultdict

try:
    from numpy import convolve as _np_convolve
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

_INT64_SAFE = 2**62


def _polymul_python(a, b):
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _polymul(a, b):
    if not a or not b:
        return []
    if _HAS_NUMPY:
        max_a = max(abs(x) for x in a)
        max_b = max(abs(x) for x in b)
        max_terms = min(len(a), len(b))
        if max_a > 0 and max_b > 0 and max_terms * max_a * max_b < _INT64_SAFE:
            return _np_convolve(a, b).astype(int).tolist()
    return _polymul_python(a, b)


def _polyadd(a, b):
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


def _coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0


def dp_all(n, adj, root):
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
            pS = [1]
            pE = [1]
            for c in children[v]:
                Ic = _polyadd(dp0[c], [0] + dp1s[c])
                pS = _polymul(pS, Ic)
                pE = _polymul(pE, dp0[c])
            dp0[v] = pS
            dp1s[v] = pE

    return children, dp0, dp1s


def find_mode(I_poly):
    m = 0
    for k in range(1, len(I_poly)):
        if I_poly[k] > I_poly[m]:
            m = k
    return m


def subtree_size(children, v):
    """BFS to count subtree size."""
    count = 0
    queue = [v]
    head = 0
    while head < len(queue):
        u = queue[head]; head += 1
        count += 1
        for c in children[u]:
            queue.append(c)
    return count


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode('ascii')]
    if data[0] <= 62:
        n = data[0]
        off = 1
    else:
        n = (data[1] << 12) | (data[2] << 6) | data[3]
        off = 4
    adj = [[] for _ in range(n)]
    bits = []
    for d in data[off:]:
        for shift in range(5, -1, -1):
            bits.append((d >> shift) & 1)
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                adj[i].append(j)
                adj[j].append(i)
            idx += 1
    return n, adj


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=16)
    parser.add_argument('--min-n', type=int, default=5)
    parser.add_argument('--geng', default='/opt/homebrew/bin/geng')
    parser.add_argument('--top-k', type=int, default=20,
                        help='Number of tightest w_k cases to report')
    args = parser.parse_args()

    import heapq

    # Collect tightest w_k cases (smallest w_k across all orderings)
    # Use a max-heap of size top_k (store negative for max-heap behavior)
    tightest = []  # (w_k, n, g6, root, s, ell, k, perm, step, term1, t23)

    for n in range(args.min_n, args.max_n + 1):
        cmd = [args.geng, '-q', str(n), f'{n-1}:{n-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode('ascii').strip()
            if not g6:
                continue
            nn, adj = parse_graph6(g6)
            deg = [len(adj[v]) for v in range(n)]

            for root in range(n):
                leaf_nbrs = [v for v in adj[root] if deg[v] == 1]
                if not leaf_nbrs:
                    continue

                children_all, dp0, dp1s = dp_all(n, adj, root)
                children_r = children_all[root]

                nonleaf_children = []
                ell = 0
                for c in children_r:
                    if len(adj[c]) == 1:
                        ell += 1
                    else:
                        nonleaf_children.append(c)

                s = len(nonleaf_children)
                if s < 2:
                    continue

                factors = [(dp0[c], dp1s[c]) for c in nonleaf_children]
                child_sizes = [subtree_size(children_all, c) for c in nonleaf_children]

                E_root = dp0[root]
                J_root = dp1s[root]
                I_poly = _polyadd(E_root, [0] + J_root)
                mode_I = find_mode(I_poly)
                if mode_I <= 1:
                    continue

                # Try natural ordering only (since w_k ≥ 0 at ALL orderings)
                perm = tuple(range(s))

                E_acc = [1]
                for _ in range(ell):
                    E_acc = _polymul(E_acc, [1, 1])
                J_acc = [1]

                for step_idx, ci in enumerate(perm):
                    g_c, h_c = factors[ci]

                    if step_idx == 0:
                        I_c = _polyadd(g_c, [0] + h_c)
                        E_acc = _polymul(E_acc, I_c)
                        J_acc = _polymul(J_acc, g_c)
                        continue

                    A = _polymul(E_acc, g_c)
                    B = _polymul(E_acc, h_c)
                    C = _polymul(J_acc, g_c)

                    for k in range(1, mode_I):
                        Ak = _coeff(A, k+1)
                        Akm = _coeff(A, k)
                        Ck = _coeff(C, k)
                        Ckm = _coeff(C, k-1)
                        Ckp = _coeff(C, k+1)
                        Bk = _coeff(B, k)
                        Bkm = _coeff(B, k-1)

                        dk_AC = Ak * Ck - Akm * Ckp
                        dkm1_BC = Bk * Ckm - Bkm * Ck
                        ck_C = Ck * Ck - Ckm * Ckp

                        Term1 = Ck * dk_AC
                        Term2 = Ckp * dkm1_BC
                        Term3 = Bk * ck_C
                        total = Term1 + Term2 + Term3

                        if Ck > 0:
                            w_k = total // Ck
                        else:
                            w_k = 0

                        star_margin = Bk * Ck - Bkm * Ckp

                        record = (w_k, n, g6, root, s, ell, k, step_idx,
                                  child_sizes, mode_I, Term1, Term2+Term3,
                                  star_margin)

                        if len(tightest) < args.top_k:
                            heapq.heappush(tightest, (-w_k, record))
                        elif -w_k > tightest[0][0]:
                            heapq.heapreplace(tightest, (-w_k, record))

                    I_c = _polyadd(g_c, [0] + h_c)
                    E_acc = _polymul(E_acc, I_c)
                    J_acc = _polymul(J_acc, g_c)

        proc.wait()
        print(f"n={n}: done")

    # Sort by w_k ascending
    results = sorted([r for _, r in tightest], key=lambda x: x[0])

    print()
    print("=" * 80)
    print(f"TOP {args.top_k} TIGHTEST w_k CASES")
    print("=" * 80)
    print(f"{'w_k':>6} | {'n':>3} | {'root':>4} | {'s':>2} | {'ell':>3} | "
          f"{'k':>3} | {'mode':>4} | {'step':>4} | {'T1':>10} | {'T2+T3':>10} | "
          f"{'⋆_margin':>10} | child_sizes")
    print("-" * 120)
    for rec in results:
        w_k, n, g6, root, s, ell, k, step, csizes, mode, t1, t23, sm = rec
        print(f"{w_k:>6} | {n:>3} | {root:>4} | {s:>2} | {ell:>3} | "
              f"{k:>3} | {mode:>4} | {step:>4} | {t1:>10} | {t23:>10} | "
              f"{sm:>10} | {csizes}")


if __name__ == '__main__':
    main()
