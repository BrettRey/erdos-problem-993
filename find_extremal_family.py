"""
Identify the extremal tree family for the Karlin rescue ratio.

From prior scans: extremal is ALWAYS s=2, ell=1, k=2, step=1 (first child
already processed, testing second child). The root has 1 leaf + 2 non-leaf
children. The two non-leaf subtrees have some structure.

Strategy: For each n, generate all trees and for each support vertex with
s=2, ell=1, find the minimum rescue ratio. Track the tree structure.

Uses a SINGLE geng call at a time (not competing with other processes).
"""

import subprocess
import sys
from fractions import Fraction

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


def subtree_size(children_all, v):
    count = 0
    queue = [v]
    head = 0
    while head < len(queue):
        u = queue[head]; head += 1
        count += 1
        for c in children_all[u]:
            queue.append(c)
    return count


def describe_subtree(children_all, v, depth=0):
    """Return a concise tree description."""
    kids = children_all[v]
    if not kids:
        return "L"
    kid_strs = sorted([describe_subtree(children_all, c, depth+1) for c in kids])
    return "(" + ",".join(kid_strs) + ")"


def scan_n(n):
    """Scan all trees at given n, find min rescue ratio at s=2, ell=1, k=2."""
    cmd = ['/opt/homebrew/bin/geng', '-q', str(n), f'{n-1}:{n-1}', '-c']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    min_ratio = None
    best = None
    tree_count = 0

    for line in proc.stdout:
        g6 = line.decode('ascii').strip()
        if not g6:
            continue
        tree_count += 1
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
            if s != 2 or ell != 1:
                continue

            # Process first non-leaf child
            g0, h0 = dp0[nonleaf_children[0]], dp1s[nonleaf_children[0]]
            I_c0 = _polyadd(g0, [0] + h0)

            E_acc = _polymul([1, 1], I_c0)  # (1+x) * I_c0
            J_acc = g0[:]  # E_c0

            # Check second non-leaf child at k=2
            g1, h1 = dp0[nonleaf_children[1]], dp1s[nonleaf_children[1]]
            A = _polymul(E_acc, g1)
            B = _polymul(E_acc, h1)
            C = _polymul(J_acc, g1)

            E_root = dp0[root]
            J_root = dp1s[root]
            I_poly = _polyadd(E_root, [0] + J_root)
            mode_I = find_mode(I_poly)
            if mode_I <= 2:
                continue

            k = 2
            Ck = _coeff(C, k)
            Ckm = _coeff(C, k-1)
            Ckp = _coeff(C, k+1)
            Bk = _coeff(B, k)
            Bkm = _coeff(B, k-1)

            star_margin = Bk * Ck - Bkm * Ckp
            if star_margin >= 0:
                continue

            Ak = _coeff(A, k+1)
            Akm = _coeff(A, k)
            dk_AC = Ak * Ck - Akm * Ckp
            Term1 = Ck * dk_AC
            t23 = Ckp * (Bk * Ckm - Bkm * Ck) + Bk * (Ck * Ck - Ckm * Ckp)

            if t23 >= 0 or Term1 <= 0:
                continue

            ratio = Fraction(Term1, abs(t23))
            if min_ratio is None or ratio < min_ratio:
                min_ratio = ratio
                sz0 = subtree_size(children_all, nonleaf_children[0])
                sz1 = subtree_size(children_all, nonleaf_children[1])
                desc0 = describe_subtree(children_all, nonleaf_children[0])
                desc1 = describe_subtree(children_all, nonleaf_children[1])
                best = {
                    'g6': g6,
                    'root': root,
                    'ratio': ratio,
                    'T1': Term1,
                    't23': t23,
                    'mode': mode_I,
                    'child_sizes': (sz0, sz1),
                    'child_types': (desc0, desc1),
                    'g0': dp0[nonleaf_children[0]],
                    'h0': dp1s[nonleaf_children[0]],
                    'g1': dp0[nonleaf_children[1]],
                    'h1': dp1s[nonleaf_children[1]],
                }

    proc.wait()
    return min_ratio, best, tree_count


def main():
    for n in range(9, 23):
        sys.stdout.write(f"n={n}: scanning... ")
        sys.stdout.flush()
        ratio, best, tc = scan_n(n)
        if ratio is None:
            print(f"no ⋆ failures (s=2,ell=1,k=2) in {tc} trees")
            continue
        print(f"min_ratio = {float(ratio):.6f} = {ratio}")
        print(f"  g6={best['g6']}, root={best['root']}, mode={best['mode']}")
        print(f"  child sizes: {best['child_sizes']}")
        print(f"  child structures: {best['child_types']}")
        print(f"  T1={best['T1']}, t23={best['t23']}")
        print(f"  g0 (E_child0) = {best['g0']}")
        print(f"  h0 (J_child0) = {best['h0']}")
        print(f"  g1 (E_child1) = {best['g1']}")
        print(f"  h1 (J_child1) = {best['h1']}")
        print()
        sys.stdout.flush()


if __name__ == '__main__':
    main()
