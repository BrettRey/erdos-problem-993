"""
Find the exact n=18 tree where the Karlin rescue ratio = 10/3,
and extract its structure.
"""

import subprocess
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


def describe_tree(children_all, v, depth=0):
    """Recursive tree description."""
    kids = children_all[v]
    if not kids:
        return f"leaf({v})"
    kid_strs = [describe_tree(children_all, c, depth+1) for c in kids]
    return f"v{v}[" + ", ".join(kid_strs) + "]"


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
    n = 18
    target_ratio = Fraction(10, 3)

    cmd = ['/opt/homebrew/bin/geng', '-q', str(n), f'{n-1}:{n-1}', '-c']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    found = False
    for line in proc.stdout:
        g6 = line.decode('ascii').strip()
        if not g6:
            continue
        nn, adj = parse_graph6(g6)
        deg = [len(adj[v]) for v in range(n)]

        root = 0
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

        factors = [(dp0[c], dp1s[c]) for c in nonleaf_children]

        E_root = dp0[root]
        J_root = dp1s[root]
        I_poly = _polyadd(E_root, [0] + J_root)
        mode_I = find_mode(I_poly)

        E_acc = _polymul([1], [1, 1])  # (1+x)^1
        J_acc = [1]

        # Step 0: process first non-leaf child
        g0, h0 = factors[0]
        I_c0 = _polyadd(g0, [0] + h0)
        E_acc = _polymul(E_acc, I_c0)
        J_acc = _polymul(J_acc, g0)

        # Step 1: check second non-leaf child
        g1, h1 = factors[1]
        A = _polymul(E_acc, g1)
        B = _polymul(E_acc, h1)
        C = _polymul(J_acc, g1)

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
        dkm1_BC = Bk * Ckm - Bkm * Ck
        ck_C = Ck * Ck - Ckm * Ckp

        Term1 = Ck * dk_AC
        t23 = Ckp * dkm1_BC + Bk * ck_C

        if t23 >= 0 or Term1 <= 0:
            continue

        ratio = Fraction(Term1, abs(t23))
        if ratio == target_ratio:
            found = True
            print(f"FOUND at g6={g6}")
            print(f"  root=0, deg={deg[0]}, mode={mode_I}")
            print(f"  ell={ell}, s={s}")
            print(f"  Non-leaf children: {nonleaf_children}")
            for i, c in enumerate(nonleaf_children):
                sz = subtree_size(children_all, c)
                print(f"    child {i} (v{c}): size={sz}, deg_in_tree={deg[c]}")
                print(f"      E_c = {dp0[c]}")
                print(f"      J_c = {dp1s[c]}")
            print(f"  Tree structure: {describe_tree(children_all, 0)}")
            print()
            print(f"  E_acc = {E_acc[:8]}...")
            print(f"  J_acc = {J_acc[:8]}...")
            print(f"  g (E_child2) = {g1}")
            print(f"  h (J_child2) = {h1}")
            print()
            print(f"  At k={k}:")
            print(f"    A = {A[:8]}...")
            print(f"    B = {B[:8]}...")
            print(f"    C = {C[:8]}...")
            print(f"    Term1 = {Term1}")
            print(f"    Term2+Term3 = {t23}")
            print(f"    Ratio = {ratio} = {float(ratio):.6f}")
            print(f"    star_margin = {star_margin}")
            print(f"    w_k = {(Term1 + t23) // Ck}")
            print()
            print(f"  Adjacency list:")
            for v in range(n):
                print(f"    {v}: {adj[v]}")
            break

    proc.wait()
    if not found:
        print("Not found!")


if __name__ == '__main__':
    main()
