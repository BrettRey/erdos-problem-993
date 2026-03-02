"""Find the exact tree achieving min (QQ+RR)/|cross| when cross < 0.

Quick scan to identify the tree for algebraic analysis.
"""

import subprocess
import sys

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
    if poly is None or k < 0 or k >= len(poly):
        return 0
    return poly[k]


def check_tree(g6):
    n, adj = parse_g6(g6)
    best = None

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
            Ic = _polyadd(dp0[c], [0] + dp1s[c])
            Ec = dp0[c]
            Rc = dp1s[c]
            factors.append((Ic, Ec, Rc, c))

        E_acc = [1]
        J_acc = [1]

        for stage, (P, Q, R, cv) in enumerate(factors, 1):
            I_old = _polyadd(E_acc, [0] + J_acc)
            e_old = _polyadd(I_old, [0] + I_old)

            eQ = _polymul(e_old, Q)
            EQ = _polymul(E_acc, Q)
            ER = _polymul(E_acc, R)
            x1xER = _polyadd([0] + ER, [0, 0] + ER)
            xER = [0] + list(ER)

            E_new = _polymul(E_acc, P)
            J_new = _polymul(J_acc, Q)
            I_new = _polyadd(E_new, [0] + J_new)
            e_new = _polyadd(I_new, [0] + I_new)

            max_k = max(len(e_new), len(E_new)) - 1

            for k in range(max_k):
                QQ = coeff(eQ, k+1) * coeff(EQ, k) - coeff(eQ, k) * coeff(EQ, k+1)
                QR = coeff(eQ, k+1) * coeff(xER, k) - coeff(eQ, k) * coeff(xER, k+1)
                RQ = coeff(x1xER, k+1) * coeff(EQ, k) - coeff(x1xER, k) * coeff(EQ, k+1)
                RR = coeff(x1xER, k+1) * coeff(xER, k) - coeff(x1xER, k) * coeff(xER, k+1)

                cross = QR + RQ
                diag = QQ + RR

                if cross < 0:
                    ratio = diag / abs(cross)
                    if best is None or ratio < best[0]:
                        best = (ratio, k, stage, s, root, cv,
                                QQ, QR, RQ, RR, diag + cross,
                                list(P), list(Q), list(R),
                                list(E_acc), list(J_acc))

            E_acc = E_new
            J_acc = J_new

    return best


# Scan n=17..22 to find the global minimum
for nn in range(17, 23):
    cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    best_n = None
    best_g6 = None

    for line in proc.stdout:
        line = line.strip()
        if not line:
            continue
        r = check_tree(line)
        if r is not None:
            if best_n is None or r[0] < best_n[0]:
                best_n = r
                best_g6 = line

    proc.wait()

    if best_n:
        ratio, k, stage, nfac, root, cv, qq, qr, rq, rr, scc, P, Q, R, E_acc, J_acc = best_n
        print(f"\nn={nn}: min diag/|cross| = {ratio:.4f}")
        print(f"  k={k}, stage={stage}, #fac={nfac}, root={root}")
        print(f"  QQ={qq}, QR={qr}, RQ={rq}, RR={rr}, SCC={scc}")
        print(f"  g6: {best_g6}")
        print(f"  Factor P (I_c): {P[:8]}...")
        print(f"  Factor Q (E_c): {Q[:8]}...")
        print(f"  Factor R (J_c): {R[:8]}...")
        print(f"  E_acc before:   {E_acc[:8]}...")
        print(f"  J_acc before:   {J_acc[:8]}...")
