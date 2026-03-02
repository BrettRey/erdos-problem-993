"""Find the tightest tree quickly (single n, with parallel)."""

import subprocess
import sys
from multiprocessing import Pool

sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul, _polyadd

GENG = '/opt/homebrew/bin/geng'

def parse_g6(g6):
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
            factors.append((Ic, dp0[c], dp1s[c]))

        E_acc = [1]
        J_acc = [1]

        for stage, (P, Q, R) in enumerate(factors, 1):
            I_old = _polyadd(E_acc, [0] + J_acc)
            e_old = _polyadd(I_old, [0] + I_old)
            eQ = _polymul(e_old, Q)
            EQ = _polymul(E_acc, Q)
            ER = _polymul(E_acc, R)
            x1xER = _polyadd([0] + ER, [0, 0] + ER)
            xER = [0] + list(ER)
            E_new = _polymul(E_acc, P)
            J_new = _polymul(J_acc, Q)

            max_k = max(len(_polyadd(_polyadd(E_new, [0] + J_new),
                                      [0] + _polyadd(E_new, [0] + J_new))),
                        len(E_new)) - 1

            for k in range(max_k):
                QQ = coeff(eQ, k+1)*coeff(EQ, k) - coeff(eQ, k)*coeff(EQ, k+1)
                QR = coeff(eQ, k+1)*coeff(xER, k) - coeff(eQ, k)*coeff(xER, k+1)
                RQ = coeff(x1xER, k+1)*coeff(EQ, k) - coeff(x1xER, k)*coeff(EQ, k+1)
                RR = coeff(x1xER, k+1)*coeff(xER, k) - coeff(x1xER, k)*coeff(xER, k+1)
                cross = QR + RQ
                if cross < 0:
                    diag = QQ + RR
                    ratio = diag / abs(cross)
                    if best is None or ratio < best[0]:
                        # Describe the factor
                        best = (ratio, k, stage, s, root, n,
                                QQ, QR, RQ, RR, diag + cross, g6,
                                P[:12], Q[:12], R[:12],
                                E_acc[:12], J_acc[:12])
            E_acc = E_new
            J_acc = J_new
    return best

def process_batch(g6_lines):
    best = None
    for g6 in g6_lines:
        if not g6.strip():
            continue
        r = check_tree(g6)
        if r is not None:
            if best is None or r[0] < best[0]:
                best = r
    return best

if __name__ == '__main__':
  for nn in [20, 21, 22]:
    cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    best = None
    batch = []
    with Pool(8) as pool:
        futures = []
        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            batch.append(line)
            if len(batch) >= 500:
                futures.append(pool.apply_async(process_batch, (batch,)))
                batch = []
        if batch:
            futures.append(pool.apply_async(process_batch, (batch,)))
        for fut in futures:
            r = fut.get()
            if r is not None:
                if best is None or r[0] < best[0]:
                    best = r
    proc.wait()

    if best:
        ratio, k, stage, nfac, root, n, qq, qr, rq, rr, scc, g6, P, Q, R, Ea, Ja = best
        print(f"\nn={nn}: min diag/|cross| = {ratio:.6f}")
        print(f"  k={k}, stage={stage}, #fac={nfac}, root={root}")
        print(f"  QQ={qq}, QR={qr}, RQ={rq}, RR={rr}, SCC={scc}")
        print(f"  g6: {g6}")
        print(f"  Factor I_c (P): {P}")
        print(f"  Factor E_c (Q): {Q}")
        print(f"  Factor J_c (R): {R}")
        print(f"  E_acc:          {Ea}")
        print(f"  J_acc:          {Ja}")
        # Degree sequence
        deg = [len(adj[v]) for _, adj in [parse_g6(g6)] for v in range(len(adj))]
        print(f"  Degree seq: {sorted(deg, reverse=True)}")
