#!/usr/bin/env python3
"""Find the s=1 minimizer tree for the one-step decomposition."""
import sys
import subprocess
from indpoly import _polymul, _polyadd


def coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0

def delta_k(F, G, k):
    return coeff(F, k+1)*coeff(G, k) - coeff(F, k)*coeff(G, k+1)

def xshift(p):
    return [0] + list(p)

def parse_g6(g6):
    s = g6.strip()
    n = ord(s[0]) - 63
    adj = [[] for _ in range(n)]
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        for sh in range(5, -1, -1):
            bits.append((v >> sh) & 1)
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
    dp0, dp1s = {}, {}
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            pS, pE = [1], [1]
            for c in children[v]:
                sc = _polyadd(dp0[c], xshift(dp1s[c]))
                pS = _polymul(pS, sc)
                pE = _polymul(pE, dp0[c])
            dp0[v] = pS
            dp1s[v] = pE
    return children[root], dp0, dp1s


for nn in [14, 16, 18, 20]:
    cmd = ["/opt/homebrew/bin/geng", "-q", str(nn),
           "%d:%d" % (nn-1, nn-1), "-c"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    min_ratio = float("inf")
    min_info = None

    for line in proc.stdout:
        line = line.strip()
        if not line:
            continue
        n, adj = parse_g6(line)
        degree = [len(adj[v]) for v in range(n)]

        for r in range(n):
            leaf_ch = [u for u in adj[r] if degree[u] == 1]
            if not leaf_ch:
                continue
            children_r, dp0, dp1s = dp_rooted(n, adj, r)
            rl = []
            rnl = []
            for c in children_r:
                nbrs = [u for u in adj[c] if u != r]
                if len(nbrs) == 0:
                    rl.append(c)
                else:
                    rnl.append(c)
            if len(rnl) != 1:
                continue

            ell = len(rl)
            E_old = [1]
            for _ in range(ell):
                E_old = _polymul(E_old, [1, 1])
            J_old = [1]

            c = rnl[0]
            Ec, Jc = dp0[c], dp1s[c]
            A = _polymul(E_old, Ec)
            B = _polymul(E_old, Jc)
            C = _polymul(J_old, Ec)
            E_new = _polyadd(A, xshift(B))
            J_new = C
            md = max(len(E_new), len(J_new), len(A), len(B), len(C))

            for k in range(md):
                Dk = delta_k(A, C, k)
                Tk = coeff(B, k)*coeff(C, k) - coeff(B, k-1)*coeff(C, k+1)
                if Tk < 0:
                    ratio = Dk / abs(Tk)
                    if ratio < min_ratio:
                        min_ratio = ratio
                        min_info = {
                            "g6": line,
                            "root": r,
                            "ell": ell,
                            "k": k,
                            "D": Dk,
                            "T": Tk,
                            "ratio": ratio,
                            "Ec": list(Ec[:8]),
                            "Jc": list(Jc[:8]),
                            "degree": sorted([degree[v] for v in range(n)], reverse=True),
                        }

    proc.wait()
    if min_info:
        info = min_info
        print("n=%d: s=1 minimizer" % nn)
        print("  ratio = %.6f = (%d-2)/(%d-6) = %.6f" % (
            info["ratio"], nn, nn, (nn-2)/(nn-6)))
        print("  root=%d, ell=%d, k=%d" % (info["root"], info["ell"], info["k"]))
        print("  D=%d, T=%d, total=%d" % (info["D"], info["T"], info["D"]+info["T"]))
        print("  E_c = %s" % info["Ec"])
        print("  J_c = %s" % info["Jc"])
        print("  degree_seq = %s" % info["degree"])
        print("  g6 = %s" % info["g6"])
        print()
    else:
        print("n=%d: no s=1 negative T_k found" % nn)
