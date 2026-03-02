#!/usr/bin/env python3
"""Quick verification of one-step decomposition for n=17,18."""
import sys
import subprocess
from collections import defaultdict
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


for nn in [17, 18]:
    cmd = ["/opt/homebrew/bin/geng", "-q", str(nn),
           f"{nn-1}:{nn-1}", "-c"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    n_trees = 0
    n_T_neg = 0
    n_T = 0
    min_ratio_n = float("inf")
    min_ratio_s = defaultdict(lambda: float("inf"))
    min_info = None

    for line in proc.stdout:
        line = line.strip()
        if not line:
            continue
        n, adj = parse_g6(line)
        n_trees += 1
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

            ell = len(rl)
            E_old = [1]
            for _ in range(ell):
                E_old = _polymul(E_old, [1, 1])
            J_old = [1]
            sc = 0

            for c in rnl:
                Ec, Jc = dp0[c], dp1s[c]
                A = _polymul(E_old, Ec)
                B = _polymul(E_old, Jc)
                C = _polymul(J_old, Ec)
                E_new = _polyadd(A, xshift(B))
                J_new = C
                md = max(len(E_new), len(J_new), len(A), len(B), len(C))
                sc += 1

                for k in range(md):
                    Dk = delta_k(A, C, k)
                    Tk = coeff(B, k)*coeff(C, k) - coeff(B, k-1)*coeff(C, k+1)
                    n_T += 1
                    if Tk < 0:
                        n_T_neg += 1
                        ratio = Dk / abs(Tk)
                        if ratio < min_ratio_s[sc]:
                            min_ratio_s[sc] = ratio
                        if ratio < min_ratio_n:
                            min_ratio_n = ratio
                            min_info = (n, r, sc, k, Dk, Tk, ratio)

                E_old, J_old = E_new, J_new

    proc.wait()
    ms = f"{min_ratio_n:.6f}" if min_ratio_n < float("inf") else "N/A"
    print(f"n={nn}: {n_trees:>8d} trees | T<0: {n_T_neg:>8d}/{n_T:>9d} | min D/|T|={ms}", flush=True)
    svals = sorted(min_ratio_s.keys())
    print("  by s: " + ", ".join(f"s={s}:{min_ratio_s[s]:.4f}" for s in svals), flush=True)
    if min_info:
        print(f"  global min: n={min_info[0]}, r={min_info[1]}, s={min_info[2]}, k={min_info[3]}, "
              f"D={min_info[4]}, T={min_info[5]}, ratio={min_info[6]:.6f}", flush=True)
