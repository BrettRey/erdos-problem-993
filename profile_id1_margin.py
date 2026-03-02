#!/usr/bin/env python3
"""Profile the margin in the W form: term1/|term2| when d_{k-1}(B,J) < 0.

W: J_k * Delta_k(A,J) + J_{k+1} * d_{k-1}(B,J) >= 0
"""
import sys
import subprocess
from collections import defaultdict
from indpoly import _polymul, _polyadd


def coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0

def delta_k(F, G, k):
    return coeff(F, k+1) * coeff(G, k) - coeff(F, k) * coeff(G, k+1)

def d_km1(U, V, k):
    return coeff(U, k) * coeff(V, k-1) - coeff(U, k-1) * coeff(V, k)

def c_k_val(V, k):
    return coeff(V, k)**2 - coeff(V, k-1) * coeff(V, k+1)

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


def main():
    max_n = 20

    # Track min W margin by s and by n
    min_W_ratio_by_s = defaultdict(lambda: float('inf'))
    min_W_ratio_by_n = defaultdict(lambda: float('inf'))
    global_min_W = float('inf')
    global_min_info = None

    # Also track: when curvature=0, what's the Karlin/|bad| ratio?
    min_karlin_only_ratio = float('inf')

    d_neg_by_s = defaultdict(int)
    total_by_s = defaultdict(int)

    for nn in range(3, max_n + 1):
        cmd = ["/opt/homebrew/bin/geng", "-q", str(nn),
               "%d:%d" % (nn-1, nn-1), "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0

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
                    J = _polymul(J_old, Ec)
                    md = max(len(A), len(B), len(J)) + 1
                    sc += 1

                    for k in range(md):
                        total_by_s[sc] += 1
                        dkm1_BJ = d_km1(B, J, k)
                        if dkm1_BJ >= 0:
                            continue

                        d_neg_by_s[sc] += 1
                        Dk_AJ = delta_k(A, J, k)
                        Jk = coeff(J, k)
                        Jkp1 = coeff(J, k+1)
                        Bk = coeff(B, k)
                        ckJ = c_k_val(J, k)

                        term1 = Jk * Dk_AJ      # Karlin, >= 0
                        term2 = Jkp1 * dkm1_BJ   # bad, < 0
                        term3 = Bk * ckJ          # curvature, >= 0

                        # W ratio: term1/|term2|
                        if term2 != 0:
                            w_ratio = term1 / abs(term2)
                            if w_ratio < min_W_ratio_by_s[sc]:
                                min_W_ratio_by_s[sc] = w_ratio
                            if w_ratio < min_W_ratio_by_n[nn]:
                                min_W_ratio_by_n[nn] = w_ratio
                            if w_ratio < global_min_W:
                                global_min_W = w_ratio
                                global_min_info = (n, r, sc, k, term1, term2, term3, w_ratio)

                            # When curvature = 0
                            if ckJ == 0:
                                kr = term1 / abs(term2)
                                if kr < min_karlin_only_ratio:
                                    min_karlin_only_ratio = kr

                    E_old = _polyadd(A, xshift(B))
                    J_old = J

        proc.wait()
        mr = min_W_ratio_by_n.get(nn, float('inf'))
        ms = "%.4f" % mr if mr < float('inf') else "N/A"
        print("n=%2d: %8d trees | min W ratio=%s" % (nn, n_trees, ms), flush=True)

    print("\n" + "="*70)
    print("W margin: term1/|term2| when d_{k-1}(B,J) < 0")
    print("="*70)
    print("\nBy s-value:")
    for s in sorted(min_W_ratio_by_s.keys()):
        print("  s=%d: min ratio=%.6f, d<0 rate=%.2f%%" % (
            s, min_W_ratio_by_s[s],
            100 * d_neg_by_s[s] / max(1, total_by_s[s])))

    print("\nBy n:")
    for n in sorted(min_W_ratio_by_n.keys()):
        print("  n=%2d: %.6f" % (n, min_W_ratio_by_n[n]))

    if global_min_info:
        info = global_min_info
        print("\nGlobal min W ratio: %.6f" % info[7])
        print("  n=%d, root=%d, s=%d, k=%d" % (info[0], info[1], info[2], info[3]))
        print("  term1=%d (Karlin), term2=%d (bad LR), term3=%d (curvature)" % (
            info[4], info[5], info[6]))
        print("  IF total = %d" % (info[4] + info[5] + info[6]))

    print("\nMin Karlin/|bad| ratio when curvature=0: %.6f" % min_karlin_only_ratio)


if __name__ == '__main__':
    main()
