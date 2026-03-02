#!/usr/bin/env python3
"""Verify the leaf-augmentation LR dominance conjecture from GPT 5.2 Pro.

Conjecture: For every tree-realizable (I, E, J) at a support vertex,
define I_tilde = I + xE. Then I_tilde ratio-dominates E, i.e.,
  d_k(I_tilde, E) = (a_{k+1}+b_k)*b_k - (a_k+b_{k-1})*b_{k+1} >= 0
for all k.

Equivalently: g_k = c_k(E) + d_k(I,E) >= 0.

Tree interpretation: I + xE is the IS poly of the tree with one extra leaf
at the root. So "leaf-augmented tree ratio-dominates root-deleted forest."

Also test the binomial smoothing lemma: (1+x)B ratio-dominates B when B is LC.
"""
import sys
import subprocess
from collections import defaultdict
from indpoly import _polymul, _polyadd


def coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0

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
    max_n = 22

    # Leaf-augmentation: d_k(I+xE, E) >= 0
    la_checks = 0
    la_fails = 0
    min_gk = float('inf')  # min g_k = c_k(E) + d_k(I,E)
    min_gk_info = None

    # Binomial smoothing: (1+x)J ratio-dominates J (when J is LC)
    bs_checks = 0
    bs_fails = 0

    # Also: I+xE ratio-dominates (1+x)E? (stronger)
    stronger_checks = 0
    stronger_fails = 0

    for nn in range(3, max_n + 1):
        cmd = ["/opt/homebrew/bin/geng", "-q", str(nn),
               "%d:%d" % (nn-1, nn-1), "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0
        n_la_fail = 0
        n_bs_fail = 0
        n_stronger_fail = 0

        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            n, adj = parse_g6(line)
            n_trees += 1
            degree = [len(adj[v]) for v in range(n)]

            # Check at ALL rootings, not just support vertices
            for r in range(n):
                children_r, dp0, dp1s = dp_rooted(n, adj, r)

                # I = dp0[r] + x*dp1s[r] = E + xJ at root
                E = dp0[r]    # exclude root
                J = dp1s[r]   # include root / x
                I = _polyadd(E, xshift(J))

                # I_tilde = I + xE
                I_tilde = _polyadd(I, xshift(E))

                # (1+x)E
                one_plus_x_E = _polyadd(E, xshift(E))

                md = max(len(I_tilde), len(E), len(J), len(one_plus_x_E))

                for k in range(md):
                    # Leaf-augmentation: d_k(I_tilde, E) >= 0
                    # = I_tilde_{k+1} * E_k - I_tilde_k * E_{k+1}
                    gk = (coeff(I_tilde, k+1) * coeff(E, k) -
                          coeff(I_tilde, k) * coeff(E, k+1))
                    la_checks += 1
                    if gk < 0:
                        la_fails += 1
                        n_la_fail += 1
                    if gk < min_gk:
                        min_gk = gk
                        min_gk_info = (n, r, k, gk)

                    # Binomial smoothing: (1+x)J ≽ J
                    # d_k((1+x)J, J) = J_k^2 - J_{k-1}*J_{k+1} = c_k(J)
                    ck_J = coeff(J, k)**2 - coeff(J, k-1)*coeff(J, k+1)
                    bs_checks += 1
                    if ck_J < 0:
                        bs_fails += 1
                        n_bs_fail += 1

                    # Stronger: I+xE ≽ (1+x)E ?
                    dk_strong = (coeff(I_tilde, k+1) * coeff(one_plus_x_E, k) -
                                 coeff(I_tilde, k) * coeff(one_plus_x_E, k+1))
                    stronger_checks += 1
                    if dk_strong < 0:
                        stronger_fails += 1
                        n_stronger_fail += 1

        proc.wait()
        print("n=%2d: %8d trees | LA_fail=%d BS_fail=%d Strong_fail=%d" % (
            nn, n_trees, n_la_fail, n_bs_fail, n_stronger_fail), flush=True)

    print("\n" + "="*70)
    print("SUMMARY n=3..%d" % max_n)
    print("="*70)
    print("Leaf-augmentation (I+xE ≽ E):  %d checks, %d fails" % (la_checks, la_fails))
    print("Binomial smoothing (J is LC):   %d checks, %d fails" % (bs_checks, bs_fails))
    print("Stronger (I+xE ≽ (1+x)E):      %d checks, %d fails" % (stronger_checks, stronger_fails))
    if min_gk_info:
        info = min_gk_info
        print("\nMin g_k = %d at n=%d, root=%d, k=%d" % (info[3], info[0], info[1], info[2]))


if __name__ == '__main__':
    main()
