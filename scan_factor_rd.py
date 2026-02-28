"""Check factor-level ratio dominance: does I_c ratio-dominate E_c at each factor?

At a support vertex r with leaf child u:
  A = prod_{c: non-leaf children} I_c, B = prod_{c: non-leaf children} E_c
where I_c = dp0[c] + x*dp1s[c], E_c = dp0[c].

Factor-level ratio dominance: d_k^(c) = (I_c)_{k+1}*(E_c)_k - (I_c)_k*(E_c)_{k+1} >= 0 for all k.
"""

import subprocess
import sys
from indpoly import _polymul, _polyadd


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


def coeff(poly, k):
    return poly[k] if 0 <= k < len(poly) else 0


geng = '/opt/homebrew/bin/geng'

total_factors = 0
negative_factors = 0  # factors where any d_k^(c) < 0
negative_examples = []

for nn in range(3, 19):
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)

    for line in proc.stdout:
        g6 = line.strip()
        n, adj = parse_g6(g6)
        if n <= 2:
            continue

        leaf_count = [0] * n
        for v in range(n):
            for u in adj[v]:
                if len(adj[u]) == 1:
                    leaf_count[v] += 1

        for r in range(n):
            if leaf_count[r] == 0:
                continue

            # Build tree rooted at r
            children = [[] for _ in range(n)]
            visited = [False] * n
            visited[r] = True
            queue = [r]
            head = 0
            while head < len(queue):
                v = queue[head]; head += 1
                for u in adj[v]:
                    if not visited[u]:
                        visited[u] = True
                        children[v].append(u)
                        queue.append(u)

            order = []
            stack = [(r, False)]
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
                    prod_S = [1]
                    prod_E = [1]
                    for c in children[v]:
                        s_c = _polyadd(dp0[c], [0] + dp1s[c])
                        prod_S = _polymul(prod_S, s_c)
                        prod_E = _polymul(prod_E, dp0[c])
                    dp0[v] = prod_S
                    dp1s[v] = prod_E

            # Check each non-leaf child factor
            nonleaf_children = [c for c in children[r] if children[c]]
            for c in nonleaf_children:
                I_c = _polyadd(dp0[c], [0] + dp1s[c])  # I_c = E_c + x*J_c
                E_c = dp0[c]

                total_factors += 1
                has_neg = False
                for k in range(max(len(I_c), len(E_c))):
                    dk = coeff(I_c, k+1)*coeff(E_c, k) - coeff(I_c, k)*coeff(E_c, k+1)
                    if dk < 0:
                        has_neg = True
                        break

                if has_neg:
                    negative_factors += 1
                    if len(negative_examples) < 5:
                        negative_examples.append((g6, n, r, c, I_c, E_c))

    proc.wait()
    print(f"n={nn:2d}: factors checked so far: {total_factors:,d}, negative: {negative_factors:,d}")

print(f"\nTotal factors: {total_factors:,d}")
print(f"Factors with d_k^(c) < 0: {negative_factors:,d}")

if negative_examples:
    print("\nExamples:")
    for g6, n, r, c, I_c, E_c in negative_examples:
        print(f"  g6={g6} n={n} root={r} child={c}")
        print(f"    I_c = {I_c}")
        print(f"    E_c = {E_c}")
        dk_seq = []
        for k in range(max(len(I_c), len(E_c))):
            dk = coeff(I_c, k+1)*coeff(E_c, k) - coeff(I_c, k)*coeff(E_c, k+1)
            dk_seq.append(dk)
        print(f"    d_k = {dk_seq}")
