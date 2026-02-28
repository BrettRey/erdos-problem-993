"""Check Condition C at the factor level.

For each factor pair (I_c, E_c) at non-leaf children of support vertices:
  d_k = (I_c)_{k+1}*(E_c)_k - (I_c)_k*(E_c)_{k+1}
  c_k = (E_c)_k^2 - (E_c)_{k-1}*(E_c)_{k+1}

Condition C: for each k >= 1 with (E_c)_{k-1} > 0:
  (E_c)_{k-1}*d_k + (E_c)_k*d_{k-1} + (I_c)_{k-1}*c_k >= 0

Also check: does factor-level LC of E_c always hold?
"""

import subprocess
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
lc_failures = 0
condC_failures = 0
condC_examples = []

for nn in range(3, 21):
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)
    n_condC_fail_this = 0

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

            nonleaf_children = [c for c in children[r] if children[c]]
            for c in nonleaf_children:
                I_c = _polyadd(dp0[c], [0] + dp1s[c])
                E_c = dp0[c]
                total_factors += 1

                # Check LC of E_c
                deg = len(E_c) - 1
                for k in range(1, deg):
                    if E_c[k]**2 < E_c[k-1]*E_c[k+1]:
                        lc_failures += 1
                        break

                # Check Condition C at factor level
                d_seq = []
                for k in range(max(len(I_c), len(E_c)) + 1):
                    dk = coeff(I_c, k+1)*coeff(E_c, k) - coeff(I_c, k)*coeff(E_c, k+1)
                    d_seq.append(dk)

                factor_fail = False
                for k in range(1, len(d_seq)):
                    ekm1 = coeff(E_c, k-1)
                    if ekm1 == 0:
                        continue
                    ek = coeff(E_c, k)
                    ekp1 = coeff(E_c, k+1)
                    ikm1 = coeff(I_c, k-1)
                    ck = ek*ek - ekm1*ekp1
                    val = ekm1*d_seq[k] + ek*d_seq[k-1] + ikm1*ck
                    if val < 0:
                        factor_fail = True
                        break

                if factor_fail:
                    condC_failures += 1
                    n_condC_fail_this += 1
                    if len(condC_examples) < 5:
                        condC_examples.append((g6, n, r, c, I_c, E_c, d_seq))

    proc.wait()
    print(f"n={nn:2d}: factors={total_factors:>10,d} LC_fail={lc_failures} condC_fail={condC_failures} (this_n={n_condC_fail_this})")

print(f"\nTotal factors: {total_factors:,d}")
print(f"Factor LC failures: {lc_failures}")
print(f"Factor Condition C failures: {condC_failures}")

if condC_examples:
    print("\nCondition C failure examples:")
    for g6, n, r, c, I_c, E_c, d_seq in condC_examples:
        print(f"  g6={g6} n={n} root={r} child={c}")
        print(f"    I_c = {I_c}")
        print(f"    E_c = {E_c}")
        print(f"    d_k = {d_seq}")
