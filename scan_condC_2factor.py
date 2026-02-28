"""Check: does the product of two Condition-C pairs satisfy Condition C?

Sample random pairs (I_1, E_1) and (I_2, E_2) from the factor pool,
form (A, B) = (I_1*I_2, E_1*E_2), and check Condition C on the product.

This tests product closure directly.
"""

import subprocess
import random
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


def check_condC(I_poly, E_poly):
    """Check Condition C for pair (I, E). Returns True if passes."""
    d_seq = []
    L = max(len(I_poly), len(E_poly)) + 1
    for k in range(L):
        dk = coeff(I_poly, k+1)*coeff(E_poly, k) - coeff(I_poly, k)*coeff(E_poly, k+1)
        d_seq.append(dk)

    for k in range(1, L):
        ekm1 = coeff(E_poly, k-1)
        if ekm1 == 0:
            continue
        ek = coeff(E_poly, k)
        ekp1 = coeff(E_poly, k+1)
        ikm1 = coeff(I_poly, k-1)
        ck = ek*ek - ekm1*ekp1
        val = ekm1*d_seq[k] + ek*d_seq[k-1] + ikm1*ck
        if val < 0:
            return False
    return True


# Collect factor pairs from trees n=3..15
geng = '/opt/homebrew/bin/geng'
factor_pool = []

for nn in range(3, 16):
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
                factor_pool.append((I_c, E_c))
    proc.wait()

print(f"Collected {len(factor_pool)} factor pairs from n=3..15")

# Deduplicate
unique_factors = list(set((tuple(I), tuple(E)) for I, E in factor_pool))
unique_factors = [(list(I), list(E)) for I, E in unique_factors]
print(f"Unique factor pairs: {len(unique_factors)}")

# Check all pairwise products (for small pool) or sample
if len(unique_factors) <= 500:
    total = 0
    fails = 0
    for i in range(len(unique_factors)):
        for j in range(i, len(unique_factors)):
            I1, E1 = unique_factors[i]
            I2, E2 = unique_factors[j]
            A = _polymul(I1, I2)
            B = _polymul(E1, E2)
            total += 1
            if not check_condC(A, B):
                fails += 1
                print(f"  FAIL: factor1=({I1},{E1}) * factor2=({I2},{E2})")
                print(f"    A={A}")
                print(f"    B={B}")

    print(f"\nPairwise products checked: {total}")
    print(f"Condition C failures: {fails}")
else:
    # Sample
    random.seed(42)
    total = 0
    fails = 0
    for _ in range(100000):
        i, j = random.sample(range(len(unique_factors)), 2)
        I1, E1 = unique_factors[i]
        I2, E2 = unique_factors[j]
        A = _polymul(I1, I2)
        B = _polymul(E1, E2)
        total += 1
        if not check_condC(A, B):
            fails += 1
            if fails <= 5:
                print(f"  FAIL: factor1=({I1},{E1}) * factor2=({I2},{E2})")

    print(f"\nSampled products checked: {total}")
    print(f"Condition C failures: {fails}")
