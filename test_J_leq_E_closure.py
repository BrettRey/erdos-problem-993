"""
Test whether J' <= E' = E1*E2 coefficientwise, where:
  I_i = E_i + x*J_i  with  J_i <= E_i coefficientwise
  J' = J1*E2 + E1*J2 + x*J1*J2
  E' = E1*E2

Part 1: Random nonneg integer polynomials (100,000 trials)
Part 2: Actual tree DP factors from indpoly.py
Part 3: Verify J <= E at single-subtree level
"""

import random
import sys
import os
import subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from indpoly import _polymul, _polyadd


def poly_shift(p):
    """Multiply polynomial by x (prepend a 0)."""
    return [0] + list(p)


def poly_sub(a, b):
    """a - b coefficientwise."""
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] -= b[i]
    return out


def check_leq(a, b):
    """Check a <= b coefficientwise."""
    la, lb = len(a), len(b)
    for i in range(max(la, lb)):
        ai = a[i] if i < la else 0
        bi = b[i] if i < lb else 0
        if ai > bi:
            return False
    return True


def random_pair(max_deg=6, max_coeff=10):
    """Generate (E, J) with J <= E coefficientwise, nonneg integers."""
    deg = random.randint(0, max_deg)
    E = [random.randint(0, max_coeff) for _ in range(deg + 1)]
    E[0] = max(E[0], 1)
    J = [random.randint(0, E[i]) for i in range(deg + 1)]
    return E, J


def test_closure(E1, J1, E2, J2):
    """
    Test J' <= E' where:
      E' = E1*E2
      J' = J1*E2 + E1*J2 + x*J1*J2
    Returns (passes, E', J', diff).
    """
    Eprime = _polymul(E1, E2)

    term1 = _polymul(J1, E2)
    term2 = _polymul(E1, J2)
    term3 = poly_shift(_polymul(J1, J2))

    Jprime = _polyadd(_polyadd(term1, term2), term3)

    diff = poly_sub(Eprime, Jprime)
    passes = all(d >= 0 for d in diff)

    return passes, Eprime, Jprime, diff


def parse_graph6(s):
    """Parse a graph6 string into (n, adjacency list)."""
    s = s.strip()
    if not s:
        return None, None
    data = [ord(c) - 63 for c in s]
    n = data[0]
    bits = []
    for byte_val in data[1:]:
        for bit in range(5, -1, -1):
            bits.append((byte_val >> bit) & 1)
    adj = [[] for _ in range(n)]
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                adj[i].append(j)
                adj[j].append(i)
            idx += 1
    return n, adj


def tree_dp(nn, adj, root):
    """Run DP on tree rooted at root. Returns (dp0, dp1, children)."""
    children_r = [[] for _ in range(nn)]
    visited = [False] * nn
    visited[root] = True
    queue = [root]
    head = 0
    while head < len(queue):
        v = queue[head]; head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                children_r[v].append(u)
                queue.append(u)

    order = []
    stack = [(root, False)]
    while stack:
        v, processed = stack.pop()
        if processed:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children_r[v]:
            stack.append((c, False))

    dp0 = [[] for _ in range(nn)]
    dp1 = [[] for _ in range(nn)]

    for v in order:
        if not children_r[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod = [1]
            for c in children_r[v]:
                summand = _polyadd(dp0[c], dp1[c])
                prod = _polymul(prod, summand)
            dp0[v] = prod

            prod = [1]
            for c in children_r[v]:
                prod = _polymul(prod, dp0[c])
            dp1[v] = [0] + prod

    return dp0, dp1, children_r


# ===== Part 1: Random polynomial pairs =====
print("=" * 60)
print("Part 1: Random polynomial pairs (100,000 trials)")
print("=" * 60)

N_TRIALS = 100_000
failures = 0
worst_violation = 0
worst_example = None

random.seed(42)

for trial in range(N_TRIALS):
    E1, J1 = random_pair()
    E2, J2 = random_pair()

    passes, Ep, Jp, diff = test_closure(E1, J1, E2, J2)

    if not passes:
        failures += 1
        min_diff = min(diff)
        if min_diff < worst_violation:
            worst_violation = min_diff
            worst_example = (E1, J1, E2, J2, Ep, Jp, diff)
        if failures <= 5:
            print(f"\nCounterexample #{failures}:")
            print(f"  E1 = {E1}, J1 = {J1}")
            print(f"  E2 = {E2}, J2 = {J2}")
            print(f"  E' = {Ep}")
            print(f"  J' = {Jp}")
            print(f"  E'-J' = {diff}")

print(f"\nRandom trials: {N_TRIALS}")
print(f"Failures: {failures}")
if failures > 0:
    print(f"Failure rate: {failures/N_TRIALS:.4%}")
    print(f"Worst violation: {worst_violation}")
else:
    print("J' <= E' holds in ALL random trials.")


# ===== Part 2: Tree DP factor product closure =====
print("\n" + "=" * 60)
print("Part 2: Tree DP factor product closure (n <= 18)")
print("=" * 60)

tree_failures = 0
total_pair_tests = 0

for n in range(4, 19):
    try:
        result = subprocess.run(
            ["/opt/homebrew/bin/geng", str(n), f"{n-1}:{n-1}", "-c", "-q"],
            capture_output=True, text=True, timeout=120
        )
        trees = result.stdout.strip().split('\n')
    except Exception as e:
        print(f"  n={n}: geng failed ({e}), skipping")
        continue

    n_failures = 0
    n_tests = 0

    for tree_str in trees:
        if not tree_str.strip():
            continue
        nn, adj = parse_graph6(tree_str)
        if nn is None:
            continue

        for root in range(nn):
            dp0, dp1, children_r = tree_dp(nn, adj, root)

            child_pairs = []
            for c in children_r[root]:
                E_c = dp0[c]
                J_c = dp1[c][1:]  # divide by x
                child_pairs.append((E_c, J_c))

            if len(child_pairs) < 2:
                continue

            E_acc = child_pairs[0][0]
            J_acc = child_pairs[0][1]

            for idx in range(1, len(child_pairs)):
                E2, J2 = child_pairs[idx]
                passes, Ep, Jp, diff = test_closure(E_acc, J_acc, E2, J2)
                n_tests += 1
                total_pair_tests += 1

                if not passes:
                    n_failures += 1
                    tree_failures += 1
                    if tree_failures <= 5:
                        print(f"\n  Tree CE (n={nn}, root={root}):")
                        print(f"    E_acc={E_acc}, J_acc={J_acc}")
                        print(f"    E2={E2}, J2={J2}")
                        print(f"    E'-J' = {diff}")

                E_acc = Ep
                J_acc = Jp

    status = "PASS" if n_failures == 0 else f"FAIL ({n_failures})"
    print(f"  n={n}: {len(trees)} trees, {n_tests} pair tests, {status}")

print(f"\nTree DP total: {total_pair_tests} pair tests, {tree_failures} failures")


# ===== Part 3: Verify J <= E at single-subtree level =====
print("\n" + "=" * 60)
print("Part 3: Verify J_c <= E_c for individual subtrees (n <= 18)")
print("=" * 60)

subtree_violations = 0
subtree_checks = 0

for n in range(2, 19):
    try:
        result = subprocess.run(
            ["/opt/homebrew/bin/geng", str(n), f"{n-1}:{n-1}", "-c", "-q"],
            capture_output=True, text=True, timeout=120
        )
        trees = result.stdout.strip().split('\n')
    except Exception:
        continue

    n_violations = 0

    for tree_str in trees:
        if not tree_str.strip():
            continue
        nn, adj = parse_graph6(tree_str)
        if nn is None:
            continue

        for root in range(nn):
            dp0, dp1, children_r = tree_dp(nn, adj, root)

            for v in range(nn):
                E_v = dp0[v]
                J_v = dp1[v][1:]
                subtree_checks += 1
                if not check_leq(J_v, E_v):
                    subtree_violations += 1
                    n_violations += 1
                    if subtree_violations <= 3:
                        print(f"  Violation: n={nn}, root={root}, v={v}")
                        print(f"    E={E_v}, J={J_v}")
                        print(f"    E-J={poly_sub(E_v, J_v)}")

    status = "PASS" if n_violations == 0 else f"FAIL ({n_violations})"
    print(f"  n={n}: {status}")

print(f"\nSubtree J<=E total: {subtree_checks} checks, {subtree_violations} violations")
