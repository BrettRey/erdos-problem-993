"""Profile the ratio e_k/b_k at each incremental stage.

SCC ≥ 0 ⟺ the ratio e_k/b_k is nondecreasing (when b_k > 0).

This script profiles:
1. The ratio e_k/b_k at each stage
2. How the correction (x(1+x)γ, xγ) changes the ratio profile
3. The "correction ratio" γ_{k-2}/γ_{k-1} (should be related to LC of γ)
4. Whether the (1+x) amplification is the key mechanism

Usage:
    python3 profile_ratio_monotonicity.py --max-n 18 --workers 1
"""

import argparse
import subprocess
import sys
import time

sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul, _polyadd

GENG = '/opt/homebrew/bin/geng'


def parse_g6(g6: str):
    s = g6.strip()
    idx = 0
    n = ord(s[idx]) - 63
    idx += 1
    adj = [[] for _ in range(n)]
    bits = []
    for ch in s[idx:]:
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=16)
    parser.add_argument('--min-n', type=int, default=10)
    args = parser.parse_args()

    # For each tight case, profile the ratio e_k/b_k before and after correction
    print("Looking for cases where the ratio change from correction is maximal...\n")

    t0 = time.time()
    worst_drop = 0.0
    worst_info = None

    for nn in range(args.min_n, args.max_n + 1):
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0

        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue

            n, adj = parse_g6(line)
            n_trees += 1

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
                    factors.append((Ic, Ec, Rc))

                E_acc = [1]
                J_acc = [1]

                for stage, (P, Q, R) in enumerate(factors, 1):
                    I_old = _polyadd(E_acc, [0] + J_acc)
                    e_old = _polyadd(I_old, [0] + I_old)

                    # After Step 1 (convolution with Q): (e_old·Q, E_old·Q)
                    eQ = _polymul(e_old, Q)
                    EQ = _polymul(E_acc, Q)

                    # After Step 2 (add correction): (e_new, b_new)
                    E_new = _polymul(E_acc, P)
                    J_new = _polymul(J_acc, Q)
                    I_new = _polyadd(E_new, [0] + J_new)
                    e_new = _polyadd(I_new, [0] + I_new)

                    # The correction
                    gamma = _polymul(E_acc, R)

                    max_k = max(len(e_new), len(E_new)) - 1

                    for k in range(1, max_k):
                        bQ_k = coeff(EQ, k)
                        bQ_kp1 = coeff(EQ, k+1)
                        bnew_k = coeff(E_new, k)
                        bnew_kp1 = coeff(E_new, k+1)

                        if bQ_k > 0 and bQ_kp1 > 0 and bnew_k > 0 and bnew_kp1 > 0:
                            # Ratio before correction
                            r_before_k = coeff(eQ, k) / bQ_k
                            r_before_kp1 = coeff(eQ, k+1) / bQ_kp1

                            # Ratio after correction
                            r_after_k = coeff(e_new, k) / bnew_k
                            r_after_kp1 = coeff(e_new, k+1) / bnew_kp1

                            # SCC = r_after_{k+1} * b_k - r_after_k * b_{k+1}
                            # which ≥ 0 iff r_after is nondecreasing

                            # Check if correction makes the ratio LESS nondecreasing
                            mono_before = r_before_kp1 - r_before_k
                            mono_after = r_after_kp1 - r_after_k

                            if mono_before > 0 and mono_after < mono_before:
                                drop = mono_before - mono_after
                                frac = mono_after / mono_before if mono_before > 0 else 0

                                if drop > worst_drop:
                                    worst_drop = drop
                                    worst_info = {
                                        'n': n, 'k': k, 'stage': stage, 'nfac': s,
                                        'g6': line,
                                        'mono_before': mono_before,
                                        'mono_after': mono_after,
                                        'frac': frac,
                                        'r_before': (r_before_k, r_before_kp1),
                                        'r_after': (r_after_k, r_after_kp1),
                                        'corr_ratio_k': (coeff(gamma, k-2) / coeff(gamma, k-1)) if coeff(gamma, k-1) > 0 else None,
                                        'corr_ratio_kp1': (coeff(gamma, k-1) / coeff(gamma, k)) if coeff(gamma, k) > 0 else None,
                                    }

                    E_acc = E_new
                    J_acc = J_new

        proc.wait()
        print(f"n={nn}: {n_trees} trees processed")

    elapsed = time.time() - t0

    print(f"\nWorst ratio monotonicity drop from correction:")
    if worst_info:
        w = worst_info
        print(f"  n={w['n']}, k={w['k']}, stage={w['stage']}, #fac={w['nfac']}")
        print(f"  Mono before correction: {w['mono_before']:.6f}")
        print(f"  Mono after correction:  {w['mono_after']:.6f}")
        print(f"  Fraction remaining:     {w['frac']:.6f}")
        print(f"  Ratio before: r_k = {w['r_before'][0]:.4f}, r_{k+1} = {w['r_before'][1]:.4f}")
        print(f"  Ratio after:  r_k = {w['r_after'][0]:.4f}, r_{k+1} = {w['r_after'][1]:.4f}")
        if w['corr_ratio_k'] is not None:
            print(f"  Correction ratio at k:   γ_{k-2}/γ_{k-1} = {w['corr_ratio_k']:.4f}")
        if w['corr_ratio_kp1'] is not None:
            print(f"  Correction ratio at k+1: γ_{k-1}/γ_k = {w['corr_ratio_kp1']:.4f}")
    else:
        print("  No drops found!")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
