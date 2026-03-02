"""Profile the E≽J correction ratio by s-value (number of non-leaf children).

For each support vertex, computes:
- s = number of non-leaf children
- The incremental Karlin main / |correction| ratio at each stage
- Min ratio per s-value

Key question: is the pendant-star (s=1) the tightest case globally,
and do s=2 cases have better ratios?
"""

import subprocess
import sys
import time

sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul, _polyadd

GENG = '/opt/homebrew/bin/geng'


def parse_g6(g6: str):
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
    if k < 0 or k >= len(poly):
        return 0
    return poly[k]


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    args = parser.parse_args()

    # Track per s-value
    # For each s: count of support vertices, count of neg corrections,
    #             min main/|corr| ratio, tightest examples
    from collections import defaultdict
    s_stats = defaultdict(lambda: {
        'vertices': 0,
        'neg_corr_checks': 0,
        'total_checks': 0,
        'min_ratio': float('inf'),
        'ej_fails': 0,
        'tightest': [],  # (ratio, n, k, stage, g6)
    })

    t0 = time.time()

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

                s_stats[s]['vertices'] += 1

                if s == 0:
                    continue  # star, trivial

                # Build incrementally
                factors = []
                for c in non_leaf_children:
                    Ic = _polyadd(dp0[c], [0] + dp1s[c])
                    Ec = dp0[c]
                    Jc = dp1s[c]
                    factors.append((Ic, Ec, Jc))

                # Include leaves: (1+x)^ell
                ell = leaf_count[root]
                E_acc = [1]
                for _ in range(ell):
                    E_acc = _polymul(E_acc, [1, 1])
                J_acc = [1]

                for stage, (Ic, Ec, Jc) in enumerate(factors, 1):
                    # Save E_old, J_old before update
                    E_old = E_acc[:]
                    J_old = J_acc[:]

                    E_acc = _polymul(E_old, Ic)
                    J_acc = _polymul(J_old, Ec)

                    # Karlin main part: A = E_old * Ec
                    A = _polymul(E_old, Ec)
                    # Correction: xB = x * E_old * Jc
                    B = _polymul(E_old, Jc)  # B without the x shift

                    max_k = max(len(E_acc), len(J_acc)) + 1

                    for k in range(max_k):
                        s_stats[s]['total_checks'] += 1

                        ek = coeff(E_acc, k)
                        ekp1 = coeff(E_acc, k + 1)
                        jk = coeff(J_acc, k)
                        jkp1 = coeff(J_acc, k + 1)

                        # E≽J minor
                        delta = ekp1 * jk - ek * jkp1
                        if delta < 0:
                            s_stats[s]['ej_fails'] += 1

                        # Karlin main part minor: Δ_k(A, J_acc)
                        ak = coeff(A, k)
                        akp1 = coeff(A, k + 1)
                        main = akp1 * jk - ak * jkp1

                        # Correction minor: Δ_k(xB, J_acc)
                        # xB has coeff (xB)_k = B_{k-1}
                        xbk = coeff(B, k - 1)
                        xbkp1 = coeff(B, k)
                        corr = xbkp1 * jk - xbk * jkp1

                        if corr < 0 and main > 0:
                            s_stats[s]['neg_corr_checks'] += 1
                            ratio = main / abs(corr)
                            if ratio < s_stats[s]['min_ratio']:
                                s_stats[s]['min_ratio'] = ratio
                            s_stats[s]['tightest'].append(
                                (ratio, n, k, stage, line.strip())
                            )

        proc.wait()
        elapsed = time.time() - t0
        print(f"n={nn}: {n_trees} trees, {elapsed:.1f}s", flush=True)

    # Sort and trim tightest lists
    for s in s_stats:
        s_stats[s]['tightest'].sort(key=lambda x: x[0])
        s_stats[s]['tightest'] = s_stats[s]['tightest'][:10]

    print(f"\n{'='*80}")
    print(f"CORRECTION RATIO BY s-VALUE (number of non-leaf children)")
    print(f"{'='*80}")
    print(f"{'s':>3s}  {'vertices':>10s}  {'checks':>12s}  {'neg_corr':>10s}  {'EJ_fails':>10s}  {'min_ratio':>10s}")
    for s in sorted(s_stats.keys()):
        st = s_stats[s]
        mr = f"{st['min_ratio']:.4f}" if st['min_ratio'] < float('inf') else "inf"
        print(f"{s:3d}  {st['vertices']:10,}  {st['total_checks']:12,}  "
              f"{st['neg_corr_checks']:10,}  {st['ej_fails']:10,}  {mr:>10s}")

    print(f"\nTightest cases per s-value:")
    for s in sorted(s_stats.keys()):
        if not s_stats[s]['tightest']:
            continue
        print(f"\n  s={s}:")
        for rank, (ratio, n, k, stage, g6) in enumerate(s_stats[s]['tightest'][:5]):
            print(f"    #{rank+1}: ratio={ratio:.4f}, n={n}, k={k}, stage={stage}, g6={g6}")

    print(f"\nElapsed: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
