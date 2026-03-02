"""Check whether E ratio-dominates J at support vertices for ALL k.

P2 is verified in the prefix (k <= mode-1). Does E_{k+1}·J_k >= E_k·J_{k+1}
hold for ALL k (including tail)?

If yes: by Karlin's TP2 theorem, (1+x) preserves ratio dominance, so
E ratio-dominates (1+x)J, giving LR_k >= 0 and SCC follows trivially.

If no: need to understand where it fails and by how much.
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
    if k < 0 or poly is None or k >= len(poly):
        return 0
    return poly[k]


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=20)
    parser.add_argument('--min-n', type=int, default=3)
    args = parser.parse_args()

    t0 = time.time()
    total_checks = 0
    total_trees = 0
    total_support = 0
    p2_fails_all_k = 0  # E_{k+1}·J_k < E_k·J_{k+1} at ANY k
    p2_fails_tail = 0   # fails at k >= mode
    p2_fails_prefix = 0  # fails at k < mode

    # Also check: (1+x)J ratio dominance
    # Does E_{k+1}·[(1+x)J]_k >= E_k·[(1+x)J]_{k+1} for all k?
    x1j_fails = 0

    min_margin_all = None
    min_margin_info = None
    min_margin_tail = None
    min_margin_tail_info = None

    # Track fails by k-mode
    fail_by_krel = {}

    for nn in range(args.min_n, args.max_n + 1):
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0
        nn_fails = 0

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
                total_support += 1

                dp0, dp1s, children = dp_rooted(n, adj, root)

                E = dp0[root]
                J = dp1s[root]
                I_poly = _polyadd(E, [0] + J)

                # Mode of I
                mode = 0
                for k in range(1, len(I_poly)):
                    if I_poly[k] > I_poly[mode]:
                        mode = k

                max_k = max(len(E), len(J)) + 1

                for k in range(max_k):
                    total_checks += 1
                    ek = coeff(E, k)
                    ekp1 = coeff(E, k+1)
                    jk = coeff(J, k)
                    jkp1 = coeff(J, k+1)

                    # P2: E_{k+1}·J_k >= E_k·J_{k+1}
                    margin = ekp1 * jk - ek * jkp1
                    if margin < 0:
                        p2_fails_all_k += 1
                        nn_fails += 1
                        krel = k - mode
                        fail_by_krel[krel] = fail_by_krel.get(krel, 0) + 1
                        if k >= mode:
                            p2_fails_tail += 1
                        else:
                            p2_fails_prefix += 1

                        if min_margin_all is None or margin < min_margin_all:
                            min_margin_all = margin
                            min_margin_info = (n, root, k, mode,
                                               ek, ekp1, jk, jkp1, line)
                        if k >= mode and (min_margin_tail is None or margin < min_margin_tail):
                            min_margin_tail = margin
                            min_margin_tail_info = (n, root, k, mode,
                                                    ek, ekp1, jk, jkp1, line)

                    # Also check (1+x)J vs E
                    x1j_k = jk + coeff(J, k-1)
                    x1j_kp1 = jkp1 + jk
                    m2 = ekp1 * x1j_k - ek * x1j_kp1
                    if m2 < 0:
                        x1j_fails += 1

        proc.wait()
        total_trees += n_trees
        elapsed = time.time() - t0
        print(f"n={nn}: {n_trees} trees, P2_fail_all_k={nn_fails}, "
              f"x1j_fails={x1j_fails}, {elapsed:.1f}s", flush=True)

    print(f"\n{'='*70}")
    print(f"RATIO DOMINANCE CHECK: E ≽ J and E ≽ (1+x)J")
    print(f"{'='*70}")
    print(f"Range: n = {args.min_n}..{args.max_n}")
    print(f"Trees:           {total_trees}")
    print(f"Support vertices:{total_support}")
    print(f"Total checks:    {total_checks}")
    print(f"\nP2 (E ≽ J) fails at ANY k:     {p2_fails_all_k}")
    print(f"  in prefix (k < mode):         {p2_fails_prefix}")
    print(f"  in tail (k >= mode):          {p2_fails_tail}")
    print(f"\nE ≽ (1+x)J fails:              {x1j_fails}")

    if min_margin_info:
        print(f"\nTightest P2 fail (any k):")
        n, root, k, mode, ek, ekp1, jk, jkp1, g6 = min_margin_info
        print(f"  n={n}, root={root}, k={k} (mode={mode}, k-mode={k-mode})")
        print(f"  E_k={ek}, E_{{k+1}}={ekp1}, J_k={jk}, J_{{k+1}}={jkp1}")
        print(f"  margin={ekp1*jk-ek*jkp1}")
        print(f"  g6={g6}")

    if min_margin_tail_info and p2_fails_tail > 0:
        print(f"\nTightest P2 fail (tail only):")
        n, root, k, mode, ek, ekp1, jk, jkp1, g6 = min_margin_tail_info
        print(f"  n={n}, root={root}, k={k} (mode={mode}, k-mode={k-mode})")
        print(f"  E_k={ek}, E_{{k+1}}={ekp1}, J_k={jk}, J_{{k+1}}={jkp1}")
        print(f"  margin={ekp1*jk-ek*jkp1}")

    if fail_by_krel:
        print(f"\nP2 fails by (k - mode):")
        for kr in sorted(fail_by_krel.keys()):
            print(f"  k-mode={kr:+3d}: {fail_by_krel[kr]:>10,}")

    print(f"\nElapsed: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
