"""Check the 2-term SCC decomposition at the FINAL (support vertex) level.

Identity: Δ_k = c_k + LR_k

where:
  Δ_k = ((1+x)I)_{k+1}·E_k - ((1+x)I)_k·E_{k+1}  (SCC)
  c_k = E_k² - E_{k-1}·E_{k+1}                      (LC gap of E)
  LR_k = E_k·(J_k+J_{k-1}) - E_{k+1}·(J_{k-1}+J_{k-2})  (LR minor of E vs (1+x)J)

Proof: (1+x)I = (1+x)E + x(1+x)J. By linearity of minors:
  Δ_k = Δ_k((1+x)E, E) + Δ_k(x(1+x)J, E) = c_k + LR_k.

Key questions:
1. How often is LR_k < 0?
2. When LR_k < 0, what fraction of c_k does |LR_k| consume?
3. What's the tightest c_k / |LR_k| ratio?
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
    id_fails = 0
    lr_neg = 0
    lr_total_neg = 0  # sum of |LR_k| when negative
    ck_total_when_lr_neg = 0  # sum of c_k when LR_k < 0
    scc_neg = 0
    ck_neg = 0
    total_checks = 0
    total_trees = 0

    min_utilization = None  # min (c_k + LR_k) / c_k when LR_k < 0 and c_k > 0
    min_util_info = None
    min_abs_margin = None  # min c_k + LR_k when LR_k < 0
    min_abs_info = None

    # Track by k relative to mode
    k_rel_counts = {}  # k - mode -> count of LR < 0

    for nn in range(args.min_n, args.max_n + 1):
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0
        nn_lr_neg = 0
        nn_checks = 0

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

                E = dp0[root]
                J = dp1s[root]
                I_poly = _polyadd(E, [0] + J)
                e_poly = _polyadd(I_poly, [0] + I_poly)  # (1+x)I

                # Find mode of I
                mode = 0
                for k in range(1, len(I_poly)):
                    if I_poly[k] > I_poly[mode]:
                        mode = k

                max_k = max(len(e_poly), len(E))

                for k in range(max_k):
                    total_checks += 1
                    nn_checks += 1

                    ek = coeff(E, k)
                    ekm1 = coeff(E, k-1)
                    ekp1 = coeff(E, k+1)
                    jk = coeff(J, k)
                    jkm1 = coeff(J, k-1)
                    jkm2 = coeff(J, k-2)

                    c_k = ek * ek - ekm1 * ekp1
                    lr_k = ek * (jk + jkm1) - ekp1 * (jkm1 + jkm2)
                    delta_k = coeff(e_poly, k+1) * ek - coeff(e_poly, k) * ekp1

                    # Verify identity
                    if delta_k != c_k + lr_k:
                        id_fails += 1

                    if c_k < 0:
                        ck_neg += 1

                    if delta_k < 0:
                        scc_neg += 1

                    if lr_k < 0:
                        lr_neg += 1
                        nn_lr_neg += 1
                        lr_total_neg += abs(lr_k)
                        ck_total_when_lr_neg += c_k

                        k_rel = k - mode
                        k_rel_counts[k_rel] = k_rel_counts.get(k_rel, 0) + 1

                        margin = c_k + lr_k  # = delta_k
                        if min_abs_margin is None or margin < min_abs_margin:
                            min_abs_margin = margin
                            min_abs_info = (n, root, k, mode, c_k, lr_k, delta_k,
                                            ek, ekm1, ekp1, jk, jkm1, jkm2, line)

                        if c_k > 0:
                            util = (c_k + lr_k) / c_k  # fraction of c_k remaining
                            if min_utilization is None or util < min_utilization:
                                min_utilization = util
                                min_util_info = (n, root, k, mode, c_k, lr_k, delta_k,
                                                 ek, ekm1, ekp1, jk, jkm1, jkm2, line)

        proc.wait()
        total_trees += n_trees
        elapsed = time.time() - t0
        pct = 100 * nn_lr_neg / max(1, nn_checks)
        print(f"n={nn}: {n_trees} trees, {nn_checks} checks, LR<0={nn_lr_neg} ({pct:.1f}%), "
              f"id_fails={id_fails}, ck<0={ck_neg}, scc<0={scc_neg}, {elapsed:.1f}s", flush=True)

    print(f"\n{'='*70}")
    print(f"2-TERM SCC DECOMPOSITION: Δ_k = c_k + LR_k")
    print(f"{'='*70}")
    print(f"Range: n = {args.min_n}..{args.max_n}")
    print(f"Trees:           {total_trees}")
    print(f"Total checks:    {total_checks}")
    print(f"Identity fails:  {id_fails}")
    print(f"c_k < 0:         {ck_neg}")
    print(f"SCC < 0:         {scc_neg}")
    print(f"LR_k < 0:        {lr_neg} ({100*lr_neg/max(1,total_checks):.2f}%)")

    if lr_neg > 0:
        print(f"\nWhen LR_k < 0:")
        print(f"  Avg |LR_k|:    {lr_total_neg/lr_neg:.1f}")
        print(f"  Avg c_k:       {ck_total_when_lr_neg/lr_neg:.1f}")
        print(f"  Avg c_k/|LR|:  {ck_total_when_lr_neg/lr_total_neg:.4f}")

    print(f"\nMin utilization ratio (c_k+LR_k)/c_k when LR<0: {min_utilization}")
    if min_util_info:
        n, root, k, mode, ck, lrk, dk, ek, ekm1, ekp1, jk, jkm1, jkm2, g6 = min_util_info
        print(f"  at n={n}, root={root}, k={k} (mode={mode}, k-mode={k-mode})")
        print(f"  c_k={ck}, LR_k={lrk}, Δ_k={dk}")
        print(f"  E: [{ekm1}, {ek}, {ekp1}]")
        print(f"  J: [{jkm2}, {jkm1}, {jk}]")
        print(f"  g6={g6.strip()}")

    print(f"\nMin absolute margin (c_k+LR_k) when LR<0: {min_abs_margin}")
    if min_abs_info:
        n, root, k, mode, ck, lrk, dk, ek, ekm1, ekp1, jk, jkm1, jkm2, g6 = min_abs_info
        print(f"  at n={n}, root={root}, k={k} (mode={mode}, k-mode={k-mode})")
        print(f"  c_k={ck}, LR_k={lrk}, Δ_k={dk}")
        print(f"  E: [{ekm1}, {ek}, {ekp1}]")
        print(f"  J: [{jkm2}, {jkm1}, {jk}]")

    # LR_k < 0 by k relative to mode
    if k_rel_counts:
        print(f"\nLR_k < 0 distribution by (k - mode):")
        for kr in sorted(k_rel_counts.keys()):
            print(f"  k-mode={kr:+3d}: {k_rel_counts[kr]:>10,}")

    print(f"\nElapsed: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
