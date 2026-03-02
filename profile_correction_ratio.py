"""Profile the correction-to-main ratio at each incremental E≽J stage.

At stage t:
  E^{(t)} = A + xB where A = E^{(t-1)}·E_t, B = E^{(t-1)}·J_t
  J^{(t)} = J^{(t-1)}·E_t

  Δ_k(E^{(t)}, J^{(t)}) = Δ_k(A, J^{(t)}) + [B_k·J^{(t)}_k - B_{k-1}·J^{(t)}_{k+1}]
                         = main_k + corr_k

  main_k ≥ 0 by Karlin (PROVED). corr_k can be negative.
  When corr_k < 0: what is main_k/|corr_k|?

Also track: min margin per stage, and which factor structure gives tightest ratios.
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
    total_stages = 0
    total_k_checks = 0
    corr_neg_stages = 0  # stages where corr < 0 at any k
    corr_neg_checks = 0  # individual (stage, k) where corr < 0

    min_ratio = None  # min main_k/|corr_k| when corr_k < 0
    min_ratio_info = None

    min_total_margin = None  # min main_k + corr_k (= Δ_k)
    min_total_info = None

    # Histogram of ratio bins
    ratio_bins = {'>10': 0, '5-10': 0, '3-5': 0, '2-3': 0, '1.5-2': 0,
                  '1-1.5': 0, '0.5-1': 0, '<0.5': 0}

    # Track by n
    per_n_min_ratio = {}

    for nn in range(args.min_n, args.max_n + 1):
        cmd = [GENG, '-q', str(nn), f'{nn-1}:{nn-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        n_trees = 0
        nn_min_ratio = None

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
                ell = leaf_count[root]

                if len(non_leaf_children) == 0:
                    continue

                # Stage 0: E^(0) = (1+x)^ℓ, J^(0) = [1]
                E_acc = [1]
                for _ in range(ell):
                    E_acc = _polymul(E_acc, [1, 1])
                J_acc = [1]

                for t_idx, c in enumerate(non_leaf_children):
                    total_stages += 1
                    Ic = _polyadd(dp0[c], [0] + dp1s[c])
                    Ec = dp0[c]
                    Jc = dp1s[c]

                    # Karlin main part: A = E_acc · E_c
                    A = _polymul(E_acc, Ec)
                    # Correction: B = E_acc · J_c
                    B = _polymul(E_acc, Jc)
                    # New products
                    E_new = _polymul(E_acc, Ic)
                    J_new = _polymul(J_acc, Ec)

                    max_k = max(len(E_new), len(J_new))
                    stage_has_neg = False

                    for k in range(max_k):
                        total_k_checks += 1

                        # main_k = Δ_k(A, J_new) = A_{k+1}·J_new_k - A_k·J_new_{k+1}
                        main_k = coeff(A, k+1) * coeff(J_new, k) - coeff(A, k) * coeff(J_new, k+1)

                        # corr_k = B_k·J_new_k - B_{k-1}·J_new_{k+1}
                        corr_k = coeff(B, k) * coeff(J_new, k) - coeff(B, k-1) * coeff(J_new, k+1)

                        total_k_margin = main_k + corr_k

                        if min_total_margin is None or total_k_margin < min_total_margin:
                            min_total_margin = total_k_margin
                            min_total_info = (n, root, t_idx, k, main_k, corr_k, line)

                        if corr_k < 0:
                            if not stage_has_neg:
                                corr_neg_stages += 1
                                stage_has_neg = True
                            corr_neg_checks += 1
                            abs_corr = -corr_k
                            if abs_corr > 0:
                                ratio = main_k / abs_corr
                            else:
                                ratio = float('inf')

                            # Bin the ratio
                            if ratio > 10:
                                ratio_bins['>10'] += 1
                            elif ratio > 5:
                                ratio_bins['5-10'] += 1
                            elif ratio > 3:
                                ratio_bins['3-5'] += 1
                            elif ratio > 2:
                                ratio_bins['2-3'] += 1
                            elif ratio > 1.5:
                                ratio_bins['1.5-2'] += 1
                            elif ratio > 1:
                                ratio_bins['1-1.5'] += 1
                            elif ratio > 0.5:
                                ratio_bins['0.5-1'] += 1
                            else:
                                ratio_bins['<0.5'] += 1

                            if min_ratio is None or ratio < min_ratio:
                                min_ratio = ratio
                                min_ratio_info = (n, root, t_idx, k, main_k, corr_k,
                                                  total_k_margin, line)

                            if nn_min_ratio is None or ratio < nn_min_ratio:
                                nn_min_ratio = ratio

                    E_acc = E_new
                    J_acc = J_new

        proc.wait()
        per_n_min_ratio[nn] = nn_min_ratio
        elapsed = time.time() - t0
        print(f"n={nn}: {n_trees} trees, stages={total_stages}, "
              f"corr_neg_stages={corr_neg_stages}, corr_neg_checks={corr_neg_checks}, "
              f"min_ratio={nn_min_ratio:.4f}" if nn_min_ratio is not None else
              f"n={nn}: {n_trees} trees, stages={total_stages}, "
              f"corr_neg_stages={corr_neg_stages}, corr_neg_checks={corr_neg_checks}, "
              f"min_ratio=N/A",
              flush=True)

    print(f"\n{'='*70}")
    print(f"CORRECTION-TO-MAIN RATIO PROFILE")
    print(f"{'='*70}")
    print(f"Range: n = {args.min_n}..{args.max_n}")
    print(f"Total stages:       {total_stages}")
    print(f"Total k-checks:     {total_k_checks}")
    print(f"Corr<0 stages:      {corr_neg_stages} ({100*corr_neg_stages/max(1,total_stages):.1f}%)")
    print(f"Corr<0 checks:      {corr_neg_checks} ({100*corr_neg_checks/max(1,total_k_checks):.1f}%)")

    print(f"\nMin main/|corr| ratio when corr<0: {min_ratio}")
    if min_ratio_info:
        n, root, tidx, k, mk, ck, tm, g6 = min_ratio_info
        print(f"  at n={n}, root={root}, factor={tidx}, k={k}")
        print(f"  main={mk}, corr={ck}, total={tm}")
        print(f"  g6={g6}")

    print(f"\nMin total margin (main+corr): {min_total_margin}")
    if min_total_info:
        n, root, tidx, k, mk, ck, g6 = min_total_info
        print(f"  at n={n}, root={root}, factor={tidx}, k={k}")
        print(f"  main={mk}, corr={ck}")
        print(f"  g6={g6}")

    print(f"\nRatio histogram (main/|corr| when corr<0):")
    for label in ['>10', '5-10', '3-5', '2-3', '1.5-2', '1-1.5', '0.5-1', '<0.5']:
        ct = ratio_bins[label]
        if ct > 0:
            print(f"  {label:>8}: {ct:>12,} ({100*ct/max(1,corr_neg_checks):.1f}%)")

    print(f"\nMin ratio by n:")
    for nn in sorted(per_n_min_ratio.keys()):
        r = per_n_min_ratio[nn]
        if r is not None:
            print(f"  n={nn:3d}: {r:.4f}")
        else:
            print(f"  n={nn:3d}: N/A")

    print(f"\nElapsed: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
