"""Check E ≽ J at each INCREMENTAL stage of the product.

At support vertex r with ℓ leaf children and non-leaf children c_1,...,c_s:
  Stage 0: E^(0) = (1+x)^ℓ, J^(0) = [1]
  Stage t: E^(t) = E^(t-1) · I_t, J^(t) = J^(t-1) · E_t

Does E^(t) ≽ J^(t) hold at every stage?

If E^(0) ≽ J^(0) (trivial) and the step E^(t-1)·I_t ≽ J^(t-1)·E_t
is preserved, we get E ≽ J by induction.

The step decomposes as:
  (E^(t-1)·I_t)_{k+1} · (J^(t-1)·E_t)_k - (E^(t-1)·I_t)_k · (J^(t-1)·E_t)_{k+1}
  = (E^(t-1)·E_t)_{k+1}·(J^(t-1)·E_t)_k - (E^(t-1)·E_t)_k·(J^(t-1)·E_t)_{k+1}
    + x·[(E^(t-1)·J_t)]_{shifted vs (J^(t-1)·E_t)]

The first part is Δ_k(A, J^(t)) where A = E^(t-1)·E_t. By Karlin:
if E^(t-1) ≽ J^(t-1) and E_t is PF2, then A ≽ J^(t). So first part ≥ 0.

This script checks whether this Karlin decomposition holds.
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


def check_ratio_dom(A, B):
    """Check if A ≽ B (A ratio-dominates B): A_{k+1}·B_k ≥ A_k·B_{k+1} for all k."""
    max_k = max(len(A) if A else 0, len(B) if B else 0)
    fails = 0
    min_margin = None
    min_k = None
    for k in range(max_k):
        ak = coeff(A, k)
        akp1 = coeff(A, k+1)
        bk = coeff(B, k)
        bkp1 = coeff(B, k+1)
        margin = akp1 * bk - ak * bkp1
        if margin < 0:
            fails += 1
        if min_margin is None or margin < min_margin:
            min_margin = margin
            min_k = k
    return fails, min_margin, min_k


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=3)
    args = parser.parse_args()

    t0 = time.time()
    total_trees = 0
    total_stages = 0
    stage_rd_fails = 0  # E^(t) ≽ J^(t) fails at any stage
    karlin_part_fails = 0  # Karlin main part A ≽ J^(t) fails
    correction_neg = 0  # xB correction term negative
    factor_rd_fails = 0  # I_t ≽ E_t fails
    factor_rd_rev_fails = 0  # E_t ≽ J_t fails (at factor level)

    min_stage_margin = None
    min_stage_info = None

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
                ell = leaf_count[root]

                if len(non_leaf_children) == 0:
                    continue

                # Stage 0: E^(0) = (1+x)^ℓ, J^(0) = [1]
                E_acc = [1]
                for _ in range(ell):
                    E_acc = _polymul(E_acc, [1, 1])
                J_acc = [1]

                for c in non_leaf_children:
                    total_stages += 1
                    Ic = _polyadd(dp0[c], [0] + dp1s[c])
                    Ec = dp0[c]
                    Jc = dp1s[c]

                    # New products
                    E_new = _polymul(E_acc, Ic)
                    J_new = _polymul(J_acc, Ec)

                    # Check E^(t) ≽ J^(t)
                    fails, margin, mk = check_ratio_dom(E_new, J_new)
                    if fails > 0:
                        stage_rd_fails += 1

                    if min_stage_margin is None or (margin is not None and margin < min_stage_margin):
                        min_stage_margin = margin
                        min_stage_info = (n, root, len(non_leaf_children), mk, margin, line)

                    # Check Karlin main part: A = E^(t-1)·E_t ≽ J^(t) = J^(t-1)·E_t
                    A = _polymul(E_acc, Ec)
                    Af, _, _ = check_ratio_dom(A, J_new)
                    if Af > 0:
                        karlin_part_fails += 1

                    # Check factor-level: I_t ≽ E_t? and E_t ≽ J_t?
                    f1, _, _ = check_ratio_dom(Ic, Ec)
                    if f1 > 0:
                        factor_rd_fails += 1

                    f2, _, _ = check_ratio_dom(Ec, Jc)
                    if f2 > 0:
                        factor_rd_rev_fails += 1

                    # Check correction term: xB vs J^(t) where B = E^(t-1)·J_t
                    B = _polymul(E_acc, Jc)
                    max_k = max(len(E_new), len(J_new))
                    for k in range(max_k):
                        # correction = B_k · J_new_k - B_{k-1} · J_new_{k+1}
                        corr = coeff(B, k) * coeff(J_new, k) - coeff(B, k-1) * coeff(J_new, k+1)
                        if corr < 0:
                            correction_neg += 1
                            break  # just count once per stage

                    E_acc = E_new
                    J_acc = J_new

        proc.wait()
        total_trees += n_trees
        elapsed = time.time() - t0
        print(f"n={nn}: {n_trees} trees, stages={total_stages}, "
              f"stage_rd_fails={stage_rd_fails}, karlin_fails={karlin_part_fails}, "
              f"factor_It≽Et_fails={factor_rd_fails}, "
              f"factor_Et≽Jt_fails={factor_rd_rev_fails}, "
              f"{elapsed:.1f}s", flush=True)

    print(f"\n{'='*70}")
    print(f"INCREMENTAL RATIO DOMINANCE E^(t) ≽ J^(t)")
    print(f"{'='*70}")
    print(f"Range: n = {args.min_n}..{args.max_n}")
    print(f"Trees:              {total_trees}")
    print(f"Total stages:       {total_stages}")
    print(f"\nE^(t) ≽ J^(t) fails: {stage_rd_fails}")
    print(f"Karlin A ≽ J^(t) fails: {karlin_part_fails}")
    print(f"Correction xB neg:   {correction_neg}")
    print(f"\nFactor I_t ≽ E_t fails: {factor_rd_fails}")
    print(f"Factor E_t ≽ J_t fails: {factor_rd_rev_fails}")
    print(f"\nMin stage margin:    {min_stage_margin}")
    if min_stage_info:
        n, root, nf, mk, margin, g6 = min_stage_info
        print(f"  at n={n}, root={root}, #factors={nf}, k={mk}")
        print(f"  g6={g6}")
    print(f"\nElapsed: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
