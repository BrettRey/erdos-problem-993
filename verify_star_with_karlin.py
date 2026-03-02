"""
Follow-up to verify_star_condition.py: at SVs where ⋆ fails,
check whether the full curvature-augmented identity still gives w_k ≥ 0.

Three terms:
  C_k · w_k = Term1 + Term2 + Term3

  Term1 = C_k · d_k(A, C)     [Karlin, ≥0]
  Term2 = C_{k+1} · d_{k-1}(B, C)
  Term3 = B_k · c_k(C)        [curvature, ≥0]

where A = E_acc·g, B = E_acc·h, C = J_acc·g.

⋆ says Term2 + Term3 ≥ 0. When ⋆ fails, we check:
  (1) Does Term1 rescue? (i.e., Term1 + Term2 + Term3 ≥ 0)
  (2) What's the ratio Term1 / |Term2 + Term3| when Term2+Term3 < 0?
"""

import subprocess
from itertools import permutations
from collections import defaultdict

try:
    from numpy import convolve as _np_convolve
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

_INT64_SAFE = 2**62


def _polymul_python(a, b):
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _polymul(a, b):
    if not a or not b:
        return []
    if _HAS_NUMPY:
        max_a = max(abs(x) for x in a)
        max_b = max(abs(x) for x in b)
        max_terms = min(len(a), len(b))
        if max_a > 0 and max_b > 0 and max_terms * max_a * max_b < _INT64_SAFE:
            return _np_convolve(a, b).astype(int).tolist()
    return _polymul_python(a, b)


def _polyadd(a, b):
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


def _coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0


def dp_all(n, adj, root):
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
            pS = [1]
            pE = [1]
            for c in children[v]:
                Ic = _polyadd(dp0[c], [0] + dp1s[c])
                pS = _polymul(pS, Ic)
                pE = _polymul(pE, dp0[c])
            dp0[v] = pS
            dp1s[v] = pE

    return children, dp0, dp1s


def find_mode(I_poly):
    m = 0
    for k in range(1, len(I_poly)):
        if I_poly[k] > I_poly[m]:
            m = k
    return m


def analyze_step(E_acc, J_acc, g, h, mode_I):
    """Full 3-term analysis at one incremental step.

    Returns dict with per-k data for k in prefix.
    """
    A = _polymul(E_acc, g)
    B = _polymul(E_acc, h)
    C = _polymul(J_acc, g)

    results = []
    star_fails = 0
    w_fails = 0
    min_star_margin = None
    min_w = None
    min_ratio = None  # Term1 / |Term2+Term3| when Term2+Term3 < 0

    for k in range(1, mode_I):
        Ak = _coeff(A, k+1)
        Akm = _coeff(A, k)
        Ck = _coeff(C, k)
        Ckm = _coeff(C, k-1)
        Ckp = _coeff(C, k+1)
        Bk = _coeff(B, k)
        Bkm = _coeff(B, k-1)

        # d_k(A,C) = A_{k+1}*C_k - A_k*C_{k+1}
        dk_AC = Ak * Ck - Akm * Ckp
        # d_{k-1}(B,C) = B_k*C_{k-1} - B_{k-1}*C_k
        dkm1_BC = Bk * Ckm - Bkm * Ck
        # c_k(C) = C_k^2 - C_{k-1}*C_{k+1}
        ck_C = Ck * Ck - Ckm * Ckp

        Term1 = Ck * dk_AC         # Karlin (≥0)
        Term2 = Ckp * dkm1_BC      # can be negative
        Term3 = Bk * ck_C          # curvature (≥0)

        star_margin = Bk * Ck - Bkm * Ckp  # ⋆ condition
        total = Term1 + Term2 + Term3       # = C_k * w_k

        # w_k = total / C_k (C_k > 0 for k < mode in prefix)
        if Ck > 0:
            w_k = total // Ck  # integer division (should be exact)
            # verify
            if total != w_k * Ck:
                w_k = total / Ck  # fallback to float
        else:
            w_k = 0

        if star_margin < 0:
            star_fails += 1
            if min_star_margin is None or star_margin < min_star_margin:
                min_star_margin = star_margin

        t23 = Term2 + Term3
        if t23 < 0 and Term1 > 0:
            ratio = Term1 / abs(t23)
            if min_ratio is None or ratio < min_ratio:
                min_ratio = ratio

        if isinstance(w_k, (int, float)) and w_k < 0:
            w_fails += 1
        if min_w is None or (isinstance(w_k, (int, float)) and w_k < min_w):
            min_w = w_k

    return {
        'star_fails': star_fails,
        'w_fails': w_fails,
        'min_star_margin': min_star_margin,
        'min_w': min_w,
        'min_ratio': min_ratio,
    }


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode('ascii')]
    if data[0] <= 62:
        n = data[0]
        off = 1
    else:
        n = (data[1] << 12) | (data[2] << 6) | data[3]
        off = 4
    adj = [[] for _ in range(n)]
    bits = []
    for d in data[off:]:
        for shift in range(5, -1, -1):
            bits.append((d >> shift) & 1)
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                adj[i].append(j)
                adj[j].append(i)
            idx += 1
    return n, adj


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=18)
    parser.add_argument('--min-n', type=int, default=5)
    parser.add_argument('--geng', default='/opt/homebrew/bin/geng')
    args = parser.parse_args()

    total_trees = 0
    total_sv = 0
    total_steps = 0

    # Step-level counts
    star_fail_steps = 0     # steps where ⋆ fails for at least one k
    w_fail_steps = 0        # steps where w_k < 0 for at least one k
    karlin_rescues = 0      # steps where ⋆ fails but w_k ≥ 0
    karlin_not_enough = 0   # steps where both ⋆ and w_k fail

    # SV-level: best ordering
    sv_no_good = 0          # SVs where no ordering has w_k ≥ 0 at all steps
    sv_star_fails_w_ok = 0  # SVs where best ordering has ⋆ fails but w≥0

    global_min_ratio = None
    global_min_w = None
    global_min_w_info = None

    for n in range(args.min_n, args.max_n + 1):
        cmd = [args.geng, '-q', str(n), f'{n-1}:{n-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        tree_count = 0
        n_sv = 0
        n_no_good = 0

        for line in proc.stdout:
            g6 = line.decode('ascii').strip()
            if not g6:
                continue
            nn, adj = parse_graph6(g6)
            tree_count += 1
            deg = [len(adj[v]) for v in range(n)]

            for root in range(n):
                leaf_nbrs = [v for v in adj[root] if deg[v] == 1]
                if not leaf_nbrs:
                    continue

                children_all, dp0, dp1s = dp_all(n, adj, root)
                children_r = children_all[root]

                nonleaf_children = []
                ell = 0
                for c in children_r:
                    if len(adj[c]) == 1:
                        ell += 1
                    else:
                        nonleaf_children.append(c)

                s = len(nonleaf_children)
                if s < 2:
                    continue

                n_sv += 1

                factors = [(dp0[c], dp1s[c]) for c in nonleaf_children]

                E_root = dp0[root]
                J_root = dp1s[root]
                I_poly = _polyadd(E_root, [0] + J_root)
                mode_I = find_mode(I_poly)
                if mode_I <= 1:
                    continue

                # Try all orderings (for s≤4) or sample
                if s <= 4:
                    orderings = list(permutations(range(s)))
                else:
                    orderings = [tuple(range(s)), tuple(range(s-1, -1, -1))]
                    import random
                    rng = random.Random(42)
                    for _ in range(min(20, s*s)):
                        perm = list(range(s))
                        rng.shuffle(perm)
                        orderings.append(tuple(perm))
                    orderings = list(set(orderings))

                best_ord_min_w = None
                best_ord_had_star_fail = False
                any_ord_all_w_ok = False

                for perm in orderings:
                    E_acc = [1]
                    for _ in range(ell):
                        E_acc = _polymul(E_acc, [1, 1])
                    J_acc = [1]

                    ord_all_w_ok = True
                    ord_any_star_fail = False
                    ord_min_w = None
                    ord_min_ratio = None

                    for step_idx, ci in enumerate(perm):
                        g_c, h_c = factors[ci]

                        if step_idx == 0:
                            I_c = _polyadd(g_c, [0] + h_c)
                            E_acc = _polymul(E_acc, I_c)
                            J_acc = _polymul(J_acc, g_c)
                            continue

                        total_steps += 1
                        res = analyze_step(E_acc, J_acc, g_c, h_c, mode_I)

                        if res['star_fails'] > 0:
                            star_fail_steps += 1
                            ord_any_star_fail = True
                            if res['w_fails'] > 0:
                                w_fail_steps += 1
                                karlin_not_enough += 1
                                ord_all_w_ok = False
                            else:
                                karlin_rescues += 1

                        if res['w_fails'] > 0:
                            ord_all_w_ok = False

                        if res['min_w'] is not None:
                            if ord_min_w is None or res['min_w'] < ord_min_w:
                                ord_min_w = res['min_w']
                        if res['min_ratio'] is not None:
                            if ord_min_ratio is None or res['min_ratio'] < ord_min_ratio:
                                ord_min_ratio = res['min_ratio']

                        I_c = _polyadd(g_c, [0] + h_c)
                        E_acc = _polymul(E_acc, I_c)
                        J_acc = _polymul(J_acc, g_c)

                    if ord_all_w_ok:
                        any_ord_all_w_ok = True
                        if ord_any_star_fail:
                            best_ord_had_star_fail = True

                    if ord_min_w is not None:
                        if best_ord_min_w is None or ord_min_w > best_ord_min_w:
                            best_ord_min_w = ord_min_w

                    if ord_min_ratio is not None:
                        if global_min_ratio is None or ord_min_ratio < global_min_ratio:
                            global_min_ratio = ord_min_ratio

                if not any_ord_all_w_ok:
                    sv_no_good += 1
                    n_no_good += 1
                elif best_ord_had_star_fail:
                    sv_star_fails_w_ok += 1

                if best_ord_min_w is not None:
                    if global_min_w is None or best_ord_min_w < global_min_w:
                        global_min_w = best_ord_min_w
                        global_min_w_info = (n, root, s, ell, g6)

        proc.wait()
        total_trees += tree_count
        total_sv += n_sv
        print(f"n={n:2d}: {tree_count:>8d} trees | s≥2 SVs={n_sv:>7d} | "
              f"no_good_ord(w)={n_no_good}")

    print()
    print("=" * 70)
    print(f"SUMMARY n={args.min_n}..{args.max_n}")
    print("=" * 70)
    print(f"Trees:                     {total_trees:>12d}")
    print(f"Support vertices (s≥2):    {total_sv:>12d}")
    print(f"Incremental steps checked: {total_steps:>12d}")
    print()
    print(f"Step-level ⋆ failures:     {star_fail_steps:>12d}")
    print(f"Step-level w_k < 0:        {w_fail_steps:>12d}")
    print(f"Karlin rescues (⋆ fail, w≥0): {karlin_rescues:>9d}")
    print(f"Karlin not enough:         {karlin_not_enough:>12d}")
    print()
    print(f"SVs with no good ordering (w_k<0 everywhere): {sv_no_good}")
    print(f"SVs where ⋆ fails but Karlin rescues: {sv_star_fails_w_ok}")
    print()
    print(f"Global min w_k (best ordering): {global_min_w}")
    if global_min_w_info:
        print(f"  at n={global_min_w_info[0]}, root={global_min_w_info[1]}, "
              f"s={global_min_w_info[2]}, ell={global_min_w_info[3]}")
    print(f"Global min Karlin/|T2+T3| ratio: {global_min_ratio}")


if __name__ == '__main__':
    main()
