"""
Test the shifted ratio condition ⋆ at s≥2 support vertices.

At each incremental step processing non-leaf child t:
  A = E_acc · g,  B = E_acc · h,  C = J_acc · g
  where g = E_t, h = J_t (child's DP factors).

⋆:  B_k · C_k  ≥  B_{k-1} · C_{k+1}   for k = 1, ..., mode-1

The curvature-augmented identity gives:
  C_k · w_k = Term1(Karlin,≥0) + C_k(B_k·C_k - B_{k-1}·C_{k+1})
So ⋆ ⟹ w_k ≥ 0 (with the Karlin term as free bonus).

We check: for each s≥2 support vertex, over ALL child orderings,
does at least one ordering satisfy ⋆ on the prefix (k < mode)?
"""

import subprocess
import sys
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


def dp_all(n, adj, root):
    """Rooted DP returning dp0[v] and dp1s[v] for all v, plus children."""
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

    dp0 = [None] * n   # E_v (exclude-root)
    dp1s = [None] * n   # J_v (include-root / x)

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
    """Return the smallest mode (argmax) of polynomial I."""
    m = 0
    for k in range(1, len(I_poly)):
        if I_poly[k] > I_poly[m]:
            m = k
    return m


def check_star_at_step(E_acc, J_acc, g, h, mode_I):
    """Check ⋆: B_k · C_k ≥ B_{k-1} · C_{k+1} for k=1..mode-1.

    B = E_acc · h, C = J_acc · g.
    Returns (passes, min_margin, worst_k).
    """
    B = _polymul(E_acc, h)
    C = _polymul(J_acc, g)

    passes = True
    min_margin = None
    worst_k = -1

    for k in range(1, mode_I):
        Bk = B[k] if k < len(B) else 0
        Bkm1 = B[k-1] if k-1 < len(B) else 0
        Ck = C[k] if k < len(C) else 0
        Ckp1 = C[k+1] if k+1 < len(C) else 0

        margin = Bk * Ck - Bkm1 * Ckp1
        if min_margin is None or margin < min_margin:
            min_margin = margin
            worst_k = k
        if margin < 0:
            passes = False

    return passes, min_margin, worst_k


def check_star_full_range(E_acc, J_acc, g, h, deg_I):
    """Check ⋆ for ALL k (not just prefix). Returns (passes, min_margin, worst_k)."""
    B = _polymul(E_acc, h)
    C = _polymul(J_acc, g)

    passes = True
    min_margin = None
    worst_k = -1

    max_k = max(len(B), len(C)) - 1
    for k in range(1, max_k):
        Bk = B[k] if k < len(B) else 0
        Bkm1 = B[k-1] if k-1 < len(B) else 0
        Ck = C[k] if k < len(C) else 0
        Ckp1 = C[k+1] if k+1 < len(C) else 0

        margin = Bk * Ck - Bkm1 * Ckp1
        if min_margin is None or margin < min_margin:
            min_margin = margin
            worst_k = k
        if margin < 0:
            passes = False

    return passes, min_margin, worst_k


def parse_graph6(s):
    """Parse a graph6 string into (n, adj)."""
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
    parser.add_argument('--max-n', type=int, default=16)
    parser.add_argument('--min-n', type=int, default=5)
    parser.add_argument('--geng', default='/opt/homebrew/bin/geng')
    parser.add_argument('--all-orderings', action='store_true',
                        help='Test all child orderings (expensive for s≥4)')
    args = parser.parse_args()

    # Global counters
    total_trees = 0
    total_sv = 0            # support vertices with s≥2
    total_steps = 0         # incremental steps at s≥2 vertices
    prefix_star_fails = 0   # steps where ⋆ fails on prefix
    full_star_fails = 0     # steps where ⋆ fails on full range
    sv_no_good_ordering = 0  # vertices where NO ordering satisfies ⋆ prefix at all steps
    sv_all_orderings_good = 0
    min_prefix_margin = None
    min_prefix_info = None
    min_full_margin = None
    min_full_info = None

    # Per-s statistics
    s_counts = defaultdict(int)
    s_prefix_pass = defaultdict(int)  # SVs where ≥1 ordering passes prefix
    s_prefix_fail = defaultdict(int)  # SVs where NO ordering passes prefix
    s_min_margin = {}  # per-s tightest prefix margin

    for n in range(args.min_n, args.max_n + 1):
        cmd = [args.geng, '-q', str(n), f'{n-1}:{n-1}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        tree_count = 0
        n_sv = 0
        n_pfails = 0

        for line in proc.stdout:
            g6 = line.decode('ascii').strip()
            if not g6:
                continue
            nn, adj = parse_graph6(g6)
            assert nn == n
            tree_count += 1
            deg = [len(adj[v]) for v in range(n)]

            # For each support vertex
            for root in range(n):
                # Must have at least one leaf neighbor
                leaf_nbrs = [v for v in adj[root] if deg[v] == 1]
                if not leaf_nbrs:
                    continue

                children_all, dp0, dp1s = dp_all(n, adj, root)
                children_r = children_all[root]

                # Separate leaf vs non-leaf children
                nonleaf_children = []
                ell = 0  # number of leaf children
                for c in children_r:
                    if not dp0[c] or (len(dp0[c]) == 1 and len(adj[c]) == 1):
                        # Leaf child: E_c = [1], J_c = [1]
                        ell += 1
                    else:
                        # Non-leaf child
                        nonleaf_children.append(c)

                s = len(nonleaf_children)
                if s < 2:
                    continue  # s=1 already proved

                n_sv += 1
                s_counts[s] += 1

                # Get factor data for each non-leaf child
                # g_c = E_c = dp0[c],  h_c = J_c = dp1s[c]
                factors = []
                for c in nonleaf_children:
                    g_c = dp0[c]
                    h_c = dp1s[c]
                    factors.append((g_c, h_c))

                # Compute I(T) for mode
                E_root = dp0[root]
                J_root = dp1s[root]
                I_poly = _polyadd(E_root, [0] + J_root)
                mode_I = find_mode(I_poly)

                if mode_I <= 1:
                    continue  # trivial prefix

                # Test orderings
                if args.all_orderings or s <= 4:
                    orderings = list(permutations(range(s)))
                else:
                    # For large s, test natural + reverse + a few random
                    orderings = [tuple(range(s)), tuple(range(s-1, -1, -1))]
                    import random
                    rng = random.Random(42)
                    for _ in range(min(20, s * s)):
                        perm = list(range(s))
                        rng.shuffle(perm)
                        orderings.append(tuple(perm))
                    orderings = list(set(orderings))

                any_ordering_passes = False
                all_orderings_pass = True
                best_min_margin_this_sv = None

                for perm in orderings:
                    # Process children in this order
                    # Initial: E_acc = (1+x)^ell, J_acc = [1]
                    E_acc = [1]
                    for _ in range(ell):
                        E_acc = _polymul(E_acc, [1, 1])
                    J_acc = [1]

                    ordering_passes = True
                    ord_min_margin = None

                    for step_idx, ci in enumerate(perm):
                        g_c, h_c = factors[ci]

                        if step_idx == 0:
                            # First child: this is the s=1 step
                            # Update E_acc, J_acc and continue
                            I_c = _polyadd(g_c, [0] + h_c)
                            E_acc = _polymul(E_acc, I_c)
                            J_acc = _polymul(J_acc, g_c)
                            continue

                        # Step ≥ 2: check ⋆
                        total_steps += 1
                        passes, margin, wk = check_star_at_step(
                            E_acc, J_acc, g_c, h_c, mode_I
                        )
                        if not passes:
                            prefix_star_fails += 1
                            ordering_passes = False

                        if margin is not None:
                            if ord_min_margin is None or margin < ord_min_margin:
                                ord_min_margin = margin

                        # Also check full range
                        f_passes, f_margin, f_wk = check_star_full_range(
                            E_acc, J_acc, g_c, h_c, len(I_poly) - 1
                        )
                        if not f_passes:
                            full_star_fails += 1

                        if f_margin is not None:
                            if min_full_margin is None or f_margin < min_full_margin:
                                min_full_margin = f_margin
                                min_full_info = (n, root, perm, step_idx, ci, f_wk)

                        # Update accumulators for next step
                        I_c = _polyadd(g_c, [0] + h_c)
                        E_acc = _polymul(E_acc, I_c)
                        J_acc = _polymul(J_acc, g_c)

                    if ordering_passes:
                        any_ordering_passes = True
                    else:
                        all_orderings_pass = False

                    if ord_min_margin is not None:
                        if best_min_margin_this_sv is None or ord_min_margin > best_min_margin_this_sv:
                            best_min_margin_this_sv = ord_min_margin

                if any_ordering_passes:
                    s_prefix_pass[s] += 1
                else:
                    sv_no_good_ordering += 1
                    s_prefix_fail[s] += 1
                    n_pfails += 1

                if all_orderings_pass:
                    sv_all_orderings_good += 1

                if best_min_margin_this_sv is not None:
                    if min_prefix_margin is None or best_min_margin_this_sv < min_prefix_margin:
                        min_prefix_margin = best_min_margin_this_sv
                        min_prefix_info = (n, root, s, ell)

                    if s not in s_min_margin or best_min_margin_this_sv < s_min_margin[s]:
                        s_min_margin[s] = best_min_margin_this_sv

        proc.wait()
        total_trees += tree_count
        total_sv += n_sv

        print(f"n={n:2d}: {tree_count:>8d} trees | s≥2 SVs={n_sv:>7d} | "
              f"no_good_ord={n_pfails}")

    print()
    print("=" * 70)
    print(f"SUMMARY n={args.min_n}..{args.max_n}")
    print("=" * 70)
    print(f"Trees:                     {total_trees:>12d}")
    print(f"Support vertices (s≥2):    {total_sv:>12d}")
    print(f"Incremental steps checked: {total_steps:>12d}")
    print()
    print(f"PREFIX ⋆ failures (step-level): {prefix_star_fails}")
    print(f"FULL   ⋆ failures (step-level): {full_star_fails}")
    print()
    print(f"SVs where NO ordering passes prefix ⋆: {sv_no_good_ordering}")
    print(f"SVs where ALL orderings pass prefix ⋆: {sv_all_orderings_good}")
    print()
    print(f"Best (highest) min prefix margin across SVs: {min_prefix_margin}")
    if min_prefix_info:
        print(f"  at n={min_prefix_info[0]}, root={min_prefix_info[1]}, "
              f"s={min_prefix_info[2]}, ell={min_prefix_info[3]}")
    print(f"Min full-range margin: {min_full_margin}")
    if min_full_info:
        print(f"  at n={min_full_info[0]}, root={min_full_info[1]}, "
              f"perm={min_full_info[2]}, step={min_full_info[3]}")

    print()
    print("Per-s breakdown:")
    for s in sorted(s_counts.keys()):
        total_s = s_counts[s]
        pass_s = s_prefix_pass.get(s, 0)
        fail_s = s_prefix_fail.get(s, 0)
        margin_s = s_min_margin.get(s, "N/A")
        print(f"  s={s}: {total_s:>7d} SVs | pass={pass_s:>7d} fail={fail_s} | "
              f"min_margin={margin_s}")


if __name__ == '__main__':
    main()
