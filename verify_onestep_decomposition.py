#!/usr/bin/env python3
"""Verify the one-step decomposition formula for E >= J ratio dominance.

At a support vertex r with ell leaf neighbors and non-leaf subtrees T_1,...,T_s,
the incremental DP builds:
  E^{(t)} = E^{(t-1)} * I_t  where I_t = E_t + x*J_t
  J^{(t)} = J^{(t-1)} * E_t

Define at each step:
  A = E_old * E_t   (convolution)
  B = E_old * J_t
  C = J_old * E_t

Then E_new = A + x*B, J_new = C.

CLAIM:
  Delta_k(A + xB, C) = Delta_k(A, C) + (B_k * C_k - B_{k-1} * C_{k+1})

where Delta_k(F, G) = F_{k+1} * G_k - F_k * G_{k+1}.

This script verifies the identity, checks Delta_k(A,C) >= 0 (Karlin main term),
and profiles the correction term T_k = B_k*C_k - B_{k-1}*C_{k+1}.
"""

import subprocess
import sys
import time
from collections import defaultdict

from indpoly import _polymul, _polyadd


# ---------------------------------------------------------------------------
# Polynomial helpers (integer arithmetic, lists)
# ---------------------------------------------------------------------------

def coeff(poly, k):
    """Safely get coefficient at index k."""
    if 0 <= k < len(poly):
        return poly[k]
    return 0


def delta_k(F, G, k):
    """Compute Delta_k(F, G) = F_{k+1} * G_k - F_k * G_{k+1}."""
    return coeff(F, k + 1) * coeff(G, k) - coeff(F, k) * coeff(G, k + 1)


def xshift(poly):
    """Multiply polynomial by x: prepend a zero."""
    return [0] + list(poly)


# ---------------------------------------------------------------------------
# DP at a root, returning per-child E_c, J_c data
# ---------------------------------------------------------------------------

def dp_rooted_with_children(n, adj, root):
    """Root tree at 'root', compute DP, return children info.

    Returns
    -------
    children_of_root : list of int
    dp0 : dict  vertex -> list[int]  (E polynomial for subtree)
    dp1s : dict  vertex -> list[int]  (J = dp1/x for subtree)
    """
    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    visited[root] = True
    queue = [root]
    head = 0
    while head < len(queue):
        v = queue[head]
        head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    # Post-order
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

    dp0 = {}   # dp0[v] = E_v
    dp1s = {}  # dp1s[v] = J_v (dp1/x)

    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            prod_S = [1]
            prod_E = [1]
            for c in children[v]:
                s_c = _polyadd(dp0[c], xshift(dp1s[c]))
                prod_S = _polymul(prod_S, s_c)
                prod_E = _polymul(prod_E, dp0[c])
            dp0[v] = prod_S
            dp1s[v] = prod_E

    return children[root], dp0, dp1s


# ---------------------------------------------------------------------------
# Parse graph6 (text mode)
# ---------------------------------------------------------------------------

def parse_g6(g6):
    """Parse graph6 string to (n, adjacency list)."""
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


# ---------------------------------------------------------------------------
# Main verification
# ---------------------------------------------------------------------------

def main():
    max_n = 20
    geng = '/opt/homebrew/bin/geng'

    # Counters
    total_identity_checks = 0
    total_identity_fails = 0
    total_karlin_checks = 0
    total_karlin_fails = 0
    total_T_neg = 0
    total_T_checks = 0

    # Min D/|T| ratio by s-value (number of non-leaf children processed so far)
    min_ratio_by_s = defaultdict(lambda: float('inf'))
    global_min_ratio = float('inf')
    global_min_info = None

    # Track by n
    stats_by_n = {}

    t0 = time.time()

    for nn in range(3, max_n + 1):
        edges = nn - 1
        cmd = [geng, '-q', str(nn), f'{edges}:{edges}', '-c']
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        n_trees = 0
        n_support = 0
        n_identity = 0
        n_id_fail = 0
        n_karlin = 0
        n_karlin_fail = 0
        n_T_checks_local = 0
        n_T_neg_local = 0
        min_ratio_n = float('inf')

        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            n, adj = parse_g6(line)
            n_trees += 1

            # Find support vertices: vertices adjacent to at least one leaf
            degree = [len(adj[v]) for v in range(n)]
            for r in range(n):
                # Check if r is a support vertex
                leaf_children = []
                nonleaf_children = []
                for u in adj[r]:
                    if degree[u] == 1:
                        leaf_children.append(u)
                    else:
                        nonleaf_children.append(u)

                if not leaf_children:
                    continue  # not a support vertex

                n_support += 1
                ell = len(leaf_children)

                # Root at r, get DP data
                children_of_r, dp0, dp1s = dp_rooted_with_children(n, adj, r)

                # Identify which children of r are leaves vs non-leaves
                # (after rooting, the classification may differ from adjacency-based)
                root_leaf_children = []
                root_nonleaf_children = []
                for c in children_of_r:
                    if not [u for u in adj[c] if u != r]:
                        # c is a leaf in the rooted tree
                        root_leaf_children.append(c)
                    else:
                        root_nonleaf_children.append(c)

                ell_rooted = len(root_leaf_children)

                # Build E_old from leaf children: (1+x)^ell
                # J_old starts at [1]
                E_old = [1]
                for _ in range(ell_rooted):
                    E_old = _polymul(E_old, [1, 1])
                J_old = [1]

                # Process non-leaf children one at a time
                s_count = 0
                for c in root_nonleaf_children:
                    E_c = dp0[c]   # subtree E
                    J_c = dp1s[c]  # subtree J (dp1/x)
                    # I_c = E_c + x*J_c

                    # Compute A, B, C
                    A = _polymul(E_old, E_c)
                    B = _polymul(E_old, J_c)
                    C = _polymul(J_old, E_c)

                    # E_new = A + x*B
                    E_new = _polyadd(A, xshift(B))
                    # J_new = C
                    J_new = C

                    # Degree of relevant polynomials
                    max_deg = max(len(E_new), len(J_new), len(A), len(B), len(C))

                    s_count += 1

                    for k in range(max_deg):
                        # IDENTITY CHECK:
                        # Delta_k(E_new, J_new) = Delta_k(A, C) + (B_k*C_k - B_{k-1}*C_{k+1})
                        lhs = delta_k(E_new, J_new, k)
                        D_k = delta_k(A, C, k)
                        Bk = coeff(B, k)
                        Ck = coeff(C, k)
                        Bkm1 = coeff(B, k - 1)
                        Ckp1 = coeff(C, k + 1)
                        T_k = Bk * Ck - Bkm1 * Ckp1
                        rhs = D_k + T_k

                        n_identity += 1
                        total_identity_checks += 1
                        if lhs != rhs:
                            n_id_fail += 1
                            total_identity_fails += 1

                        # KARLIN CHECK: Delta_k(A, C) >= 0
                        n_karlin += 1
                        total_karlin_checks += 1
                        if D_k < 0:
                            n_karlin_fail += 1
                            total_karlin_fails += 1

                        # T_k sign tracking
                        n_T_checks_local += 1
                        total_T_checks += 1
                        if T_k < 0:
                            n_T_neg_local += 1
                            total_T_neg += 1
                            # Compute D/|T| ratio
                            if T_k != 0:
                                ratio = D_k / abs(T_k)
                                if ratio < min_ratio_by_s[s_count]:
                                    min_ratio_by_s[s_count] = ratio
                                if ratio < global_min_ratio:
                                    global_min_ratio = ratio
                                    global_min_info = {
                                        'n': n,
                                        'g6': line,
                                        'root': r,
                                        's': s_count,
                                        'k': k,
                                        'D_k': D_k,
                                        'T_k': T_k,
                                        'ratio': ratio,
                                    }
                                if ratio < min_ratio_n:
                                    min_ratio_n = ratio

                    # Update E_old, J_old for next stage
                    E_old = E_new
                    J_old = J_new

        proc.wait()
        elapsed = time.time() - t0

        stats_by_n[nn] = {
            'trees': n_trees,
            'support_rootings': n_support,
            'identity_checks': n_identity,
            'identity_fails': n_id_fail,
            'karlin_checks': n_karlin,
            'karlin_fails': n_karlin_fail,
            'T_checks': n_T_checks_local,
            'T_neg': n_T_neg_local,
            'min_D_over_T': min_ratio_n if min_ratio_n < float('inf') else None,
        }

        min_str = f"{min_ratio_n:.6f}" if min_ratio_n < float('inf') else "N/A"
        print(f"n={nn:2d}: {n_trees:>10,d} trees | "
              f"support={n_support:>10,d} | "
              f"id_fail={n_id_fail} karlin_fail={n_karlin_fail} | "
              f"T<0: {n_T_neg_local:>8,d}/{n_T_checks_local:>10,d} | "
              f"min D/|T|={min_str} | "
              f"{elapsed:.1f}s")

    elapsed = time.time() - t0
    print(f"\n{'='*80}")
    print(f"SUMMARY (n=3..{max_n}, {elapsed:.1f}s)")
    print(f"{'='*80}")
    print(f"Total identity checks: {total_identity_checks:,d}")
    print(f"Identity failures:     {total_identity_fails:,d}")
    print(f"Total Karlin checks:   {total_karlin_checks:,d}")
    print(f"Karlin failures:       {total_karlin_fails:,d}")
    print(f"Total T_k checks:      {total_T_checks:,d}")
    print(f"T_k < 0 count:         {total_T_neg:,d} "
          f"({100*total_T_neg/max(1,total_T_checks):.2f}%)")

    print(f"\nMin D/|T| ratio by s-value:")
    for s in sorted(min_ratio_by_s.keys()):
        print(f"  s={s}: {min_ratio_by_s[s]:.6f}")

    if global_min_ratio < float('inf'):
        print(f"\nGlobal minimum D/|T| ratio: {global_min_ratio:.6f}")
        if global_min_info:
            info = global_min_info
            print(f"  Achieved at: n={info['n']}, root={info['root']}, "
                  f"s={info['s']}, k={info['k']}")
            print(f"  D_k={info['D_k']}, T_k={info['T_k']}")
            print(f"  g6={info['g6'][:40]}")
    else:
        print(f"\nNo T_k < 0 cases found (correction term always non-negative).")


if __name__ == '__main__':
    main()
