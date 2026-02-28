#!/usr/bin/env python3
"""Verify GPT Instance 1's Condition C reformulation claims."""
import sys, subprocess
sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul, _polyadd
import networkx as nx
import numpy as np


def rooted_dp(T, root):
    """Return dict {v: (dp0, dp1)} for tree T rooted at root.
    dp0[v] = poly for subtree(v), v excluded
    dp1[v] = poly for subtree(v), v included (leading 0 = x * ...)
    """
    n = T.number_of_nodes()
    parent = {root: -1}
    children = {v: [] for v in T.nodes()}
    visited = {root}
    queue = [root]
    head = 0
    while head < len(queue):
        v = queue[head]; head += 1
        for u in T.neighbors(v):
            if u not in visited:
                visited.add(u)
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

    dp0 = {}
    dp1 = {}
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            prod = [1]
            for c in children[v]:
                prod = _polymul(prod, _polyadd(dp0[c], dp1[c]))
            dp0[v] = prod
            prod = [1]
            for c in children[v]:
                prod = _polymul(prod, dp0[c])
            dp1[v] = [0] + prod
    return dp0, dp1, children


def to_np(poly, length):
    """Pad polynomial list to numpy array of given length."""
    a = np.zeros(length, dtype=float)
    a[:len(poly)] = poly
    return a


def test_gpt_claims():
    """Test on P_4 (path on 4 vertices) at support vertex 1."""
    # Build P_4: 0-1-2-3
    T = nx.path_graph(4)

    A = np.array([1, 2], dtype=float)  # I_2 = 1+2x
    B = np.array([1, 1], dtype=float)  # E_2 = 1+x

    print("=== P_4 at support vertex 1 (leaf child 0) ===")
    print(f"A = {A}  (I_2 = 1+2x)")
    print(f"B = {B}  (E_2 = 1+x)")

    maxdeg = max(len(A), len(B))
    a = np.zeros(maxdeg + 2)
    b = np.zeros(maxdeg + 2)
    a[:len(A)] = A
    b[:len(B)] = B

    print("\n--- Coefficient check ---")
    for i in range(maxdeg + 1):
        print(f"  a_{i} = {a[i]},  b_{i} = {b[i]}")

    e = np.zeros(maxdeg + 2)
    for k in range(maxdeg + 2):
        e[k] = a[k] + (a[k-1] if k > 0 else 0)

    print(f"\ne (coeffs of (1+x)*A): {e[:maxdeg+2]}")

    print("\n--- k=1 test ---")
    k = 1

    d_k_ours = a[k+1]*b[k] - a[k]*b[k+1]
    d_km1_ours = a[k]*b[k-1] - a[k-1]*b[k]

    d_k_gpt = a[k]*b[k+1] - a[k+1]*b[k]
    d_km1_gpt = a[k-1]*b[k] - a[k]*b[k-1]

    c_k = b[k]**2 - b[k-1]*b[k+1]

    print(f"d_1 (ours) = {d_k_ours}")
    print(f"d_0 (ours) = {d_km1_ours}")
    print(f"d_1 (GPT)  = {d_k_gpt}")
    print(f"d_0 (GPT)  = {d_km1_gpt}")
    print(f"c_1 = {c_k}")

    SCC_ours = b[k-1]*d_k_ours + b[k]*d_km1_ours + a[k-1]*c_k
    print(f"\nSCC (our sign) = b_0*d_1 + b_1*d_0 + a_0*c_1 = {b[k-1]}*{d_k_ours} + {b[k]}*{d_km1_ours} + {a[k-1]}*{c_k} = {SCC_ours}")

    SCC_gpt = b[k-1]*d_k_gpt + b[k]*d_km1_gpt + a[k-1]*c_k
    print(f"S_k (GPT sign) = b_0*d_1 + b_1*d_0 + a_0*c_1 = {b[k-1]}*{d_k_gpt} + {b[k]}*{d_km1_gpt} + {a[k-1]}*{c_k} = {SCC_gpt}")

    Delta_k = e[k+1]*b[k] - e[k]*b[k+1]
    print(f"\nDelta_1 = e_2*b_1 - e_1*b_2 = {e[k+1]}*{b[k]} - {e[k]}*{b[k+1]} = {Delta_k}")
    print(f"b_0 * Delta_1 = {b[k-1] * Delta_k}")

    print(f"\nSCC (our sign) == b_0*Delta_1? {SCC_ours == b[k-1]*Delta_k}")
    print(f"S_k (GPT sign) == b_0*Delta_1? {SCC_gpt == b[k-1]*Delta_k}")

    # --- Systematic test ---
    print("\n\n=== Systematic test on all trees n <= 12 ===")

    sign_error_count = 0
    our_identity_fails = 0
    total_checks = 0

    for n in range(3, 13):
        result = subprocess.run(
            ['/opt/homebrew/bin/geng', '-q', str(n), f'{n-1}:{n-1}', '-c'],
            capture_output=True
        )
        for line in result.stdout.strip().split(b'\n'):
            if not line:
                continue
            T = nx.from_graph6_bytes(line)

            for r in T.nodes():
                leaves = [v for v in T.neighbors(r) if T.degree(v) == 1]
                if not leaves:
                    continue
                non_leaf_children = [v for v in T.neighbors(r) if T.degree(v) > 1]
                if not non_leaf_children:
                    continue

                dp0, dp1, children = rooted_dp(T, r)

                # A = prod_{c non-leaf child of r} I_c,  B = prod_{c} E_c
                A_poly = np.array([1.0])
                B_poly = np.array([1.0])
                for c in non_leaf_children:
                    Ic = np.array(_polyadd(dp0[c], dp1[c]), dtype=float)
                    Ec = np.array(dp0[c], dtype=float)
                    A_poly = np.convolve(A_poly, Ic)
                    B_poly = np.convolve(B_poly, Ec)

                maxd = max(len(A_poly), len(B_poly))
                aa = np.zeros(maxd + 2)
                bb = np.zeros(maxd + 2)
                aa[:len(A_poly)] = A_poly
                bb[:len(B_poly)] = B_poly

                ee = np.zeros(maxd + 2)
                for kk in range(maxd + 2):
                    ee[kk] = aa[kk] + (aa[kk-1] if kk > 0 else 0)

                for k in range(1, maxd + 1):
                    if bb[k-1] <= 0:
                        continue
                    total_checks += 1

                    dk_ours = aa[k+1]*bb[k] - aa[k]*bb[k+1]
                    dkm1_ours = aa[k]*bb[k-1] - aa[k-1]*bb[k]
                    ck = bb[k]**2 - bb[k-1]*bb[k+1]

                    SCC = bb[k-1]*dk_ours + bb[k]*dkm1_ours + aa[k-1]*ck
                    Delta = ee[k+1]*bb[k] - ee[k]*bb[k+1]
                    bDelta = bb[k-1]*Delta

                    if abs(SCC - bDelta) > 1e-6:
                        our_identity_fails += 1

                    dk_gpt = aa[k]*bb[k+1] - aa[k+1]*bb[k]
                    dkm1_gpt = aa[k-1]*bb[k] - aa[k]*bb[k-1]
                    Sk_gpt = bb[k-1]*dk_gpt + bb[k]*dkm1_gpt + aa[k-1]*ck

                    if abs(Sk_gpt - bDelta) > 1e-6:
                        sign_error_count += 1

    print(f"Total checks: {total_checks}")
    print(f"Our identity (SCC = b_{{k-1}}*Delta_k) failures: {our_identity_fails}")
    print(f"GPT's S_k (flipped d) = b_{{k-1}}*Delta_k failures: {sign_error_count}")
    print(f"\nConclusion: GPT's sign for d_k is {'WRONG' if sign_error_count > 0 else 'correct'}")
    print(f"Our identity is {'CORRECT' if our_identity_fails == 0 else 'WRONG'}")

    # What does GPT's S_k actually equal?
    print("\n\n=== What does GPT's S_k equal? ===")
    k = 1
    aa = np.zeros(4)
    bb = np.zeros(4)
    aa[:2] = [1, 2]
    bb[:2] = [1, 1]
    ee = np.zeros(4)
    for kk in range(4):
        ee[kk] = aa[kk] + (aa[kk-1] if kk > 0 else 0)

    dk_gpt = aa[k]*bb[k+1] - aa[k+1]*bb[k]
    dkm1_gpt = aa[k-1]*bb[k] - aa[k]*bb[k-1]
    ck = bb[k]**2 - bb[k-1]*bb[k+1]
    Sk_gpt = bb[k-1]*dk_gpt + bb[k]*dkm1_gpt + aa[k-1]*ck

    two_ac_minus_bDelta = 2*aa[k-1]*ck - bb[k-1]*(ee[k+1]*bb[k] - ee[k]*bb[k+1])
    print(f"P_4, k=1: GPT's S_k = {Sk_gpt}")
    print(f"P_4, k=1: 2*a_0*c_1 - b_0*Delta_1 = {two_ac_minus_bDelta}")


test_gpt_claims()
