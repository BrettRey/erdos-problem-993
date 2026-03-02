#!/usr/bin/env python3
"""Verify GPT 5.2 Pro's claimed counterexample to leaf-augmentation.

Tree: T_{3,4} "broom" — root v with 3 children w_i,
each w_i has 4 children x_i^j, each x_i^j has 1 leaf y_i^j.
n = 1 + 3 + 12 + 12 = 28.

Claim: rooting at r = x_1^1 (degree-2 vertex) gives g_14 < 0,
i.e., I_tilde_{15} * E_{14} - I_tilde_{14} * E_{15} = -334 < 0.
"""
import sys
sys.path.insert(0, ".")
from indpoly import _polymul, _polyadd


def xshift(p):
    return [0] + list(p)


def coeff(p, k):
    return p[k] if 0 <= k < len(p) else 0


def delta_k(F, G, k):
    return coeff(F, k+1) * coeff(G, k) - coeff(F, k) * coeff(G, k+1)


def build_broom(m, t):
    """Build T_{m,t} broom tree as adjacency list.

    Returns (n, adj, vertex_labels) where vertex_labels maps index to role.
    """
    n = 1 + m + m*t + m*t
    adj = [[] for _ in range(n)]
    labels = {}

    v = 0  # root
    labels[v] = 'v'

    next_id = 1
    w_ids = []
    for i in range(m):
        w = next_id; next_id += 1
        w_ids.append(w)
        labels[w] = f'w_{i}'
        adj[v].append(w)
        adj[w].append(v)

    x_ids = []
    for i in range(m):
        for j in range(t):
            x = next_id; next_id += 1
            x_ids.append(x)
            labels[x] = f'x_{i}_{j}'
            adj[w_ids[i]].append(x)
            adj[x].append(w_ids[i])

    y_ids = []
    for i in range(m):
        for j in range(t):
            y = next_id; next_id += 1
            y_ids.append(y)
            labels[y] = f'y_{i}_{j}'
            x_idx = i * t + j
            adj[x_ids[x_idx]].append(y)
            adj[y].append(x_ids[x_idx])

    assert next_id == n
    return n, adj, labels, x_ids, w_ids, y_ids


def dp_rooted(n, adj, root):
    """Standard tree DP. Returns dp0[v], dp1s[v] for all v."""
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

    dp0, dp1s = {}, {}
    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1s[v] = [1]
        else:
            pS, pE = [1], [1]
            for c in children[v]:
                sc = _polyadd(dp0[c], xshift(dp1s[c]))
                pS = _polymul(pS, sc)
                pE = _polymul(pE, dp0[c])
            dp0[v] = pS
            dp1s[v] = pE

    return dp0, dp1s


def main():
    m, t = 3, 4
    n, adj, labels, x_ids, w_ids, y_ids = build_broom(m, t)
    print(f"Tree T_{{{m},{t}}}: n={n}")
    print(f"Degrees: {sorted([len(adj[v]) for v in range(n)], reverse=True)}")

    # Root at x_1^1 (= x_ids[0], which is the first x vertex under w_0)
    r = x_ids[0]
    print(f"\nRooting at r={r} (label={labels[r]}), degree={len(adj[r])}")

    dp0, dp1s = dp_rooted(n, adj, r)

    E = dp0[r]  # exclude-root
    J = dp1s[r]  # include-root / x
    I_poly = _polyadd(E, xshift(J))  # full IS poly
    I_tilde = _polyadd(I_poly, xshift(E))  # I + xE

    print(f"\nE coefficients (len={len(E)}):")
    print(E)
    print(f"\nI_tilde coefficients (len={len(I_tilde)}):")
    print(I_tilde)

    # Compare with GPT's claimed values
    gpt_E = [1, 27, 326, 2340, 11182, 37693, 92484, 167669, 225112, 221520,
             155643, 74182, 21662, 3030, 52, 1]
    gpt_Itilde = [1, 28, 353, 2666, 13522, 48875, 130177, 260153, 392781, 446632,
                  377163, 229825, 95844, 24692, 3090, 53, 1]

    print(f"\nGPT E matches: {E == gpt_E}")
    print(f"GPT I_tilde matches: {I_tilde == gpt_Itilde}")

    # Check all g_k = d_k(I_tilde, E)
    print(f"\nLR minors g_k = d_k(I_tilde, E):")
    max_deg = max(len(I_tilde), len(E))
    any_negative = False
    for k in range(max_deg):
        gk = delta_k(I_tilde, E, k)
        marker = " <-- NEGATIVE!" if gk < 0 else ""
        if gk < 0:
            any_negative = True
        print(f"  k={k:2d}: g_k = {gk}{marker}")

    print(f"\nAny negative g_k: {any_negative}")

    # Also check at OTHER rootings of this tree
    print(f"\n{'='*60}")
    print("Checking ALL rootings of T_{3,4}:")
    print(f"{'='*60}")

    for r2 in range(n):
        dp0_2, dp1s_2 = dp_rooted(n, adj, r2)
        E2 = dp0_2[r2]
        J2 = dp1s_2[r2]
        I2 = _polyadd(E2, xshift(J2))
        It2 = _polyadd(I2, xshift(E2))

        max_d = max(len(It2), len(E2))
        fails = []
        for k in range(max_d):
            gk = delta_k(It2, E2, k)
            if gk < 0:
                fails.append((k, gk))

        if fails:
            print(f"  r={r2:2d} ({labels[r2]:>8s}), deg={len(adj[r2])}: "
                  f"FAILS at {len(fails)} indices: {fails[:5]}")

    # Check: does it fail at support vertices?
    print(f"\nSupport vertices (adjacent to a leaf):")
    for r2 in range(n):
        is_support = any(len(adj[u]) == 1 for u in adj[r2])
        if not is_support:
            continue
        dp0_2, dp1s_2 = dp_rooted(n, adj, r2)
        E2 = dp0_2[r2]
        J2 = dp1s_2[r2]
        I2 = _polyadd(E2, xshift(J2))
        It2 = _polyadd(I2, xshift(E2))

        max_d = max(len(It2), len(E2))
        fails = []
        for k in range(max_d):
            gk = delta_k(It2, E2, k)
            if gk < 0:
                fails.append((k, gk))

        status = "FAIL" if fails else "OK"
        print(f"  r={r2:2d} ({labels[r2]:>8s}), deg={len(adj[r2])}: {status}"
              + (f" at {fails[:3]}" if fails else ""))


if __name__ == '__main__':
    main()
