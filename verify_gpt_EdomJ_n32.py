"""
Verify GPT's claim that E ≽ J (ratio dominance) fails at n=32.

Construction: support vertex r with:
  - 1 pendant leaf
  - 1 P_2 subtree (r -- a -- b)
  - 1 T_{3,4} broom subtree (v has 3 children w_i, each w_i has 4 children x_{i,j},
    each x_{i,j} has 1 leaf y_{i,j})

Total: r(1) + leaf(1) + a,b(2) + v(1) + 3 w's + 12 x's + 12 y's = 32 vertices.

E = dp[r][0] (exclude root), J = dp[r][1] (include root, before x-shift).
E ≽ J means d_k = E_{k+1}·J_k − E_k·J_{k+1} ≥ 0 for all k.
"""

import sys
sys.path.insert(0, '/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993')
from indpoly import _polymul
import networkx as nx
from collections import defaultdict


def build_tree():
    """Build the n=32 tree described by GPT."""
    G = nx.Graph()
    r = 'r'
    # Pendant leaf
    G.add_edge(r, 'leaf')
    # P_2 subtree
    G.add_edge(r, 'a')
    G.add_edge('a', 'b')
    # T_{3,4} broom subtree
    G.add_edge(r, 'v')
    for i in range(3):
        wi = f'w{i}'
        G.add_edge('v', wi)
        for j in range(4):
            xij = f'x{i}{j}'
            G.add_edge(wi, xij)
            yij = f'y{i}{j}'
            G.add_edge(xij, yij)
    return G, r


def compute_dp(G, root):
    """Compute dp[v] = (exc, inc) for tree G rooted at root.

    dp[v][0] = polynomial for IS in subtree(v), v excluded
    dp[v][1] = polynomial for IS in subtree(v), v included (WITHOUT the x factor)

    So the full IS polynomial is dp[root][0] + x * dp[root][1].
    """
    # BFS to get traversal order and children
    parent = {root: None}
    order = [root]
    queue = [root]
    while queue:
        u = queue.pop(0)
        for v in G.neighbors(u):
            if v not in parent:
                parent[v] = u
                order.append(v)
                queue.append(v)

    children = defaultdict(list)
    for v in order:
        if parent[v] is not None:
            children[parent[v]].append(v)

    # DP bottom-up
    dp = {}

    def polyadd(a, b):
        result = list(a) + [0] * max(0, len(b) - len(a))
        for i, val in enumerate(b):
            if i < len(result):
                result[i] += val
            else:
                result.append(val)
        return result

    for v in reversed(order):
        if not children[v]:  # leaf
            dp[v] = ([1], [1])  # exc = 1, inc = 1 (coeff of x^0 in the "included" factor)
        else:
            exc = [1]  # product of (dp[c][0] + x * dp[c][1]) for children
            inc = [1]  # product of dp[c][0] for children
            for c in children[v]:
                c_exc, c_inc = dp[c]
                # dp[c][0] + x * dp[c][1]: shift c_inc by 1, add to c_exc
                c_total = list(c_exc)
                while len(c_total) < len(c_inc) + 1:
                    c_total.append(0)
                for i, val in enumerate(c_inc):
                    c_total[i + 1] += val
                exc = _polymul(exc, c_total)
                inc = _polymul(inc, c_exc)
            dp[v] = (exc, inc)

    return dp


def main():
    G, r = build_tree()
    n = G.number_of_nodes()
    m = G.number_of_edges()
    print(f"Tree: {n} nodes, {m} edges")
    print(f"Is tree: {nx.is_tree(G)}")
    print(f"Degree of r: {G.degree(r)}")
    print(f"Neighbors of r: {list(G.neighbors(r))}")

    # Verify r is a support vertex (has at least one leaf child)
    leaf_children = [v for v in G.neighbors(r) if G.degree(v) == 1]
    print(f"Leaf children of r: {leaf_children} (count = {len(leaf_children)})")

    dp = compute_dp(G, r)
    E = dp[r][0]  # exclude root
    J = dp[r][1]  # include root (before x-shift)

    print(f"\nE (exclude root): degree {len(E)-1}, {len(E)} coefficients")
    print(f"J (include root):  degree {len(J)-1}, {len(J)} coefficients")
    print(f"\nE = {E}")
    print(f"J = {J}")

    # Cross-check: IS polynomial = E + x*J
    I_poly = list(E)
    while len(I_poly) < len(J) + 1:
        I_poly.append(0)
    for i, val in enumerate(J):
        I_poly[i + 1] += val

    print(f"\nIS polynomial (first 20 coeffs): {I_poly[:20]}")

    # Also cross-check with indpoly.independence_poly via networkx -> adjacency list
    from indpoly import independence_poly
    # Build adjacency list with integer labels
    node_list = list(G.nodes())
    node_to_int = {v: i for i, v in enumerate(node_list)}
    adj = [[] for _ in range(n)]
    for u, v in G.edges():
        adj[node_to_int[u]].append(node_to_int[v])
        adj[node_to_int[v]].append(node_to_int[u])
    I_check = independence_poly(n, adj)
    print(f"indpoly check    (first 20 coeffs): {I_check[:20]}")
    if I_poly == I_check:
        print("Cross-check: MATCH")
    else:
        print("Cross-check: MISMATCH!")
        print(f"  DP-based:  {I_poly}")
        print(f"  indpoly:   {I_check}")

    # Check unimodality
    from indpoly import is_unimodal
    print(f"\nIS polynomial unimodal: {is_unimodal(I_poly)}")

    # Mode
    mode_idx = max(range(len(I_poly)), key=lambda i: I_poly[i])
    print(f"Mode index: {mode_idx}, I[mode] = {I_poly[mode_idx]}")

    # Check E ≽ J: d_k = E_{k+1}*J_k - E_k*J_{k+1} >= 0 for all valid k
    print("\n" + "=" * 60)
    print("LR minors d_k = E_{k+1}*J_k - E_k*J_{k+1}:")
    print("=" * 60)
    max_k = max(len(E), len(J))
    failures = []
    for k in range(max_k):
        ek = E[k] if k < len(E) else 0
        ek1 = E[k + 1] if k + 1 < len(E) else 0
        jk = J[k] if k < len(J) else 0
        jk1 = J[k + 1] if k + 1 < len(J) else 0
        dk = ek1 * jk - ek * jk1
        marker = ""
        if dk < 0:
            failures.append((k, dk))
            marker = "  *** NEGATIVE ***"
        print(f"  k={k:2d}: E[{k+1}]={ek1:>15d}  J[{k}]={jk:>15d}  "
              f"E[{k}]={ek:>15d}  J[{k+1}]={jk1:>15d}  d_k = {dk}{marker}")

    print()
    if failures:
        print(f"*** E ≽ J FAILS at {len(failures)} index/indices: ***")
        for k, dk in failures:
            print(f"    k={k}: d_k = {dk}")
    else:
        print("*** E ≽ J HOLDS at all indices ***")

    # Also check if J is prefix-LC (helpful context)
    from indpoly import is_log_concave
    print(f"\nJ log-concave: {is_log_concave(J)}")
    print(f"E log-concave: {is_log_concave(E)}")


if __name__ == '__main__':
    main()
