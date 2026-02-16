"""Check log-concavity of A(x) and I(T') for subdivisions of all trees.

For each tree T and each edge uv:
  T' = subdivision of T at edge uv (insert new vertex w between u,v)
  I(T') = I(T) + A(x)
  where A(x) = Q_u * Q_v + x * P_u * P_v
  (Q_u = x * R_u, i.e., R_u shifted up by 1)

We check:
  1. Is A(x) log-concave?
  2. Is I(T') log-concave?
  3. Is I(T') unimodal?
"""

import subprocess
import sys
import time
from collections import deque

from indpoly import (
    _polyadd,
    _polymul,
    independence_poly,
    is_log_concave,
    is_unimodal,
)

GENG = "/opt/homebrew/bin/geng"


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode("ascii")]
    n = data[0]
    idx = 1
    adj = [[] for _ in range(n)]
    bit = 5
    word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5
                idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj


def get_edges(n, adj):
    """Return list of edges (u, v) with u < v."""
    edges = []
    for u in range(n):
        for v in adj[u]:
            if u < v:
                edges.append((u, v))
    return edges


def split_at_edge(n, adj, u, v):
    """Split tree at edge uv into two vertex sets (BFS from u avoiding v)."""
    A = set()
    queue = deque([u])
    A.add(u)
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y not in A and y != v:
                A.add(y)
                queue.append(y)
    return A, set(range(n)) - A


def forest_poly(n, adj, vertices):
    """Compute independence polynomial of the induced subgraph on vertices.

    For a forest (subgraph of a tree minus some vertices), this works
    by building the restricted adjacency list and calling independence_poly.
    """
    if not vertices:
        return [1]
    vlist = sorted(vertices)
    k = len(vlist)
    mapping = {old: new for new, old in enumerate(vlist)}
    sub_adj = [[] for _ in range(k)]
    for u in vlist:
        for w in adj[u]:
            if w in vertices:
                sub_adj[mapping[u]].append(mapping[w])
    return independence_poly(k, sub_adj)


def compute_A(n, adj, u, v):
    """Compute A(x) for subdividing edge uv.

    Split at edge uv:
      Side_u contains u, Side_v contains v.
      P_u = IS poly of Side_u with u included (rooted at u, u in set)
      R_u = IS poly of Side_u with u excluded
      Similarly P_v, R_v.

    Then A(x) = Q_u * Q_v + x * P_u * P_v
    where Q_u = x * R_u (shift R_u by 1).
    """
    side_u, side_v = split_at_edge(n, adj, u, v)

    # Compute rooted IS polys for side_u rooted at u
    P_u, R_u = rooted_is_poly(adj, side_u, u)
    P_v, R_v = rooted_is_poly(adj, side_v, v)

    # Q_u = x * R_u = [0] + R_u
    Q_u = [0] + R_u
    Q_v = [0] + R_v

    # A = Q_u * Q_v + x * P_u * P_v
    term1 = _polymul(Q_u, Q_v)
    term2 = [0] + _polymul(P_u, P_v)  # x * P_u * P_v
    A = _polyadd(term1, term2)
    return A


def rooted_is_poly(adj, vertices, root):
    """Compute IS polynomial of forest on `vertices`, rooted at `root`.

    Returns (P, R) where:
      P = poly counting IS that include root
      R = poly counting IS that exclude root
    """
    if len(vertices) == 1:
        return ([1], [1])  # P=[1] (just root, size 1 -> coeff of x^1... wait)

    # Actually: P counts IS including root, R counts IS excluding root.
    # For a single vertex: P = x (one IS of size 1), R = 1 (empty set).
    # But we store as coefficient lists: P = [0, 1], R = [1].
    # Hmm, let's be careful. The IS polynomial is sum over IS S of x^|S|.
    # P(x) = sum_{S: root in S} x^|S|, R(x) = sum_{S: root not in S} x^|S|.
    # Single vertex: P(x) = x -> [0, 1]; R(x) = 1 -> [1].

    # Build subtree adj restricted to vertices
    vset = set(vertices)

    # BFS to get parent-child structure
    parent = {root: -1}
    order = []
    queue = deque([root])
    while queue:
        x = queue.popleft()
        order.append(x)
        for y in adj[x]:
            if y in vset and y not in parent:
                parent[y] = x
                queue.append(y)

    # DP bottom-up
    # dp_in[v] = poly for subtree(v) with v included
    # dp_out[v] = poly for subtree(v) with v excluded
    dp_in = {}
    dp_out = {}

    for v in reversed(order):
        children = [y for y in adj[v] if y in vset and parent.get(y) == v]
        if not children:
            dp_in[v] = [0, 1]  # x
            dp_out[v] = [1]  # 1
        else:
            # dp_in[v] = x * prod(dp_out[c] for c in children)
            prod_out = [1]
            for c in children:
                prod_out = _polymul(prod_out, dp_out[c])
            dp_in[v] = [0] + prod_out  # multiply by x

            # dp_out[v] = prod(dp_in[c] + dp_out[c] for c in children)
            prod_both = [1]
            for c in children:
                prod_both = _polymul(prod_both, _polyadd(dp_in[c], dp_out[c]))
            dp_out[v] = prod_both

    return dp_in[root], dp_out[root]


def enumerate_trees(n):
    """Enumerate all non-isomorphic trees on n vertices using geng."""
    cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for line in proc.stdout:
        yield line.decode("ascii").strip()
    proc.wait()


def check_tree(g6):
    """Check all edge subdivisions of a tree.

    Returns (n_edges, a_lc_fails, it_lc_fails, it_uni_fails, details).
    """
    n, adj = parse_graph6(g6)
    edges = get_edges(n, adj)

    I_T = independence_poly(n, adj)

    a_lc_fails = 0
    it_lc_fails = 0
    it_uni_fails = 0
    details = []

    for u, v in edges:
        A = compute_A(n, adj, u, v)
        I_Tp = _polyadd(I_T, A)

        a_lc = is_log_concave(A)
        it_lc = is_log_concave(I_Tp)
        it_uni = is_unimodal(I_Tp)

        if not a_lc:
            a_lc_fails += 1
        if not it_lc:
            it_lc_fails += 1
        if not it_uni:
            it_uni_fails += 1
            details.append((g6, u, v, I_T, A, I_Tp))

    return len(edges), a_lc_fails, it_lc_fails, it_uni_fails, details


def main():
    max_n = 19
    if len(sys.argv) > 1:
        max_n = int(sys.argv[1])

    print(f"Subdivision A(x) properties check, n up to {max_n}", flush=True)
    print("=" * 72, flush=True)

    t0 = time.time()

    for n in range(3, max_n + 1):
        tn = time.time()
        total_trees = 0
        total_edges = 0
        total_a_lc_fails = 0
        total_it_lc_fails = 0
        total_it_uni_fails = 0
        all_details = []

        for g6 in enumerate_trees(n):
            total_trees += 1
            n_edges, a_lc, it_lc, it_uni, details = check_tree(g6)
            total_edges += n_edges
            total_a_lc_fails += a_lc
            total_it_lc_fails += it_lc
            total_it_uni_fails += it_uni
            all_details.extend(details)

        elapsed = time.time() - tn
        print(
            f"n={n:2d} | trees={total_trees:>10,} | edges checked={total_edges:>12,} "
            f"| A lc_fail={total_a_lc_fails} | I(T') lc_fail={total_it_lc_fails} "
            f"| I(T') uni_fail={total_it_uni_fails} | {elapsed:.1f}s",
            flush=True,
        )

        if all_details:
            print(f"  *** FAILURES at n={n}:", flush=True)
            for g6, u, v, IT, A, ITp in all_details[:5]:
                print(f"      tree={g6} edge=({u},{v})", flush=True)
                print(f"      I(T)  = {IT}", flush=True)
                print(f"      A(x)  = {A}", flush=True)
                print(f"      I(T') = {ITp}", flush=True)

    total_elapsed = time.time() - t0
    print("=" * 72, flush=True)
    print(f"Total time: {total_elapsed:.1f}s", flush=True)


if __name__ == "__main__":
    main()
