#!/usr/bin/env python3
"""Shared helpers for attack4 structural-induction verification scripts."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from typing import Iterator

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import _polyadd, _polymul, independence_poly


@dataclass
class BridgeDecomposition:
    """Canonical degree-2 bridge decomposition data for one tree."""

    n: int
    g6: str
    adj: list[list[int]]

    poly_t: list[int]
    m_t: int

    leaf: int
    support: int
    u: int

    b_adj: list[list[int]]
    b_poly: list[int]
    u_in_b: int

    children_of_u: list[int]
    child_sizes: list[int]
    f_list: list[list[int]]
    g_list: list[list[int]]
    h_list: list[list[int]]

    p_poly: list[int]
    q_poly: list[int]
    g_poly: list[int]


def get_coeff(poly: list[int], k: int) -> int:
    """Return coefficient [x^k]poly, with 0 outside range."""
    if 0 <= k < len(poly):
        return poly[k]
    return 0


def poly_prod(polys: list[list[int]]) -> list[int]:
    """Multiply a list of coefficient polynomials."""
    out = [1]
    for poly in polys:
        out = _polymul(out, poly)
    return out


def iter_tree_g6(min_n: int, max_n: int, geng: str) -> Iterator[tuple[int, list[list[int]], str]]:
    """Yield (n, adjacency, graph6) for all trees in the size range."""
    for n in range(min_n, max_n + 1):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert proc.stdout is not None
        for raw in proc.stdout:
            nn, adj = parse_graph6(raw)
            yield nn, adj, raw.decode("ascii").strip()


def choose_canonical_deg2_leaf(adj: list[list[int]]) -> tuple[int, int, int] | None:
    """Return canonical (leaf, support, u) with deg(support)=2, or None.

    Canonical choice matches prior scans: pick a leaf whose parent has minimum
    degree; break ties by leaf index.
    """
    n = len(adj)
    if n < 2:
        return None

    deg = [len(nei) for nei in adj]
    leaves = [v for v in range(n) if deg[v] == 1]
    if not leaves:
        return None

    min_parent_deg = min(deg[adj[l][0]] for l in leaves)
    leaf = min(l for l in leaves if deg[adj[l][0]] == min_parent_deg)
    support = adj[leaf][0]
    if deg[support] != 2:
        return None

    u = adj[support][0] if adj[support][1] == leaf else adj[support][1]
    return leaf, support, u


def _rooted_child_factors(adj: list[list[int]], root: int) -> tuple[
    list[int],
    list[list[int]],
    list[list[int]],
    list[int],
]:
    """Compute rooted DP and subtree sizes.

    Returns:
      children_of_root,
      dp0,
      dp1,
      subtree_size
    """
    n = len(adj)
    if n == 0:
        return [], [], [], []

    children = [[] for _ in range(n)]
    seen = [False] * n
    seen[root] = True
    bfs_queue = [root]
    head = 0
    while head < len(bfs_queue):
        v = bfs_queue[head]
        head += 1
        for w in adj[v]:
            if not seen[w]:
                seen[w] = True
                children[v].append(w)
                bfs_queue.append(w)

    order: list[int] = []
    stack: list[tuple[int, bool]] = [(root, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    dp0: list[list[int]] = [[] for _ in range(n)]
    dp1: list[list[int]] = [[] for _ in range(n)]
    subtree_size = [0] * n

    for v in order:
        if not children[v]:
            dp0[v] = [1]
            dp1[v] = [0, 1]
            subtree_size[v] = 1
            continue

        prod0 = [1]
        prod1 = [1]
        size_v = 1
        for c in children[v]:
            prod0 = _polymul(prod0, _polyadd(dp0[c], dp1[c]))
            prod1 = _polymul(prod1, dp0[c])
            size_v += subtree_size[c]

        dp0[v] = prod0
        dp1[v] = [0] + prod1
        subtree_size[v] = size_v

    return children[root], dp0, dp1, subtree_size


def bridge_decomposition(
    n: int,
    adj: list[list[int]],
    g6: str,
    require_dleaf: bool = True,
) -> BridgeDecomposition | None:
    """Build canonical bridge decomposition data for a tree.

    Returns None when the tree is filtered out (`d_leaf<=1` requested but not
    satisfied) or when canonical leaf does not have degree-2 support.
    """
    if require_dleaf and not is_dleaf_le_1(n, adj):
        return None

    leaf_triplet = choose_canonical_deg2_leaf(adj)
    if leaf_triplet is None:
        return None

    leaf, support, u = leaf_triplet

    poly_t = independence_poly(n, adj)
    m_t = mode_index_leftmost(poly_t)

    b_adj = remove_vertices(adj, {leaf, support})
    if not b_adj:
        return None

    keep = [v for v in range(n) if v not in {leaf, support}]
    idx_map = {old: new for new, old in enumerate(keep)}
    u_in_b = idx_map[u]

    children_of_u, dp0, dp1, subtree_size = _rooted_child_factors(b_adj, u_in_b)

    p_poly = dp0[u_in_b]
    q_poly = dp1[u_in_b]
    b_poly = _polyadd(p_poly, q_poly)

    f_list: list[list[int]] = []
    g_list: list[list[int]] = []
    h_list: list[list[int]] = []
    child_sizes: list[int] = []

    for c in children_of_u:
        g = dp0[c]
        h = dp1[c]
        f = _polyadd(g, h)
        f_list.append(f)
        g_list.append(g)
        h_list.append(h)
        child_sizes.append(subtree_size[c])

    g_poly = poly_prod(g_list) if g_list else [1]

    # Sanity checks for the canonical decomposition.
    if q_poly != [0] + g_poly:
        raise RuntimeError("Identity failure: Q != x*G")
    if p_poly != (poly_prod(f_list) if f_list else [1]):
        raise RuntimeError("Identity failure: P != product_i f_i")

    return BridgeDecomposition(
        n=n,
        g6=g6,
        adj=adj,
        poly_t=poly_t,
        m_t=m_t,
        leaf=leaf,
        support=support,
        u=u,
        b_adj=b_adj,
        b_poly=b_poly,
        u_in_b=u_in_b,
        children_of_u=children_of_u,
        child_sizes=child_sizes,
        f_list=f_list,
        g_list=g_list,
        h_list=h_list,
        p_poly=p_poly,
        q_poly=q_poly,
        g_poly=g_poly,
    )
