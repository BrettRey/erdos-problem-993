#!/usr/bin/env python3
"""Check d_leaf ≤ 1 trees at large n for high mode.

Construct specific families of d_leaf ≤ 1 trees and check if any are
high-mode (mode > floor(n/3)+1).

Families:
1. Caterpillars: spine with whiskers at each/some vertices
2. Lobsters: spine with branches of length 2
3. Trees with many degree-3 branch points connected by paths
"""

import sys
sys.path.insert(0, ".")
from indpoly import independence_poly


def make_caterpillar(spine_length, whisker_positions=None):
    """Caterpillar: path with optional whiskers (leaf at each position).

    Returns (n, adj).
    """
    s = spine_length
    if whisker_positions is None:
        whisker_positions = list(range(s))  # whisker at every vertex

    n = s + len(whisker_positions)
    adj = [[] for _ in range(n)]

    # Spine edges: 0-1-2-..-(s-1)
    for i in range(s - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)

    # Whiskers: vertex s, s+1, ... attached to positions
    leaf_idx = s
    for pos in whisker_positions:
        adj[pos].append(leaf_idx)
        adj[leaf_idx].append(pos)
        leaf_idx += 1

    return n, adj


def make_lobster(spine_length, branch_positions=None):
    """Lobster: spine with branches of length 2 at some positions.

    Each branch is: spine_v -- intermediate -- leaf.
    d_leaf(spine_v) = 0 (intermediate is not a leaf).
    d_leaf(intermediate) = 1 (one leaf-child).
    So d_leaf ≤ 1 everywhere.
    """
    s = spine_length
    if branch_positions is None:
        branch_positions = list(range(s))

    n = s + 2 * len(branch_positions)
    adj = [[] for _ in range(n)]

    for i in range(s - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)

    next_v = s
    for pos in branch_positions:
        # intermediate vertex
        adj[pos].append(next_v)
        adj[next_v].append(pos)
        # leaf
        adj[next_v].append(next_v + 1)
        adj[next_v + 1].append(next_v)
        next_v += 2

    return n, adj


def make_binary_tree(depth):
    """Complete binary tree of given depth. d_leaf = 0 or 2 for internal.
    BUT: this has d_leaf = 2 for internal nodes! Not d_leaf ≤ 1.
    Skip this.
    """
    pass


def make_branchy_path(segment_length, num_segments):
    """Path of segments joined by degree-3 vertices, each with 1 whisker.

    Structure: [path_seg]-[branch_v(+whisker)]-[path_seg]-...

    All d_leaf ≤ 1: branch vertices have 1 whisker (d_leaf = 1),
    path endpoints might become leaves.
    """
    # Simpler: just a path of length L with branches (whiskers) every k steps
    L = segment_length * num_segments
    positions = list(range(0, L, segment_length))
    return make_caterpillar(L, positions)


def mode_of(poly):
    return max(range(len(poly)), key=lambda i: poly[i])


def check_dleaf(n, adj):
    """Verify all d_leaf ≤ 1."""
    leaves = set()
    for v in range(n):
        if len(adj[v]) == 1:
            leaves.add(v)
    for v in range(n):
        if v in leaves:
            continue
        d_leaf_v = sum(1 for u in adj[v] if u in leaves)
        if d_leaf_v >= 2:
            return False
    return True


def main():
    print("LARGE d_leaf ≤ 1 TREES: HIGH-MODE CHECK", flush=True)
    print("=" * 70, flush=True)
    print(flush=True)

    # 1. Full caterpillars (whisker at every spine vertex)
    print("FULL CATERPILLARS (whisker at every spine vertex):", flush=True)
    for s in [10, 20, 50, 100, 200, 500]:
        n, adj = make_caterpillar(s)
        assert check_dleaf(n, adj), f"d_leaf check failed for caterpillar s={s}"
        ell = sum(1 for v in range(n) if len(adj[v]) == 1)
        poly = independence_poly(n, adj)
        m = mode_of(poly)
        thr = n // 3 + 1
        high = "HIGH!" if m > thr else ""
        print(f"  s={s:4d}: n={n:5d}, ℓ={ell:4d}, mode={m:4d}, "
              f"thr={thr:4d}, m/thr={m/thr:.3f} {high}", flush=True)

    # 2. Half caterpillars (whisker at every other spine vertex)
    # Skip positions 0, 1, s-2, s-1 to avoid spine-endpoint d_leaf conflicts
    print(flush=True)
    print("HALF CATERPILLARS (whisker at every other interior vertex):", flush=True)
    for s in [10, 20, 50, 100, 200, 500]:
        positions = [i for i in range(2, s - 2, 2)]
        n, adj = make_caterpillar(s, positions)
        assert check_dleaf(n, adj), f"d_leaf fail s={s}, pos={positions[:5]}"
        ell = sum(1 for v in range(n) if len(adj[v]) == 1)
        poly = independence_poly(n, adj)
        m = mode_of(poly)
        thr = n // 3 + 1
        high = "HIGH!" if m > thr else ""
        print(f"  s={s:4d}: n={n:5d}, ℓ={ell:4d}, mode={m:4d}, "
              f"thr={thr:4d}, m/thr={m/thr:.3f} {high}", flush=True)

    # 3. Sparse caterpillars (whisker every 3rd vertex, skip endpoints)
    print(flush=True)
    print("SPARSE CATERPILLARS (whisker every 3rd interior vertex):", flush=True)
    for s in [9, 21, 51, 99, 201, 501]:
        positions = [i for i in range(2, s - 2, 3)]
        n, adj = make_caterpillar(s, positions)
        assert check_dleaf(n, adj), f"d_leaf fail s={s}"
        ell = sum(1 for v in range(n) if len(adj[v]) == 1)
        poly = independence_poly(n, adj)
        m = mode_of(poly)
        thr = n // 3 + 1
        high = "HIGH!" if m > thr else ""
        print(f"  s={s:4d}: n={n:5d}, ℓ={ell:4d}, mode={m:4d}, "
              f"thr={thr:4d}, m/thr={m/thr:.3f} {high}", flush=True)

    # 4. Lobsters (branches of length 2)
    print(flush=True)
    print("FULL LOBSTERS (branch at every spine vertex):", flush=True)
    for s in [10, 20, 50, 100, 200]:
        n, adj = make_lobster(s)
        assert check_dleaf(n, adj), f"d_leaf check failed for lobster s={s}"
        ell = sum(1 for v in range(n) if len(adj[v]) == 1)
        poly = independence_poly(n, adj)
        m = mode_of(poly)
        thr = n // 3 + 1
        high = "HIGH!" if m > thr else ""
        print(f"  s={s:4d}: n={n:5d}, ℓ={ell:4d}, mode={m:4d}, "
              f"thr={thr:4d}, m/thr={m/thr:.3f} {high}", flush=True)

    # 5. Asymmetric: caterpillar with whiskers only on interior of first half
    print(flush=True)
    print("ASYMMETRIC: whiskers on interior of first half:", flush=True)
    for s in [20, 50, 100, 200, 500]:
        positions = list(range(2, s // 2))
        n, adj = make_caterpillar(s, positions)
        assert check_dleaf(n, adj), f"d_leaf fail s={s}"
        ell = sum(1 for v in range(n) if len(adj[v]) == 1)
        poly = independence_poly(n, adj)
        m = mode_of(poly)
        thr = n // 3 + 1
        high = "HIGH!" if m > thr else ""
        print(f"  s={s:4d}: n={n:5d}, ℓ={ell:4d}, mode={m:4d}, "
              f"thr={thr:4d}, m/thr={m/thr:.3f} {high}", flush=True)

    # 6. "Combs": path with a single long branch
    # This creates a tree where one vertex has degree 3 but d_leaf = 1
    print(flush=True)
    print("COMBS: path of length L with one whisker at midpoint:", flush=True)
    for L in [10, 20, 50, 100, 200, 500]:
        positions = [L // 2]
        n, adj = make_caterpillar(L, positions)
        assert check_dleaf(n, adj)
        ell = sum(1 for v in range(n) if len(adj[v]) == 1)
        poly = independence_poly(n, adj)
        m = mode_of(poly)
        thr = n // 3 + 1
        high = "HIGH!" if m > thr else ""
        print(f"  L={L:4d}: n={n:5d}, ℓ={ell:4d}, mode={m:4d}, "
              f"thr={thr:4d}, m/thr={m/thr:.3f} {high}", flush=True)

    print(flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
