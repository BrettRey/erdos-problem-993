#!/usr/bin/env python3
"""
Explore augmented injection: containment + swap edges.

Key insight from explore_injection.py: the simple containment bipartite graph
(add one vertex) fails for ~20% of trees due to maximal independent sets
below the mode.

Here we augment the bipartite graph with "swap" edges:
  S ~ T if T = (S \ {u}) ∪ {v, w} for some u ∈ S, v,w ∉ S
  (remove one vertex, add two — requires T to be independent)

Question: Does the augmented graph ALWAYS have a matching saturating
level k for k < mode(I(T))?

If yes, this opens a path to an injection proof of unimodality.
"""

import subprocess
import sys
import time
from collections import defaultdict

sys.path.insert(0, ".")
from indpoly import independence_poly


def parse_graph6(s):
    s = s.strip()
    if not s:
        return 0, []
    idx = 0
    if ord(s[0]) - 63 < 63:
        n = ord(s[0]) - 63
        idx = 1
    else:
        idx = 1
        n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63)
            idx += 1
    adj = [[] for _ in range(n)]
    bits = []
    for c in s[idx:]:
        val = ord(c) - 63
        for b in range(5, -1, -1):
            bits.append((val >> b) & 1)
    bit_idx = 0
    for j in range(n):
        for i in range(j):
            if bit_idx < len(bits) and bits[bit_idx]:
                adj[i].append(j)
                adj[j].append(i)
            bit_idx += 1
    return n, adj


def enumerate_independent_sets(n, adj):
    nbr = [set(adj[v]) for v in range(n)]
    levels = defaultdict(list)

    def backtrack(v, current, forbidden):
        if v == n:
            levels[len(current)].append(frozenset(current))
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return levels


def hopcroft_karp(left_adj, n_left, n_right):
    match_l = [-1] * n_left
    match_r = [-1] * n_right

    def bfs():
        dist = [0] * n_left
        queue = []
        for i in range(n_left):
            if match_l[i] == -1:
                dist[i] = 0
                queue.append(i)
            else:
                dist[i] = float("inf")
        found = False
        qi = 0
        while qi < len(queue):
            i = queue[qi]
            qi += 1
            for j in left_adj[i]:
                ni = match_r[j]
                if ni == -1:
                    found = True
                elif dist[ni] == float("inf"):
                    dist[ni] = dist[i] + 1
                    queue.append(ni)
        return found, dist

    def dfs(i, dist):
        for j in left_adj[i]:
            ni = match_r[j]
            if ni == -1 or (dist[ni] == dist[i] + 1 and dfs(ni, dist)):
                match_l[i] = j
                match_r[j] = i
                return True
        dist[i] = float("inf")
        return False

    while True:
        found, dist = bfs()
        if not found:
            break
        for i in range(n_left):
            if match_l[i] == -1:
                dfs(i, dist)

    return sum(1 for m in match_l if m != -1), match_l, match_r


def check_augmented_matching(levels, k, n, adj):
    """Check matching in augmented bipartite graph (containment + swap edges).

    Returns (saturates, matching_size, left_size, right_size,
             n_containment_edges, n_swap_edges,
             n_maximal_left, min_left_deg, max_left_deg).
    """
    left = levels.get(k, [])
    right = levels.get(k + 1, [])
    if not left:
        return True, 0, 0, len(right), 0, 0, 0, 0, 0
    if not right:
        return False, 0, len(left), 0, 0, 0, 0, 0, 0

    nbr = [set(adj[v]) for v in range(n)]
    left_idx = {s: i for i, s in enumerate(left)}
    right_idx = {s: j for j, s in enumerate(right)}
    n_left = len(left)
    n_right = len(right)

    # Build adjacency with both containment and swap edges
    left_neighbors = [set() for _ in range(n_left)]  # use sets to deduplicate
    n_contain = 0
    n_swap = 0

    # Containment edges: T = S ∪ {v} for some v
    for j, t in enumerate(right):
        for v in t:
            s = t - {v}
            if s in left_idx:
                i = left_idx[s]
                if j not in left_neighbors[i]:
                    left_neighbors[i].add(j)
                    n_contain += 1

    # Swap edges: T = (S \ {u}) ∪ {v, w} for some u ∈ S, v,w ∉ S, T independent
    for i, s in enumerate(left):
        for u in s:
            # Remove u from S
            s_minus_u = s - {u}
            # Find vertices that can be added: not in s, not adjacent to any vertex in s_minus_u
            forbidden_by_remaining = set()
            for x in s_minus_u:
                forbidden_by_remaining |= nbr[x]
            candidates = [v for v in range(n)
                          if v not in s and v not in forbidden_by_remaining]
            # Try all pairs of candidates that form an independent pair
            for ci, v in enumerate(candidates):
                for w in candidates[ci + 1:]:
                    if w not in nbr[v]:
                        t = s_minus_u | {v, w}
                        if t in right_idx:
                            j = right_idx[t]
                            if j not in left_neighbors[i]:
                                left_neighbors[i].add(j)
                                n_swap += 1

    # Count maximal ISes (degree 0 in containment-only graph)
    # and compute degrees in augmented graph
    n_maximal = 0
    degs = []
    for i in range(n_left):
        d = len(left_neighbors[i])
        degs.append(d)
        # Check if this IS is maximal (no containment edges)
        s = left[i]
        can_extend = any(v not in s and not (nbr[v] & s) for v in range(n))
        if not can_extend:
            n_maximal += 1

    left_adj_lists = [sorted(left_neighbors[i]) for i in range(n_left)]
    min_d = min(degs) if degs else 0
    max_d = max(degs) if degs else 0

    msize, _, _ = hopcroft_karp(left_adj_lists, n_left, n_right)
    return (msize == n_left, msize, n_left, n_right,
            n_contain, n_swap, n_maximal, min_d, max_d)


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def main():
    max_n = 16  # Swap edge enumeration is expensive; keep reasonable

    for n in range(3, max_n + 1):
        t0 = time.time()
        print(f"\n{'='*60}", flush=True)
        print(f"n = {n}", flush=True)
        print(f"{'='*60}", flush=True)

        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]
        n_trees = len(lines)

        aug_perfect = 0
        aug_failures = []
        total_levels = 0
        total_contain_edges = 0
        total_swap_edges = 0
        total_maximal = 0
        zero_deg_in_augmented = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)
            levels = enumerate_independent_sets(n, adj_data)

            tree_ok = True
            for k in range(mode):
                lk = levels.get(k, [])
                if not lk:
                    continue
                total_levels += 1

                (ok, msize, lsz, rsz,
                 nc, ns, nmax, min_d, max_d) = check_augmented_matching(
                    levels, k, n, adj_data)

                total_contain_edges += nc
                total_swap_edges += ns
                total_maximal += nmax
                if min_d == 0:
                    zero_deg_in_augmented += 1

                if not ok:
                    tree_ok = False
                    if len(aug_failures) < 20:
                        aug_failures.append({
                            "n": n, "tree": tidx, "g6": line.strip(),
                            "k": k, "match": msize, "left": lsz, "right": rsz,
                            "contain_edges": nc, "swap_edges": ns,
                            "maximal": nmax, "min_deg": min_d,
                            "poly": poly,
                        })

            if tree_ok:
                aug_perfect += 1

        elapsed = time.time() - t0
        print(f"  Trees: {n_trees}  ({elapsed:.1f}s)", flush=True)
        print(f"  Augmented matching saturates all levels: {aug_perfect}/{n_trees} "
              f"({100*aug_perfect/max(n_trees,1):.1f}%)")
        print(f"  Total levels checked: {total_levels}")
        print(f"  Total containment edges: {total_contain_edges}")
        print(f"  Total swap edges: {total_swap_edges}")
        print(f"  Total maximal IS below mode: {total_maximal}")
        print(f"  Levels with min-degree 0 in augmented graph: {zero_deg_in_augmented}")

        if aug_failures:
            print(f"\n  *** AUGMENTED MATCHING FAILURES ({len(aug_failures)}): ***")
            for f in aug_failures[:10]:
                print(f"    tree={f['tree']} k={f['k']} match={f['match']}/{f['left']} "
                      f"contain={f['contain_edges']} swap={f['swap_edges']} "
                      f"maximal={f['maximal']} min_deg={f['min_deg']}")
                print(f"    poly={f['poly']}")
        else:
            print(f"  ** NO augmented matching failures! **")

    print(f"\n\n{'='*60}")
    print("CONCLUSION")
    print(f"{'='*60}")
    print("If augmented matching (containment + swap) always saturates,")
    print("then a two-phase injection proof is feasible:")
    print("  Phase 1: Extend non-maximal IS via containment (add vertex)")
    print("  Phase 2: Transform maximal IS via swap (remove 1, add 2)")


if __name__ == "__main__":
    main()
