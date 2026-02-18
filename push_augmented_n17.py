#!/usr/bin/env python3
"""
Push augmented matching check to n=17 (and 18 if feasible).
Focused: only check the augmented matching, skip edge counting.
Use early termination on any failure.
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

    return sum(1 for m in match_l if m != -1)


def check_augmented_matching(levels, k, n, adj):
    """Check if augmented bipartite graph has matching saturating level k."""
    left = levels.get(k, [])
    right = levels.get(k + 1, [])
    if not left:
        return True
    if not right:
        return False

    nbr = [set(adj[v]) for v in range(n)]
    left_idx = {s: i for i, s in enumerate(left)}
    right_idx = {s: j for j, s in enumerate(right)}
    n_left = len(left)
    n_right = len(right)

    left_neighbors = [set() for _ in range(n_left)]

    # Containment edges
    for j, t in enumerate(right):
        for v in t:
            s = t - {v}
            if s in left_idx:
                left_neighbors[left_idx[s]].add(j)

    # Swap edges
    for i, s in enumerate(left):
        for u in s:
            s_minus_u = s - {u}
            forbidden_by_remaining = set()
            for x in s_minus_u:
                forbidden_by_remaining |= nbr[x]
            candidates = [v for v in range(n)
                          if v not in s and v not in forbidden_by_remaining]
            for ci, v in enumerate(candidates):
                for w in candidates[ci + 1:]:
                    if w not in nbr[v]:
                        t = s_minus_u | {v, w}
                        if t in right_idx:
                            left_neighbors[i].add(right_idx[t])

    left_adj_lists = [sorted(left_neighbors[i]) for i in range(n_left)]
    msize = hopcroft_karp(left_adj_lists, n_left, n_right)
    return msize == n_left


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def main():
    for n in [17, 18]:
        t0 = time.time()
        print(f"\n{'='*60}", flush=True)
        print(f"n = {n}", flush=True)
        print(f"{'='*60}", flush=True)

        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]
        n_trees = len(lines)
        print(f"Trees: {n_trees}", flush=True)

        failures = []
        checked = 0
        levels_checked = 0

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
                levels_checked += 1

                ok = check_augmented_matching(levels, k, n, adj_data)
                if not ok:
                    tree_ok = False
                    failures.append({
                        "tree": tidx, "g6": line.strip(), "k": k,
                        "poly": poly, "mode": mode,
                    })
                    print(f"\n*** FAILURE at tree={tidx} k={k} ***", flush=True)
                    print(f"    poly={poly}", flush=True)
                    break  # skip remaining levels for this tree

            checked += 1
            if checked % 2000 == 0:
                elapsed = time.time() - t0
                rate = checked / elapsed
                eta = (n_trees - checked) / rate if rate > 0 else 0
                print(f"  {checked}/{n_trees} trees ({elapsed:.0f}s, "
                      f"~{eta:.0f}s remaining, {levels_checked} levels, "
                      f"{len(failures)} failures)", flush=True)

        elapsed = time.time() - t0
        print(f"\nn={n} DONE: {checked}/{n_trees} trees in {elapsed:.1f}s", flush=True)
        print(f"  Levels checked: {levels_checked}", flush=True)
        print(f"  Failures: {len(failures)}", flush=True)

        if failures:
            print("  FAILURE DETAILS:")
            for f in failures[:10]:
                print(f"    tree={f['tree']} k={f['k']} poly={f['poly']}")
            break  # Stop on any failure
        else:
            print(f"  ** ALL TREES PASS **", flush=True)


if __name__ == "__main__":
    main()
