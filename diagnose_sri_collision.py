#!/usr/bin/env python3
"""Diagnose SRI collisions: check if containment and swap ever target the same IS.

The original explore_sri.py test_containment_first checked each (s, t) pair
individually using the known forward rule, but didn't check whether two different
left elements map to the same right element (collision). This script checks.
"""

import subprocess
import sys
sys.path.insert(0, ".")
from indpoly import independence_poly
from collections import defaultdict


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


def is_maximal_is(s, n, nbr):
    for v in range(n):
        if v not in s and not (nbr[v] & s):
            return False
    return True


def compute_mode(poly):
    return max(range(len(poly)), key=lambda i: poly[i])


def main():
    print("COLLISION DIAGNOSIS: does containment-first SRI have injective forward map?")
    print("=" * 75)

    for n in range(5, 14):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_trees = len(lines)
        total_levels = 0
        collision_levels = 0
        example_shown = False

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            deg = [len(adj_data[v]) for v in range(tn)]
            poly = independence_poly(tn, adj_data)
            mode = compute_mode(poly)
            levels = enumerate_independent_sets(tn, adj_data)
            order = list(range(tn))

            for k in range(mode):
                left = levels.get(k, [])
                right = levels.get(k + 1, [])
                if not left:
                    continue

                total_levels += 1
                left_set = set(left)
                right_set = set(right)

                # Build forward map (same as original containment-first)
                forward_map = {}
                for s in left:
                    if is_maximal_is(s, tn, nbr):
                        # Swap: find canonical leaf-swap triple
                        candidates = []
                        for u in sorted(s):
                            leaf_nbrs = sorted(
                                v for v in nbr[u] if v not in s and deg[v] == 1)
                            for i, v in enumerate(leaf_nbrs):
                                for w in leaf_nbrs[i + 1:]:
                                    t = (s - {u}) | {v, w}
                                    if t in right_set:
                                        candidates.append((u, v, w, t))
                        if candidates:
                            candidates.sort(
                                key=lambda x: (-deg[x[0]], x[0], x[1], x[2]))
                            u, v, w, t = candidates[0]
                            forward_map[s] = (t, "swap")
                    else:
                        # Containment: first free vertex
                        for v in order:
                            if v not in s and not (nbr[v] & s):
                                t = s | {v}
                                forward_map[s] = (t, "contain")
                                break

                # Check for collisions
                target_sources = defaultdict(list)
                for s, (t, rule) in forward_map.items():
                    target_sources[t].append((s, rule))

                has_collision = any(
                    len(sources) > 1 for sources in target_sources.values())

                if has_collision:
                    collision_levels += 1
                    if not example_shown and n <= 8:
                        example_shown = True
                        for t, sources in target_sources.items():
                            if len(sources) > 1:
                                print(f"\n  COLLISION at n={tn} tree={tidx} "
                                      f"k={k} mode={mode}")
                                print(f"  Target T = {sorted(t)}")
                                for s, rule in sources:
                                    mx = "maximal" if is_maximal_is(s, tn, nbr) \
                                        else "non-max"
                                    print(f"    Source S={sorted(s)} "
                                          f"rule={rule} ({mx})")

        pct = 100 * collision_levels / total_levels if total_levels else 0
        print(f"n={n:2d}: {n_trees:5d} trees, {total_levels:5d} levels, "
              f"{collision_levels:4d} collisions ({pct:.1f}%)")


if __name__ == "__main__":
    main()
