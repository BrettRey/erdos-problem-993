#!/usr/bin/env python3
"""
Explore injection/matching structure in independence complexes of trees.

Three untouched approaches from Stanley's taxonomy:
1. Direct injection (containment matching between levels)
2. Chain decomposition structure
3. Maximal IS obstruction analysis

Key questions:
- Does the containment bipartite graph between levels k and k+1 always
  have a matching saturating level k (for k < mode)?
- Are there maximal independent sets of size < mode? (These block the
  simple "add a vertex" injection.)
- Do canonical vertex-selection rules yield injective maps?
"""

import subprocess
import sys
import json
import time
from collections import defaultdict

sys.path.insert(0, ".")
from indpoly import independence_poly


def parse_graph6(s):
    """Parse graph6 format to adjacency list."""
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
    """Enumerate all independent sets organized by size. Returns dict: size -> list of frozensets."""
    nbr = [set(adj[v]) for v in range(n)]
    levels = defaultdict(list)

    def backtrack(v, current, forbidden):
        if v == n:
            levels[len(current)].append(frozenset(current))
            return
        # Skip v
        backtrack(v + 1, current, forbidden)
        # Include v if not forbidden
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return levels


def hopcroft_karp(left_adj, n_left, n_right):
    """Hopcroft-Karp maximum bipartite matching.
    left_adj[i] = list of right neighbors of left node i.
    Returns matching size.
    """
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


def check_level_matching(levels, k):
    """Check containment matching between levels k and k+1.
    Returns (matching_saturates_left, matching_size, left_size, right_size,
             min_left_degree, max_left_degree).
    """
    left = levels.get(k, [])
    right = levels.get(k + 1, [])
    if not left:
        return True, 0, 0, len(right), 0, 0
    if not right:
        return False, 0, len(left), 0, 0, 0

    left_idx = {s: i for i, s in enumerate(left)}
    n_left = len(left)
    n_right = len(right)

    # Build adjacency: for each right set T, T \ {v} for v in T gives left neighbors
    left_adj = [[] for _ in range(n_left)]
    for j, t in enumerate(right):
        for v in t:
            s = t - {v}
            if s in left_idx:
                left_adj[left_idx[s]].append(j)

    degs = [len(left_adj[i]) for i in range(n_left)]
    min_d = min(degs) if degs else 0
    max_d = max(degs) if degs else 0

    msize, _, _ = hopcroft_karp(left_adj, n_left, n_right)
    return msize == n_left, msize, n_left, n_right, min_d, max_d


def try_canonical_injections(levels, k, n, adj):
    """Try canonical vertex-selection rules for injection level k -> k+1."""
    left = levels.get(k, [])
    if not left:
        return {}

    nbr = [set(adj[v]) for v in range(n)]
    deg = [len(adj[v]) for v in range(n)]

    # BFS from vertex 0 for depth-based rules
    depth = [0] * n
    visited = [False] * n
    queue = [0]
    visited[0] = True
    qi = 0
    while qi < len(queue):
        v = queue[qi]
        qi += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                depth[u] = depth[v] + 1
                queue.append(u)

    def available(s):
        return [v for v in range(n) if v not in s and not (nbr[v] & s)]

    rules = {
        "smallest": lambda avail, s: min(avail),
        "largest": lambda avail, s: max(avail),
        "near_root": lambda avail, s: min(avail, key=lambda v: (depth[v], v)),
        "far_root": lambda avail, s: max(avail, key=lambda v: (depth[v], -v)),
        "min_degree": lambda avail, s: min(avail, key=lambda v: (deg[v], v)),
        "max_degree": lambda avail, s: max(avail, key=lambda v: (deg[v], -v)),
    }

    results = {}
    for rule_name, pick_fn in rules.items():
        images = set()
        collisions = 0
        unmappable = 0
        for s in left:
            avail = available(s)
            if not avail:
                unmappable += 1
                continue
            v = pick_fn(avail, s)
            image = s | frozenset({v})
            if image in images:
                collisions += 1
            images.add(image)
        is_inj = collisions == 0 and unmappable == 0
        results[rule_name] = {
            "injective": is_inj,
            "collisions": collisions,
            "unmappable": unmappable,
            "distinct_images": len(images),
            "total": len(left),
        }
    return results


def analyze_maximal_is(n, adj, levels, mode):
    """Count maximal independent sets at each level below the mode."""
    nbr = [set(adj[v]) for v in range(n)]
    counts = {}  # level -> (num_maximal, num_total)
    for size in range(mode):
        sets = levels.get(size, [])
        if not sets:
            continue
        n_maximal = 0
        for s in sets:
            can_extend = any(v not in s and not (nbr[v] & s) for v in range(n))
            if not can_extend:
                n_maximal += 1
        if n_maximal > 0:
            counts[size] = (n_maximal, len(sets))
    return counts


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def main():
    max_n_matching = 14       # Check matching through this n
    max_n_injection = 12      # Check canonical injections through this n
    max_n_maximal = 14        # Check maximal IS through this n

    all_results = {}

    for n in range(3, max_n_matching + 1):
        t0 = time.time()
        print(f"\n{'='*60}")
        print(f"n = {n}")
        print(f"{'='*60}")

        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]
        n_trees = len(lines)

        # Counters
        matching_perfect_count = 0
        matching_fail_trees = []
        maximal_below_mode_count = 0
        maximal_examples = []
        total_level_checks = 0

        # Injection rule stats (only for n <= max_n_injection)
        rule_totals = defaultdict(lambda: {"levels": 0, "injective": 0})

        # Degree distribution in containment graph
        min_deg_hist = defaultdict(int)  # min_left_degree -> count of (tree, level) pairs

        for tidx, line in enumerate(lines):
            tn, adj = parse_graph6(line)
            poly = independence_poly(n, adj)
            mode = compute_mode(poly)
            levels = enumerate_independent_sets(n, adj)

            # --- Maximal IS below mode ---
            if n <= max_n_maximal:
                max_is = analyze_maximal_is(n, adj, levels, mode)
                if max_is:
                    maximal_below_mode_count += 1
                    if len(maximal_examples) < 20:
                        maximal_examples.append({
                            "tree_idx": tidx, "g6": line.strip(),
                            "poly": poly, "mode": mode,
                            "maximal_levels": {str(k): v for k, v in max_is.items()},
                        })

            # --- Containment matching ---
            tree_ok = True
            for k in range(mode):
                lk = levels.get(k, [])
                if not lk:
                    continue
                total_level_checks += 1
                ok, msize, lsz, rsz, min_d, max_d = check_level_matching(levels, k)
                min_deg_hist[min_d] += 1
                if not ok:
                    tree_ok = False
                    if len(matching_fail_trees) < 20:
                        matching_fail_trees.append({
                            "tree_idx": tidx, "g6": line.strip(), "k": k,
                            "matching": msize, "left": lsz, "right": rsz,
                            "min_deg": min_d,
                        })

            if tree_ok:
                matching_perfect_count += 1

            # --- Canonical injections ---
            if n <= max_n_injection:
                for k in range(mode):
                    if not levels.get(k, []):
                        continue
                    inj = try_canonical_injections(levels, k, n, adj)
                    for rname, res in inj.items():
                        rule_totals[rname]["levels"] += 1
                        if res["injective"]:
                            rule_totals[rname]["injective"] += 1

        elapsed = time.time() - t0

        # --- Report ---
        print(f"  Trees: {n_trees}  ({elapsed:.1f}s)")
        print(f"  Matching saturates left (all levels): {matching_perfect_count}/{n_trees} "
              f"({100*matching_perfect_count/max(n_trees,1):.1f}%)")
        if matching_fail_trees:
            print(f"  *** MATCHING FAILURES: {len(matching_fail_trees)} ***")
            for mf in matching_fail_trees[:5]:
                print(f"      tree={mf['tree_idx']} k={mf['k']} match={mf['matching']}/{mf['left']} min_deg={mf['min_deg']}")

        print(f"  Level checks: {total_level_checks}")
        print(f"  Min-left-degree histogram: {dict(sorted(min_deg_hist.items()))}")

        if n <= max_n_maximal:
            print(f"  Trees with maximal IS below mode: {maximal_below_mode_count}/{n_trees} "
                  f"({100*maximal_below_mode_count/max(n_trees,1):.1f}%)")
            if maximal_examples:
                for ex in maximal_examples[:3]:
                    print(f"    Example: poly={ex['poly']} mode={ex['mode']} maximal_at={ex['maximal_levels']}")

        if n <= max_n_injection:
            print(f"  Canonical injection results:")
            for rname in sorted(rule_totals.keys()):
                rt = rule_totals[rname]
                pct = 100 * rt["injective"] / max(rt["levels"], 1)
                print(f"    {rname}: {rt['injective']}/{rt['levels']} levels injective ({pct:.1f}%)")

        all_results[n] = {
            "n_trees": n_trees,
            "matching_perfect": matching_perfect_count,
            "matching_failures": len(matching_fail_trees),
            "maximal_below_mode": maximal_below_mode_count,
            "total_level_checks": total_level_checks,
            "min_deg_hist": dict(min_deg_hist),
            "injection_stats": {k: dict(v) for k, v in rule_totals.items()},
            "elapsed": elapsed,
        }

    # Save
    with open("results/injection_exploration.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    # --- Grand summary ---
    print(f"\n\n{'='*60}")
    print("GRAND SUMMARY")
    print(f"{'='*60}")
    print("\nMatching (Hall's condition via Hopcroft-Karp):")
    for nv in sorted(all_results.keys()):
        r = all_results[nv]
        print(f"  n={nv:2d}: {r['matching_perfect']}/{r['n_trees']} trees fully matchable, "
              f"{r['matching_failures']} failures")

    print("\nMaximal IS below mode:")
    for nv in sorted(all_results.keys()):
        r = all_results[nv]
        if "maximal_below_mode" in r:
            print(f"  n={nv:2d}: {r['maximal_below_mode']}/{r['n_trees']} trees "
                  f"({100*r['maximal_below_mode']/max(r['n_trees'],1):.1f}%)")

    print("\nCanonical injection success rates (levels where rule is injective):")
    # Aggregate across all n
    agg = defaultdict(lambda: {"levels": 0, "injective": 0})
    for nv in sorted(all_results.keys()):
        r = all_results[nv]
        for rname, stats in r.get("injection_stats", {}).items():
            agg[rname]["levels"] += stats["levels"]
            agg[rname]["injective"] += stats["injective"]
    for rname in sorted(agg.keys()):
        a = agg[rname]
        pct = 100 * a["injective"] / max(a["levels"], 1)
        print(f"  {rname}: {a['injective']}/{a['levels']} ({pct:.1f}%)")


if __name__ == "__main__":
    main()
