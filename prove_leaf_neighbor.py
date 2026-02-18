#!/usr/bin/env python3
"""
Comprehensive analysis of proof routes for the Leaf-Neighbor Property.

For any tree T and any maximal IS S of size k < mode(I(T)), some u in S
has >= 2 leaf-neighbors (degree-1 vertices adjacent to u, outside S).

We test four proof routes:
  A: Pigeonhole via |L_out| > k
  B: Degree bound: max deg in S >= k+1
  C: Support vertex: S contains a support vertex with >= 2 leaf-children
  D: Combined / structural characterization of residual cases

Goal: find a route (or combination) that covers 100% of cases, enabling
a clean proof.
"""

import subprocess
import sys
import time
from collections import Counter, defaultdict

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


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def find_maximal_is_below_mode(n, adj, mode):
    """Find all maximal independent sets of size < mode."""
    nbr = [set(adj[v]) for v in range(n)]
    results = []

    def backtrack(v, current, forbidden):
        if v == n:
            k = len(current)
            if k >= mode:
                return
            s = frozenset(current)
            can_extend = any(u not in s and not (nbr[u] & s) for u in range(n))
            if not can_extend:
                results.append((s, k))
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])

    backtrack(0, [], set())
    return results


def analyze_case(n, adj, s, k, mode):
    """Analyze a single (tree, maximal IS) pair for all proof routes."""
    deg = [len(adj[v]) for v in range(n)]
    nbr = [set(adj[v]) for v in range(n)]
    leaves = {v for v in range(n) if deg[v] == 1}
    n_leaves = len(leaves)

    L_in = s & leaves
    L_out = leaves - s
    n_lout = len(L_out)
    n_lin = len(L_in)

    # Per-vertex leaf-neighbor counts (leaf-neighbors outside S)
    leaf_nbrs = {}
    for u in s:
        leaf_nbrs[u] = sum(1 for v in nbr[u] if v in L_out)

    max_leaf_nbrs = max(leaf_nbrs.values()) if leaf_nbrs else 0
    n_with_2plus = sum(1 for u in s if leaf_nbrs[u] >= 2)

    # Route A: pigeonhole via |L_out| > k
    route_a = n_lout > k  # strictly, need |L_out| >= k+1

    # Route B: degree bound max deg(u) >= k+1 for u in S
    max_deg_in_s = max(deg[u] for u in s)
    route_b = max_deg_in_s >= k + 1

    # Route C: support vertex with >= 2 leaf-children in S
    # A support vertex has >= 1 leaf-child. If u in S is support with >= 2
    # leaf-children, those children are in L_out (S is independent), giving
    # >= 2 leaf-neighbors.
    # Actually, leaf-children of u in the tree = leaf-neighbors of u in the
    # tree. If u in S, all its neighbors are outside S, so leaf-children
    # that are leaves are in L_out.
    route_c = False
    for u in s:
        # Count leaf-children: neighbors of u that are leaves
        # (same as leaf_nbrs[u] since u in S => neighbors outside S)
        if leaf_nbrs[u] >= 2:
            route_c = True
            break

    # Note: Route C is equivalent to "some u in S has >= 2 leaf-neighbors"
    # which is exactly the LNP. So Route C = LNP directly.
    # What we really want for Route C is a structural sufficient condition:
    # S contains a support vertex (vertex adjacent to a leaf) with >= 2
    # leaf-children in the TREE (regardless of S).
    # This is slightly different: a vertex might have 3 leaf-children in
    # the tree, but 1 could be in S (as a leaf of the tree that's in S).
    # Wait: S is independent, so if u in S and v is a child of u, v not in S.
    # So for u in S, leaf-children of u in the tree = leaf-neighbors of u
    # outside S. They're the same.

    # Route C': Does the tree have a support vertex with >= 2 leaf-children,
    # AND is that vertex in S?
    support_verts_with_2plus = set()
    for v in range(n):
        leaf_children = sum(1 for w in adj[v] if deg[w] == 1)
        if leaf_children >= 2:
            support_verts_with_2plus.add(v)
    route_c_prime = bool(s & support_verts_with_2plus)

    # Route D: the "concentrated" route — even without a support vertex
    # with 2 leaf-children, can we still guarantee 2 leaf-neighbors from
    # non-leaf support vertices? This happens when u in S has leaf-neighbors
    # from different parts of the tree.

    # Private neighbor analysis
    priv = {}
    shared = {}
    for u in s:
        p = 0
        sh = 0
        for v in nbr[u]:
            if v not in s:
                s_nbrs_of_v = nbr[v] & s
                if len(s_nbrs_of_v) == 1:
                    p += 1
                else:
                    sh += 1
        priv[u] = p
        shared[u] = sh

    max_priv = max(priv.values()) if priv else 0

    return {
        "n_leaves": n_leaves,
        "n_lout": n_lout,
        "n_lin": n_lin,
        "max_leaf_nbrs": max_leaf_nbrs,
        "n_with_2plus": n_with_2plus,
        "route_a": route_a,
        "route_b": route_b,
        "route_c": route_c,  # same as LNP
        "route_c_prime": route_c_prime,
        "max_deg_in_s": max_deg_in_s,
        "max_priv": max_priv,
        "leaf_nbrs": leaf_nbrs,
        "priv": priv,
        "shared": shared,
        "support_2plus": support_verts_with_2plus,
    }


def classify_tree(n, adj):
    """Classify tree structure: path, star, caterpillar, spider, etc."""
    deg = [len(adj[v]) for v in range(n)]
    max_deg = max(deg)
    leaves = sum(1 for v in range(n) if deg[v] == 1)
    internal = [v for v in range(n) if deg[v] >= 2]

    if max_deg <= 2:
        return "path"
    if max_deg == n - 1:
        return "star"

    # Double star: exactly 2 internal vertices, both with high degree
    if len(internal) == 2:
        return f"double-star"

    # Caterpillar: removing all leaves gives a path
    if internal:
        internal_set = set(internal)
        internal_degs = [sum(1 for w in adj[v] if w in internal_set) for v in internal]
        if all(d <= 2 for d in internal_degs):
            return "caterpillar"

    # Spider: one high-degree vertex, rest are paths
    high_deg = [v for v in range(n) if deg[v] >= 3]
    if len(high_deg) == 1:
        return "spider"

    return "other"


def main():
    max_n = 16

    print("=" * 75)
    print("LEAF-NEIGHBOR PROPERTY: Comprehensive Proof Route Analysis")
    print("=" * 75)

    total = 0
    route_a_ok = 0
    route_b_ok = 0
    route_c_ok = 0  # LNP itself
    route_c_prime_ok = 0  # support vertex in S with >= 2 leaf-children
    both_ab_ok = 0
    neither_ab = 0  # not covered by A or B

    # Pigeonhole failure analysis
    pigeonhole_fails = []
    # Degree bound failure analysis
    degree_fails = []
    # Neither A nor B
    neither_cases = []
    # Tight LNP cases (max_leaf_nbrs == 2)
    tight_lnp = []
    # Route C' failure analysis
    c_prime_fails = []

    # Per-n stats
    per_n = defaultdict(lambda: {"total": 0, "a": 0, "b": 0, "c_prime": 0,
                                  "neither_ab": 0})

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)

            max_is_list = find_maximal_is_below_mode(n, adj_data, mode)
            tree_type = classify_tree(n, adj_data) if max_is_list else None

            for s, k in max_is_list:
                total += 1
                info = analyze_case(n, adj_data, s, k, mode)

                per_n[n]["total"] += 1

                if info["route_a"]:
                    route_a_ok += 1
                    per_n[n]["a"] += 1
                if info["route_b"]:
                    route_b_ok += 1
                    per_n[n]["b"] += 1
                if info["route_c"]:
                    route_c_ok += 1
                if info["route_c_prime"]:
                    route_c_prime_ok += 1
                    per_n[n]["c_prime"] += 1
                if info["route_a"] and info["route_b"]:
                    both_ab_ok += 1

                if not info["route_a"] and not info["route_b"]:
                    neither_ab += 1
                    per_n[n]["neither_ab"] += 1
                    neither_cases.append({
                        "n": n, "tree": tidx, "k": k, "mode": mode,
                        "g6": line.strip(), "tree_type": tree_type,
                        "info": info, "s": sorted(s),
                    })

                if not info["route_a"]:
                    pigeonhole_fails.append({
                        "n": n, "k": k, "mode": mode,
                        "n_leaves": info["n_leaves"],
                        "n_lout": info["n_lout"],
                        "tree_type": tree_type,
                        "route_b": info["route_b"],
                        "route_c_prime": info["route_c_prime"],
                        "max_leaf_nbrs": info["max_leaf_nbrs"],
                    })

                if not info["route_b"]:
                    degree_fails.append({
                        "n": n, "k": k, "mode": mode,
                        "max_deg_in_s": info["max_deg_in_s"],
                        "tree_type": tree_type,
                        "route_a": info["route_a"],
                        "max_leaf_nbrs": info["max_leaf_nbrs"],
                    })

                if not info["route_c_prime"]:
                    c_prime_fails.append({
                        "n": n, "k": k, "mode": mode,
                        "tree_type": tree_type,
                        "max_leaf_nbrs": info["max_leaf_nbrs"],
                        "n_leaves": info["n_leaves"],
                        "support_2plus": info["support_2plus"],
                        "s": sorted(s),
                    })

                if info["max_leaf_nbrs"] == 2:
                    tight_lnp.append({
                        "n": n, "tree": tidx, "k": k, "mode": mode,
                        "g6": line.strip(), "tree_type": tree_type,
                        "n_leaves": info["n_leaves"],
                        "n_lout": info["n_lout"],
                        "leaf_nbrs": {u: info["leaf_nbrs"][u] for u in sorted(s)},
                        "priv": {u: info["priv"][u] for u in sorted(s)},
                        "shared": {u: info["shared"][u] for u in sorted(s)},
                        "s": sorted(s),
                        "route_a": info["route_a"],
                        "route_b": info["route_b"],
                        "route_c_prime": info["route_c_prime"],
                    })

        elapsed = time.time() - t0
        pn = per_n[n]
        print(f"n={n:2d} ({len(lines):5d} trees, {elapsed:5.1f}s): "
              f"cases={pn['total']:5d}  "
              f"A={pn['a']:5d}  B={pn['b']:5d}  "
              f"C'={pn['c_prime']:5d}  "
              f"neither_AB={pn['neither_ab']:3d}",
              flush=True)

    # ===== SUMMARY =====
    print(f"\n\n{'='*75}")
    print("ROUTE COVERAGE SUMMARY")
    print(f"{'='*75}")
    print(f"Total maximal IS below mode: {total}")
    print(f"  Route A (|L_out| > k):         {route_a_ok:6d}/{total} "
          f"({100*route_a_ok/total:.2f}%)")
    print(f"  Route B (max deg >= k+1):      {route_b_ok:6d}/{total} "
          f"({100*route_b_ok/total:.2f}%)")
    print(f"  Route C' (support w/ 2+ in S): {route_c_prime_ok:6d}/{total} "
          f"({100*route_c_prime_ok/total:.2f}%)")
    print(f"  Route A or B:                  {route_a_ok + route_b_ok - both_ab_ok:6d}/{total} "
          f"({100*(route_a_ok + route_b_ok - both_ab_ok)/total:.2f}%)")
    print(f"  LNP (the property itself):     {route_c_ok:6d}/{total} "
          f"({100*route_c_ok/total:.2f}%)")

    # ===== PIGEONHOLE FAILURES =====
    print(f"\n{'='*75}")
    print(f"ROUTE A FAILURES (pigeonhole, {len(pigeonhole_fails)} cases)")
    print(f"{'='*75}")
    if pigeonhole_fails:
        # How many are covered by Route B?
        b_covers = sum(1 for pf in pigeonhole_fails if pf["route_b"])
        c_covers = sum(1 for pf in pigeonhole_fails if pf["route_c_prime"])
        print(f"  Covered by Route B: {b_covers}/{len(pigeonhole_fails)}")
        print(f"  Covered by Route C': {c_covers}/{len(pigeonhole_fails)}")
        print(f"\n  Tree types among A-failures:")
        type_counts = Counter(pf["tree_type"] for pf in pigeonhole_fails)
        for tt, cnt in type_counts.most_common():
            print(f"    {tt}: {cnt}")
        print(f"\n  |L_out| - k distribution:")
        diff_counts = Counter(pf["n_lout"] - pf["k"] for pf in pigeonhole_fails)
        for d in sorted(diff_counts.keys()):
            print(f"    |L_out|-k = {d}: {diff_counts[d]}")
        print(f"\n  Sample A-failures (first 10):")
        for pf in pigeonhole_fails[:10]:
            print(f"    n={pf['n']} k={pf['k']} mode={pf['mode']} "
                  f"leaves={pf['n_leaves']} L_out={pf['n_lout']} "
                  f"type={pf['tree_type']} max_leaf_nbrs={pf['max_leaf_nbrs']} "
                  f"B={'Y' if pf['route_b'] else 'N'} "
                  f"C'={'Y' if pf['route_c_prime'] else 'N'}")

    # ===== NEITHER A NOR B =====
    print(f"\n{'='*75}")
    print(f"CASES COVERED BY NEITHER A NOR B ({neither_ab} cases)")
    print(f"{'='*75}")
    if neither_cases:
        # All of these should still have the LNP (max_leaf_nbrs >= 2)
        lnp_holds = sum(1 for nc in neither_cases
                        if nc["info"]["max_leaf_nbrs"] >= 2)
        print(f"  LNP holds: {lnp_holds}/{neither_ab}")
        c_prime_covers = sum(1 for nc in neither_cases
                             if nc["info"]["route_c_prime"])
        print(f"  Route C' covers: {c_prime_covers}/{neither_ab}")
        print(f"\n  Tree types:")
        type_counts = Counter(nc["tree_type"] for nc in neither_cases)
        for tt, cnt in type_counts.most_common():
            print(f"    {tt}: {cnt}")
        print(f"\n  All neither-AB cases:")
        for nc in neither_cases[:30]:
            info = nc["info"]
            print(f"    n={nc['n']} k={nc['k']} mode={nc['mode']} "
                  f"type={nc['tree_type']} "
                  f"leaves={info['n_leaves']} L_out={info['n_lout']} "
                  f"max_deg_in_S={info['max_deg_in_s']} "
                  f"max_leaf_nbrs={info['max_leaf_nbrs']} "
                  f"C'={'Y' if info['route_c_prime'] else 'N'}")
            print(f"      S={nc['s']}  leaf_nbrs={info['leaf_nbrs']}")
    else:
        print("  *** Route A or B covers ALL cases! ***")

    # ===== ROUTE C' FAILURES =====
    print(f"\n{'='*75}")
    print(f"ROUTE C' FAILURES ({len(c_prime_fails)} cases)")
    print(f"{'='*75}")
    if c_prime_fails:
        # Cases where S doesn't contain a support vertex with >= 2 leaf-children
        # But LNP still holds — how?
        lnp_holds = sum(1 for cf in c_prime_fails
                        if cf["max_leaf_nbrs"] >= 2)
        print(f"  LNP still holds: {lnp_holds}/{len(c_prime_fails)}")
        print(f"\n  These cases have max_leaf_nbrs >= 2 via vertices whose")
        print(f"  leaf-neighbors come from the tree structure, not from being")
        print(f"  a support vertex with 2+ leaf-children in the tree.")
        print(f"\n  Sample C'-failures (first 15):")
        for cf in c_prime_fails[:15]:
            print(f"    n={cf['n']} k={cf['k']} mode={cf['mode']} "
                  f"type={cf['tree_type']} "
                  f"leaves={cf['n_leaves']} max_leaf_nbrs={cf['max_leaf_nbrs']} "
                  f"support_2plus_in_tree={sorted(cf['support_2plus'])}")
    else:
        print("  *** Route C' covers ALL cases! ***")

    # ===== TIGHT LNP CASES =====
    print(f"\n{'='*75}")
    print(f"TIGHT LNP CASES (max_leaf_nbrs = 2, {len(tight_lnp)} cases)")
    print(f"{'='*75}")
    if tight_lnp:
        print(f"\n  Tree types:")
        type_counts = Counter(tl["tree_type"] for tl in tight_lnp)
        for tt, cnt in type_counts.most_common():
            print(f"    {tt}: {cnt}")
        print(f"\n  Route coverage of tight cases:")
        ta = sum(1 for tl in tight_lnp if tl["route_a"])
        tb = sum(1 for tl in tight_lnp if tl["route_b"])
        tcp = sum(1 for tl in tight_lnp if tl["route_c_prime"])
        print(f"    Route A: {ta}/{len(tight_lnp)}")
        print(f"    Route B: {tb}/{len(tight_lnp)}")
        print(f"    Route C': {tcp}/{len(tight_lnp)}")
        print(f"\n  All tight cases (max 30):")
        for tl in tight_lnp[:30]:
            print(f"    n={tl['n']} k={tl['k']} mode={tl['mode']} "
                  f"type={tl['tree_type']} "
                  f"leaves={tl['n_leaves']} L_out={tl['n_lout']} "
                  f"A={'Y' if tl['route_a'] else 'N'} "
                  f"B={'Y' if tl['route_b'] else 'N'} "
                  f"C'={'Y' if tl['route_c_prime'] else 'N'}")
            print(f"      S={tl['s']}  leaf_nbrs={tl['leaf_nbrs']}")
    else:
        print("  *** No tight cases (min max_leaf_nbrs >= 3)! ***")

    # ===== KEY QUESTION: Does |L_out| > k always hold? =====
    # Actually the data says NO (156 failures). So let's check:
    # Does |L_out| >= k always hold? (weaker: pigeonhole gives >= 1 each)
    # The data shows |L_out|-k = -1 in 5 cases. So |L_out| < k is possible!
    print(f"\n{'='*75}")
    print("KEY STRUCTURAL RELATIONSHIPS")
    print(f"{'='*75}")

    # Relationship: n_leaves vs k
    print("\n  Leaf count vs k (all cases):")
    leaf_vs_k = Counter()
    for n in range(5, max_n + 1):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]
        for line in lines:
            tn, adj_data = parse_graph6(line)
            deg = [len(adj_data[v]) for v in range(n)]
            leaves = sum(1 for v in range(n) if deg[v] == 1)
            poly = independence_poly(n, adj_data)
            mode = compute_mode(poly)
            for s, k in find_maximal_is_below_mode(n, adj_data, mode):
                leaf_vs_k[(leaves >= 2 * k + 1, leaves >= k + 2)] = \
                    leaf_vs_k.get((leaves >= 2*k+1, leaves >= k+2), 0) + 1

    print(f"    l >= 2k+1 (strong pigeonhole): "
          f"{sum(v for (a,b),v in leaf_vs_k.items() if a)}")
    print(f"    l >= k+2: "
          f"{sum(v for (a,b),v in leaf_vs_k.items() if b)}")

    # ===== PROOF ROUTE ASSESSMENT =====
    print(f"\n\n{'='*75}")
    print("PROOF ROUTE ASSESSMENT")
    print(f"{'='*75}")

    if route_c_prime_ok == total:
        print("\n  *** Route C' covers 100%! ***")
        print("  Proof path: every maximal IS below mode contains a support")
        print("  vertex with >= 2 leaf-children in the tree.")
    elif route_a_ok == total:
        print("\n  *** Route A covers 100%! ***")
        print("  Proof path: |L_out| > k always holds when k < mode.")
    elif route_a_ok + route_b_ok - both_ab_ok == total:
        print("\n  *** Routes A+B combined cover 100%! ***")
        print("  Proof path: either |L_out| > k or max deg in S >= k+1.")
    else:
        a_or_b = route_a_ok + route_b_ok - both_ab_ok
        gap = total - a_or_b
        print(f"\n  Routes A+B leave {gap} uncovered cases.")
        if route_c_prime_ok + a_or_b >= total:
            print("  Adding Route C' covers the rest.")
        else:
            print("  Need additional structural argument for remaining cases.")

    print(f"\n  Route C' alone: {route_c_prime_ok}/{total} "
          f"({100*route_c_prime_ok/total:.2f}%)")
    print(f"  Gap from C': {total - route_c_prime_ok} cases")
    if total - route_c_prime_ok > 0:
        print("  The C'-gap cases have max_leaf_nbrs >= 2 via a vertex whose")
        print("  leaf-neighbors are NOT all children in a rooting — they come")
        print("  from multiple directions in the tree.")


if __name__ == "__main__":
    main()
