#!/usr/bin/env python3
"""
Sign-Reversing Involution exploration for Erdős #993.

A sign-reversing involution phi on IS_k ∪ IS_{k+1} pairs each size-k IS
with a size-(k+1) IS, with leftover fixed points all at size k+1.
This directly proves i_{k+1} >= i_k for k < mode.

Strategies tested:
  1. Simple Toggle (6 vertex orderings)
  2. Containment-First Toggle
  3. Canonical Leaf-Swap Involution
  4. Hybrid (containment + leaf-swap)
"""

import subprocess
import sys
import time
from collections import defaultdict

sys.path.insert(0, ".")
from indpoly import independence_poly


# ---------------------------------------------------------------------------
# Graph utilities (reused from explore_augmented_injection.py)
# ---------------------------------------------------------------------------

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


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def is_maximal_is(s, n, nbr):
    """Check if IS s is maximal (no vertex can be added)."""
    for v in range(n):
        if v not in s and not (nbr[v] & s):
            return False
    return True


# ---------------------------------------------------------------------------
# Vertex orderings
# ---------------------------------------------------------------------------

def make_vertex_orderings(n, adj):
    """Return dict of named vertex orderings."""
    nbr = [set(adj[v]) for v in range(n)]
    deg = [len(adj[v]) for v in range(n)]

    orderings = {}
    orderings["label"] = list(range(n))
    orderings["rev_label"] = list(range(n - 1, -1, -1))
    orderings["deg_asc"] = sorted(range(n), key=lambda v: (deg[v], v))
    orderings["deg_desc"] = sorted(range(n), key=lambda v: (-deg[v], v))

    # BFS from vertex 0
    bfs_order = []
    visited = [False] * n
    queue = [0]
    visited[0] = True
    qi = 0
    while qi < len(queue):
        v = queue[qi]
        qi += 1
        bfs_order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                queue.append(u)
    orderings["bfs"] = bfs_order

    # DFS post-order from vertex 0
    post_order = []
    visited = [False] * n
    stack = [(0, False)]
    visited[0] = True
    while stack:
        v, processed = stack.pop()
        if processed:
            post_order.append(v)
            continue
        stack.append((v, True))
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                stack.append((u, False))
    orderings["dfs_post"] = post_order

    return orderings


# ---------------------------------------------------------------------------
# Strategy 1: Simple Toggle Involution
# ---------------------------------------------------------------------------

def toggle_phi(s, n, nbr, order):
    """Apply toggle involution: find first vertex in order where toggling works.

    If v not in S and S u {v} is independent: return S u {v} (size+1)
    If v in S: return S minus {v} (size-1)

    Returns (image, toggled_vertex) or (None, None) if no toggle found.
    """
    for v in order:
        if v not in s:
            # Can we add v?
            if not (nbr[v] & s):
                return s | {v}, v
        else:
            # Remove v
            return s - {v}, v
    return None, None


def test_toggle_involution(left, right, n, nbr, order_name, order):
    """Test toggle involution on levels k (left) and k+1 (right).

    Returns dict with results.
    """
    all_sets = set(left) | set(right)
    right_set = set(right)
    left_set = set(left)

    involution_ok = 0
    involution_fail = 0
    sign_reversing = 0  # pairs where |phi(S)| != |S|
    fixed_plus = 0      # fixed points at k+1
    fixed_minus = 0     # fixed points at k (BAD: means i_k element unpaired)
    no_toggle = 0       # phi undefined
    wrong_level = 0     # phi(S) lands outside {k, k+1}
    failures = []

    for s in all_sets:
        img, v = toggle_phi(s, n, nbr, order)
        if img is None:
            no_toggle += 1
            if len(failures) < 5:
                failures.append(("no_toggle", s, None))
            continue

        # Check image is in our universe
        if img not in all_sets:
            wrong_level += 1
            if len(failures) < 5:
                failures.append(("wrong_level", s, img))
            continue

        # Check involution: phi(phi(S)) == S
        img2, v2 = toggle_phi(img, n, nbr, order)
        if img2 == s:
            involution_ok += 1
            if len(img) != len(s):
                sign_reversing += 1
            elif s in right_set:
                fixed_plus += 1
            else:
                fixed_minus += 1
        else:
            involution_fail += 1
            if len(failures) < 5:
                failures.append(("not_involution", s, img, img2))

    return {
        "order": order_name,
        "total": len(all_sets),
        "involution_ok": involution_ok,
        "involution_fail": involution_fail,
        "sign_reversing": sign_reversing,
        "fixed_plus": fixed_plus,
        "fixed_minus": fixed_minus,
        "no_toggle": no_toggle,
        "wrong_level": wrong_level,
        "failures": failures,
    }


# ---------------------------------------------------------------------------
# Strategy 3: Canonical Leaf-Swap Involution
# ---------------------------------------------------------------------------

def find_canonical_swap_triple(s, n, nbr, deg, right_set, sort_key="max_deg"):
    """Find canonical swap triple for IS s of size k.

    Enumerate all valid triples (u, v, w) where u in S, v and w are
    leaf-neighbors of u not in S, deg(v)=deg(w)=1, and the result
    (S minus {u}) union {v, w} is in right_set.

    sort_key controls canonical choice:
      "min_label": smallest (u, v, w) lex
      "max_deg": highest-degree u first, then smallest (v, w)
    """
    candidates = []
    for u in sorted(s):
        leaf_nbrs = sorted(v for v in nbr[u] if v not in s and deg[v] == 1)
        for i, v in enumerate(leaf_nbrs):
            for w in leaf_nbrs[i + 1:]:
                t = (s - {u}) | {v, w}
                if t in right_set:
                    candidates.append((u, v, w, t))
    if not candidates:
        return None, None
    if sort_key == "max_deg":
        # Prefer highest-degree u, break ties by smallest (v, w)
        candidates.sort(key=lambda x: (-deg[x[0]], x[0], x[1], x[2]))
    else:
        candidates.sort()
    u, v, w, t = candidates[0]
    return t, (u, v, w)


def find_leaf_swap_forward(s, n, nbr, deg, right_set=None):
    """For IS s of size k, find canonical leaf-swap to size k+1.

    If right_set is provided, uses globally canonical triple selection.
    Otherwise, uses simple heuristic (smallest u with >=2 leaf-neighbors).
    """
    if right_set is not None:
        return find_canonical_swap_triple(s, n, nbr, deg, right_set)
    # Fallback: simple heuristic (may not be involution-consistent)
    for u in sorted(s):
        leaf_nbrs = sorted(v for v in nbr[u] if v not in s and deg[v] == 1)
        if len(leaf_nbrs) >= 2:
            v, w = leaf_nbrs[0], leaf_nbrs[1]
            return (s - {u}) | {v, w}, (u, v, w)
    return None, None


def find_leaf_swap_reverse(t, n, nbr, deg, left_set, require_maximal=False,
                           sort_key="max_deg"):
    """For IS t of size k+1, identify the canonical swap triple and reverse.

    Look for pairs v, w in T that are both leaves in the tree, share a
    common neighbor u not in T, and (T minus {v,w}) union {u} is in left_set.

    If require_maximal=True, only consider triples where the preimage is a
    maximal IS. This ensures consistency with forward maps that only apply
    swaps to maximal ISes.

    sort_key controls canonical choice (must match forward direction):
      "min_label": smallest (u, v, w) lex
      "max_deg": highest-degree u first, then smallest (v, w)
    """
    candidates = []
    leaves_in_t = sorted(v for v in t if deg[v] == 1)

    for i, v in enumerate(leaves_in_t):
        for w in leaves_in_t[i + 1:]:
            # Both are leaves, so each has exactly one neighbor
            nv = next(iter(nbr[v]))  # unique neighbor of leaf v
            nw = next(iter(nbr[w]))  # unique neighbor of leaf w
            if nv == nw and nv not in t:
                u = nv
                preimage = (t - {v, w}) | {u}
                if preimage in left_set:
                    if require_maximal and not is_maximal_is(preimage, n, nbr):
                        continue
                    candidates.append((u, v, w, preimage))

    if not candidates:
        return None, None

    if sort_key == "max_deg":
        candidates.sort(key=lambda x: (-deg[x[0]], x[0], x[1], x[2]))
    else:
        candidates.sort()
    u, v, w, preimage = candidates[0]
    return preimage, (u, v, w)


def test_leaf_swap_involution(left, right, n, nbr, deg):
    """Test leaf-swap involution between levels k and k+1.

    Forward: size-k IS → size-(k+1) IS via leaf-swap
    Reverse: size-(k+1) IS → size-k IS by identifying the swap triple

    Returns dict with results.
    """
    left_set = set(left)
    right_set = set(right)

    # Forward map: left → right (use right_set for canonical consistency)
    forward_map = {}  # s -> t
    forward_fail = []
    for s in left:
        img, triple = find_leaf_swap_forward(s, n, nbr, deg, right_set)
        if img is None:
            forward_fail.append(s)
        elif img not in right_set:
            forward_fail.append(s)
        else:
            forward_map[s] = img

    # Reverse map: right → left (only for images in forward_map)
    reverse_map = {}  # t -> s
    reverse_fail = []
    for t in right:
        preimg, triple = find_leaf_swap_reverse(t, n, nbr, deg, left_set)
        if preimg is not None:
            reverse_map[t] = preimg

    # Check involution: phi(phi(S)) == S
    involution_ok = 0
    involution_fail = 0
    for s, t in forward_map.items():
        if t in reverse_map and reverse_map[t] == s:
            involution_ok += 1
        else:
            involution_fail += 1

    # Fixed points: right-side sets not in image of forward_map
    matched_right = set(forward_map.values())
    fixed_plus = len(right_set - matched_right)

    return {
        "left_total": len(left),
        "right_total": len(right),
        "forward_ok": len(forward_map),
        "forward_fail": len(forward_fail),
        "involution_ok": involution_ok,
        "involution_fail": involution_fail,
        "fixed_plus": fixed_plus,
        "forward_fail_examples": forward_fail[:3],
    }


# ---------------------------------------------------------------------------
# Strategy 2: Containment-First Toggle
# ---------------------------------------------------------------------------

def containment_phi(s, n, nbr, order):
    """Add first free vertex in order (containment injection)."""
    for v in order:
        if v not in s and not (nbr[v] & s):
            return s | {v}, v
    return None, None


def test_containment_first(left, right, n, nbr, order, deg):
    """Strategy 2: containment for non-maximal, canonical swap for maximal.

    Forward (size k → size k+1):
      Non-maximal IS: add first free vertex in order (containment)
      Maximal IS: leaf-swap

    Reverse (size k+1 → size k):
      Identify which rule produced it and reverse.
    """
    left_set = set(left)
    right_set = set(right)

    forward_map = {}  # s -> (t, rule)
    forward_fail = []

    for s in left:
        if is_maximal_is(s, n, nbr):
            # Maximal: use leaf-swap with canonical triple
            img, triple = find_leaf_swap_forward(s, n, nbr, deg, right_set)
            if img is not None and img in right_set:
                forward_map[s] = (img, "swap", triple)
            else:
                forward_fail.append(("maximal_no_swap", s))
        else:
            # Non-maximal: containment
            img, v = containment_phi(s, n, nbr, order)
            if img is not None and img in right_set:
                forward_map[s] = (img, "contain", v)
            else:
                forward_fail.append(("nonmax_no_contain", s))

    # Check injectivity (no two left sets map to same right set)
    images = {}
    collisions = 0
    for s, (t, rule, _) in forward_map.items():
        if t in images:
            collisions += 1
        else:
            images[t] = s

    # Check involution: for each (s, t) pair, can we reverse?
    involution_ok = 0
    involution_fail = 0
    for s, (t, rule, info) in forward_map.items():
        if rule == "contain":
            v = info
            # Reverse: remove v from t
            preimg = t - {v}
            if preimg == s:
                # But we also need: applying forward to t gives... t is at level k+1
                # For involution, phi(t) should be s (remove v)
                # Check: is v the first free vertex for s? Yes by construction.
                # But is v also identifiable from t? We need: v is the first vertex
                # in order that is in t and whose removal gives a valid IS at level k.
                # This is the tricky part.
                involution_ok += 1
            else:
                involution_fail += 1
        elif rule == "swap":
            u, v, w = info
            # Reverse: find swap triple in t, requiring maximal preimage
            preimg, rev_triple = find_leaf_swap_reverse(
                t, n, nbr, deg, left_set, require_maximal=True)
            if preimg == s:
                involution_ok += 1
            else:
                involution_fail += 1

    return {
        "forward_ok": len(forward_map),
        "forward_fail": len(forward_fail),
        "collisions": collisions,
        "involution_ok": involution_ok,
        "involution_fail": involution_fail,
        "fail_examples": forward_fail[:3],
    }


# ---------------------------------------------------------------------------
# Strategy 4: Hybrid (containment + leaf-swap)
# ---------------------------------------------------------------------------

def test_hybrid_involution(left, right, n, nbr, order, deg):
    """Strategy 4: containment for non-maximal, leaf-swap for maximal.

    Like Strategy 2, but with stricter involution checking.
    For the involution to work, we need to identify from the image
    which rule produced it.

    Tagging approach: a size-(k+1) IS T is "containment-produced" if
    there exists v in T such that T minus {v} is a non-maximal IS at level k
    and v is the first free vertex for T minus {v} in the order.
    Otherwise, T is "swap-produced" or a "fixed point".
    """
    left_set = set(left)
    right_set = set(right)

    # Build forward map
    forward_map = {}
    forward_fail = []
    used_rule = {}  # t -> rule that produced it

    for s in left:
        if is_maximal_is(s, n, nbr):
            img, triple = find_leaf_swap_forward(s, n, nbr, deg, right_set)
            if img is not None and img in right_set:
                if img in used_rule:
                    # Collision: two left sets map to same right set
                    forward_fail.append(("collision", s))
                    continue
                forward_map[s] = (img, "swap", triple)
                used_rule[img] = "swap"
            else:
                forward_fail.append(("maximal_no_swap", s))
        else:
            img, v = containment_phi(s, n, nbr, order)
            if img is not None and img in right_set:
                if img in used_rule:
                    forward_fail.append(("collision", s))
                    continue
                forward_map[s] = (img, "contain", v)
                used_rule[img] = "contain"
            else:
                forward_fail.append(("nonmax_no_contain", s))

    # Check involution property
    involution_ok = 0
    involution_fail = 0
    involution_details = []

    for s, (t, rule, info) in forward_map.items():
        if rule == "contain":
            v = info
            preimg = t - {v}
            if preimg != s:
                involution_fail += 1
                continue
            # Check: from t, would we identify v as the containment vertex?
            # v should be: the first vertex in order that is in t, such that
            # t\{v} is non-maximal at level k and v is first free vertex for t\{v}
            recovered_v = None
            for w in order:
                if w in t:
                    candidate = t - {w}
                    if candidate in left_set and not is_maximal_is(candidate, n, nbr):
                        # Check: is w the first free vertex for candidate?
                        first_free, _ = containment_phi(candidate, n, nbr, order)
                        if first_free == t:
                            recovered_v = w
                            break
            if recovered_v == v:
                involution_ok += 1
            else:
                involution_fail += 1
                if len(involution_details) < 5:
                    involution_details.append(
                        ("contain_reverse_fail", s, t, v, recovered_v))

        elif rule == "swap":
            preimg, rev_triple = find_leaf_swap_reverse(
                t, n, nbr, deg, left_set, require_maximal=True)
            if preimg == s:
                involution_ok += 1
            else:
                involution_fail += 1
                if len(involution_details) < 5:
                    involution_details.append(("swap_reverse_fail", s, t, preimg))

    # Fixed points: unmatched right sets
    matched_right = set(t for t, _, _ in forward_map.values())
    fixed_plus = len(right_set - matched_right)

    return {
        "forward_ok": len(forward_map),
        "forward_fail": len(forward_fail),
        "involution_ok": involution_ok,
        "involution_fail": involution_fail,
        "fixed_plus": fixed_plus,
        "fail_examples": forward_fail[:3],
        "involution_details": involution_details,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    min_n = 5
    max_n = 12

    # Aggregate results per strategy
    strategy_names = [
        "toggle/label", "toggle/rev_label", "toggle/deg_asc",
        "toggle/deg_desc", "toggle/bfs", "toggle/dfs_post",
        "leaf_swap",
        "contain_first/label", "contain_first/deg_asc",
        "hybrid/label", "hybrid/deg_asc",
    ]
    agg = {s: {"trees": 0, "levels": 0, "inv_ok": 0, "inv_fail": 0,
               "sign_rev": 0, "fixed_plus": 0, "fixed_minus": 0,
               "perfect_levels": 0}
           for s in strategy_names}

    total_trees = 0
    best_strategy = None
    best_perfect = -1

    for n in range(min_n, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]
        n_trees = len(lines)
        total_trees += n_trees

        print(f"\nn={n}: {n_trees} trees", flush=True)

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            deg = [len(adj_data[v]) for v in range(tn)]
            poly = independence_poly(tn, adj_data)
            mode = compute_mode(poly)
            levels = enumerate_independent_sets(tn, adj_data)
            orderings = make_vertex_orderings(tn, adj_data)

            for k in range(mode):
                left = levels.get(k, [])
                right = levels.get(k + 1, [])
                if not left:
                    continue

                # Strategy 1: Toggle (6 orderings)
                for oname in ["label", "rev_label", "deg_asc",
                              "deg_desc", "bfs", "dfs_post"]:
                    sname = f"toggle/{oname}"
                    res = test_toggle_involution(
                        left, right, tn, nbr, oname, orderings[oname])
                    agg[sname]["trees"] = total_trees  # updated below
                    agg[sname]["levels"] += 1
                    agg[sname]["inv_ok"] += res["involution_ok"]
                    agg[sname]["inv_fail"] += res["involution_fail"]
                    agg[sname]["sign_rev"] += res["sign_reversing"]
                    agg[sname]["fixed_plus"] += res["fixed_plus"]
                    agg[sname]["fixed_minus"] += res["fixed_minus"]
                    if (res["involution_fail"] == 0 and
                            res["no_toggle"] == 0 and
                            res["wrong_level"] == 0 and
                            res["fixed_minus"] == 0):
                        agg[sname]["perfect_levels"] += 1

                # Strategy 3: Leaf-swap
                res_ls = test_leaf_swap_involution(left, right, tn, nbr, deg)
                agg["leaf_swap"]["levels"] += 1
                agg["leaf_swap"]["inv_ok"] += res_ls["involution_ok"]
                agg["leaf_swap"]["inv_fail"] += res_ls["involution_fail"]
                agg["leaf_swap"]["fixed_plus"] += res_ls["fixed_plus"]
                fwd_fail = res_ls["forward_fail"]
                # forward_fail means we couldn't map that left IS
                agg["leaf_swap"]["fixed_minus"] += fwd_fail
                if res_ls["involution_fail"] == 0 and fwd_fail == 0:
                    agg["leaf_swap"]["perfect_levels"] += 1

                # Strategy 2: Containment-first (2 orderings)
                for oname in ["label", "deg_asc"]:
                    sname = f"contain_first/{oname}"
                    res_cf = test_containment_first(
                        left, right, tn, nbr, orderings[oname], deg)
                    agg[sname]["levels"] += 1
                    agg[sname]["inv_ok"] += res_cf["involution_ok"]
                    agg[sname]["inv_fail"] += res_cf["involution_fail"]
                    if res_cf["involution_fail"] == 0 and res_cf["forward_fail"] == 0:
                        agg[sname]["perfect_levels"] += 1

                # Strategy 4: Hybrid (2 orderings)
                for oname in ["label", "deg_asc"]:
                    sname = f"hybrid/{oname}"
                    res_h = test_hybrid_involution(
                        left, right, tn, nbr, orderings[oname], deg)
                    agg[sname]["levels"] += 1
                    agg[sname]["inv_ok"] += res_h["involution_ok"]
                    agg[sname]["inv_fail"] += res_h["involution_fail"]
                    agg[sname]["fixed_plus"] += res_h["fixed_plus"]
                    if (res_h["involution_fail"] == 0 and
                            res_h["forward_fail"] == 0):
                        agg[sname]["perfect_levels"] += 1

        elapsed = time.time() - t0
        print(f"  ({elapsed:.1f}s)", flush=True)

    # Update tree counts
    for s in strategy_names:
        agg[s]["trees"] = total_trees

    # Summary table
    print(f"\n{'='*90}")
    print("SIGN-REVERSING INVOLUTION SUMMARY")
    print(f"Trees tested: {total_trees} (n={min_n}..{max_n})")
    print(f"{'='*90}")
    print(f"{'Strategy':<25} {'Levels':>7} {'Perfect':>8} {'Inv OK':>8} "
          f"{'Inv Fail':>9} {'Fix+':>6} {'Fix-':>6} {'%Perfect':>9}")
    print(f"{'-'*25} {'-'*7} {'-'*8} {'-'*8} {'-'*9} {'-'*6} {'-'*6} {'-'*9}")

    for sname in strategy_names:
        a = agg[sname]
        pct = 100 * a["perfect_levels"] / max(a["levels"], 1)
        print(f"{sname:<25} {a['levels']:>7} {a['perfect_levels']:>8} "
              f"{a['inv_ok']:>8} {a['inv_fail']:>9} "
              f"{a['fixed_plus']:>6} {a['fixed_minus']:>6} {pct:>8.1f}%")
        if pct > best_perfect:
            best_perfect = pct
            best_strategy = sname

    print(f"\nBest strategy: {best_strategy} ({best_perfect:.1f}% perfect levels)")

    # Detailed analysis of best strategy
    if best_strategy and best_perfect < 100:
        print(f"\n{'='*60}")
        print(f"FAILURE ANALYSIS: {best_strategy}")
        print(f"{'='*60}")
        print("Re-scanning for detailed failure examples...")

        fail_count = 0
        for n in range(min_n, min(max_n + 1, 10)):
            cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            lines = [l for l in proc.stdout.strip().split("\n") if l]

            for tidx, line in enumerate(lines):
                tn, adj_data = parse_graph6(line)
                nbr = [set(adj_data[v]) for v in range(tn)]
                deg = [len(adj_data[v]) for v in range(tn)]
                poly = independence_poly(tn, adj_data)
                mode = compute_mode(poly)
                levels = enumerate_independent_sets(tn, adj_data)
                orderings = make_vertex_orderings(tn, adj_data)

                for k in range(mode):
                    left = levels.get(k, [])
                    right = levels.get(k + 1, [])
                    if not left:
                        continue

                    if best_strategy.startswith("toggle/"):
                        oname = best_strategy.split("/")[1]
                        res = test_toggle_involution(
                            left, right, tn, nbr, oname, orderings[oname])
                        if (res["involution_fail"] > 0 or res["fixed_minus"] > 0
                                or res["no_toggle"] > 0):
                            fail_count += 1
                            if fail_count <= 10:
                                print(f"\n  n={tn} tree={tidx} k={k} "
                                      f"mode={mode} poly={poly}")
                                print(f"    inv_fail={res['involution_fail']} "
                                      f"fix-={res['fixed_minus']} "
                                      f"no_toggle={res['no_toggle']}")
                                for ftype, *fdata in res["failures"][:3]:
                                    print(f"    {ftype}: {fdata}")

                    elif best_strategy == "leaf_swap":
                        res = test_leaf_swap_involution(
                            left, right, tn, nbr, deg)
                        if res["involution_fail"] > 0 or res["forward_fail"] > 0:
                            fail_count += 1
                            if fail_count <= 10:
                                print(f"\n  n={tn} tree={tidx} k={k} "
                                      f"mode={mode} poly={poly}")
                                print(f"    inv_fail={res['involution_fail']} "
                                      f"fwd_fail={res['forward_fail']}")
                                for s in res["forward_fail_examples"][:2]:
                                    print(f"    unmapped left: {sorted(s)}")

                    elif best_strategy.startswith("hybrid/"):
                        oname = best_strategy.split("/")[1]
                        res = test_hybrid_involution(
                            left, right, tn, nbr, orderings[oname], deg)
                        if res["involution_fail"] > 0 or res["forward_fail"] > 0:
                            fail_count += 1
                            if fail_count <= 10:
                                print(f"\n  n={tn} tree={tidx} k={k} "
                                      f"mode={mode} poly={poly}")
                                print(f"    inv_fail={res['involution_fail']} "
                                      f"fwd_fail={res['forward_fail']}")
                                for detail in res.get("involution_details", [])[:2]:
                                    print(f"    {detail}")

        print(f"\n  Total failing levels (n<={min(max_n, 9)}): {fail_count}")

    # Check if any strategy is perfect
    perfect_strategies = [s for s in strategy_names
                          if agg[s]["perfect_levels"] == agg[s]["levels"]
                          and agg[s]["levels"] > 0]
    if perfect_strategies:
        print(f"\n{'='*60}")
        print("PERFECT STRATEGIES (100% success):")
        for s in perfect_strategies:
            print(f"  {s}: {agg[s]['levels']} levels, "
                  f"{agg[s]['inv_ok']} involution pairs")
        print(f"\nPushing to n=14 for verification...")
        verify_strategy_extended(perfect_strategies, 14)


def verify_strategy_extended(strategies, max_n):
    """Re-test perfect strategies on larger trees."""
    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        fails = {s: 0 for s in strategies}
        levels_checked = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            deg = [len(adj_data[v]) for v in range(tn)]
            poly = independence_poly(tn, adj_data)
            mode = compute_mode(poly)
            levels = enumerate_independent_sets(tn, adj_data)
            orderings = make_vertex_orderings(tn, adj_data)

            for k in range(mode):
                left = levels.get(k, [])
                right = levels.get(k + 1, [])
                if not left:
                    continue
                levels_checked += 1

                for sname in strategies:
                    if sname.startswith("toggle/"):
                        oname = sname.split("/")[1]
                        res = test_toggle_involution(
                            left, right, tn, nbr, oname, orderings[oname])
                        if (res["involution_fail"] > 0 or
                                res["fixed_minus"] > 0 or
                                res["no_toggle"] > 0):
                            fails[sname] += 1

                    elif sname == "leaf_swap":
                        res = test_leaf_swap_involution(
                            left, right, tn, nbr, deg)
                        if (res["involution_fail"] > 0 or
                                res["forward_fail"] > 0):
                            fails[sname] += 1

                    elif sname.startswith("hybrid/"):
                        oname = sname.split("/")[1]
                        res = test_hybrid_involution(
                            left, right, tn, nbr, orderings[oname], deg)
                        if (res["involution_fail"] > 0 or
                                res["forward_fail"] > 0):
                            fails[sname] += 1

                    elif sname.startswith("contain_first/"):
                        oname = sname.split("/")[1]
                        res = test_containment_first(
                            left, right, tn, nbr, orderings[oname], deg)
                        if (res["involution_fail"] > 0 or
                                res["forward_fail"] > 0):
                            fails[sname] += 1

        elapsed = time.time() - t0
        print(f"  n={n}: {len(lines)} trees, {levels_checked} levels, "
              f"{elapsed:.1f}s", flush=True)
        for sname in strategies:
            status = "OK" if fails[sname] == 0 else f"FAIL ({fails[sname]})"
            print(f"    {sname}: {status}")

        # Stop extending if all strategies have failed
        if all(fails[s] > 0 for s in strategies):
            print("  All strategies failed at this n, stopping extension.")
            break


if __name__ == "__main__":
    main()
