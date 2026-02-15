#!/usr/bin/env python3
"""
Sign-Reversing Involution with private-neighbor swaps.

Modification of explore_sri.py: uses private neighbors (not just leaf-
neighbors) for the swap step. This handles the n=18 LNP counterexample
where no vertex has >= 2 leaf-neighbors but max_priv = 5.

Strategy: Containment-First with private-neighbor swap.

Forward (size k -> size k+1):
  Non-maximal IS: add first free vertex in label order (containment)
  Maximal IS: canonical private-neighbor swap
    - u = max-degree vertex in S with priv(u) >= 2
    - v, w = two smallest-label private neighbors of u
    - T = (S \ {u}) ∪ {v, w}

Reverse (size k+1 -> size k):
  For each u ∉ T with exactly 2 T-neighbors {v, w}:
    S' = (T \ {v,w}) ∪ {u}
    If S' is maximal IS, and canonical forward of S' picks (u,v,w):
      phi(T) = S'
  If no swap reverse found, try containment reverse.
  Otherwise: T is a fixed point.
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


def compute_mode(poly):
    if not poly:
        return 0
    return max(range(len(poly)), key=lambda i: poly[i])


def is_maximal_is(s, n, nbr):
    for v in range(n):
        if v not in s and not (nbr[v] & s):
            return False
    return True


def compute_private_neighbors(s, n, nbr):
    """For each u in S, return list of private neighbors (sorted by label)."""
    privs = {}
    for u in s:
        priv_list = []
        for v in sorted(nbr[u]):
            if v not in s:
                s_nbrs = nbr[v] & s
                if len(s_nbrs) == 1:  # only u
                    priv_list.append(v)
        privs[u] = priv_list
    return privs


def canonical_swap_triple(s, n, nbr, deg):
    """Find canonical swap triple for maximal IS s.

    Returns (u, v, w) where:
      u = max-degree vertex in S with >= 2 private neighbors
          (break ties by smallest label)
      v, w = two smallest-label private neighbors of u
    Returns None if no vertex has >= 2 private neighbors.
    """
    privs = compute_private_neighbors(s, n, nbr)

    # Find u: max degree among those with priv >= 2, break ties smallest label
    candidates = [(u, privs[u]) for u in s if len(privs[u]) >= 2]
    if not candidates:
        return None

    candidates.sort(key=lambda x: (-deg[x[0]], x[0]))
    u, priv_list = candidates[0]
    v, w = priv_list[0], priv_list[1]
    return (u, v, w)


def phi_forward(s, k, n, nbr, deg, right_set, order):
    """Apply containment-first forward map.

    Returns (image, rule, info) or (None, "fail", None).
    """
    if is_maximal_is(s, n, nbr):
        # Maximal: private-neighbor swap
        triple = canonical_swap_triple(s, n, nbr, deg)
        if triple is None:
            return None, "maximal_no_swap", None
        u, v, w = triple
        t = (s - {u}) | {v, w}
        if t in right_set:
            return t, "swap", (u, v, w)
        else:
            return None, "swap_not_in_right", (u, v, w)
    else:
        # Non-maximal: containment (add first free vertex)
        for v in order:
            if v not in s and not (nbr[v] & s):
                t = s | {v}
                if t in right_set:
                    return t, "contain", v
                else:
                    return None, "contain_not_in_right", v
        return None, "nonmax_no_free", None


def phi_reverse(t, k_plus_1, n, nbr, deg, left_set, order):
    """Identify the preimage of t under phi.

    CRITICAL: try containment reverse FIRST. Non-maximal preimages always
    use containment in the forward direction, so they must be identified
    first. Only try swap reverse if containment doesn't match.

    Returns (preimage, rule, info) or (None, "fixed", None).
    """
    # Try containment reverse FIRST: find v ∈ T such that T \ {v} is
    # non-maximal IS at level k, and v is the first free vertex for T \ {v}
    for v in order:
        if v not in t:
            continue
        preimage = t - {v}
        if preimage not in left_set:
            continue
        if is_maximal_is(preimage, n, nbr):
            continue  # maximal preimages use swap, not containment
        # Check: is v the first free vertex for preimage?
        for w in order:
            if w not in preimage and not (nbr[w] & preimage):
                if w == v:
                    return preimage, "contain", v
                else:
                    break  # first free vertex is w ≠ v

    # Try swap reverse: find u ∉ T with exactly 2 T-neighbors
    swap_candidates = []
    for u in range(n):
        if u in t:
            continue
        t_nbrs = nbr[u] & t
        if len(t_nbrs) == 2:
            v, w = sorted(t_nbrs)
            preimage = (t - {v, w}) | {u}
            # preimage is automatically independent (u has no neighbor in T\{v,w})
            if preimage not in left_set:
                continue
            if not is_maximal_is(preimage, n, nbr):
                continue
            # Check canonical: would forward map of preimage pick (u, v, w)?
            canon = canonical_swap_triple(preimage, n, nbr, deg)
            if canon == (u, v, w):
                swap_candidates.append((u, v, w, preimage))

    if swap_candidates:
        # Should be unique if canonical selection is consistent
        if len(swap_candidates) == 1:
            u, v, w, preimage = swap_candidates[0]
            return preimage, "swap", (u, v, w)
        else:
            # Multiple swap preimages: ambiguity
            return None, "swap_ambiguous", swap_candidates

    return None, "fixed", None


def test_involution(left, right, n, nbr, deg, order):
    """Test the containment-first private-neighbor involution.

    Returns dict with detailed results.
    """
    left_set = set(left)
    right_set = set(right)
    k = len(left[0]) if left else 0

    # Forward map
    forward_map = {}
    forward_fails = []
    for s in left:
        t, rule, info = phi_forward(s, k, n, nbr, deg, right_set, order)
        if t is not None:
            forward_map[s] = (t, rule, info)
        else:
            forward_fails.append((rule, s, info))

    # Check injectivity
    image_counts = defaultdict(int)
    for s, (t, _, _) in forward_map.items():
        image_counts[t] += 1
    collisions = sum(1 for c in image_counts.values() if c > 1)

    # Full involution check: for each element in left ∪ right, apply phi
    # and check phi(phi(x)) == x
    involution_ok = 0
    involution_fail = 0
    sign_reversing = 0
    fixed_plus = 0  # fixed points at k+1 (good)
    fixed_minus = 0  # fixed points at k (bad)

    # Check forward then reverse
    for s, (t, rule, info) in forward_map.items():
        preimg, rev_rule, rev_info = phi_reverse(
            t, k + 1, n, nbr, deg, left_set, order)
        if preimg == s:
            involution_ok += 1
            sign_reversing += 1
        else:
            involution_fail += 1

    # Check unmatched right sets (fixed points at k+1)
    matched_right = set(t for t, _, _ in forward_map.values())
    for t in right:
        if t not in matched_right:
            # Verify it's truly a fixed point (reverse gives None)
            preimg, rev_rule, _ = phi_reverse(
                t, k + 1, n, nbr, deg, left_set, order)
            if preimg is None:
                fixed_plus += 1
            else:
                # This right set maps back to something but wasn't matched forward
                involution_fail += 1

    # Check unmatched left sets (should not happen for involution)
    for s in left:
        if s not in forward_map:
            fixed_minus += 1

    perfect = (involution_fail == 0 and fixed_minus == 0 and
               len(forward_fails) == 0)

    return {
        "perfect": perfect,
        "forward_ok": len(forward_map),
        "forward_fail": len(forward_fails),
        "collisions": collisions,
        "involution_ok": involution_ok,
        "involution_fail": involution_fail,
        "sign_reversing": sign_reversing,
        "fixed_plus": fixed_plus,
        "fixed_minus": fixed_minus,
        "forward_fail_examples": forward_fails[:3],
    }


def main():
    max_n = 18

    print("=" * 75)
    print("CONTAINMENT-FIRST SRI with PRIVATE-NEIGHBOR SWAPS")
    print("=" * 75)

    total_trees = 0
    total_levels = 0
    perfect_levels = 0
    failed_levels = 0
    fail_examples = []

    for n in range(5, max_n + 1):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_trees = len(lines)
        total_trees += n_trees
        n_levels = 0
        n_perfect = 0
        n_fail = 0

        for tidx, line in enumerate(lines):
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            deg = [len(adj_data[v]) for v in range(tn)]
            poly = independence_poly(tn, adj_data)
            mode = compute_mode(poly)
            levels = enumerate_independent_sets(tn, adj_data)
            order = list(range(tn))  # label order

            for k in range(mode):
                left = levels.get(k, [])
                right = levels.get(k + 1, [])
                if not left:
                    continue

                total_levels += 1
                n_levels += 1

                res = test_involution(left, right, tn, nbr, deg, order)

                if res["perfect"]:
                    perfect_levels += 1
                    n_perfect += 1
                else:
                    failed_levels += 1
                    n_fail += 1
                    if len(fail_examples) < 10:
                        fail_examples.append({
                            "n": tn, "tree": tidx, "k": k, "mode": mode,
                            "g6": line.strip(),
                            "res": res,
                        })

        elapsed = time.time() - t0
        print(f"n={n:2d} ({n_trees:6d} trees, {elapsed:7.1f}s): "
              f"levels={n_levels:6d}  "
              f"perfect={n_perfect:6d}  "
              f"fail={n_fail:3d}",
              flush=True)

    print(f"\n{'='*75}")
    print("SUMMARY")
    print(f"{'='*75}")
    print(f"Total trees: {total_trees}")
    print(f"Total levels: {total_levels}")
    print(f"Perfect levels: {perfect_levels}/{total_levels} "
          f"({100*perfect_levels/total_levels:.4f}%)")
    print(f"Failed levels: {failed_levels}")

    if failed_levels == 0:
        print(f"\n*** SRI with private-neighbor swaps: 100% success through n={max_n}! ***")
    else:
        print(f"\nFailure examples:")
        for fe in fail_examples:
            res = fe["res"]
            print(f"  n={fe['n']} tree={fe['tree']} k={fe['k']} mode={fe['mode']}")
            print(f"    fwd_ok={res['forward_ok']} fwd_fail={res['forward_fail']} "
                  f"inv_ok={res['involution_ok']} inv_fail={res['involution_fail']} "
                  f"collisions={res['collisions']}")
            for rule, s, info in res["forward_fail_examples"]:
                print(f"    {rule}: S={sorted(s)} info={info}")


if __name__ == "__main__":
    main()
