#!/usr/bin/env python3
"""Verify the pendant tree mode bound for larger n.

For each tree T with a maximal IS S where all priv(u) <= 1:
- Check if mode(I(T)) <= k = |S|
- Record the tightest cases (mode closest to k)
- Characterize the structure

This is the remaining gap (Part 3) in the 1-Private Mode Conjecture.
"""

import subprocess
import sys
import time
import json

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


def find_all_maximal_is(n, adj):
    nbr = [set(adj[v]) for v in range(n)]
    result = []
    def backtrack(v, current, forbidden):
        if v == n:
            s = frozenset(current)
            for w in range(n):
                if w not in s and not (nbr[w] & s):
                    return
            result.append(s)
            return
        backtrack(v + 1, current, forbidden)
        if v not in forbidden:
            backtrack(v + 1, current + [v], forbidden | nbr[v])
    backtrack(0, [], set())
    return result


def main():
    print("PENDANT TREE MODE BOUND VERIFICATION")
    print("=" * 70)
    print()
    print("For each tree: find maximal IS with all priv <= 1.")
    print("Check: mode(I(T)) <= k for these IS.")
    print("Also track: ratio i_{k+1}/i_k (how close mode is to k+1)")
    print()

    results = {}

    for n in range(5, 21):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_pendant_trees = 0
        n_pendant_is = 0
        n_mode_leq_k = 0
        n_mode_gt_k = 0
        worst_ratio = 0.0  # max i_{k+1}/i_k
        worst_tree_info = None
        tight_cases = []  # cases where mode == k

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode_idx = max(range(len(poly)), key=lambda i: poly[i])
            nbr = [set(adj_data[v]) for v in range(tn)]

            all_mis = find_all_maximal_is(tn, adj_data)
            tree_has_pendant_is = False

            for s in all_mis:
                k = len(s)
                if k >= mode_idx:
                    continue

                # Compute priv for each u in S
                max_priv = 0
                all_priv_leq_1 = True
                priv_counts = {}
                for u in s:
                    priv_u = 0
                    for v in adj_data[u]:
                        if v not in s:
                            s_nbrs = len(nbr[v] & s)
                            if s_nbrs == 1:
                                priv_u += 1
                    priv_counts[u] = priv_u
                    max_priv = max(max_priv, priv_u)
                    if priv_u > 1:
                        all_priv_leq_1 = False

                if not all_priv_leq_1:
                    continue

                # This IS has all priv <= 1: it's a "1-Private" or "sub-1-Private" IS
                n_pendant_is += 1
                if not tree_has_pendant_is:
                    tree_has_pendant_is = True
                    n_pendant_trees += 1

                # Check mode bound
                if mode_idx <= k:
                    n_mode_leq_k += 1
                else:
                    n_mode_gt_k += 1
                    print(f"  *** VIOLATION at n={n}: k={k}, mode={mode_idx}, poly={poly}")

                # Track ratio i_{k+1}/i_k
                if k + 1 < len(poly) and poly[k] > 0:
                    ratio = poly[k + 1] / poly[k]
                    if ratio > worst_ratio:
                        worst_ratio = ratio
                        worst_tree_info = {
                            "n": n,
                            "k": k,
                            "mode": mode_idx,
                            "ratio": ratio,
                            "i_k": poly[k],
                            "i_k1": poly[k + 1],
                            "graph6": line.strip(),
                        }

                # Track tight cases (mode == k)
                if mode_idx == k and n <= 18:
                    # Compute structure
                    m_out = 0
                    n_priv0 = sum(1 for u in s if priv_counts[u] == 0)
                    n_priv1 = sum(1 for u in s if priv_counts[u] == 1)
                    for v in range(tn):
                        if v not in s:
                            s_count = len(nbr[v] & s)
                            if s_count == 0:
                                m_out += 1  # not dominated (shouldn't happen for maximal IS)
                    tight_cases.append({
                        "n": n, "k": k, "mode": mode_idx,
                        "priv0": n_priv0, "priv1": n_priv1,
                        "m_out": m_out,
                        "ratio": poly[k + 1] / poly[k] if k + 1 < len(poly) and poly[k] > 0 else 0,
                    })

        elapsed = time.time() - t0

        if n_pendant_is > 0:
            print(f"n={n:2d}: {n_pendant_trees:5d} trees, {n_pendant_is:6d} IS, "
                  f"mode<=k: {n_mode_leq_k}, violations: {n_mode_gt_k}, "
                  f"worst ratio: {worst_ratio:.4f} ({elapsed:.1f}s)")
        else:
            print(f"n={n:2d}: no 1-Private IS below mode ({elapsed:.1f}s)")

        results[n] = {
            "pendant_trees": n_pendant_trees,
            "pendant_is": n_pendant_is,
            "mode_leq_k": n_mode_leq_k,
            "violations": n_mode_gt_k,
            "worst_ratio": worst_ratio,
            "worst_tree": worst_tree_info,
            "tight_cases": tight_cases[:5] if tight_cases else [],
        }

        if n_mode_gt_k > 0:
            print("*** COUNTEREXAMPLE FOUND! Stopping. ***")
            break

    print()
    print("=" * 70)
    if all(r["violations"] == 0 for r in results.values()):
        print("No violations found. Pendant tree mode bound holds.")
    else:
        print("VIOLATIONS FOUND!")

    # Save results
    with open("results/pendant_mode_bound.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print("Results saved to results/pendant_mode_bound.json")


if __name__ == "__main__":
    main()
