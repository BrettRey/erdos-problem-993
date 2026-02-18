#!/usr/bin/env python3
"""Check whether every vertex in a maximal IS below mode has ≥ 1 external private neighbor.

The Bollobás-Cockayne lemma (as cited in our notes) claims priv(u) ≥ 1 for all u in any
maximal IS. This is NOT true in general (K_{1,4} all-leaves IS is a counterexample).
But it might hold for IS below the mode.

If it fails for some IS below mode, then:
- PNP failure does NOT imply 1-Private (could have mixed priv = 0 / priv = 1)
- The 1-Private Mode Conjecture check is INCOMPLETE (only checked 1-Private, not "all priv ≤ 1")
- We need to verify the broader claim for the proof

Three checks:
1. For ALL maximal IS (any size): how many have some vertex with priv = 0?
2. For maximal IS BELOW MODE: how many have some vertex with priv = 0?
3. For maximal IS below mode with all priv ≤ 1 (PNP failures): are they all 1-Private?
"""

import subprocess
import sys
import time

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


def compute_priv(u, s, nbr):
    """Count external private neighbors of u in S."""
    count = 0
    for v in nbr[u]:
        if v not in s:
            if len(nbr[v] & s) == 1:
                count += 1
    return count


def main():
    print("EXTERNAL PRIVATE NEIGHBOR VERIFICATION")
    print("=" * 70)
    print()
    print("Checking whether priv(u) ≥ 1 for all u in maximal IS below mode.")
    print("If any vertex has priv = 0 in IS below mode, the 1-Private reduction has a gap.")
    print()

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_trees = len(lines)
        # Counters for maximal IS below mode
        n_mis_below_mode = 0
        n_mis_with_priv0 = 0     # IS where some vertex has priv = 0
        n_pnp_failures = 0       # IS where ALL vertices have priv ≤ 1
        n_pnp_fail_not_1priv = 0 # PNP failures that are NOT 1-Private (mixed 0/1)

        # Also track: among ALL maximal IS (any size)
        n_mis_total = 0
        n_mis_any_priv0 = 0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            nbr = [set(adj_data[v]) for v in range(tn)]
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            all_mis = find_all_maximal_is(tn, adj_data)

            for s in all_mis:
                k = len(s)
                privs = [compute_priv(u, s, nbr) for u in s]
                min_priv = min(privs)
                max_priv = max(privs)

                n_mis_total += 1
                if min_priv == 0:
                    n_mis_any_priv0 += 1

                if k < mode:
                    n_mis_below_mode += 1
                    if min_priv == 0:
                        n_mis_with_priv0 += 1
                        print(f"  PRIV=0 BELOW MODE: n={tn} k={k} mode={mode} "
                              f"privs={sorted(privs)} S={sorted(s)}")

                    if max_priv <= 1:
                        # PNP failure
                        n_pnp_failures += 1
                        has_zero = any(p == 0 for p in privs)
                        if has_zero:
                            n_pnp_fail_not_1priv += 1
                            print(f"  PNP FAIL NOT 1-PRIV: n={tn} k={k} mode={mode} "
                                  f"privs={sorted(privs)} S={sorted(s)}")

        elapsed = time.time() - t0
        print(f"n={n:2d}: {n_trees:6d} trees, "
              f"total_MIS={n_mis_total:7d}, any_priv0={n_mis_any_priv0:6d} "
              f"({100*n_mis_any_priv0/max(1,n_mis_total):.1f}%), "
              f"below_mode={n_mis_below_mode:6d}, "
              f"bm_priv0={n_mis_with_priv0}, "
              f"pnp_fail={n_pnp_failures}, "
              f"pnp_fail_not_1priv={n_pnp_fail_not_1priv} "
              f"({elapsed:.1f}s)")

    print()
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print()
    print("If bm_priv0 = 0: every maximal IS below mode has priv ≥ 1 for all vertices.")
    print("  → Bollobás-Cockayne (external) holds below mode; 1-Private reduction is sound.")
    print()
    print("If bm_priv0 > 0 but pnp_fail_not_1priv = 0: priv=0 exists below mode,")
    print("  but PNP still holds (some other vertex has priv ≥ 2). No proof gap.")
    print()
    print("If pnp_fail_not_1priv > 0: there exist PNP-failing IS below mode that are")
    print("  NOT 1-Private. The 1-Private check is INCOMPLETE for the proof.")


if __name__ == "__main__":
    main()
