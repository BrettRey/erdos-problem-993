#!/usr/bin/env python3
"""Analyze whether PNP can be proved recursively through subtrees.

For each maximal IS S below mode with max_priv >= 2:
  - Find the vertex u with max_priv
  - Is u the hub? Or a subtree vertex?
  - For the subtree containing u: does pigeonhole (n_w > 3*|S_w|) explain priv >= 2?

If the subtree pigeonhole always explains priv >= 2, then a recursive
argument works: apply PNB within each subtree, and show that at least
one subtree is "undersaturated" (n_w > 3*gamma_w).
"""

import subprocess
import sys
import time
from collections import deque

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
    count = 0
    for v in nbr[u]:
        if v not in s:
            if len(nbr[v] & s) == 1:
                count += 1
    return count


def main():
    print("RECURSIVE PNP: WHERE DOES priv >= 2 COME FROM?")
    print("=" * 70)
    print()

    # For each high-mode tree, for each IS below mode:
    # - Which vertex provides max_priv?
    # - Is it the hub or a subtree vertex?
    # - Does pigeonhole (n_local > 3*k_local) explain it?

    for n in range(8, 17):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        threshold = n // 3 + 2
        third_bound = (n + 1 + 2) // 3  # ceil((n+1)/3)

        n_high = 0
        n_mis_below = 0
        n_hub_provides = 0  # max_priv vertex is the hub
        n_subtree_provides = 0  # max_priv vertex is NOT the hub
        n_pigeonhole_global = 0  # global pigeonhole (n >= 3k) explains
        n_pigeonhole_local = 0  # local pigeonhole (in some subtree) explains

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])

            if mode < threshold:
                continue

            n_high += 1
            nbr = [set(adj_data[v]) for v in range(tn)]

            # Find hub
            hub = max(range(tn), key=lambda v: len(adj_data[v]))

            all_mis = find_all_maximal_is(tn, adj_data)
            for s in all_mis:
                k = len(s)
                if k >= mode:
                    continue

                n_mis_below += 1

                # Find vertex with max priv
                privs = {u: compute_priv(u, s, nbr) for u in s}
                max_priv_vertex = max(privs, key=privs.get)
                max_priv_val = privs[max_priv_vertex]

                if max_priv_vertex == hub:
                    n_hub_provides += 1
                else:
                    n_subtree_provides += 1

                # Global pigeonhole: n >= 3k?
                if tn >= 3 * k:
                    n_pigeonhole_global += 1

                # Local pigeonhole: does pigeonhole work in max_priv_vertex's
                # "local neighborhood"?
                # Compute the "effective subtree" around max_priv_vertex
                # For simplicity: check if any vertex in S has
                # deg(u) > 2*(k-1) (degree-based bound: priv(u) >= deg(u) - (k-1))
                # This gives priv >= 2 when deg(u) >= k+1
                max_deg_in_s = max(len(nbr[u]) for u in s)
                if max_deg_in_s >= k + 1:
                    n_pigeonhole_local += 1

        elapsed = time.time() - t0
        if n_mis_below > 0:
            print(f"n={n:2d}: {n_high:4d} hm trees, {n_mis_below:6d} IS below mode, "
                  f"hub_provides={n_hub_provides} ({100*n_hub_provides/n_mis_below:.0f}%), "
                  f"deg_bound={n_pigeonhole_local} ({100*n_pigeonhole_local/n_mis_below:.0f}%), "
                  f"global_ph={n_pigeonhole_global} ({100*n_pigeonhole_global/n_mis_below:.0f}%) "
                  f"({elapsed:.1f}s)")
        else:
            print(f"n={n:2d}: {n_high:4d} hm trees, no IS below mode ({elapsed:.1f}s)")

    print()
    print("=" * 70)
    print("DEGREE BOUND: priv(u) >= deg(u) - (k-1). If max_deg_in_S >= k+1,")
    print("then some vertex has priv >= 2. 'deg_bound' shows % covered.")
    print("=" * 70)


if __name__ == "__main__":
    main()
