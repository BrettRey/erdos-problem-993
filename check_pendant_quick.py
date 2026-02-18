#!/usr/bin/env python3
"""Quick check: pendant tree mode bound + R-independence, n <= 16.

For each tree, find maximal IS below mode with all priv <= 1.
Check: (1) mode <= k, (2) is R independent?
"""

import subprocess, sys, time
sys.path.insert(0, ".")
from indpoly import independence_poly


def parse_graph6(s):
    s = s.strip()
    if not s:
        return 0, []
    idx = 0
    if ord(s[0]) - 63 < 63:
        n = ord(s[0]) - 63; idx = 1
    else:
        idx = 1; n = 0
        for _ in range(3):
            n = n * 64 + (ord(s[idx]) - 63); idx += 1
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
                adj[i].append(j); adj[j].append(i)
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
    import sys
    print("PENDANT TREE: MODE BOUND + R-INDEPENDENCE", flush=True)
    print("=" * 55, flush=True)

    totals = {"pendant_is": 0, "mode_ok": 0, "mode_fail": 0,
              "R_indep": 0, "R_not_indep": 0}

    for n in range(5, 19):
        t0 = time.time()
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = [l for l in proc.stdout.strip().split("\n") if l]

        n_pis, n_ok, n_fail, n_ri, n_rni, max_er, max_ratio = 0, 0, 0, 0, 0, 0, 0.0

        for line in lines:
            tn, adj_data = parse_graph6(line)
            poly = independence_poly(tn, adj_data)
            mode = max(range(len(poly)), key=lambda i: poly[i])
            nbr = [set(adj_data[v]) for v in range(tn)]
            all_mis = find_all_maximal_is(tn, adj_data)

            for s in all_mis:
                k = len(s)
                if k >= mode:
                    continue
                # Check all priv <= 1
                all_leq1 = True
                privates = set()
                for u in s:
                    priv_u = 0
                    for v in adj_data[u]:
                        if v not in s and len(nbr[v] & s) == 1:
                            priv_u += 1
                            privates.add(v)
                    if priv_u > 1:
                        all_leq1 = False
                        break
                if not all_leq1:
                    continue

                n_pis += 1

                # Mode bound
                if mode <= k:
                    n_ok += 1
                else:
                    n_fail += 1
                    print(f"  MODE VIOLATION n={tn} k={k} mode={mode}", flush=True)

                # R-independence
                R = set(range(tn)) - s - privates
                eR = sum(1 for r in R for r2 in adj_data[r] if r2 in R and r2 > r)
                max_er = max(max_er, eR)
                if eR == 0:
                    n_ri += 1
                else:
                    n_rni += 1

                # Ratio
                if k + 1 < len(poly) and poly[k] > 0:
                    ratio = poly[k+1] / poly[k]
                    max_ratio = max(max_ratio, ratio)

        elapsed = time.time() - t0
        totals["pendant_is"] += n_pis
        totals["mode_ok"] += n_ok; totals["mode_fail"] += n_fail
        totals["R_indep"] += n_ri; totals["R_not_indep"] += n_rni

        if n_pis > 0:
            print(f"n={n:2d}: {n_pis:5d} pIS, mode_ok={n_ok}, mode_fail={n_fail}, "
                  f"R_ind={n_ri}, R_not={n_rni}, max_eR={max_er}, "
                  f"max_ratio={max_ratio:.4f} ({elapsed:.1f}s)", flush=True)
        else:
            print(f"n={n:2d}: 0 pendant IS below mode ({elapsed:.1f}s)", flush=True)

    print(flush=True)
    print(f"TOTALS: {totals}", flush=True)


if __name__ == "__main__":
    main()
