#!/usr/bin/env python3
"""Verify Conjecture A for n=23: d_leaf <= 1 implies mode <= floor(n/3)+1 = 8.

Uses geng with res/mod partitioning across 8 workers for speed.
"""

import subprocess
import sys
import time
from multiprocessing import Pool

from indpoly import independence_poly

N = 23
BOUND = N // 3 + 1  # floor(23/3) + 1 = 8
NUM_WORKERS = 8
GENG = "/opt/homebrew/bin/geng"


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


def check_dleaf(n, adj):
    """Return True if all vertices have d_leaf <= 1."""
    leaves = set()
    for v in range(n):
        if len(adj[v]) == 1:
            leaves.add(v)
    for v in range(n):
        if v in leaves:
            continue
        d_leaf = sum(1 for u in adj[v] if u in leaves)
        if d_leaf >= 2:
            return False
    return True


def worker(res):
    """Process partition res out of NUM_WORKERS."""
    cmd = [GENG, str(N), f"{N-1}:{N-1}", "-c", "-q",
           f"{res}/{NUM_WORKERS}"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    total = 0
    dleaf_count = 0
    violations = []
    tightest_mode = 0  # track highest mode among d_leaf<=1 trees

    for raw_line in proc.stdout:
        line = raw_line.decode("ascii").strip()
        if not line:
            continue
        n, adj = parse_graph6(line)
        total += 1

        if not check_dleaf(n, adj):
            continue

        dleaf_count += 1
        poly = independence_poly(n, adj)
        mode = max(range(len(poly)), key=lambda k: poly[k])

        if mode > tightest_mode:
            tightest_mode = mode

        if mode > BOUND:
            violations.append((line, mode, poly))

        if total % 500_000 == 0:
            print(f"  worker {res}: {total} trees, "
                  f"{dleaf_count} d_leaf<=1, "
                  f"{len(violations)} violations so far",
                  flush=True)

    proc.wait()
    return total, dleaf_count, violations, tightest_mode


def main():
    print(f"Conjecture A verification for n={N}", flush=True)
    print(f"Bound: mode <= floor({N}/3)+1 = {BOUND}", flush=True)
    print(f"Workers: {NUM_WORKERS}", flush=True)
    print(flush=True)

    t0 = time.time()

    with Pool(NUM_WORKERS) as pool:
        results = pool.map(worker, range(NUM_WORKERS))

    elapsed = time.time() - t0

    grand_total = sum(r[0] for r in results)
    grand_dleaf = sum(r[1] for r in results)
    all_violations = []
    best_mode = 0
    for r in results:
        all_violations.extend(r[2])
        if r[3] > best_mode:
            best_mode = r[3]

    print(f"\n{'='*60}", flush=True)
    print(f"RESULTS for n={N}", flush=True)
    print(f"{'='*60}", flush=True)
    print(f"Total trees:       {grand_total:>12,}", flush=True)
    print(f"d_leaf <= 1 trees: {grand_dleaf:>12,}", flush=True)
    print(f"Violations:        {len(all_violations):>12,}", flush=True)
    print(f"Tightest mode:     {best_mode:>12} (bound = {BOUND})", flush=True)
    print(f"Elapsed:           {elapsed:>12.1f}s", flush=True)

    if all_violations:
        print(f"\nVIOLATIONS FOUND:", flush=True)
        for g6, mode, poly in all_violations[:20]:
            print(f"  graph6={g6}  mode={mode}  poly={poly}", flush=True)
    else:
        print(f"\nConjecture A HOLDS for n={N}: "
              f"all {grand_dleaf:,} d_leaf<=1 trees have mode <= {BOUND}",
              flush=True)


if __name__ == "__main__":
    main()
