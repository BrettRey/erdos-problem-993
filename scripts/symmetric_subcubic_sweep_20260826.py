#!/usr/bin/env python3
"""Exhaustive sweep of 3-fold symmetric subcubic trees (2026-08-26).

Motivation: the first subcubic LC failure found by this census (n=34) is a
centre with three pairwise-isomorphic 11-vertex arms. That structure supplies a
cheap targeted slice; it does not imply that symmetry is generically enriched
among failures. The exhaustive ladder reaches n=36 with 2.6e10 trees, whereas
the 3-fold symmetric sweep through arm size 22 examines about 5.4 million arms
and reaches n=67.

Construction. A 3-fold symmetric subcubic tree is a centre vertex c joined to
three copies of one rooted arm. The centre then has degree 3; each arm root is
joined to c so it may have at most 2 children, and every other arm vertex may
have at most 2 children, so arms are exactly the unordered rooted trees with at
most two children per node (Wedderburn-Etherington).

Independence polynomial, exactly and fast: if the arm has F = polynomial with
arm-root EXCLUDED and Gp = arm-root INCLUDED, then

    I(T) = (F + Gp)^3 + x * F^3

since the centre is either excluded (each arm free) or included (each arm root
excluded). All arithmetic is exact Python integers.

Also sweeps the 2-fold symmetric case (an edge joining two isomorphic rooted
arms, arms allowed up to 2 children at the root => centre degree <= 3), which is
the other natural symmetric shape.

Reports every log-concavity failure with its break positions and depth
dist = k - ceil((2*alpha-1)/3), and ALARMS loudly on any non-unimodal sequence
(that would be a counterexample to Erdos Problem 993).

Positive control: the n=34 witness must be re-found by the 3-fold sweep at arm
size 11. The run aborts if it is not.
"""
import json
import sys
from functools import lru_cache

sys.setrecursionlimit(100000)

WITNESS_POLY = [1, 34, 528, 4967, 31651, 144723, 490680, 1256998, 2456187,
                3669088, 4172235, 3571401, 2256341, 1018880, 311700, 58485,
                5325, 72, 1]


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in enumerate(b):
            if y:
                out[i + j] += x * y
    return out


def add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(n)]


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


# ---------------------------------------------------------------- arm shapes
# An arm is an unordered rooted tree with <= 2 children per node, encoded
# canonically as a sorted tuple of child encodings.

@lru_cache(maxsize=None)
def arms_of_size(m):
    """All canonical unordered rooted trees on m nodes with <= 2 children."""
    if m == 1:
        return (((),),)[0:1] if False else ((),)  # the leaf: empty child tuple
    out = set()
    # one child subtree of size m-1
    for c in arms_of_size(m - 1):
        out.add((c,))
    # two child subtrees of sizes i and m-1-i, unordered
    for i in range(1, m - 1):
        j = m - 1 - i
        if i > j:
            continue
        for a in arms_of_size(i):
            for b in arms_of_size(j):
                pair = tuple(sorted((a, b)))
                out.add(pair)
    return tuple(sorted(out))


@lru_cache(maxsize=None)
def arm_polys(arm):
    """(F, G) for a rooted arm: F = root excluded, G = root included."""
    F, G = [1], [0, 1]           # leaf: F = 1, G = x
    for child in arm:
        cf, cg = arm_polys(child)
        F = mul(F, add(cf, cg))  # root excluded: child free
        G = mul(G, cf)           # root included: child root excluded
    return (tuple(trim(list(F))), tuple(trim(list(G))))


def arm_size(arm):
    return 1 + sum(arm_size(c) for c in arm)


def is_unimodal(seq):
    rising = True
    for i in range(1, len(seq)):
        if rising and seq[i] < seq[i - 1]:
            rising = False
        elif not rising and seq[i] > seq[i - 1]:
            return False
    return True


def lc_breaks(seq):
    return [k for k in range(1, len(seq) - 1)
            if seq[k] * seq[k] < seq[k - 1] * seq[k + 1]]


def validate(poly, n, label):
    """Every independence polynomial of a graph on n vertices satisfies these.

    Added 2026-08-26 after a sign/shift bug in the 2-fold construction produced
    sequences ending `..., 174, 0, 1` and three spurious NON-UNIMODAL alarms.
    An independent set of size k+1 always contains one of size k, so internal
    zeros are impossible; i_0 = 1 and i_1 = n are free checks on the whole
    construction. Any of these firing means the CODE is wrong, not the tree.
    """
    if poly[0] != 1:
        raise AssertionError(f"{label}: i_0 = {poly[0]} != 1")
    if len(poly) > 1 and poly[1] != n:
        raise AssertionError(f"{label}: i_1 = {poly[1]} != n = {n}")
    if any(c < 0 for c in poly):
        raise AssertionError(f"{label}: negative coefficient in {poly}")
    nz = [i for i, c in enumerate(poly) if c != 0]
    if nz and nz != list(range(nz[0], nz[-1] + 1)):
        raise AssertionError(
            f"{label}: INTERNAL ZERO in {poly} — impossible for an "
            f"independence polynomial; the construction is wrong")


def analyse(poly, label, results, alarms, n=None):
    poly = trim(list(poly))
    if n is not None:
        validate(poly, n, label)
    alpha = len(poly) - 1
    if not is_unimodal(poly):
        alarms.append({"label": label, "poly": poly, "alpha": alpha})
        print("!" * 70)
        print(f"ALARM: NON-UNIMODAL — {label}")
        print(f"  poly = {poly}")
        print("!" * 70, flush=True)
        return
    br = lc_breaks(poly)
    if br:
        thr = -(-(2 * alpha - 1) // 3)
        results.append({
            "label": label, "n": None, "alpha": alpha, "poly": poly,
            "breaks": br, "dists": [k - thr for k in br],
        })


def main():
    max_arm = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    results, alarms = [], []
    found_witness = False
    counts = {}

    print(f"=== 3-fold symmetric subcubic trees, arm size 1..{max_arm} "
          f"(n = 3m+1 up to {3*max_arm+1}) ===", flush=True)
    for m in range(1, max_arm + 1):
        arms = arms_of_size(m)
        counts[m] = len(arms)
        n = 3 * m + 1
        nfail = 0
        for arm in arms:
            F, G = arm_polys(arm)
            s = add(list(F), list(G))
            poly = add(mul(mul(s, s), s), [0] + mul(mul(list(F), list(F)), list(F)))
            poly = trim(poly)
            if poly == WITNESS_POLY:
                found_witness = True
            before = len(results)
            analyse(poly, f"sym3 m={m} n={n} arm={arm}", results, alarms, n=n)
            nfail += len(results) - before
        for r in results:
            if r["n"] is None:
                r["n"] = n
        print(f"  m={m:2d} n={n:3d}: {len(arms):>7} arms, {nfail} LC failure(s)",
              flush=True)

    print(f"\n=== 2-fold symmetric (edge between two isomorphic arms) ===",
          flush=True)
    two_fold = []
    for m in range(1, max_arm + 1):
        arms = arms_of_size(m)
        n = 2 * m
        nfail = 0
        for arm in arms:
            F, G = arm_polys(arm)
            s = add(list(F), list(G))
            # Edge between the two arm roots: all pairs of independent sets,
            # minus those including BOTH roots. G already carries the root's x,
            # so the subtracted term is G*G with NO extra shift. (The shift was
            # the 2026-08-26 bug; see validate().)
            poly = add(mul(s, s), [-c for c in mul(list(G), list(G))])
            poly = trim(poly)
            before = len(two_fold)
            analyse(poly, f"sym2 m={m} n={n} arm={arm}", two_fold, alarms, n=n)
            nfail += len(two_fold) - before
        for r in two_fold:
            if r["n"] is None:
                r["n"] = n
        if nfail:
            print(f"  m={m:2d} n={n:3d}: {len(arms):>7} arms, {nfail} LC failure(s)",
                  flush=True)

    print()
    if not found_witness:
        print("POSITIVE CONTROL FAILED: the n=34 witness was not re-found.")
        raise SystemExit(1)
    print("positive control PASSED: n=34 witness re-found by the 3-fold sweep")
    print(f"3-fold LC failures: {len(results)}")
    print(f"2-fold LC failures: {len(two_fold)}")
    print(f"UNIMODALITY ALARMS: {len(alarms)}")

    if results:
        print("\n--- 3-fold failures ---")
        for r in sorted(results, key=lambda r: (r["n"], r["breaks"])):
            print(f"  n={r['n']} alpha={r['alpha']} breaks={r['breaks']} "
                  f"dist={r['dists']}")
            print(f"     arm={r['label'].split('arm=')[1]}")

    out = {"max_arm": max_arm, "arm_counts": counts,
           "sym3_failures": results, "sym2_failures": two_fold,
           "alarms": alarms}
    with open("results/symmetric_subcubic_sweep_20260826.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nsaved results/symmetric_subcubic_sweep_20260826.json")


if __name__ == "__main__":
    main()
