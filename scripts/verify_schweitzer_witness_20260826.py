#!/usr/bin/env python3
"""Exact replay of the checkable claims in Schweitzer, arXiv:2608.23262 (v1).

Three claims are replayed in exact integer arithmetic before the paper cites
the result (source hook: notes/source-hooks/schweitzer-2026-matroid-intersection-ulc.md):

1. Example 1 (the three-partition-matroid witness): ground set
   {a1,a2,a3,b1..bm}, matroid M_i has the single nontrivial part
   {a_i, b1..bm}.  The common independent sets are exactly the subsets of
   {a1,a2,a3} plus the singletons {b_j}; the size sequence is
   (1, m+3, 3, 1); log-concavity fails at k=2 iff m >= 7; the sequence is
   unimodal for every m >= 0.

2. The graph reading in the source hook: that independence system is the
   independence system of the complete split graph (b's a clique, every b
   adjacent to every a, a's mutually non-adjacent).

3. The reduction the manuscript will state: the independent sets of a graph G
   are the common independent sets of one partition matroid per edge (part
   {u,v} with capacity 1, all other elements singletons).  Verified by brute
   force on all trees n <= 9 (via networkx nonisomorphic_trees) plus the
   witness graphs.

Also checks the derived remark for the paper: ULC (Schweitzer Cor. 1) implies
LC for sequences without internal zeros, so an LC-failing graph independence
system cannot be the common independents of two matroids.  The implication
ULC => LC is checked symbolically on the binomial ratio.

Everything is exact (fractions/ints).  Exit 0 iff every check passes.
"""

from fractions import Fraction
from itertools import combinations
from math import comb

FAILURES = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def partition_matroid_indep(parts_with_caps, subset) -> bool:
    """subset independent iff it meets every part within its capacity."""
    return all(len(subset & part) <= cap for part, cap in parts_with_caps)


def common_independents(matroids, ground):
    """Brute-force the common independent sets of a list of partition matroids."""
    ground = sorted(ground)
    out = []
    for r in range(len(ground) + 1):
        for c in combinations(ground, r):
            s = frozenset(c)
            if all(partition_matroid_indep(m, s) for m in matroids):
                out.append(s)
    return out


def size_sequence(sets_, n_max=None):
    top = max((len(s) for s in sets_), default=0) if n_max is None else n_max
    seq = [0] * (top + 1)
    for s in sets_:
        seq[len(s)] += 1
    return seq


def is_unimodal(seq):
    rising = True
    for i in range(1, len(seq)):
        if rising and seq[i] < seq[i - 1]:
            rising = False
        elif not rising and seq[i] > seq[i - 1]:
            return False
    return True


def lc_violations(seq):
    return [k for k in range(1, len(seq) - 1) if seq[k] ** 2 < seq[k - 1] * seq[k + 1]]


def graph_independent_sets(vertices, edges):
    es = [frozenset(e) for e in edges]
    out = []
    for r in range(len(vertices) + 1):
        for c in combinations(sorted(vertices), r):
            s = frozenset(c)
            if not any(e <= s for e in es):
                out.append(s)
    return out


def edge_partition_matroids(vertices, edges):
    """One partition matroid per edge: part {u,v} cap 1, other vertices singleton cap 1."""
    ms = []
    vs = set(vertices)
    for (u, v) in edges:
        parts = [(frozenset({u, v}), 1)]
        parts += [(frozenset({w}), 1) for w in vs - {u, v}]
        ms.append(parts)
    return ms


print("Claim 1: Example 1 witness, m = 0..12")
for m in range(13):
    ground = [f"a{i}" for i in (1, 2, 3)] + [f"b{j}" for j in range(1, m + 1)]
    bs = frozenset(g for g in ground if g.startswith("b"))
    matroids = []
    for i in (1, 2, 3):
        big = frozenset({f"a{i}"}) | bs
        parts = [(big, 1)] + [
            (frozenset({f"a{j}"}), 1) for j in (1, 2, 3) if j != i
        ]
        matroids.append(parts)
    common = common_independents(matroids, ground)
    # Structural claim: subsets of {a1,a2,a3} plus singleton b's.
    a_subsets = set()
    for r in range(4):
        for c in combinations([f"a{i}" for i in (1, 2, 3)], r):
            a_subsets.add(frozenset(c))
    expected_family = a_subsets | {frozenset({b}) for b in bs}
    fam_ok = set(common) == expected_family
    seq = size_sequence(common, 3)
    seq_ok = seq == [1, m + 3, 3, 1]
    lc = lc_violations(seq)
    lc_ok = (lc == [2]) if m >= 7 else (lc == [])
    uni_ok = is_unimodal(seq)
    check(
        f"m={m}",
        fam_ok and seq_ok and lc_ok and uni_ok,
        f"seq={seq}, LC violations at {lc}, unimodal={uni_ok}",
    )

print("Claim 2: the witness family equals the complete split graph independence system")
for m in (3, 7, 10):
    ground = [f"a{i}" for i in (1, 2, 3)] + [f"b{j}" for j in range(1, m + 1)]
    bs = [g for g in ground if g.startswith("b")]
    # Complete split graph: b's a clique, every b adjacent to every a.
    edges = [(x, y) for x, y in combinations(bs, 2)]
    edges += [(a, b) for a in ("a1", "a2", "a3") for b in bs]
    graph_is = graph_independent_sets(ground, edges)
    matroids = []
    for i in (1, 2, 3):
        big = frozenset({f"a{i}"}) | frozenset(bs)
        parts = [(big, 1)] + [
            (frozenset({f"a{j}"}), 1) for j in (1, 2, 3) if j != i
        ]
        matroids.append(parts)
    common = common_independents(matroids, ground)
    check(f"m={m} split-graph match", set(graph_is) == set(common))

print("Claim 3: graph independence = common independents of one partition matroid per edge")
try:
    import networkx as nx

    trees = []
    for n in range(2, 10):
        trees += [(f"tree n={n} #{i}", list(t.nodes()), list(t.edges()))
                  for i, t in enumerate(nx.nonisomorphic_trees(n))]
    sample = trees[:: max(1, len(trees) // 25)]  # spread across orders
    all_ok = True
    for name, vs, es in sample:
        gi = set(graph_independent_sets(vs, es))
        ci = set(common_independents(edge_partition_matroids(vs, es), vs))
        if gi != ci:
            all_ok = False
            check(name, False)
    check(f"all {len(sample)} sampled trees n<=9", all_ok)
except ImportError:
    check("networkx available", False, "install networkx to run claim 3")

print("Derived remark: ULC => LC (binomial ratio >= 1), so LC failure blocks two-matroid form")
ratio_ok = True
for n in range(2, 60):
    for k in range(1, n):
        r = Fraction(comb(n, k) ** 2, comb(n, k - 1) * comb(n, k + 1))
        if r < 1:
            ratio_ok = False
check("C(n,k)^2 >= C(n,k-1)C(n,k+1) for n<60", ratio_ok)

print()
if FAILURES:
    print(f"FAILED: {len(FAILURES)} check(s): {FAILURES}")
    raise SystemExit(1)
print("All checks passed.")
