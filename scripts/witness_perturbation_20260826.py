#!/usr/bin/env python3
"""Targeted mixed-arm follow-up to the n=34/n=36 subcubic failures.

The witness is a centre with three isomorphic 11-vertex arms. Two questions the
full-symmetry sweep cannot answer:

  Q1. Does the failure need all three arms equal, or does it survive perturbing
      one or two of them? (If it survives, symmetry is incidental and the real
      driver is some arm statistic.)
  Q2. Is the witness arm special, or does any arm of the same "shape class"
      work?

The default run extracts the 17 distinct rooted arms occurring in the 15 n=36
certificates, adds the n=34 witness arm, and tests every unordered triple from
that small seed pool. This directly covers the observed arm-size range 2..25;
it replaces the abandoned all-arms 9..12 sweep, which required roughly 900
million Python triples and was not a sensible default.

The older {W,W,X} and {W,X,Y} broad perturbations remain available behind
``--broad``. Exact integer arithmetic throughout; every polynomial is validated
(i_0 = 1, i_1 = n, nonnegative coefficients, no internal zeros).

Reports LC failures with break depth, and ALARMS on any non-unimodal sequence.
"""
import argparse
import json
import re
import sys
from functools import lru_cache
from pathlib import Path

import networkx as nx

sys.setrecursionlimit(100000)

# The n=34 witness arm, as re-derived by the symmetric sweep's positive control.
WITNESS_ARM = (((((),),),), ((((),), ((),)),))
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


def add(*ps):
    n = max(len(p) for p in ps)
    return [sum(p[i] if i < len(p) else 0 for p in ps) for i in range(n)]


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


@lru_cache(maxsize=None)
def arms_of_size(m):
    if m == 1:
        return ((),)
    out = set()
    for c in arms_of_size(m - 1):
        out.add((c,))
    for i in range(1, m - 1):
        j = m - 1 - i
        if i > j:
            continue
        for a in arms_of_size(i):
            for b in arms_of_size(j):
                out.add(tuple(sorted((a, b))))
    return tuple(sorted(out))


@lru_cache(maxsize=None)
def arm_polys(arm):
    F, G = [1], [0, 1]
    for child in arm:
        cf, cg = arm_polys(child)
        F = mul(F, add(list(cf), list(cg)))
        G = mul(G, list(cf))
    return (tuple(trim(F)), tuple(trim(G)))


@lru_cache(maxsize=None)
def arm_size(arm):
    return 1 + sum(arm_size(c) for c in arm)


def rooted_arm(graph, root, parent=-1):
    """Canonical tuple encoding compatible with ``arms_of_size``."""
    return tuple(sorted(
        rooted_arm(graph, child, root)
        for child in graph[root]
        if child != parent
    ))


def graph_from_arms(arms):
    graph = nx.Graph()
    graph.add_node(0)
    next_vertex = 1

    def add_arm(arm, parent):
        nonlocal next_vertex
        vertex = next_vertex
        next_vertex += 1
        graph.add_edge(parent, vertex)
        for child in arm:
            add_arm(child, vertex)

    for arm in arms:
        add_arm(arm, 0)
    return graph


def unrooted_graph_code(graph):
    """AHU-style canonical code for an unrooted tree."""
    centres = nx.center(graph)
    if len(centres) == 1:
        return ("vertex", rooted_arm(graph, centres[0]))
    assert len(centres) == 2
    left, right = centres
    sides = sorted((
        rooted_arm(graph, left, right),
        rooted_arm(graph, right, left),
    ))
    return ("edge", tuple(sides))


def unrooted_tree_code(arms):
    """Canonical code, used to avoid duplicate arm presentations."""
    return unrooted_graph_code(graph_from_arms(arms))


def census_seed_arms(summary_path):
    """Extract the distinct rooted arms from the n=36 FAIL certificates."""
    summary = json.loads(Path(summary_path).read_text())
    seeds = set()
    for line in summary["fail_lines"]:
        parents = [
            int(value)
            for value in re.search(r" par=([^ ]+)", line).group(1).split(",")
        ]
        graph = nx.Graph()
        graph.add_nodes_from(range(len(parents)))
        graph.add_edges_from((i, parents[i] - 1) for i in range(1, len(parents)))
        centres = nx.center(graph)
        assert len(centres) == 1
        centre = centres[0]
        without_centre = graph.copy()
        without_centre.remove_node(centre)
        for component in nx.connected_components(without_centre):
            vertices = set(component)
            root = next(v for v in graph[centre] if v in vertices)
            seeds.add(rooted_arm(without_centre.subgraph(vertices), root))
    assert len(seeds) == 17, f"expected 17 distinct n=36 arms, found {len(seeds)}"
    return tuple(sorted(seeds))


def census_tree_ids(summary_path):
    """Canonical IDs of the 15 exhaustive n=36 failure certificates."""
    summary = json.loads(Path(summary_path).read_text())
    identifiers = set()
    for line in summary["fail_lines"]:
        parents = [
            int(value)
            for value in re.search(r" par=([^ ]+)", line).group(1).split(",")
        ]
        graph = nx.Graph()
        graph.add_nodes_from(range(len(parents)))
        graph.add_edges_from((i, parents[i] - 1) for i in range(1, len(parents)))
        code = unrooted_graph_code(graph)
        identifiers.add(repr(code))
    assert len(identifiers) == 15
    return identifiers


def tree_poly(arms):
    """Centre joined to the given arms: prod(F_i+G_i) + x*prod(F_i)."""
    excl, incl = [1], [1]
    for a in arms:
        F, G = arm_polys(a)
        excl = mul(excl, add(list(F), list(G)))
        incl = mul(incl, list(F))
    return trim(add(excl, [0] + incl))


def validate(poly, n, label):
    assert poly[0] == 1, f"{label}: i_0 != 1"
    assert len(poly) < 2 or poly[1] == n, f"{label}: i_1={poly[1]} != n={n}"
    assert all(c >= 0 for c in poly), f"{label}: negative coefficient"
    nz = [i for i, c in enumerate(poly) if c]
    assert not nz or nz == list(range(nz[0], nz[-1] + 1)), \
        f"{label}: INTERNAL ZERO — construction is wrong"


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


def check(arms, results, alarms, tag):
    n = 1 + sum(arm_size(a) for a in arms)
    poly = tree_poly(arms)
    label = f"{tag} n={n} arms={[arm_size(a) for a in arms]}"
    validate(poly, n, label)
    if not is_unimodal(poly):
        alarms.append({"tag": tag, "n": n, "arms": [str(a) for a in arms],
                       "poly": poly})
        print("!" * 70 + f"\nALARM NON-UNIMODAL: {label}\n  {poly}\n" + "!" * 70,
              flush=True)
        return True
    br = lc_breaks(poly)
    if br:
        alpha = len(poly) - 1
        thr = -(-(2 * alpha - 1) // 3)
        results.append({"tag": tag, "n": n, "alpha": alpha, "breaks": br,
                        "dists": [k - thr for k in br],
                        "tree_id": repr(unrooted_tree_code(arms)),
                        "arms": [str(a) for a in arms], "poly": poly})
        return True
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary",
        default="results/subcubic_census_20260826/summary_n36.json",
    )
    parser.add_argument(
        "--output",
        default="results/witness_perturbation_20260826.json",
    )
    parser.add_argument(
        "--broad",
        action="store_true",
        help="also run the older WWX and WXY all-arm perturbations",
    )
    parser.add_argument("--cap-one", type=int, default=20)
    parser.add_argument("--cap-two", type=int, default=13)
    args = parser.parse_args()
    results, alarms = [], []

    # sanity: the witness itself must reproduce
    assert tree_poly([WITNESS_ARM] * 3) == WITNESS_POLY, \
        "witness arm does not reproduce the stored n=34 polynomial"
    print("positive control: witness arm reproduces the n=34 polynomial\n")

    seeds = tuple(sorted(set(census_seed_arms(args.summary)) | {WITNESS_ARM}))
    print("=== Witness-seeded mixed arms ===", flush=True)
    print(
        f"  {len(seeds)} distinct rooted arms; sizes "
        f"{sorted({arm_size(arm) for arm in seeds})}",
        flush=True,
    )
    seeded_hits = 0
    seeded_multisets = 0
    seeded_codes = set()
    for i, A in enumerate(seeds):
        for j in range(i, len(seeds)):
            B = seeds[j]
            for k in range(j, len(seeds)):
                C = seeds[k]
                seeded_multisets += 1
                code = unrooted_tree_code([A, B, C])
                if code in seeded_codes:
                    continue
                seeded_codes.add(code)
                before = len(results)
                check([A, B, C], results, alarms, "seeded_ABC")
                seeded_hits += len(results) - before
    print(
        f"  {seeded_multisets} arm multisets, {len(seeded_codes)} distinct trees, "
        f"{seeded_hits} LC failure(s)",
        flush=True,
    )
    expected_n36 = census_tree_ids(args.summary)
    seeded_n36 = {
        result["tree_id"]
        for result in results
        if result["tag"] == "seeded_ABC" and result["n"] == 36
    }
    assert seeded_n36 == expected_n36, "seeded sweep did not reproduce the n=36 set exactly"
    print("  positive control: reproduces the exact 15-tree n=36 failure set", flush=True)

    broad_tested = {"WWX": 0, "WXY": 0}
    if args.broad:
        print(
            f"\n=== Broad WWX: |X| <= {args.cap_one} ===",
            flush=True,
        )
        for m in range(1, args.cap_one + 1):
            hits = 0
            for X in arms_of_size(m):
                broad_tested["WWX"] += 1
                before = len(results)
                check([WITNESS_ARM, WITNESS_ARM, X], results, alarms, "WWX")
                hits += len(results) - before
            print(
                f"  |X|={m:2d} n={23+m:3d}: {len(arms_of_size(m)):>7} arms, "
                f"{hits} failure(s)",
                flush=True,
            )

        print(
            f"\n=== Broad WXY: |X|,|Y| <= {args.cap_two} ===",
            flush=True,
        )
        pool = [a for m in range(1, args.cap_two + 1) for a in arms_of_size(m)]
        hits = 0
        for i, X in enumerate(pool):
            for Y in pool[i:]:
                broad_tested["WXY"] += 1
                before = len(results)
                check([WITNESS_ARM, X, Y], results, alarms, "WXY")
                hits += len(results) - before
        print(
            f"  {broad_tested['WXY']} pairs tested, {hits} failure(s)",
            flush=True,
        )

    print(f"\n=== TOTALS ===")
    print(f"LC failures: {len(results)}")
    print(f"UNIMODALITY ALARMS: {len(alarms)}")
    by_n = {}
    for r in results:
        by_n.setdefault(r["n"], []).append(r)
    for n in sorted(by_n):
        ds = sorted({d for r in by_n[n] for d in r["dists"]})
        print(f"  n={n}: {len(by_n[n])} failure(s), break depths {ds}")
    if results:
        best = min(results, key=lambda r: min(r["dists"]))
        print(f"\nshallowest break depth found: dist={min(best['dists'])} "
              f"at n={best['n']} (alpha={best['alpha']}, breaks={best['breaks']})")

    with open(args.output, "w") as fh:
        json.dump(
            {
                "summary": args.summary,
                "seed_arms": len(seeds),
                "seeded_arm_multisets": seeded_multisets,
                "seeded_distinct_trees": len(seeded_codes),
                "reproduces_exact_n36_failure_set": True,
                "broad": args.broad,
                "broad_tested": broad_tested,
                "failures": results,
                "alarms": alarms,
            },
            fh,
            indent=1,
        )
    print(f"saved {args.output}")


if __name__ == "__main__":
    main()
