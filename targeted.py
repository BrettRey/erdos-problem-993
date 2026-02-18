#!/usr/bin/env python3
"""Targeted search for unimodality failures in specific tree families.

Tests tree families known to violate log-concavity (Galvin 2025,
Kadrawi & Levit 2023, Ramos & Sun 2025) at much larger vertex counts
than exhaustive enumeration allows.

Families implemented:
  1. Spherically symmetric trees T_{m,t,1} (Galvin 2025)
  2. Spherically symmetric trees T_{m,t,d} (generalized depth)
  3. Caterpillars with varying pendant patterns
  4. "Broom" trees (path + star)
  5. Double stars and spiders
  6. Random trees biased toward degree-2 vertices adjacent to leaves
     (the structural motif Ramos & Sun 2025 found characteristic of
     log-concavity failures)
"""

import argparse
import json
import os
import random
import time

from indpoly import (
    independence_poly,
    is_log_concave,
    is_unimodal,
    log_concavity_ratio,
    near_miss_ratio,
)


# ---------------------------------------------------------------------------
# Tree family generators
# ---------------------------------------------------------------------------

def make_T_m_t_1(m: int, t: int) -> tuple[int, list[list[int]]]:
    """Galvin's spherically symmetric tree T_{m,t,1}.

    Structure (rooted):
      - root (vertex 0) has m children w_1..w_m
      - each w_i has t children x_{i,1}..x_{i,t}
      - each x_{i,j} has 1 child y_{i,j}

    Total vertices: 1 + m + mt + mt = 1 + m + 2mt = 1 + m(1 + 2t)
    Independence number: m(1 + t)
    Log-concavity failure at position mt + 2 (for t <= m <= 2^{t/16}).
    """
    n = 1 + m * (1 + 2 * t)
    adj: list[list[int]] = [[] for _ in range(n)]

    # root = 0
    # w_i = 1..m
    # x_{i,j} = m+1 .. m+mt  (in blocks of t per w_i)
    # y_{i,j} = m+mt+1 .. m+2mt

    def _add_edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    for i in range(m):
        w = 1 + i
        _add_edge(0, w)
        for j in range(t):
            x = 1 + m + i * t + j
            _add_edge(w, x)
            y = 1 + m + m * t + i * t + j
            _add_edge(x, y)

    return n, adj


def make_T_m_t_d(m: int, t: int, d: int) -> tuple[int, list[list[int]]]:
    """Generalized spherically symmetric tree with depth d.

    Structure: root has m children, each has t children, and so on
    for d levels. At the deepest level, each vertex has 1 child (leaf).

    This generalizes T_{m,t,1} to deeper trees.
    """
    adj: list[list[int]] = [[]]
    # BFS-style construction
    current_layer = [0]
    for level in range(d):
        next_layer = []
        if level == 0:
            branching_factor = m
        else:
            branching_factor = t
        for v in current_layer:
            for _ in range(branching_factor):
                u = len(adj)
                adj.append([])
                adj[v].append(u)
                adj[u].append(v)
                next_layer.append(u)
        current_layer = next_layer

    # Final layer: each gets 1 leaf child
    for v in current_layer:
        u = len(adj)
        adj.append([])
        adj[v].append(u)
        adj[u].append(v)

    return len(adj), adj


def make_caterpillar(spine_len: int, pendants: list[int]) -> tuple[int, list[list[int]]]:
    """Caterpillar tree: path of spine_len vertices, each with pendants[i] leaves.

    pendants[i] = number of pendant leaves at spine vertex i.
    len(pendants) must equal spine_len.
    """
    assert len(pendants) == spine_len
    n = spine_len + sum(pendants)
    adj: list[list[int]] = [[] for _ in range(n)]

    def _add_edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    # Spine: 0, 1, ..., spine_len-1
    for i in range(spine_len - 1):
        _add_edge(i, i + 1)

    # Pendants
    next_v = spine_len
    for i in range(spine_len):
        for _ in range(pendants[i]):
            _add_edge(i, next_v)
            next_v += 1

    return n, adj


def make_broom(path_len: int, star_size: int) -> tuple[int, list[list[int]]]:
    """Broom tree: path of path_len vertices with star_size leaves at one end."""
    n = path_len + star_size
    adj: list[list[int]] = [[] for _ in range(n)]

    def _add_edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    for i in range(path_len - 1):
        _add_edge(i, i + 1)
    for i in range(star_size):
        _add_edge(path_len - 1, path_len + i)

    return n, adj


def make_spider(legs: list[int]) -> tuple[int, list[list[int]]]:
    """Spider tree: central vertex with paths of given lengths.

    legs[i] = length of leg i (number of edges, so leg has legs[i] vertices
    beyond the centre).
    """
    n = 1 + sum(legs)
    adj: list[list[int]] = [[] for _ in range(n)]

    def _add_edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    next_v = 1
    for leg_len in legs:
        prev = 0
        for _ in range(leg_len):
            _add_edge(prev, next_v)
            prev = next_v
            next_v += 1

    return n, adj


def make_double_star(a: int, b: int) -> tuple[int, list[list[int]]]:
    """Double star: two adjacent vertices with a and b pendant leaves."""
    n = 2 + a + b
    adj: list[list[int]] = [[] for _ in range(n)]

    def _add_edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    _add_edge(0, 1)
    for i in range(a):
        _add_edge(0, 2 + i)
    for i in range(b):
        _add_edge(1, 2 + a + i)

    return n, adj


def make_ramos_sun_style(n_target: int, rng: random.Random) -> tuple[int, list[list[int]]]:
    """Generate a random tree biased toward the structural motif Ramos & Sun
    found characteristic of log-concavity failures: many degree-2 vertices
    adjacent to leaves, with independence number near n/2 + 1.

    Uses Prüfer code generation with bias toward creating the right structure.
    """
    n = n_target
    # Strategy: build a tree with ~n/2 leaves by creating degree-2 "bridges"
    # that connect to leaves, similar to subdivided stars.

    # Start with a small core (random tree on ~n/4 vertices)
    core_size = max(3, n // 4)
    adj: list[list[int]] = [[] for _ in range(n)]

    def _add_edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    # Build core via random Prüfer code
    if core_size <= 2:
        if core_size == 2:
            _add_edge(0, 1)
        next_v = core_size
    else:
        prufer = [rng.randrange(core_size) for _ in range(core_size - 2)]
        degree = [1] * core_size
        for v in prufer:
            degree[v] += 1
        for v in prufer:
            for u in range(core_size):
                if degree[u] == 1:
                    _add_edge(u, v)
                    degree[u] -= 1
                    degree[v] -= 1
                    break
        # Connect last two degree-1 vertices
        last = [u for u in range(core_size) if degree[u] == 1]
        if len(last) == 2:
            _add_edge(last[0], last[1])
        next_v = core_size

    # Attach degree-2 bridges + leaves to core vertices
    core_vertices = list(range(core_size))
    while next_v < n - 1:
        # Pick a random core vertex, attach bridge + leaf
        v = rng.choice(core_vertices)
        bridge = next_v
        _add_edge(v, bridge)
        next_v += 1
        if next_v < n:
            leaf = next_v
            _add_edge(bridge, leaf)
            next_v += 1

    # If one vertex left, attach to a random vertex
    if next_v < n:
        v = rng.randrange(next_v)
        _add_edge(v, next_v)
        next_v += 1

    return n, adj


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def analyze_tree(n: int, adj: list[list[int]], label: str) -> dict:
    """Compute all metrics for a single tree."""
    poly = independence_poly(n, adj)
    uni = is_unimodal(poly)
    lc = is_log_concave(poly)
    lc_r, lc_pos = log_concavity_ratio(poly)
    nm_r, nm_pos = near_miss_ratio(poly)

    return {
        "label": label,
        "n": n,
        "alpha": len(poly) - 1,
        "unimodal": uni,
        "log_concave": lc,
        "lc_ratio": round(lc_r, 12),
        "lc_pos": lc_pos,
        "nm_ratio": round(nm_r, 12),
        "nm_pos": nm_pos,
        "poly": poly,
    }


def _strip_poly(result: dict) -> dict:
    """Drop polynomial payload unless needed for a counterexample."""
    if result.get("unimodal", True):
        result.pop("poly", None)
    return result


def _summarize_family(results: list[dict]) -> dict:
    """Summarize a family run for reproducibility outputs."""
    total = len(results)
    lc_fail = sum(1 for r in results if not r["log_concave"])
    non_uni = sum(1 for r in results if not r["unimodal"])
    if results:
        best = max(results, key=lambda r: r["nm_ratio"])
        best_nm = best["nm_ratio"]
        best_label = best["label"]
        best_n = best["n"]
    else:
        best_nm = None
        best_label = None
        best_n = None
    return {
        "count": total,
        "lc_failures": lc_fail,
        "non_unimodal": non_uni,
        "best_nm": best_nm,
        "best_label": best_label,
        "best_n": best_n,
    }


def run_galvin_family(max_n: int) -> list[dict]:
    """Scan T_{m,t,1} for all feasible (m,t) up to max_n vertices."""
    results = []
    for t in range(1, 200):
        if 1 + t * (1 + 2 * t) > max_n:  # smallest m=t already too big
            break
        for m in range(t, 10000):
            n = 1 + m * (1 + 2 * t)
            if n > max_n:
                break
            n_tree, adj = make_T_m_t_1(m, t)
            result = analyze_tree(n_tree, adj, f"T_{{{m},{t},1}}")
            _strip_poly(result)
            results.append(result)
            if not result["log_concave"]:
                flag = "*** NON-UNIMODAL ***" if not result["unimodal"] else "LC-fail"
                print(
                    f"  {result['label']:>20s}  n={n:>5}  alpha={result['alpha']:>4}  "
                    f"lc_ratio={result['lc_ratio']:.8f}  nm={result['nm_ratio']:.8f}  "
                    f"{flag}",
                    flush=True,
                )
    return results


def _sst_size(m: int, t: int, d: int) -> int:
    """Compute vertex count of T_{m,t,d} without building it.

    Layer 0: 1 vertex (root)
    Layer 1: m vertices
    Layer 2..d: each vertex in prev layer has t children
    Layer d+1: each vertex in layer d has 1 child (leaf)
    """
    total = 1 + m  # root + layer 1
    layer_size = m
    for _ in range(d - 1):
        layer_size *= t
        total += layer_size
    total += layer_size  # leaf layer (1 child each)
    return total


def run_generalized_sst(max_n: int) -> list[dict]:
    """Scan T_{m,t,d} for various depths."""
    results = []
    for d in range(2, 6):
        for t in range(2, 20):
            for m in range(2, 200):
                n_est = _sst_size(m, t, d)
                if n_est > max_n:
                    break
                n_tree, adj = make_T_m_t_d(m, t, d)
                result = analyze_tree(n_tree, adj, f"T_{{{m},{t},{d}}}")
                _strip_poly(result)
                results.append(result)
                if not result["log_concave"]:
                    flag = "*** NON-UNIMODAL ***" if not result["unimodal"] else "LC-fail"
                    print(
                        f"  {result['label']:>20s}  n={n_tree:>5}  alpha={result['alpha']:>4}  "
                        f"lc_ratio={result['lc_ratio']:.8f}  nm={result['nm_ratio']:.8f}  "
                        f"{flag}",
                        flush=True,
                    )
    return results


def run_caterpillars(max_n: int) -> list[dict]:
    """Scan caterpillars with various pendant patterns."""
    results = []
    rng = random.Random(42)

    # Systematic: constant pendant count k, varying spine length
    for k in range(1, 10):
        for spine in range(4, max_n // (k + 1) + 1):
            n = spine + spine * k
            if n > max_n:
                break
            pendants = [k] * spine
            n_tree, adj = make_caterpillar(spine, pendants)
            result = analyze_tree(n_tree, adj, f"cat({spine},[{k}]*{spine})")
            _strip_poly(result)
            results.append(result)

    # Alternating pendant patterns
    for k1 in range(0, 6):
        for k2 in range(k1 + 1, 8):
            for spine in range(4, 200):
                pendants = [k1 if i % 2 == 0 else k2 for i in range(spine)]
                n = spine + sum(pendants)
                if n > max_n:
                    break
                n_tree, adj = make_caterpillar(spine, pendants)
                label = f"cat({spine},[{k1},{k2}]alt)"
                result = analyze_tree(n_tree, adj, label)
                _strip_poly(result)
                results.append(result)

    # Random pendant patterns
    for _ in range(1000):
        spine = rng.randint(5, min(50, max_n // 2))
        pendants = [rng.randint(0, 8) for _ in range(spine)]
        n = spine + sum(pendants)
        if n > max_n:
            continue
        n_tree, adj = make_caterpillar(spine, pendants)
        result = analyze_tree(n_tree, adj, f"cat_rand({spine})")
        _strip_poly(result)
        results.append(result)

    # Report any interesting ones
    for r in results:
        if not r["log_concave"]:
            flag = "*** NON-UNIMODAL ***" if not r["unimodal"] else "LC-fail"
            print(
                f"  {r['label']:>30s}  n={r['n']:>5}  nm={r['nm_ratio']:.8f}  {flag}"
            )
    return results


def run_spiders_and_brooms(max_n: int) -> list[dict]:
    """Scan spiders and brooms."""
    results = []

    # Spiders: many legs of length 2 or 3 (Ramos & Sun motif)
    for num_legs in range(3, 200):
        for leg_len in [2, 3]:
            n = 1 + num_legs * leg_len
            if n > max_n:
                break
            legs = [leg_len] * num_legs
            n_tree, adj = make_spider(legs)
            result = analyze_tree(n_tree, adj, f"spider({num_legs}x{leg_len})")
            _strip_poly(result)
            results.append(result)

    # Mixed-leg spiders
    for n2 in range(1, 100):
        for n3 in range(1, 100):
            n = 1 + 2 * n2 + 3 * n3
            if n > max_n:
                break
            legs = [2] * n2 + [3] * n3
            n_tree, adj = make_spider(legs)
            result = analyze_tree(n_tree, adj, f"spider({n2}x2+{n3}x3)")
            _strip_poly(result)
            results.append(result)

    # Brooms
    for path_len in range(2, max_n):
        for star_size in range(2, max_n):
            n = path_len + star_size
            if n > max_n:
                break
            n_tree, adj = make_broom(path_len, star_size)
            result = analyze_tree(n_tree, adj, f"broom({path_len},{star_size})")
            _strip_poly(result)
            results.append(result)

    for r in results:
        if not r["log_concave"]:
            flag = "*** NON-UNIMODAL ***" if not r["unimodal"] else "LC-fail"
            print(
                f"  {r['label']:>30s}  n={r['n']:>5}  nm={r['nm_ratio']:.8f}  {flag}"
            )
    return results


def run_random_ramos_sun(max_n: int, count: int = 5000) -> list[dict]:
    """Generate random trees with Ramos & Sun structural bias."""
    results = []
    rng = random.Random(42)
    for i in range(count):
        n = rng.randint(26, max_n)
        n_tree, adj = make_ramos_sun_style(n, rng)
        result = analyze_tree(n_tree, adj, f"rs_rand_{i}(n={n})")
        results.append(result)
        if not result["unimodal"]:
            print(f"  *** NON-UNIMODAL: {result['label']} ***")
            print(f"      poly = {result['poly']}")
        _strip_poly(result)

    # Report top log-concavity failures
    lc_fails = [r for r in results if not r["log_concave"]]
    lc_fails.sort(key=lambda r: -r["lc_ratio"])
    if lc_fails:
        print(f"\n  {len(lc_fails)} log-concavity failures in random trees:")
        for r in lc_fails[:10]:
            flag = "*** NON-UNIMODAL ***" if not r["unimodal"] else "LC-fail"
            print(
                f"    {r['label']:>25s}  n={r['n']:>5}  "
                f"lc={r['lc_ratio']:.8f}  nm={r['nm_ratio']:.8f}  {flag}"
            )
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Targeted search for unimodality failures in specific tree families"
    )
    parser.add_argument(
        "--max-n", type=int, default=500, help="Maximum vertex count (default: 500)"
    )
    parser.add_argument(
        "--random-count",
        type=int,
        default=5000,
        help="Number of random Ramos-Sun-style trees (default: 5000)",
    )
    args = parser.parse_args()

    print(f"Targeted search for unimodality failures (max n = {args.max_n})", flush=True)
    print(f"{'='*70}\n", flush=True)

    all_results = []
    family_summaries: dict[str, dict] = {}

    # 1. Galvin's T_{m,t,1}
    print("1. Galvin spherically symmetric trees T_{{m,t,1}}...", flush=True)
    t0 = time.time()
    galvin = run_galvin_family(args.max_n)
    all_results.extend(galvin)
    lc_fails = sum(1 for r in galvin if not r["log_concave"])
    family_summaries["galvin_T_m_t_1"] = _summarize_family(galvin)
    print(f"   {len(galvin)} trees tested, {lc_fails} LC failures, "
          f"{time.time()-t0:.1f}s\n", flush=True)

    # 2. Generalized spherically symmetric trees
    print("2. Generalized spherically symmetric trees T_{{m,t,d}}...", flush=True)
    t0 = time.time()
    gen_sst = run_generalized_sst(args.max_n)
    all_results.extend(gen_sst)
    lc_fails = sum(1 for r in gen_sst if not r["log_concave"])
    family_summaries["generalized_T_m_t_d"] = _summarize_family(gen_sst)
    print(f"   {len(gen_sst)} trees tested, {lc_fails} LC failures, "
          f"{time.time()-t0:.1f}s\n", flush=True)

    # 3. Caterpillars
    print("3. Caterpillars...", flush=True)
    t0 = time.time()
    cats = run_caterpillars(args.max_n)
    all_results.extend(cats)
    lc_fails = sum(1 for r in cats if not r["log_concave"])
    family_summaries["caterpillars"] = _summarize_family(cats)
    print(f"   {len(cats)} trees tested, {lc_fails} LC failures, "
          f"{time.time()-t0:.1f}s\n", flush=True)

    # 4. Spiders and brooms
    print("4. Spiders and brooms...", flush=True)
    t0 = time.time()
    sb = run_spiders_and_brooms(args.max_n)
    all_results.extend(sb)
    lc_fails = sum(1 for r in sb if not r["log_concave"])
    family_summaries["spiders_and_brooms"] = _summarize_family(sb)
    print(f"   {len(sb)} trees tested, {lc_fails} LC failures, "
          f"{time.time()-t0:.1f}s\n", flush=True)

    # 5. Random Ramos-Sun-style trees
    print(f"5. Random Ramos-Sun-style trees ({args.random_count})...", flush=True)
    t0 = time.time()
    rs = run_random_ramos_sun(args.max_n, args.random_count)
    all_results.extend(rs)
    lc_fails = sum(1 for r in rs if not r["log_concave"])
    family_summaries["random_ramos_sun"] = _summarize_family(rs)
    print(f"   {len(rs)} trees tested, {lc_fails} LC failures, "
          f"{time.time()-t0:.1f}s\n", flush=True)

    # Summary
    total = len(all_results)
    total_lc_fail = sum(1 for r in all_results if not r["log_concave"])
    total_non_uni = sum(1 for r in all_results if not r["unimodal"])

    print(f"{'='*70}")
    print(f"SUMMARY: {total:,} trees tested across all families")
    print(f"  Log-concavity failures: {total_lc_fail}")
    print(f"  Unimodality failures:   {total_non_uni}")

    if total_non_uni > 0:
        print(f"\n*** COUNTEREXAMPLE(S) FOUND! ***")
        for r in all_results:
            if not r["unimodal"]:
                print(f"  {r['label']}: n={r['n']}, poly={r['poly']}")

    # Top near-misses across all families
    all_results.sort(key=lambda r: -r["nm_ratio"])
    print(f"\nTop 20 near-misses (highest a_{{j+1}}/a_j in descending tail):")
    for r in all_results[:20]:
        lc_tag = "LC-fail" if not r["log_concave"] else "LC-ok"
        print(
            f"  nm={r['nm_ratio']:.10f}  {r['label']:>30s}  "
            f"n={r['n']:>5}  alpha={r['alpha']:>4}  {lc_tag}"
        )

    # Save
    os.makedirs("results", exist_ok=True)
    # Save without full polynomials (too large)
    save_results = []
    for r in all_results[:500]:  # top 500 by near-miss
        save_r = {k: v for k, v in r.items() if k != "poly"}
        save_results.append(save_r)
    path = f"results/targeted_n{args.max_n}.json"
    with open(path, "w") as f:
        json.dump({
            "max_n": args.max_n,
            "total_tested": total,
            "lc_failures": total_lc_fail,
            "non_unimodal": total_non_uni,
            "top_results": save_results,
        }, f, indent=2)
    print(f"\nResults saved to {path}")

    fam_path = "results/targeted_families.json"
    with open(fam_path, "w") as f:
        json.dump(
            {
                "description": "Per-family summary for targeted.py runs",
                "source": {
                    "script": "targeted.py",
                    "max_n": args.max_n,
                    "random_count": args.random_count,
                },
                "seed": 42,
                "families": family_summaries,
                "totals": {
                    "total_tested": total,
                    "lc_failures": total_lc_fail,
                    "non_unimodal": total_non_uni,
                },
            },
            f,
            indent=2,
        )
    print(f"Per-family summary saved to {fam_path}")


if __name__ == "__main__":
    main()
