#!/usr/bin/env python3
"""Attack 1: mode-superadditivity checks for PLC polynomial products.

Target inequality (leftmost mode):
    mode(f * g) >= mode(f) + mode(g)

This script stress-tests the claim on:
1) simple PLC families (including binomials),
2) adversarial PLC sequences,
3) products of tree independence polynomials.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import time
from typing import Any

from indpoly import _polymul, independence_poly, is_log_concave
from trees import trees, trees_geng_raw


def mode_left(poly: list[int]) -> int:
    """Leftmost index attaining the maximum coefficient."""
    best_i = 0
    best_v = poly[0]
    for i, v in enumerate(poly):
        if v > best_v:
            best_v = v
            best_i = i
    return best_i


def mode_right(poly: list[int]) -> int:
    """Rightmost index attaining the maximum coefficient."""
    best_i = 0
    best_v = poly[0]
    for i, v in enumerate(poly):
        if v >= best_v:
            best_v = v
            best_i = i
    return best_i


def mean_index(poly: list[int]) -> float:
    z = sum(poly)
    if z == 0:
        return 0.0
    return sum(i * c for i, c in enumerate(poly)) / z


def is_positive_lc(poly: list[int]) -> bool:
    return bool(poly) and all(c > 0 for c in poly) and is_log_concave(poly)


def pair_metrics(f: list[int], g: list[int]) -> dict[str, Any]:
    h = _polymul(f, g)
    mf_l = mode_left(f)
    mg_l = mode_left(g)
    mh_l = mode_left(h)
    mf_r = mode_right(f)
    mg_r = mode_right(g)
    mh_r = mode_right(h)
    return {
        "mode_left_f": mf_l,
        "mode_left_g": mg_l,
        "mode_left_fg": mh_l,
        "lhs_left": mh_l,
        "rhs_left": mf_l + mg_l,
        "deficit_left": (mf_l + mg_l) - mh_l,
        "mode_right_f": mf_r,
        "mode_right_g": mg_r,
        "mode_right_fg": mh_r,
        "lhs_right": mh_r,
        "rhs_right": mf_r + mg_r,
        "deficit_right": (mf_r + mg_r) - mh_r,
        "fg": h,
    }


def update_pair_stats(stats: dict[str, Any], f: list[int], g: list[int], tag: str) -> None:
    m = pair_metrics(f, g)
    stats["checked_pairs"] += 1

    if m["lhs_left"] < m["rhs_left"]:
        stats["fail_left"] += 1
        if stats["first_fail_left"] is None:
            stats["first_fail_left"] = {
                "tag": tag,
                "f": f,
                "g": g,
                "fg": m["fg"],
                "mode_left_f": m["mode_left_f"],
                "mode_left_g": m["mode_left_g"],
                "mode_left_fg": m["mode_left_fg"],
                "deficit_left": m["deficit_left"],
                "mean_f": mean_index(f),
                "mean_g": mean_index(g),
                "mean_fg": mean_index(m["fg"]),
            }
    if m["lhs_right"] < m["rhs_right"]:
        stats["fail_right"] += 1
        if stats["first_fail_right"] is None:
            stats["first_fail_right"] = {
                "tag": tag,
                "f": f,
                "g": g,
                "fg": m["fg"],
                "mode_right_f": m["mode_right_f"],
                "mode_right_g": m["mode_right_g"],
                "mode_right_fg": m["mode_right_fg"],
                "deficit_right": m["deficit_right"],
                "mean_f": mean_index(f),
                "mean_g": mean_index(g),
                "mean_fg": mean_index(m["fg"]),
            }

    if m["deficit_left"] > stats["max_deficit_left"]:
        stats["max_deficit_left"] = m["deficit_left"]
        stats["max_deficit_left_witness"] = {
            "tag": tag,
            "f": f,
            "g": g,
            "fg": m["fg"],
            "mode_left_f": m["mode_left_f"],
            "mode_left_g": m["mode_left_g"],
            "mode_left_fg": m["mode_left_fg"],
        }

    if m["deficit_right"] > stats["max_deficit_right"]:
        stats["max_deficit_right"] = m["deficit_right"]
        stats["max_deficit_right_witness"] = {
            "tag": tag,
            "f": f,
            "g": g,
            "fg": m["fg"],
            "mode_right_f": m["mode_right_f"],
            "mode_right_g": m["mode_right_g"],
            "mode_right_fg": m["mode_right_fg"],
        }


def init_pair_stats() -> dict[str, Any]:
    return {
        "checked_pairs": 0,
        "fail_left": 0,
        "fail_right": 0,
        "max_deficit_left": 0,
        "max_deficit_right": 0,
        "first_fail_left": None,
        "first_fail_right": None,
        "max_deficit_left_witness": None,
        "max_deficit_right_witness": None,
    }


def run_binomial_checks(max_n: int) -> dict[str, Any]:
    out = init_pair_stats()
    for n1 in range(1, max_n + 1):
        f = [math.comb(n1, k) for k in range(n1 + 1)]
        for n2 in range(1, max_n + 1):
            g = [math.comb(n2, k) for k in range(n2 + 1)]
            update_pair_stats(out, f, g, f"binom({n1})*binom({n2})")
    out["max_n"] = max_n
    return out


def run_simple_explicit_cases() -> list[dict[str, Any]]:
    cases = [
        ("arith_4_self", [1, 2, 3, 4], [1, 2, 3, 4]),
        ("edge_poly_self", [1, 2], [1, 2]),
        ("arith_5_self", [1, 2, 3, 4, 5], [1, 2, 3, 4, 5]),
        ("geom_5_self", [1, 2, 4, 8, 16], [1, 2, 4, 8, 16]),
    ]
    out = []
    for name, f, g in cases:
        m = pair_metrics(f, g)
        out.append(
            {
                "name": name,
                "f": f,
                "g": g,
                "fg": m["fg"],
                "lhs_left": m["lhs_left"],
                "rhs_left": m["rhs_left"],
                "deficit_left": m["deficit_left"],
                "lhs_right": m["lhs_right"],
                "rhs_right": m["rhs_right"],
                "deficit_right": m["deficit_right"],
                "mean_f": mean_index(f),
                "mean_g": mean_index(g),
                "mean_fg": mean_index(m["fg"]),
                "is_plc_f": is_positive_lc(f),
                "is_plc_g": is_positive_lc(g),
                "is_plc_fg": is_positive_lc(m["fg"]),
            }
        )
    return out


def random_lc_sequence(
    rng: random.Random,
    length: int,
    max_coeff: int,
    bias_rising: bool,
) -> list[int]:
    """Generate a positive integer LC sequence by local LC constraints."""
    while True:
        a0 = rng.randint(1, max_coeff)
        a1 = rng.randint(1, max_coeff)
        seq = [a0, a1]
        ok = True
        for _ in range(2, length):
            hi = min(max_coeff, (seq[-1] * seq[-1]) // seq[-2])
            if hi < 1:
                ok = False
                break
            lo = 1
            if bias_rising:
                lo = min(max_coeff, seq[-1])
                if lo > hi:
                    lo = 1
            seq.append(rng.randint(lo, hi))
        if not ok:
            continue
        if is_positive_lc(seq):
            return seq


def build_adversarial_pool(
    rng: random.Random,
    count: int,
    len_min: int,
    len_max: int,
    max_coeff: int,
) -> list[list[int]]:
    pool: set[tuple[int, ...]] = set()

    # Deterministic "high-mode" PLC families.
    for L in range(max(2, len_min), len_max + 1):
        pool.add(tuple(range(1, L + 1)))
        pool.add(tuple(2**i for i in range(L)))
        pool.add(tuple(3**i for i in range(L)))

    # Random LC sequences.
    while len(pool) < count:
        L = rng.randint(len_min, len_max)
        bias = rng.random() < 0.7
        seq = random_lc_sequence(rng, L, max_coeff=max_coeff, bias_rising=bias)
        pool.add(tuple(seq))

    seqs = [list(t) for t in pool]
    seqs.sort(key=lambda s: (mode_left(s) - mean_index(s), len(s), s), reverse=True)
    return seqs


def run_adversarial_checks(
    rng: random.Random,
    pool_size: int,
    top_k_allpairs: int,
    random_pairs: int,
    len_min: int,
    len_max: int,
    max_coeff: int,
) -> dict[str, Any]:
    pool = build_adversarial_pool(
        rng=rng,
        count=pool_size,
        len_min=len_min,
        len_max=len_max,
        max_coeff=max_coeff,
    )
    top = pool[: min(top_k_allpairs, len(pool))]

    out = init_pair_stats()
    out["pool_size"] = len(pool)
    out["top_k_allpairs"] = len(top)
    out["random_pairs"] = random_pairs

    # Exhaustive pairs on top skew-heavy candidates.
    for i, f in enumerate(top):
        for j in range(i, len(top)):
            g = top[j]
            update_pair_stats(out, f, g, f"adversarial_top[{i}]*[{j}]")

    # Random pairs from full pool.
    for t in range(random_pairs):
        f = pool[rng.randrange(len(pool))]
        g = pool[rng.randrange(len(pool))]
        update_pair_stats(out, f, g, f"adversarial_random[{t}]")

    # Summarize skew statistics.
    skew_vals = [mode_left(s) - mean_index(s) for s in pool]
    out["max_mode_minus_mean_in_pool"] = max(skew_vals) if skew_vals else None
    out["min_mode_minus_mean_in_pool"] = min(skew_vals) if skew_vals else None
    out["median_mode_minus_mean_in_pool"] = sorted(skew_vals)[len(skew_vals) // 2] if skew_vals else None
    out["top_sequence_examples"] = pool[:10]
    return out


def run_tree_exhaustive_unique(max_n: int) -> dict[str, Any]:
    """Exhaustive pair checks over unique tree IS polynomials up to max_n."""
    poly_source: dict[tuple[int, ...], dict[str, Any]] = {}
    seen_trees = 0

    for n in range(1, max_n + 1):
        for nn, adj, raw in trees_geng_raw(n):
            seen_trees += 1
            poly = tuple(independence_poly(nn, adj))
            if poly not in poly_source:
                poly_source[poly] = {"n": nn, "g6": raw.decode("ascii").strip()}

    polys = [list(t) for t in poly_source]
    out = init_pair_stats()
    out["max_n"] = max_n
    out["tree_count_seen"] = seen_trees
    out["unique_poly_count"] = len(polys)

    for i, f in enumerate(polys):
        for j in range(i, len(polys)):
            g = polys[j]
            tag = f"tree_unique[{i}]*[{j}]"
            before_fail = out["fail_left"]
            update_pair_stats(out, f, g, tag)
            if out["fail_left"] > before_fail and out["first_fail_left"] is not None:
                out["first_fail_left"]["source_f"] = poly_source[tuple(f)]
                out["first_fail_left"]["source_g"] = poly_source[tuple(g)]
                break
        if out["first_fail_left"] is not None:
            break

    # Continue full scan even after first failure for fail-rate and max deficit.
    # Restart full pass (without breaking) for aggregate totals.
    out_full = init_pair_stats()
    out_full["max_n"] = max_n
    out_full["tree_count_seen"] = seen_trees
    out_full["unique_poly_count"] = len(polys)
    for i, f in enumerate(polys):
        for j in range(i, len(polys)):
            g = polys[j]
            update_pair_stats(out_full, f, g, f"tree_unique[{i}]*[{j}]")
    out_full["first_fail_left"] = out["first_fail_left"]
    return out_full


def run_tree_random_pairs(
    max_n: int,
    random_pairs: int,
    rng: random.Random,
) -> dict[str, Any]:
    """Random pair checks over all tree IS polynomials (with multiplicity)."""
    polys: list[list[int]] = []
    for n in range(1, max_n + 1):
        for nn, adj in trees(n, backend="geng"):
            polys.append(independence_poly(nn, adj))

    out = init_pair_stats()
    out["max_n"] = max_n
    out["poly_count_with_multiplicity"] = len(polys)
    out["random_pairs"] = random_pairs

    for t in range(random_pairs):
        f = polys[rng.randrange(len(polys))]
        g = polys[rng.randrange(len(polys))]
        update_pair_stats(out, f, g, f"tree_random[{t}]")

    return out


def write_json(path: str, payload: dict[str, Any]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def main() -> None:
    ap = argparse.ArgumentParser(description="Attack 1 PLC mode-superadditivity scan.")
    ap.add_argument("--seed", type=int, default=20260219)
    ap.add_argument("--binom-max-n", type=int, default=40)
    ap.add_argument("--adv-pool-size", type=int, default=3000)
    ap.add_argument("--adv-top-k-allpairs", type=int, default=350)
    ap.add_argument("--adv-random-pairs", type=int, default=200000)
    ap.add_argument("--adv-len-min", type=int, default=3)
    ap.add_argument("--adv-len-max", type=int, default=12)
    ap.add_argument("--adv-max-coeff", type=int, default=80)
    ap.add_argument("--tree-exhaustive-max-n", type=int, default=12)
    ap.add_argument("--tree-random-max-n", type=int, default=16)
    ap.add_argument("--tree-random-pairs", type=int, default=250000)
    ap.add_argument(
        "--out",
        default="results/attack1_mode_superadditivity_plc_2026_02_19.json",
    )
    args = ap.parse_args()

    rng = random.Random(args.seed)
    t0 = time.time()

    print("Attack 1: mode(f*g) >= mode(f)+mode(g) for PLC sequences", flush=True)
    print(f"seed={args.seed}", flush=True)

    simple_cases = run_simple_explicit_cases()
    binom = run_binomial_checks(max_n=args.binom_max_n)
    adv = run_adversarial_checks(
        rng=rng,
        pool_size=args.adv_pool_size,
        top_k_allpairs=args.adv_top_k_allpairs,
        random_pairs=args.adv_random_pairs,
        len_min=args.adv_len_min,
        len_max=args.adv_len_max,
        max_coeff=args.adv_max_coeff,
    )
    tree_exhaustive = run_tree_exhaustive_unique(max_n=args.tree_exhaustive_max_n)
    tree_random = run_tree_random_pairs(
        max_n=args.tree_random_max_n,
        random_pairs=args.tree_random_pairs,
        rng=rng,
    )

    payload = {
        "params": vars(args),
        "simple_cases": simple_cases,
        "binomial_family": binom,
        "adversarial_plc": adv,
        "tree_unique_exhaustive": tree_exhaustive,
        "tree_random_pairs": tree_random,
        "wall_s": time.time() - t0,
    }
    write_json(args.out, payload)

    def fmt_stats(name: str, s: dict[str, Any]) -> str:
        return (
            f"{name}: checked={s['checked_pairs']:,}, "
            f"fail_left={s['fail_left']:,}, fail_right={s['fail_right']:,}, "
            f"max_deficit_left={s['max_deficit_left']}, "
            f"max_deficit_right={s['max_deficit_right']}"
        )

    print(fmt_stats("binomial", binom), flush=True)
    print(fmt_stats("adversarial", adv), flush=True)
    print(fmt_stats("tree_unique_exhaustive", tree_exhaustive), flush=True)
    print(fmt_stats("tree_random_pairs", tree_random), flush=True)
    print(f"Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
