#!/usr/bin/env python3
"""Adversarial kill-test of the universal C2 base facts beyond pendant weight 24.

Background: ``notes/c2_bounded_pendant_core_2026-08-08.md`` proves that, for a
fixed pendant-arm pair ``(a, b)`` on two hubs, four base facts imply
log-concavity of the two-hub tree for *every* connector length:

    B is log-concave;  C is log-concave;  C ~_p B;  xC ~_p B.

All pairs with total pendant weight <= 24 were checked exhaustively. This
script hunts for a violation of any of these facts at total pendant weight
strictly above 24. ``C`` is the independence polynomial of the spider formed
by merging the two arm multisets at one hub, so C log-concavity is already a
theorem (Li-Li-Yang-Zhang, arXiv:2501.04245: all spiders are strongly
log-concave). Any apparent C-LC violation found here is therefore a harness
bug, not a mathematical finding, and is reported as a self-check.

All arithmetic is exact (Python ints / fractions.Fraction). Floats are used
only for human-readable display of margins.

Definitions (matching the task and the note)
---------------------------------------------
P_{-2} = [0] (treated as the zero polynomial), P_{-1} = P_0 = [1],
P_1 = [1, 1], P_n = P_{n-1} + x*P_{n-2}.

For arm multisets a = (a_1..a_r), b = (b_1..b_s) of positive integers:
    Q_a = prod_i P_{a_i},   R_a = prod_i P_{a_i - 1}   (empty product = [1])
    B = Q_a Q_b + x(R_a Q_b + Q_a R_b)
    C = Q_a Q_b + x R_a R_b

Partial synchronization A ~_p D: for all m >= n >= 0,
    a_m d_n + a_n d_m >= a_{m+1} d_{n-1} + a_{n-1} d_{m+1}
with out-of-range coefficients read as 0.

Margin convention
------------------
For every inequality "required >= bound" checked here, margin is defined as
``required - bound`` (an exact integer; margin >= 0 means the fact holds at
that index), and the normalized margin is the exact Fraction
``margin / max(1, |bound|)``. For log-concavity, required = a_k^2,
bound = a_{k-1} a_{k+1}. For partial synchronization at indices
(lower <= upper), required = a_upper b_lower + a_lower b_upper,
bound = a_{upper+1} b_{lower-1} + a_{lower-1} b_{upper+1}.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Iterable

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from indpoly import independence_poly, is_log_concave
from verify_c2_bounded_pendant_core_20260808 import (
    adjacent_and_contracted,
    connector_adjacency,
    connector_formula,
    rooted_arm_state,
)
from verify_connector_partial_sync_route_20260808 import shift

Poly = list[int]
ArmVector = tuple[int, ...]
State = tuple[ArmVector, Poly, Poly]

SEED = 993  # deterministic; Erdos problem number, fixed for reproducibility.


# --------------------------------------------------------------------------
# Core polynomial construction, reusing the tested helpers from the bounded
# pendant-core certificate.
# --------------------------------------------------------------------------


def state(arms: ArmVector) -> State:
    return (arms, *rooted_arm_state(arms))


def build_B_C(a: ArmVector, b: ArmVector) -> tuple[Poly, Poly]:
    return adjacent_and_contracted(state(a), state(b))


def spider_adjacency(arms: ArmVector) -> list[list[int]]:
    """Adjacency list for the single-hub spider with the given arm lengths."""
    order = 1 + sum(arms)
    adjacency = [[] for _ in range(order)]

    def add_edge(u: int, v: int) -> None:
        adjacency[u].append(v)
        adjacency[v].append(u)

    next_vertex = 1
    for length in arms:
        previous = 0
        for _ in range(length):
            add_edge(previous, next_vertex)
            previous = next_vertex
            next_vertex += 1
    assert next_vertex == order
    return adjacency


def independent_dp_confirm(a: ArmVector, b: ArmVector, B: Poly, C: Poly) -> dict:
    """Independently rebuild B (t=0 two-hub tree) and C (merged spider) via
    direct rooted-tree DP and confirm they equal the algebraic B, C."""
    left_state = state(a)
    right_state = state(b)
    b_formula = connector_formula(left_state, right_state, 0)
    b_adjacency = connector_adjacency(a, b, 0)
    b_dp = independence_poly(len(b_adjacency), b_adjacency)

    combined = tuple(sorted(a + b))
    c_adjacency = spider_adjacency(combined)
    c_dp = independence_poly(len(c_adjacency), c_adjacency)

    return {
        "B_formula_matches_algebraic": b_formula == B,
        "B_dp_matches_algebraic": b_dp == B,
        "C_dp_matches_algebraic": c_dp == C,
    }


# --------------------------------------------------------------------------
# Exact margin computation.
# --------------------------------------------------------------------------


def lc_min_margin(poly: Poly) -> dict | None:
    """Minimum normalized log-concavity margin over all interior indices."""
    best_k = None
    best_margin = None
    best_lhs = best_rhs = None
    best_num = best_den = None  # margin / den, den > 0, compared via cross-mult
    for k in range(1, len(poly) - 1):
        lhs = poly[k] * poly[k]
        rhs = poly[k - 1] * poly[k + 1]
        margin = lhs - rhs
        den = max(1, abs(rhs))
        if best_num is None or margin * best_den < best_num * den:
            best_num, best_den = margin, den
            best_k, best_margin, best_lhs, best_rhs = k, margin, lhs, rhs
    if best_k is None:
        return None
    return {
        "k": best_k,
        "margin": best_margin,
        "lhs": best_lhs,
        "rhs": best_rhs,
        "normalized": Fraction(best_margin, best_den),
    }


def banded_pairs(degree: int, top_k: int = 60, bottom_k: int = 60, diag_width: int = 10):
    """(lower, upper) pairs near the top corner, bottom corner, and a thin
    near-diagonal stripe across the whole range. Cheap alternative to the
    full O(d^2) scan for large degrees."""
    limit = degree + 2
    pairs = set()
    top_start = max(0, limit - top_k)
    for lower in range(top_start, limit):
        for upper in range(lower, limit):
            pairs.add((lower, upper))
    bottom_end = min(limit, bottom_k)
    for lower in range(0, bottom_end):
        for upper in range(lower, min(limit, lower + bottom_k)):
            pairs.add((lower, upper))
    for lower in range(limit):
        for upper in range(lower, min(limit, lower + diag_width)):
            pairs.add((lower, upper))
    return pairs


def partial_sync_min_margin(
    a: Poly,
    b: Poly,
    pairs: Iterable[tuple[int, int]] | None = None,
) -> dict | None:
    """Minimum normalized partial-synchronization margin.

    If ``pairs`` is None, scans the full (lower, upper) grid (O(d^2)).
    Otherwise scans only the supplied pairs (e.g. from banded_pairs).
    """
    degree = max(len(a), len(b)) - 1
    pad = degree + 3
    aa = a + [0] * (pad - len(a))
    bb = b + [0] * (pad - len(b))

    def coeff(p: Poly, i: int) -> int:
        return 0 if i < 0 else p[i]

    if pairs is None:
        limit = degree + 2
        pairs = (
            (lower, upper)
            for lower in range(limit)
            for upper in range(lower, limit)
        )

    best_lower = best_upper = None
    best_margin = best_required = best_bound = None
    best_num = best_den = None
    for lower, upper in pairs:
        required = coeff(aa, upper) * coeff(bb, lower) + coeff(aa, lower) * coeff(bb, upper)
        bound = coeff(aa, upper + 1) * coeff(bb, lower - 1) + coeff(aa, lower - 1) * coeff(bb, upper + 1)
        if required == 0 and bound == 0:
            # Vacuous: both sides are forced to 0 by out-of-range or
            # shift-induced zero coefficients. Not an informative margin.
            continue
        margin = required - bound
        den = max(1, abs(bound))
        if best_num is None or margin * best_den < best_num * den:
            best_num, best_den = margin, den
            best_lower, best_upper = lower, upper
            best_margin, best_required, best_bound = margin, required, bound
    if best_lower is None:
        return None
    return {
        "lower": best_lower,
        "upper": best_upper,
        "margin": best_margin,
        "required": best_required,
        "bound": best_bound,
        "normalized": Fraction(best_margin, best_den),
    }


def evaluate_pair(a: ArmVector, b: ArmVector, full_sync: bool = True, band=None) -> dict:
    a = tuple(sorted(a))
    b = tuple(sorted(b))
    B, C = build_B_C(a, b)
    degree = len(B) - 1

    b_lc = lc_min_margin(B)
    c_lc = lc_min_margin(C)

    pairs = None if full_sync else banded_pairs(degree, **(band or {}))
    c_sync_b = partial_sync_min_margin(C, B, pairs)
    xc_sync_b = partial_sync_min_margin(shift(C), B, pairs)

    facts = {"B_lc": b_lc, "C_sync_B": c_sync_b, "xC_sync_B": xc_sync_b}
    overall = min(facts.items(), key=lambda kv: kv[1]["normalized"] if kv[1] else Fraction(1))

    return {
        "a": list(a),
        "b": list(b),
        "weight": sum(a) + sum(b),
        "degree": degree,
        "full_sync": full_sync,
        "B_lc": b_lc,
        "C_lc_self_check": c_lc,
        "C_sync_B": c_sync_b,
        "xC_sync_B": xc_sync_b,
        "min_open_fact": overall[0],
        "min_open_fact_normalized": overall[1]["normalized"] if overall[1] else None,
    }


def is_violation(result: dict) -> list[str]:
    """Which of the three open facts (B_lc, C_sync_B, xC_sync_B) actually
    fail (normalized margin < 0) in this exact result."""
    failed = []
    for key in ("B_lc", "C_sync_B", "xC_sync_B"):
        fact = result[key]
        if fact is not None and fact["normalized"] < 0:
            failed.append(key)
    return failed


def frac_str(x: Fraction | None) -> str | None:
    return None if x is None else f"{x.numerator}/{x.denominator}"


def result_to_json(result: dict) -> dict:
    out = dict(result)
    out["min_open_fact_normalized"] = frac_str(result["min_open_fact_normalized"])
    for key in ("B_lc", "C_lc_self_check", "C_sync_B", "xC_sync_B"):
        fact = result[key]
        if fact is not None:
            fact = dict(fact)
            fact["normalized"] = frac_str(fact["normalized"])
            fact["normalized_float"] = float(result[key]["normalized"])
        out[key] = fact
    return out


# --------------------------------------------------------------------------
# Search-space generators (structured grid).
# --------------------------------------------------------------------------


def family_1r_ls(r_max: int, s_max: int, l_max: int):
    for r in range(1, r_max + 1):
        for s in range(1, s_max + 1):
            for l in range(1, l_max + 1):
                yield (1,) * r, (l,) * s


def family_kk_mm(k_max: int, m_max: int):
    for k in range(1, k_max + 1):
        for m in range(1, m_max + 1):
            yield (k, k), (m, m)


def family_1r_M(r_values: Iterable[int], m_range: Iterable[int]):
    for r in r_values:
        for m in m_range:
            yield (1,) * r, (m,)


THREE_WAY_SHAPES = [
    (1, 2, 3),
    (1, 2, 3, 4),
    (2, 3, 5),
    (1, 1, 2, 3, 5, 8),
    (3, 5, 7),
    (1, 4, 9),
]


def family_three_way(m_range: Iterable[int]):
    for shape in THREE_WAY_SHAPES:
        for m in m_range:
            yield tuple(shape), (m,)
        for mult in (2, 3):
            yield tuple(shape), tuple(x * mult for x in shape)


def family_11_2h_2h(h_range: Iterable[int]):
    for h in h_range:
        if h < 1:
            continue
        yield (1, 1), (2 * h, 2 * h)


def family_11_2hp1_2hp3(h_range: Iterable[int]):
    for h in h_range:
        if h < 0:
            continue
        yield (1, 1), (2 * h + 1, 2 * h + 3)


def random_partition(total: int, rng: random.Random, style: str) -> tuple[int, ...]:
    """Random positive-integer partition of ``total`` in one of three styles."""
    if total <= 0:
        return ()
    if style == "many_ones_one_giant":
        if total == 1:
            return (1,)
        giant = rng.randint(max(1, total // 2), total - 1)
        rest = total - giant
        return tuple(sorted((giant,) + (1,) * rest))
    if style == "balanced":
        parts = []
        remaining = total
        target_count = max(1, rng.randint(2, 8))
        while remaining > 0 and len(parts) < target_count - 1:
            share = max(1, remaining // (target_count - len(parts)))
            share = rng.randint(1, max(1, min(remaining - (target_count - len(parts) - 1), share * 2)))
            share = min(share, remaining)
            parts.append(share)
            remaining -= share
        if remaining > 0:
            parts.append(remaining)
        return tuple(sorted(p for p in parts if p > 0))
    if style == "geometric":
        parts = []
        remaining = total
        length = 1
        while remaining > 0:
            take = min(remaining, length)
            parts.append(take)
            remaining -= take
            length *= 2
        return tuple(sorted(p for p in parts if p > 0))
    raise ValueError(f"unknown style: {style}")


def random_pairs(count: int, min_weight: int, max_weight: int, seed: int):
    rng = random.Random(seed)
    styles = ["many_ones_one_giant", "balanced", "geometric"]
    for _ in range(count):
        total = rng.randint(min_weight, max_weight)
        left_weight = rng.randint(1, total - 1) if total > 1 else 1
        right_weight = total - left_weight
        left_style = rng.choice(styles)
        right_style = rng.choice(styles)
        a = random_partition(left_weight, rng, left_style)
        b = random_partition(right_weight, rng, right_style)
        if not a and not b:
            continue
        yield a, b


# --------------------------------------------------------------------------
# Hill climbing.
# --------------------------------------------------------------------------


def mutate(a: ArmVector, b: ArmVector, rng: random.Random, max_weight: int) -> tuple[ArmVector, ArmVector]:
    a, b = list(a), list(b)
    move = rng.choice(["add", "remove", "tweak", "transfer"])
    sides = [a, b]

    def total_weight() -> int:
        return sum(a) + sum(b)

    if move == "add" and total_weight() < max_weight:
        side = rng.choice(sides)
        side.append(rng.randint(1, 3))
    elif move == "remove":
        side = rng.choice([s for s in sides if s] or sides)
        if side:
            side.pop(rng.randrange(len(side)))
    elif move == "tweak":
        side = rng.choice([s for s in sides if s] or sides)
        if side:
            idx = rng.randrange(len(side))
            delta = rng.choice([-1, 1])
            new_val = side[idx] + delta
            if new_val <= 0:
                side.pop(idx)
            elif total_weight() + delta <= max_weight:
                side[idx] = new_val
    elif move == "transfer":
        src, dst = (a, b) if rng.random() < 0.5 else (b, a)
        if src:
            idx = rng.randrange(len(src))
            dst.append(src.pop(idx))

    return tuple(sorted(a)), tuple(sorted(b))


def margin_of(result: dict) -> Fraction:
    return result["min_open_fact_normalized"]


def hill_climb(
    seeds: list[tuple[ArmVector, ArmVector]],
    steps: int,
    restarts: int,
    max_weight: int,
    seed: int,
    budget_s: float | None = None,
) -> tuple[list[dict], dict | None, bool]:
    """Return (trajectory summaries, best result found overall, stopped_early)."""
    rng = random.Random(seed)
    trajectories = []
    global_best = None
    start_time = time.time()
    stopped_early = False

    for seed_idx, (start_a, start_b) in enumerate(seeds):
        if stopped_early:
            break
        for restart in range(restarts):
            if budget_s is not None and (time.time() - start_time) > budget_s:
                stopped_early = True
                break
            current_a, current_b = start_a, start_b
            current_result = evaluate_pair(current_a, current_b, full_sync=True)
            best_in_run = current_result
            accepted_steps = 0
            for step in range(steps):
                cand_a, cand_b = mutate(current_a, current_b, rng, max_weight)
                if not cand_a and not cand_b:
                    continue
                cand_result = evaluate_pair(cand_a, cand_b, full_sync=True)
                if margin_of(cand_result) < margin_of(current_result):
                    current_a, current_b = cand_a, cand_b
                    current_result = cand_result
                    accepted_steps += 1
                    if margin_of(current_result) < margin_of(best_in_run):
                        best_in_run = current_result
                violations = is_violation(cand_result)
                if violations:
                    global_best = cand_result
                    trajectories.append({
                        "seed_index": seed_idx,
                        "restart": restart,
                        "start": [list(start_a), list(start_b)],
                        "steps_run": step + 1,
                        "accepted_steps": accepted_steps,
                        "end": [list(cand_a), list(cand_b)],
                        "end_min_margin": str(margin_of(cand_result)),
                        "violation_found": True,
                        "violated_facts": violations,
                    })
                    return trajectories, global_best, False

            if global_best is None or margin_of(best_in_run) < margin_of(global_best):
                global_best = best_in_run

            trajectories.append({
                "seed_index": seed_idx,
                "restart": restart,
                "start": [list(start_a), list(start_b)],
                "start_min_margin": str(margin_of(current_result)) if steps == 0 else None,
                "steps_run": steps,
                "accepted_steps": accepted_steps,
                "end": [list(current_a), list(current_b)],
                "end_min_margin": str(margin_of(current_result)),
                "end_min_fact": current_result["min_open_fact"],
                "violation_found": False,
            })

    return trajectories, global_best, stopped_early


# --------------------------------------------------------------------------
# Tracking and orchestration.
# --------------------------------------------------------------------------


class Tracker:
    def __init__(self) -> None:
        self.top_candidates: list[dict] = []
        self.violations: list[dict] = []
        self.global_min: dict[str, dict | None] = {
            "B_lc": None,
            "C_sync_B": None,
            "xC_sync_B": None,
        }
        self.phase_coverage: list[dict] = []
        self.evaluated_total = 0

    def _update_global_min(self, result: dict) -> None:
        for key in ("B_lc", "C_sync_B", "xC_sync_B"):
            fact = result[key]
            if fact is None:
                continue
            current = self.global_min[key]
            if current is None or fact["normalized"] < current["fact"]["normalized"]:
                self.global_min[key] = {
                    "fact": fact,
                    "a": result["a"],
                    "b": result["b"],
                    "weight": result["weight"],
                    "full_sync": result["full_sync"],
                }

    def record(self, result: dict) -> list[str]:
        self.evaluated_total += 1
        self._update_global_min(result)
        self.top_candidates.append(result)
        if len(self.top_candidates) > 400:
            self.top_candidates.sort(key=lambda r: r["min_open_fact_normalized"])
            self.top_candidates = self.top_candidates[:60]
        violations = is_violation(result)
        if violations:
            confirm = independent_dp_confirm(
                tuple(result["a"]), tuple(result["b"]),
                *build_B_C(tuple(result["a"]), tuple(result["b"])),
            )
            entry = {
                "a": result["a"],
                "b": result["b"],
                "violated_facts": violations,
                "result": result_to_json(result),
                "independent_dp_confirmation": confirm,
            }
            self.violations.append(entry)
        return violations

    def finalize_top(self, n: int = 20) -> list[dict]:
        self.top_candidates.sort(key=lambda r: r["min_open_fact_normalized"])
        return self.top_candidates[:n]


def run_phase(
    tracker: Tracker,
    name: str,
    pairs: Iterable[tuple[ArmVector, ArmVector]],
    full_sync: bool,
    band: dict | None,
    time_budget_s: float | None,
    stop_on_violation: bool = True,
) -> dict:
    start = time.time()
    count = 0
    stopped_early = False
    found_violation = False
    for a, b in pairs:
        if not a and not b:
            continue
        result = evaluate_pair(a, b, full_sync=full_sync, band=band)
        violations = tracker.record(result)
        count += 1
        if violations and stop_on_violation:
            found_violation = True
            break
        if time_budget_s is not None and (time.time() - start) > time_budget_s:
            stopped_early = True
            break
    elapsed = time.time() - start
    coverage = {
        "phase": name,
        "evaluated": count,
        "elapsed_s": round(elapsed, 2),
        "full_sync": full_sync,
        "band": band,
        "stopped_early_on_time_budget": stopped_early,
        "found_violation": found_violation,
    }
    tracker.phase_coverage.append(coverage)
    print(
        f"[{name}] evaluated={count} elapsed={elapsed:.1f}s "
        f"full_sync={full_sync} stopped_early={stopped_early} "
        f"violation={found_violation}",
        flush=True,
    )
    return coverage


def build_phases(args: argparse.Namespace) -> list[dict]:
    """Return the ordered list of structured-grid phase specs."""
    phases = []

    phases.append({
        "name": "grid_1r_ls",
        "pairs": family_1r_ls(args.r_max, args.s_max, args.l_max),
        "full_sync": True,
        "band": None,
        "time_budget_s": args.phase_budget_s,
    })

    phases.append({
        "name": "family_kk_mm",
        "pairs": family_kk_mm(args.kk_k_max, args.kk_m_max),
        "full_sync": True,
        "band": None,
        "time_budget_s": args.phase_budget_s,
    })

    phases.append({
        "name": "family_11_mm_dense",
        "pairs": [((1, 1), (m, m)) for m in range(1, args.mm_dense_max + 1)],
        "full_sync": True,
        "band": None,
        "time_budget_s": args.phase_budget_s,
    })

    phases.append({
        "name": "family_11_mm_sparse",
        # Largest parameter first: if the phase time budget cuts this short,
        # the highest-value (largest, least-covered) points are checked
        # before the smaller ones already covered by the dense tier.
        "pairs": [((1, 1), (m, m)) for m in sorted(args.mm_sparse_points, reverse=True)],
        "full_sync": True,
        "band": None,
        "time_budget_s": args.mm_sparse_budget_s,
    })

    phases.append({
        "name": "family_11_mm_cheap_sweep",
        "pairs": [
            ((1, 1), (m, m))
            for m in range(1, args.mm_cheap_max + 1, args.mm_cheap_step)
        ],
        "full_sync": False,
        "band": {"top_k": 80, "bottom_k": 80, "diag_width": 15},
        "time_budget_s": args.mm_cheap_budget_s,
    })

    phases.append({
        "name": "family_11_2h_2h_dense",
        "pairs": family_11_2h_2h(range(1, args.neighbor_dense_max + 1)),
        "full_sync": True,
        "band": None,
        "time_budget_s": args.phase_budget_s,
    })

    phases.append({
        "name": "family_11_2h_2h_cheap_sweep",
        "pairs": family_11_2h_2h(
            range(1, args.neighbor_cheap_max + 1, args.neighbor_cheap_step)
        ),
        "full_sync": False,
        "band": {"top_k": 80, "bottom_k": 80, "diag_width": 15},
        "time_budget_s": args.neighbor_cheap_budget_s,
    })

    phases.append({
        "name": "family_11_2hp1_2hp3_dense",
        "pairs": family_11_2hp1_2hp3(range(0, args.neighbor_dense_max + 1)),
        "full_sync": True,
        "band": None,
        "time_budget_s": args.phase_budget_s,
    })

    phases.append({
        "name": "family_11_2hp1_2hp3_cheap_sweep",
        "pairs": family_11_2hp1_2hp3(
            range(0, args.neighbor_cheap_max + 1, args.neighbor_cheap_step)
        ),
        "full_sync": False,
        "band": {"top_k": 80, "bottom_k": 80, "diag_width": 15},
        "time_budget_s": args.neighbor_cheap_budget_s,
    })

    phases.append({
        "name": "family_1r_M_dense",
        "pairs": family_1r_M(range(1, 13), range(1, args.single_arm_dense_max + 1)),
        "full_sync": True,
        "band": None,
        "time_budget_s": args.phase_budget_s,
    })

    phases.append({
        "name": "family_1r_M_cheap_sweep",
        "pairs": family_1r_M(
            args.single_arm_reps,
            range(1, args.single_arm_cheap_max + 1, args.single_arm_cheap_step),
        ),
        "full_sync": False,
        "band": {"top_k": 80, "bottom_k": 80, "diag_width": 15},
        "time_budget_s": args.single_arm_cheap_budget_s,
    })

    phases.append({
        "name": "three_way_mixes",
        "pairs": family_three_way(range(1, args.three_way_max + 1)),
        "full_sync": True,
        "band": None,
        "time_budget_s": args.phase_budget_s,
    })

    phases.append({
        "name": "random_search",
        "pairs": random_pairs(
            args.random_count, args.random_min_weight, args.random_max_weight, SEED
        ),
        "full_sync": True,
        "band": None,
        "time_budget_s": args.phase_budget_s,
    })

    return phases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-max", type=int, default=12)
    parser.add_argument("--s-max", type=int, default=12)
    parser.add_argument("--l-max", type=int, default=40)
    parser.add_argument("--kk-k-max", type=int, default=5)
    parser.add_argument("--kk-m-max", type=int, default=200)
    parser.add_argument("--mm-dense-max", type=int, default=700)
    parser.add_argument(
        "--mm-sparse-points", type=int, nargs="*",
        default=[800, 900, 1000, 1200, 1500, 1750, 2000, 2500, 3000],
    )
    parser.add_argument("--mm-sparse-budget-s", type=float, default=900.0)
    parser.add_argument("--mm-cheap-max", type=int, default=3000)
    parser.add_argument("--mm-cheap-step", type=int, default=10)
    parser.add_argument("--mm-cheap-budget-s", type=float, default=900.0)
    parser.add_argument("--neighbor-dense-max", type=int, default=300)
    parser.add_argument("--neighbor-cheap-max", type=int, default=800)
    parser.add_argument("--neighbor-cheap-step", type=int, default=5)
    parser.add_argument("--neighbor-cheap-budget-s", type=float, default=300.0)
    parser.add_argument("--single-arm-dense-max", type=int, default=300)
    parser.add_argument("--single-arm-cheap-max", type=int, default=1000)
    parser.add_argument("--single-arm-cheap-step", type=int, default=10)
    parser.add_argument("--single-arm-cheap-budget-s", type=float, default=400.0)
    parser.add_argument(
        "--single-arm-reps", type=int, nargs="*", default=[1, 6, 12]
    )
    parser.add_argument("--three-way-max", type=int, default=300)
    parser.add_argument("--random-count", type=int, default=2000)
    parser.add_argument("--random-min-weight", type=int, default=25)
    parser.add_argument("--random-max-weight", type=int, default=150)
    parser.add_argument("--phase-budget-s", type=float, default=240.0)
    parser.add_argument("--hillclimb-steps", type=int, default=200)
    parser.add_argument("--hillclimb-restarts", type=int, default=10)
    parser.add_argument("--hillclimb-seeds", type=int, default=20)
    parser.add_argument("--hillclimb-max-weight", type=int, default=220)
    parser.add_argument("--hillclimb-budget-s", type=float, default=900.0)
    parser.add_argument("--skip-hillclimb", action="store_true")
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/c2_universal_killtest_20260811.json"),
    )
    args = parser.parse_args()

    overall_start = time.time()
    tracker = Tracker()

    for spec in build_phases(args):
        coverage = run_phase(
            tracker,
            spec["name"],
            spec["pairs"],
            spec["full_sync"],
            spec["band"],
            spec["time_budget_s"],
        )
        if coverage["found_violation"]:
            print(f"VIOLATION FOUND in phase {spec['name']}; stopping further phases.")
            break

    hillclimb_trajectories = []
    hillclimb_best = None
    hillclimb_stopped_early = False
    already_violated = bool(tracker.violations)
    if not args.skip_hillclimb and not already_violated:
        seeds_pool = [
            (tuple(r["a"]), tuple(r["b"]))
            for r in tracker.finalize_top(200)
            if r["weight"] <= args.hillclimb_max_weight
        ][: args.hillclimb_seeds]
        if len(seeds_pool) < args.hillclimb_seeds:
            print(
                f"Only {len(seeds_pool)} seeds under weight cap "
                f"{args.hillclimb_max_weight} available for hill-climb "
                f"(requested {args.hillclimb_seeds})."
            )
        hc_start = time.time()
        hillclimb_trajectories, hillclimb_best, hillclimb_stopped_early = hill_climb(
            seeds_pool,
            args.hillclimb_steps,
            args.hillclimb_restarts,
            args.hillclimb_max_weight,
            SEED,
            budget_s=args.hillclimb_budget_s,
        )
        print(
            f"[hill_climb] seeds={len(seeds_pool)} elapsed={time.time()-hc_start:.1f}s "
            f"trajectories={len(hillclimb_trajectories)} "
            f"stopped_early={hillclimb_stopped_early}",
            flush=True,
        )
        # Record hill-climb evaluations into the global tracker too, by
        # re-evaluating the best-found pair from each trajectory endpoint.
        for traj in hillclimb_trajectories:
            end_a, end_b = tuple(traj["end"][0]), tuple(traj["end"][1])
            if not end_a and not end_b:
                continue
            end_result = evaluate_pair(end_a, end_b, full_sync=True)
            tracker.record(end_result)

    top20 = tracker.finalize_top(20)

    scan_two_branch_note = (
        "scan_two_branch_lc.py exhaustively enumerates C2 trees by total "
        "vertex order n (default max_n=24), which bounds pendant weight plus "
        "hubs plus connector vertices together and uses a fixed (small) "
        "connector length implied by that order; it therefore does not reach "
        "pendant weight above 24 combined with an unbounded connector, and "
        "does not cover the search space explored here. Not rerun."
    )

    report = {
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "question": (
            "Does any pair of pendant-arm multisets on two hubs, with total "
            "pendant weight strictly above 24, violate B log-concavity, "
            "C ~_p B, or xC ~_p B, where B is the adjacent-hub (t=0) "
            "polynomial and C is the contracted-hub (merged spider) "
            "polynomial? C log-concavity is checked only as a harness "
            "self-check (it is a theorem: Li-Li-Yang-Zhang, "
            "arXiv:2501.04245, all spiders are strongly log-concave)."
        ),
        "seed": SEED,
        "scan_two_branch_lc_coverage_note": scan_two_branch_note,
        "phase_coverage": tracker.phase_coverage,
        "evaluated_total": tracker.evaluated_total,
        "global_min_margins": {
            key: (
                None if value is None else {
                    "normalized_fraction": frac_str(value["fact"]["normalized"]),
                    "normalized_float": float(value["fact"]["normalized"]),
                    "a": value["a"],
                    "b": value["b"],
                    "weight": value["weight"],
                    "full_sync": value["full_sync"],
                    "location": {
                        k: v for k, v in value["fact"].items()
                        if k not in ("normalized",)
                    },
                }
            )
            for key, value in tracker.global_min.items()
        },
        "tightest_20_pairs": [result_to_json(r) for r in top20],
        "violations": tracker.violations,
        "violation_count": len(tracker.violations),
        "hillclimb": {
            "seeds_requested": args.hillclimb_seeds,
            "steps": args.hillclimb_steps,
            "restarts": args.hillclimb_restarts,
            "max_weight_cap": args.hillclimb_max_weight,
            "trajectories": hillclimb_trajectories,
            "best_found_min_margin": (
                None if hillclimb_best is None
                else str(hillclimb_best["min_open_fact_normalized"])
            ),
            "skipped": args.skip_hillclimb or already_violated,
            "stopped_early_on_time_budget": hillclimb_stopped_early,
        },
        "elapsed_total_s": round(time.time() - overall_start, 2),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "out": str(args.out),
        "evaluated_total": tracker.evaluated_total,
        "violation_count": len(tracker.violations),
        "elapsed_total_s": report["elapsed_total_s"],
    }))


if __name__ == "__main__":
    main()
