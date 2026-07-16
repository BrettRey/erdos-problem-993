#!/usr/bin/env python3
"""Certified root stress tests for two literature-defined tree families.

The first lane reconstructs the Jerrum--Patel balanced trees from
arXiv:2510.01466v2, Section 6.  It keeps separate the uniformly subdivided
complete-binary family stated in Lemma 17 and the phase-truncated periodic
branching words produced by the proof's appeal to Buys' Lemma 13.

The second lane reconstructs the Bautista-Ramos--Guillen-Galvan--
Gomez-Salgado pattern families from arXiv:2603.14204.  It keeps separate the
P2-pendant families governed by the BKW circle |z+1/3|=1/3 and the S_{2,n}
families with consecutive log-concavity failures and a different
equimodular locus.

All graph polynomials and recurrence checks use exact integers.  Root balls,
modulus ordering, threshold membership, and angles use Arb through
python-flint.  Floating-point values in the JSON are reporting midpoints and
endpoints only; no classification decision uses them.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.fspath(ROOT))

from indpoly import independence_poly, is_unimodal  # noqa: E402
from scripts.analyze_prufer_corpus import make_li_tree  # noqa: E402


TAUS = (1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0)
DEFAULT_OUTPUT = ROOT / "results" / "literature_root_stress_20260716.json"


def poly_trim(poly: Sequence[int]) -> list[int]:
    """Return a coefficient list without trailing zero coefficients."""
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out or [0]


def poly_add(left: Sequence[int], right: Sequence[int]) -> list[int]:
    """Add exact coefficient lists."""
    out = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    return poly_trim(out)


def poly_neg(poly: Sequence[int]) -> list[int]:
    """Negate an exact coefficient list."""
    return [-value for value in poly]


def poly_sub(left: Sequence[int], right: Sequence[int]) -> list[int]:
    """Subtract exact coefficient lists."""
    return poly_add(left, poly_neg(right))


def poly_mul(left: Sequence[int], right: Sequence[int]) -> list[int]:
    """Multiply exact coefficient lists using arbitrary-precision integers."""
    out = [0] * (len(left) + len(right) - 1)
    for i, a_i in enumerate(left):
        if a_i == 0:
            continue
        for j, b_j in enumerate(right):
            if b_j:
                out[i + j] += a_i * b_j
    return poly_trim(out)


def poly_pow(poly: Sequence[int], exponent: int) -> list[int]:
    """Raise an exact coefficient list to a nonnegative integer power."""
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    result = [1]
    base = list(poly)
    power = exponent
    while power:
        if power & 1:
            result = poly_mul(result, base)
        power >>= 1
        if power:
            base = poly_mul(base, base)
    return result


def poly_shift(poly: Sequence[int], amount: int = 1) -> list[int]:
    """Multiply an exact coefficient list by x**amount."""
    if amount < 0:
        raise ValueError("shift amount must be nonnegative")
    return [0] * amount + list(poly)


def lc_failure_indices(poly: Sequence[int]) -> list[int]:
    """Return every exact log-concavity failure index."""
    return [
        index
        for index in range(1, len(poly) - 1)
        if poly[index] * poly[index]
        < poly[index - 1] * poly[index + 1]
    ]


def validate_tree(adj: Sequence[Sequence[int]]) -> dict[str, Any]:
    """Validate a simple adjacency list as a connected tree."""
    n = len(adj)
    if n == 0:
        raise AssertionError("tree must be nonempty")
    if any(vertex in neighbors for vertex, neighbors in enumerate(adj)):
        raise AssertionError("self-loop in tree constructor")
    for vertex, neighbors in enumerate(adj):
        if len(set(neighbors)) != len(neighbors):
            raise AssertionError("parallel adjacency entry in tree constructor")
        for neighbor in neighbors:
            if vertex not in adj[neighbor]:
                raise AssertionError("asymmetric adjacency list")
    edge_count = sum(map(len, adj)) // 2
    seen = {0}
    stack = [0]
    while stack:
        vertex = stack.pop()
        for neighbor in adj[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    if len(seen) != n or edge_count != n - 1:
        raise AssertionError(
            f"constructor did not produce a tree: n={n}, edges={edge_count}, "
            f"connected={len(seen)}"
        )
    histogram: dict[str, int] = {}
    for neighbors in adj:
        key = str(len(neighbors))
        histogram[key] = histogram.get(key, 0) + 1
    return {
        "order": n,
        "edges": edge_count,
        "edge_list": [
            [vertex, neighbor]
            for vertex, neighbors in enumerate(adj)
            for neighbor in neighbors
            if vertex < neighbor
        ],
        "max_degree": max(map(len, adj), default=0),
        "degree_histogram": histogram,
    }


def balanced_tree_from_word(word: Sequence[int]) -> list[list[int]]:
    """Build a balanced rooted tree from bottom-to-top branching factors.

    Start with one leaf.  Reading ``word`` left to right, add a new root
    with ``branching`` identical copies of the current rooted tree below it.
    """
    adj: list[list[int]] = [[]]
    root = 0
    for branching in word:
        if branching not in (1, 2):
            raise ValueError("literature branching words use only 1 and 2")
        old_adj = adj
        old_root = root
        old_n = len(old_adj)
        adj = [[] for _ in range(1 + branching * old_n)]
        root = 0
        for copy_index in range(branching):
            offset = 1 + copy_index * old_n
            for vertex, neighbors in enumerate(old_adj):
                for neighbor in neighbors:
                    if vertex < neighbor:
                        u = offset + vertex
                        v = offset + neighbor
                        adj[u].append(v)
                        adj[v].append(u)
            child_root = offset + old_root
            adj[root].append(child_root)
            adj[child_root].append(root)
    for neighbors in adj:
        neighbors.sort()
    return adj


def compressed_state_from_word(
    word: Sequence[int],
) -> tuple[list[int], list[int]]:
    """Exact homogeneous recurrence without constructing the adjacency list."""
    excluded = [1]
    included = [0, 1]
    for branching in word:
        old_excluded = excluded
        old_included = included
        total = poly_add(old_excluded, old_included)
        excluded = poly_pow(total, branching)
        included = poly_shift(poly_pow(old_excluded, branching))
    return excluded, included


def jerrum_patel_word(subdivisions_k: int, height: int) -> tuple[int, ...]:
    """Return the exact Lemma-17 branching word, bottom to top."""
    if subdivisions_k < 0 or height < 0:
        raise ValueError("k and height must be nonnegative")
    period = (1,) * (2 * subdivisions_k) + (2,)
    return period * height


def phase_truncated_words(
    subdivisions_k: int, repetitions: int
) -> Iterable[tuple[str, tuple[int, ...]]]:
    """Enumerate all endpoint phases of a finite periodic branching word."""
    if subdivisions_k <= 0:
        return
    period = (1,) * (2 * subdivisions_k) + (2,)
    full = period * repetitions
    period_length = len(period)
    for start_trim in range(period_length):
        for end_trim in range(period_length):
            stop = len(full) - end_trim if end_trim else len(full)
            word = full[start_trim:stop]
            if word and 2 in word:
                label = (
                    f"jp-phase-k{subdivisions_k}-r{repetitions}"
                    f"-start{start_trim}-end{end_trim}"
                )
                yield label, word


def s2_tree(arms: int) -> tuple[list[list[int]], int]:
    """Return S_{2,arms}, rooted at its center."""
    if arms < 0:
        raise ValueError("number of arms must be nonnegative")
    adj = [[] for _ in range(1 + 2 * arms)]
    for arm in range(arms):
        middle = 1 + 2 * arm
        leaf = middle + 1
        adj[0].append(middle)
        adj[middle].append(0)
        adj[middle].append(leaf)
        adj[leaf].append(middle)
    return adj, 0


def t1_tree(arms: int) -> tuple[list[list[int]], int]:
    """Return T_{1,arms}, rooted at the endpoint above S_{2,arms}."""
    star, star_root = s2_tree(arms)
    adj = [[] for _ in range(len(star) + 1)]
    for vertex, neighbors in enumerate(star):
        for neighbor in neighbors:
            if vertex < neighbor:
                u = vertex + 1
                v = neighbor + 1
                adj[u].append(v)
                adj[v].append(u)
    adj[0].append(star_root + 1)
    adj[star_root + 1].append(0)
    return adj, 0


def pattern_repeat(
    base_adj: Sequence[Sequence[int]],
    base_root: int,
    pendant_adj: Sequence[Sequence[int]],
    pendant_root: int,
    width: int,
    depth: int,
) -> list[list[int]]:
    """Build (G^v:H^w)_width^(depth) from Definition 2.1."""
    if width < 0 or depth < 0:
        raise ValueError("pattern width and depth must be nonnegative")
    adj = [list(neighbors) for neighbors in base_adj]
    for _ in range(depth):
        scale_root = len(adj)
        adj.append([])
        adj[base_root].append(scale_root)
        adj[scale_root].append(base_root)
        for _copy in range(width):
            offset = len(adj)
            adj.extend([] for _ in pendant_adj)
            for vertex, neighbors in enumerate(pendant_adj):
                for neighbor in neighbors:
                    if vertex < neighbor:
                        u = offset + vertex
                        v = offset + neighbor
                        adj[u].append(v)
                        adj[v].append(u)
            copied_root = offset + pendant_root
            adj[scale_root].append(copied_root)
            adj[copied_root].append(scale_root)
    for neighbors in adj:
        neighbors.sort()
    return adj


P1 = [1, 1]
P2 = [1, 2]


def scale_poly(arms: int) -> list[int]:
    """Return I(S_{2,arms};x) = (1+2x)^arms + x(1+x)^arms."""
    return poly_add(poly_pow(P2, arms), poly_shift(poly_pow(P1, arms)))


def circle_family_poly(k: int, offset: int) -> list[int]:
    """Return the unstarred Bautista 3,k,k+offset polynomial."""
    return poly_add(
        poly_mul(
            scale_poly(3),
            poly_mul(scale_poly(k), scale_poly(k + offset)),
        ),
        poly_shift(poly_pow(P2, 3 + 2 * k + offset)),
    )


def consecutive_family_poly(
    ell: int, pendant_arms: int, width: int, depth: int
) -> list[int]:
    """Return U_m from Theorem 4.11 / Equation (12)."""
    s_ell = scale_poly(ell)
    s_n = scale_poly(pendant_arms)
    t_kn = poly_add(
        poly_pow(s_n, width),
        poly_shift(poly_pow(P2, width * pendant_arms)),
    )
    return poly_add(
        poly_mul(s_ell, poly_pow(t_kn, depth)),
        poly_shift(
            poly_mul(
                poly_pow(P2, ell),
                poly_pow(s_n, width * depth),
            )
        ),
    )


def circle_recurrence_residual(k: int, offset: int) -> list[int]:
    """Return the third-order recurrence residual for 3,k,k+offset."""
    if k < 3:
        raise ValueError("circle recurrence check requires k >= 3")
    roots = [poly_pow(P2, 2), poly_mul(P1, P2), poly_pow(P1, 2)]
    e1 = poly_add(poly_add(roots[0], roots[1]), roots[2])
    e2 = poly_add(
        poly_add(poly_mul(roots[0], roots[1]), poly_mul(roots[0], roots[2])),
        poly_mul(roots[1], roots[2]),
    )
    e3 = poly_mul(poly_mul(roots[0], roots[1]), roots[2])
    return poly_add(
        poly_sub(
            circle_family_poly(k, offset),
            poly_mul(e1, circle_family_poly(k - 1, offset)),
        ),
        poly_sub(
            poly_mul(e2, circle_family_poly(k - 2, offset)),
            poly_mul(e3, circle_family_poly(k - 3, offset)),
        ),
    )


def consecutive_recurrence_residual(
    ell: int, pendant_arms: int, width: int, depth: int
) -> list[int]:
    """Return the second-order recurrence residual for U_m."""
    if depth < 2:
        raise ValueError("consecutive recurrence check requires depth >= 2")
    s_n = scale_poly(pendant_arms)
    s_power = poly_pow(s_n, width)
    t_kn = poly_add(
        s_power,
        poly_shift(poly_pow(P2, width * pendant_arms)),
    )
    coefficient_1 = poly_add(t_kn, s_power)
    coefficient_2 = poly_mul(t_kn, s_power)
    return poly_add(
        poly_sub(
            consecutive_family_poly(ell, pendant_arms, width, depth),
            poly_mul(
                coefficient_1,
                consecutive_family_poly(
                    ell, pendant_arms, width, depth - 1
                ),
            ),
        ),
        poly_mul(
            coefficient_2,
            consecutive_family_poly(ell, pendant_arms, width, depth - 2),
        ),
    )


@dataclass(frozen=True)
class StressCase:
    """One exact tree-polynomial case and its source-specific checks."""

    label: str
    lane: str
    source: str
    params: dict[str, int]
    adj: list[list[int]]
    formula_poly: list[int]
    expected_degree: int | None = None
    expected_breaks: tuple[int, ...] | None = None
    locus: str | None = None


def make_jp_case(
    label: str,
    lane: str,
    subdivisions_k: int,
    word: Sequence[int],
) -> StressCase:
    """Construct and cross-label a Jerrum--Patel balanced tree case."""
    adj = balanced_tree_from_word(word)
    excluded, included = compressed_state_from_word(word)
    formula = poly_add(excluded, included)
    return StressCase(
        label=label,
        lane=lane,
        source="arXiv:2510.01466v2, Section 6 / Lemma 17 proof",
        params={
            "subdivisions_k": subdivisions_k,
            "word_length": len(word),
            "branch_vertices": sum(value == 2 for value in word),
        },
        adj=adj,
        formula_poly=formula,
        locus="jerrum_patel_positive_axis",
    )


def build_cases(include_phase: bool) -> list[StressCase]:
    """Build the curated, tractable literature stress corpus."""
    cases: list[StressCase] = []

    for subdivisions_k, heights in ((0, range(2, 9)), (1, range(2, 6))):
        for height in heights:
            word = jerrum_patel_word(subdivisions_k, height)
            expected_order = (
                (2 * subdivisions_k + 1) * 2 ** (height + 1)
                - (4 * subdivisions_k + 1)
            )
            case = make_jp_case(
                label=f"jp-exact-k{subdivisions_k}-h{height}",
                lane="jerrum_patel_stated_exact_family",
                subdivisions_k=subdivisions_k,
                word=word,
            )
            if len(case.adj) != expected_order:
                raise AssertionError("Jerrum--Patel order formula mismatch")
            cases.append(case)

    if include_phase:
        for label, word in phase_truncated_words(1, repetitions=5):
            cases.append(
                make_jp_case(
                    label=label,
                    lane="jerrum_patel_proof_phase_truncation",
                    subdivisions_k=1,
                    word=word,
                )
            )

    for offset in (0, 1, 2):
        for k in (4, 8, 12, 16, 24, 32):
            adj = make_li_tree(k, k + offset)
            cases.append(
                StressCase(
                    label=f"bautista-circle-3-{k}-{k + offset}",
                    lane="bautista_p2_bkw_circle",
                    source="arXiv:2603.14204, Theorem 4.1",
                    params={"k": k, "offset": offset},
                    adj=adj,
                    formula_poly=circle_family_poly(k, offset),
                    expected_degree=2 * k + offset + 6,
                    locus="bkw_circle",
                )
            )

    consecutive_specs = (
        (2, 4, 2, 9, (92,)),
        (3, 4, 2, 16, (162, 163)),
        (7, 5, 2, 13, (161, 162, 163)),
    )
    for ell, pendant_arms, width, depth, breaks in consecutive_specs:
        base_adj, base_root = t1_tree(ell)
        pendant_adj, pendant_root = s2_tree(pendant_arms)
        adj = pattern_repeat(
            base_adj,
            base_root,
            pendant_adj,
            pendant_root,
            width,
            depth,
        )
        cases.append(
            StressCase(
                label=(
                    f"bautista-consecutive-l{ell}-n{pendant_arms}"
                    f"-k{width}-m{depth}"
                ),
                lane="bautista_s2_consecutive_noncircle",
                source="arXiv:2603.14204, Theorem 4.11 and Corollaries 4.12--4.14",
                params={
                    "ell": ell,
                    "pendant_arms": pendant_arms,
                    "width": width,
                    "depth": depth,
                },
                adj=adj,
                formula_poly=consecutive_family_poly(
                    ell, pendant_arms, width, depth
                ),
                expected_degree=width * (pendant_arms + 1) * depth + ell + 1,
                expected_breaks=breaks,
                locus="consecutive_equimodular",
            )
        )

    return cases


def arb_interval(value: Any) -> dict[str, Any]:
    """Serialize an Arb interval without using it for decisions."""
    return {
        "lower": float(value.lower()),
        "upper": float(value.upper()),
        "midpoint": float(value.mid()),
        "ball": str(value),
    }


def root_record(
    z: Any,
    multiplicity: int,
    beta_modulus: Any,
    arb: Any,
) -> dict[str, Any]:
    """Serialize one upper-half-plane root ball and derived intervals."""
    modulus = abs(z)
    argument = z.arg()
    positive_angle = argument
    negative_axis_deviation = arb.pi() - argument
    return {
        "multiplicity": multiplicity,
        "real": arb_interval(z.real),
        "imag": arb_interval(z.imag),
        "modulus": arb_interval(modulus),
        "modulus_ratio": arb_interval(modulus / beta_modulus),
        "positive_axis_angle": arb_interval(positive_angle),
        "negative_axis_deviation": arb_interval(negative_axis_deviation),
    }


def certified_unique_minimum(values: Sequence[Any]) -> tuple[int, bool]:
    """Choose a midpoint minimum and certify whether its ball is separated."""
    if not values:
        raise ValueError("cannot minimize an empty interval sequence")
    index = min(range(len(values)), key=lambda item: float(values[item].mid()))
    certified = all(
        other == index or values[index].upper() < values[other].lower()
        for other in range(len(values))
    )
    return index, certified


def neutral_target_ball(subdivisions_k: int, fmpz_poly: Any) -> Any | None:
    """Return a certified positive neutral-activity ball for k=0,1,2."""
    coefficient_map = {
        0: [-4, 1],
        1: [-4, -36, -108, -120, -40, 1],
        2: [
            -4,
            -68,
            -476,
            -1776,
            -3820,
            -4776,
            -3312,
            -1120,
            -132,
            1,
        ],
    }
    coefficients = coefficient_map.get(subdivisions_k)
    if coefficients is None:
        return None
    roots = fmpz_poly(coefficients).complex_roots()
    positive = [
        z.real
        for z, _multiplicity in roots
        if z.imag.contains(0) and z.real.lower() > 0
    ]
    if len(positive) != 1:
        raise AssertionError("neutral polynomial did not isolate one positive root")
    return positive[0]


def spectrum_certificate(
    poly: Sequence[int],
    locus: str | None,
    params: dict[str, int],
    initial_precision: int,
) -> dict[str, Any]:
    """Isolate roots and certify the minimum-modulus root and root metrics."""
    try:
        from flint import arb, ctx, fmpz_poly
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise RuntimeError(
            "python-flint is required; use an external venv, e.g. "
            "python3 -m venv /tmp/erdos993-flint-venv && "
            "/tmp/erdos993-flint-venv/bin/pip install python-flint"
        ) from exc

    roots = None
    beta_index = -1
    beta_separated = False
    precision = initial_precision
    while precision <= 4 * initial_precision:
        ctx.prec = precision
        candidate_roots = fmpz_poly(list(poly)).complex_roots()
        moduli = [abs(z) for z, _multiplicity in candidate_roots]
        beta_index, beta_separated = certified_unique_minimum(moduli)
        roots = candidate_roots
        if beta_separated:
            break
        precision *= 2
    if roots is None or not beta_separated:
        raise AssertionError("minimum-modulus root was not separated by Arb balls")

    polynomial = fmpz_poly(list(poly))
    if sum(multiplicity for _z, multiplicity in roots) != polynomial.degree():
        raise AssertionError("root multiplicities do not sum to polynomial degree")
    if any(abs(z).contains(0) for z, _multiplicity in roots):
        raise AssertionError("a root enclosure contains zero")

    beta, beta_multiplicity = roots[beta_index]
    if beta_multiplicity != 1 or not beta.imag.contains(0) or beta.real.upper() >= 0:
        raise AssertionError(
            "minimum-modulus root is not certified simple negative real"
        )
    beta_modulus = abs(beta)

    upper_roots = [
        (z, multiplicity)
        for z, multiplicity in roots
        if z.imag.lower() > 0
    ]
    serialized = [
        root_record(z, multiplicity, beta_modulus, arb)
        for z, multiplicity in upper_roots
    ]
    ratios = [abs(z) / beta_modulus for z, _multiplicity in upper_roots]
    positive_angles = [z.arg() for z, _multiplicity in upper_roots]
    negative_deviations = [
        arb.pi() - z.arg() for z, _multiplicity in upper_roots
    ]

    result: dict[str, Any] = {
        "precision_bits": precision,
        "degree": polynomial.degree(),
        "squarefree": polynomial.gcd(polynomial.derivative()).degree() == 0,
        "root_ball_count": len(roots),
        "multiplicity_sum": sum(multiplicity for _z, multiplicity in roots),
        "beta": {
            "multiplicity": beta_multiplicity,
            "real": arb_interval(beta.real),
            "modulus": arb_interval(beta_modulus),
            "certified_unique_minimum": beta_separated,
        },
        "nonreal_conjugate_pairs": sum(
            multiplicity for _z, multiplicity in upper_roots
        ),
    }

    if upper_roots:
        ratio_index, ratio_unique = certified_unique_minimum(ratios)
        angle_index, angle_unique = certified_unique_minimum(positive_angles)
        result["minimum_nonreal_ratio"] = {
            "certified_unique": ratio_unique,
            "root": serialized[ratio_index],
        }
        result["minimum_positive_axis_angle"] = {
            "certified_unique": angle_unique,
            "root": serialized[angle_index],
        }

        frontier: dict[str, Any] = {}
        for tau in TAUS:
            inside = [
                index
                for index, ratio in enumerate(ratios)
                if ratio.upper() <= tau
            ]
            ambiguous = [
                index
                for index, ratio in enumerate(ratios)
                if ratio.lower() <= tau < ratio.upper()
            ]
            if inside:
                lower = max(
                    float(negative_deviations[index].lower())
                    for index in inside
                )
                upper = max(
                    float(negative_deviations[index].upper())
                    for index in inside
                )
                max_deviation = {"lower": lower, "upper": upper}
            else:
                max_deviation = None
            frontier[str(tau)] = {
                "certified_inside_pairs": len(inside),
                "threshold_ambiguous_pairs": len(ambiguous),
                "max_negative_axis_deviation": max_deviation,
            }
        result["dominance_angle_frontier"] = frontier

    if locus == "bkw_circle" and upper_roots:
        one_third = arb(1) / 3
        residuals = [
            abs(abs(z + one_third) - one_third)
            for z, _multiplicity in upper_roots
        ]
        residual_index, residual_unique = certified_unique_minimum(residuals)
        result["bkw_circle_closest_root"] = {
            "certified_unique": residual_unique,
            "residual": arb_interval(residuals[residual_index]),
            "root": serialized[residual_index],
        }

    if locus == "consecutive_equimodular" and upper_roots:
        pendant_arms = params["pendant_arms"]
        width = params["width"]
        s_n = fmpz_poly(scale_poly(pendant_arms))
        t_kn = fmpz_poly(
            poly_add(
                poly_pow(scale_poly(pendant_arms), width),
                poly_shift(poly_pow(P2, width * pendant_arms)),
            )
        )
        residuals = [
            abs(abs(t_kn(z)) - abs(s_n(z) ** width))
            for z, _multiplicity in upper_roots
        ]
        residual_index, residual_unique = certified_unique_minimum(residuals)
        result["consecutive_locus_closest_root"] = {
            "certified_unique": residual_unique,
            "residual": arb_interval(residuals[residual_index]),
            "root": serialized[residual_index],
        }

    if locus == "jerrum_patel_positive_axis" and upper_roots:
        subdivisions_k = params["subdivisions_k"]
        target = neutral_target_ball(subdivisions_k, fmpz_poly)
        if target is not None:
            distances = [
                abs(z - target) for z, _multiplicity in upper_roots
            ]
            target_index, target_unique = certified_unique_minimum(distances)
            result["neutral_target"] = arb_interval(target)
            result["neutral_target_closest_root"] = {
                "certified_unique": target_unique,
                "distance": arb_interval(distances[target_index]),
                "root": serialized[target_index],
            }

    return result


def audit_case(case: StressCase, precision: int) -> dict[str, Any]:
    """Run exact graph, formula, recurrence, shape, and root checks."""
    tree = validate_tree(case.adj)
    dp_poly = independence_poly(len(case.adj), case.adj)
    if dp_poly != case.formula_poly:
        raise AssertionError(f"formula/DP mismatch for {case.label}")
    if dp_poly[0] != 1 or dp_poly[1] != len(case.adj):
        raise AssertionError(f"constant/linear coefficient mismatch for {case.label}")
    degree = len(dp_poly) - 1
    if case.expected_degree is not None and degree != case.expected_degree:
        raise AssertionError(f"degree mismatch for {case.label}")
    failures = lc_failure_indices(dp_poly)
    if case.expected_breaks is not None and failures != list(case.expected_breaks):
        raise AssertionError(
            f"LC failures for {case.label}: expected {case.expected_breaks}, "
            f"found {failures}"
        )
    if not is_unimodal(dp_poly):
        raise AssertionError(f"unexpected non-unimodal literature tree {case.label}")

    if case.lane == "bautista_p2_bkw_circle":
        recurrence = circle_recurrence_residual(
            case.params["k"], case.params["offset"]
        )
    elif case.lane == "bautista_s2_consecutive_noncircle":
        recurrence = consecutive_recurrence_residual(
            case.params["ell"],
            case.params["pendant_arms"],
            case.params["width"],
            case.params["depth"],
        )
    else:
        recurrence = [0]
    if poly_trim(recurrence) != [0]:
        raise AssertionError(f"recurrence residual for {case.label} is nonzero")

    started = time.time()
    spectrum = spectrum_certificate(
        dp_poly,
        case.locus,
        case.params,
        initial_precision=precision,
    )
    return {
        "label": case.label,
        "lane": case.lane,
        "source": case.source,
        "params": case.params,
        "tree": tree,
        "polynomial": {
            "degree": degree,
            "constant": dp_poly[0],
            "linear": dp_poly[1],
            "coefficients": dp_poly,
            "formula_matches_tree_dp": True,
            "recurrence_residual_zero": True,
            "unimodal": True,
            "log_concavity_failure_indices": failures,
        },
        "spectrum": spectrum,
        "root_seconds": time.time() - started,
    }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"JSON result path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=192,
        help="initial Arb precision in bits (default: 192)",
    )
    parser.add_argument(
        "--skip-phase",
        action="store_true",
        help="skip the Jerrum--Patel proof-supported phase truncations",
    )
    parser.add_argument(
        "--labels",
        nargs="*",
        help="run only exact case labels listed here",
    )
    return parser.parse_args()


def main() -> None:
    """Run the curated certified literature-family audit."""
    args = parse_args()
    if args.precision < 64:
        raise ValueError("use at least 64 bits of Arb precision")
    cases = build_cases(include_phase=not args.skip_phase)
    if args.labels:
        requested = set(args.labels)
        cases = [case for case in cases if case.label in requested]
        missing = requested - {case.label for case in cases}
        if missing:
            raise ValueError(f"unknown case labels: {sorted(missing)}")

    report: dict[str, Any] = {
        "schema_version": 1,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "sources": {
            "jerrum_patel": "arXiv:2510.01466v2, Section 6 / Lemma 17",
            "bautista_et_al": (
                "arXiv:2603.14204v1, Definition 2.1, "
                "Theorems 2.3, 4.1, 4.11"
            ),
        },
        "rigor": {
            "polynomials": "exact integer DP plus exact closed-form/recurrence replay",
            "root_isolation": "python-flint fmpz_poly.complex_roots (Arb)",
            "root_ranking": (
                "Arb modulus intervals; precision doubled until beta separates"
            ),
            "thresholds": (
                "classified only when the full Arb ratio interval lies inside"
            ),
            "floats": "JSON reporting only; not used for certification decisions",
        },
        "jerrum_patel_scope_warning": (
            "Lemma 17 states the exact uniformly subdivided binary family, "
            "but the cited normal-family step directly yields phase-truncated "
            "periodic branching words. Both lanes are audited separately."
        ),
        "environment": {
            "python": sys.version,
            "python_flint": version("python-flint"),
            "platform": platform.platform(),
            "initial_precision_bits": args.precision,
        },
        "taus": list(TAUS),
        "cases": [],
    }

    started = time.time()
    for index, case in enumerate(cases, start=1):
        print(
            f"[{index:02d}/{len(cases):02d}] {case.label} "
            f"n={len(case.adj)} degree={len(case.formula_poly) - 1}",
            flush=True,
        )
        report["cases"].append(audit_case(case, args.precision))
    report["elapsed_seconds"] = time.time() - started

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"wrote {args.output} ({len(report['cases'])} cases, "
        f"{report['elapsed_seconds']:.1f}s)",
        flush=True,
    )


if __name__ == "__main__":
    main()
