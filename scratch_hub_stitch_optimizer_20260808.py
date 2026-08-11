#!/usr/bin/env python3
"""Turn archived forest near-valleys into trees by optimized hub stitching.

For rooted tree components ``(C_i, r_i)``, joining a new hub to every ``r_i``
gives the exact tree polynomial

    product_i I(C_i; x) + x product_i I(C_i-r_i; x).

If the hub also has ``ell`` pendant leaves, the first term gains a factor
``(1+x)^ell``.  The archived forest semigroup search already found products
with a shallow dip followed by a rebound.  This script chooses every
attachment vertex and the number of hub leaves to make the positive second
term favor the rebound over the dip.

Long-double arithmetic ranks coordinate-search candidates only.  Every
finalist is replayed with Python integers, and any counterexample is rebuilt
as an explicit tree and independently checked by the generic tree DP.
"""

from __future__ import annotations

import argparse
import json
import random
import re
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from indpoly import _polyadd, _polymul, independence_poly, is_unimodal
from scratch_existential_split_mode_20260807 import tree_vertex_splits
from scripts.valley_search import bouquet_adj, bouquet_poly, valley_score


@dataclass
class Component:
    label: str
    order: int
    adj: list[list[int]]
    full: list[int]
    attachments: list[tuple[int, list[int]]]
    float_attachments: list[np.ndarray]


def parse_bush_label(label: str) -> list[tuple[int, ...]]:
    """Recover the deterministic bouquet grammar from a stored seed label."""
    body = label.split("bush(", 1)[1].split(")|", 1)[0]
    gadgets: list[tuple[int, ...]] = []
    for token in body.split(","):
        match = re.fullmatch(r"([us])(\d+)p(\d+)", token)
        if match is None:
            raise ValueError(f"unrecognized bush token {token!r}")
        kind, children_text, pendant_text = match.groups()
        children = int(children_text)
        pendant = int(pendant_text)
        legs = [pendant + 1] * children
        if kind == "s":
            legs[0] += 1
        gadgets.append(tuple(sorted(legs)))
    return gadgets


def load_component(label: str) -> Component:
    gadgets = parse_bush_label(label)
    order, adj = bouquet_adj(gadgets)
    full = bouquet_poly(gadgets)
    replay, splits = tree_vertex_splits(adj)
    if list(replay) != full:
        raise AssertionError(f"component replay failed for {label}")

    # Different vertices often have identical deletion polynomials.  Retain a
    # deterministic representative of each exact option.
    distinct: dict[tuple[int, ...], int] = {}
    for split in splits:
        distinct.setdefault(split.a_poly, split.vertex)
    attachments = [
        (vertex, list(poly))
        for poly, vertex in distinct.items()
    ]
    scale = np.longdouble(sum(full))
    float_attachments = [
        np.asarray(poly, dtype=np.longdouble) / scale
        for _, poly in attachments
    ]
    return Component(
        label=label,
        order=order,
        adj=adj,
        full=full,
        attachments=attachments,
        float_attachments=float_attachments,
    )


def normalized_forest(components: list[Component]) -> np.ndarray:
    out = np.asarray([1.0], dtype=np.longdouble)
    for component in components:
        profile = np.asarray(component.full, dtype=np.longdouble)
        profile /= np.longdouble(sum(component.full))
        out = np.convolve(out, profile)
    return out


def binomial_poly(leaves: int) -> list[int]:
    out = [1]
    for _ in range(leaves):
        out = _polymul(out, [1, 1])
    return out


def float_product(
    components: list[Component], choices: list[int], skip: int | None = None,
) -> np.ndarray:
    out = np.asarray([1.0], dtype=np.longdouble)
    for index, (component, choice) in enumerate(zip(components, choices)):
        if index == skip:
            continue
        out = np.convolve(out, component.float_attachments[choice])
    return out


def float_threshold_rebound(
    profile: np.ndarray, threshold: float, min_separation: int,
) -> tuple[float, int, int]:
    """Return a later-peak/dip ratio only after a nontrivial prior drop."""
    if len(profile) < 3:
        return 0.0, -1, -1
    prefix = np.maximum.accumulate(profile)
    suffix = np.maximum.accumulate(profile[::-1])[::-1]
    best_ratio = np.longdouble(0)
    best_b = -1
    for b in range(1, len(profile) - min_separation):
        if profile[b] == 0:
            continue
        if prefix[b - 1] < (1.0 + threshold) * profile[b]:
            continue
        ratio = suffix[b + min_separation] / profile[b]
        if ratio > best_ratio:
            best_ratio = ratio
            best_b = b
    if best_b < 0:
        return 0.0, -1, -1
    target = suffix[best_b + min_separation]
    rise = next(
        c for c in range(best_b + min_separation, len(profile))
        if profile[c] == target
    )
    return float(best_ratio), best_b, rise


def exact_threshold_rebound(
    poly: list[int], threshold_numerator: int, threshold_denominator: int,
    min_separation: int,
) -> dict[str, object]:
    prefix = [0] * len(poly)
    running = poly[0]
    for index in range(1, len(poly)):
        prefix[index] = running
        running = max(running, poly[index])
    suffix = [0] * len(poly)
    suffix[-1] = poly[-1]
    for index in range(len(poly) - 2, -1, -1):
        suffix[index] = max(poly[index], suffix[index + 1])

    best_num, best_den = 0, 1
    best_b = best_c = -1
    for b in range(1, len(poly) - min_separation):
        if (
            prefix[b] * threshold_denominator
            < (threshold_denominator + threshold_numerator) * poly[b]
        ):
            continue
        numerator = suffix[b + min_separation]
        denominator = poly[b]
        if numerator * best_den > best_num * denominator:
            best_num, best_den = numerator, denominator
            best_b = b
            best_c = next(
                c for c in range(b + min_separation, len(poly))
                if poly[c] == numerator
            )
    return {
        "ratio_numerator": best_num,
        "ratio_denominator": best_den,
        "ratio": best_num / best_den,
        "dip": best_b,
        "rise": best_c,
        "witness": best_num > best_den,
    }


def stitched_float(
    forest: np.ndarray, components: list[Component], choices: list[int], leaves: int,
) -> np.ndarray:
    first = np.convolve(
        forest,
        np.asarray(binomial_poly(leaves), dtype=np.longdouble),
    )
    second = float_product(components, choices)
    size = max(len(first), len(second) + 1)
    total = np.zeros(size, dtype=np.longdouble)
    total[:len(first)] += first
    total[1:len(second) + 1] += second
    return total


def coordinate_optimize(
    forest: np.ndarray,
    components: list[Component],
    choices: list[int],
    leaves: int,
    max_sweeps: int,
    threshold: float,
    min_separation: int,
) -> tuple[list[int], tuple[float, int, int]]:
    """Coordinate ascent over exact deletion-state choices, ranked in float."""
    choices = choices.copy()
    current_profile = stitched_float(forest, components, choices, leaves)
    current = float_threshold_rebound(
        current_profile, threshold, min_separation,
    )
    for _ in range(max_sweeps):
        changed = False
        for index, component in enumerate(components):
            other = float_product(components, choices, skip=index)
            first = np.convolve(
                forest,
                np.asarray(binomial_poly(leaves), dtype=np.longdouble),
            )
            best_choice = choices[index]
            best_score = current
            for option, option_poly in enumerate(component.float_attachments):
                second = np.convolve(other, option_poly)
                size = max(len(first), len(second) + 1)
                profile = np.zeros(size, dtype=np.longdouble)
                profile[:len(first)] += first
                profile[1:len(second) + 1] += second
                score = float_threshold_rebound(
                    profile, threshold, min_separation,
                )
                if score[0] > best_score[0] + 1e-18:
                    best_choice = option
                    best_score = score
            if best_choice != choices[index]:
                choices[index] = best_choice
                current = best_score
                changed = True
            else:
                # Earlier coordinates may have changed the context.
                current = float_threshold_rebound(
                    stitched_float(forest, components, choices, leaves),
                    threshold,
                    min_separation,
                )
        if not changed:
            break
    return choices, current


def exact_stitched_poly(
    components: list[Component], choices: list[int], leaves: int,
) -> list[int]:
    first = [1]
    second = [1]
    for component, choice in zip(components, choices):
        first = _polymul(first, component.full)
        second = _polymul(second, component.attachments[choice][1])
    first = _polymul(first, binomial_poly(leaves))
    return _polyadd(first, [0] + second)


def build_stitched_tree(
    components: list[Component], choices: list[int], leaves: int,
) -> list[list[int]]:
    total_order = 1 + leaves + sum(component.order for component in components)
    adj: list[list[int]] = [[] for _ in range(total_order)]

    def edge(left: int, right: int) -> None:
        adj[left].append(right)
        adj[right].append(left)

    offset = 1
    for _ in range(leaves):
        edge(0, offset)
        offset += 1
    for component, choice in zip(components, choices):
        attachment = component.attachments[choice][0]
        edge(0, offset + attachment)
        for left, neighbors in enumerate(component.adj):
            for right in neighbors:
                if left < right:
                    edge(offset + left, offset + right)
        offset += component.order
    if offset != total_order:
        raise AssertionError("stitched-tree order mismatch")
    return adj


def encode_exact(
    source: str,
    row: dict[str, object],
    components: list[Component],
    choices: list[int],
    leaves: int,
    poly: list[int],
    include_poly: bool,
    threshold_numerator: int,
    threshold_denominator: int,
    min_separation: int,
) -> dict[str, object]:
    score = valley_score(poly)
    return {
        "source": source,
        "source_factors": row["factors"],
        "component_labels": [component.label for component in components],
        "component_orders": [component.order for component in components],
        "attachment_vertices": [
            component.attachments[choice][0]
            for component, choice in zip(components, choices)
        ],
        "hub_leaves": leaves,
        "order": 1 + leaves + sum(component.order for component in components),
        "alpha": len(poly) - 1,
        "unimodal": is_unimodal(poly),
        "valley": score,
        "objective": exact_threshold_rebound(
            poly,
            threshold_numerator,
            threshold_denominator,
            min_separation,
        ),
        "polynomial": poly if include_poly else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--max-hub-leaves", type=int, default=8)
    parser.add_argument("--max-sweeps", type=int, default=8)
    parser.add_argument("--threshold", type=float, default=0.001)
    parser.add_argument("--min-separation", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260808)
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/hub_stitch_optimizer_20260808.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)
    exact_threshold = Fraction(str(args.threshold))
    threshold_numerator = exact_threshold.numerator
    threshold_denominator = exact_threshold.denominator
    started = time.time()

    sources = [
        Path("results/forest_valley_semigroup_20260807.json"),
        Path("results/forest_valley_later_slope_sep2_20260807.json"),
        Path("results/forest_valley_later_slope_sep5_20260807.json"),
        Path("results/forest_valley_later_slope_sep10_20260807.json"),
    ]
    rows: list[tuple[str, dict[str, object]]] = []
    seen: set[tuple[str, ...]] = set()
    for source in sources:
        data = json.loads(source.read_text())
        for row in data["top_exact"]:
            key = tuple(row["labels"])
            if key in seen:
                continue
            seen.add(key)
            rows.append((str(source), row))

    cache: dict[str, Component] = {}
    best: dict[str, object] | None = None
    counterexample: dict[str, object] | None = None
    evaluated_finalists = 0
    for row_index, (source, row) in enumerate(rows):
        labels = row["labels"]
        components = []
        for label in labels:
            component = cache.get(label)
            if component is None:
                component = load_component(label)
                cache[label] = component
            components.append(component)
        forest = normalized_forest(components)
        for leaves in range(args.max_hub_leaves + 1):
            starts = [[0] * len(components)]
            for _ in range(args.restarts - 1):
                starts.append([
                    rng.randrange(len(component.attachments))
                    for component in components
                ])
            for start_choices in starts:
                choices, _ = coordinate_optimize(
                    forest,
                    components,
                    start_choices,
                    leaves,
                    args.max_sweeps,
                    args.threshold,
                    args.min_separation,
                )
                poly = exact_stitched_poly(components, choices, leaves)
                evaluated_finalists += 1
                record = encode_exact(
                    source,
                    row,
                    components,
                    choices,
                    leaves,
                    poly,
                    False,
                    threshold_numerator,
                    threshold_denominator,
                    args.min_separation,
                )
                if best is None or record["objective"]["ratio"] > best["objective"]["ratio"]:
                    best = record
                if not is_unimodal(poly):
                    adj = build_stitched_tree(components, choices, leaves)
                    replay = independence_poly(len(adj), adj)
                    if replay != poly:
                        raise AssertionError("explicit stitched-tree replay failed")
                    counterexample = encode_exact(
                        source,
                        row,
                        components,
                        choices,
                        leaves,
                        poly,
                        True,
                        threshold_numerator,
                        threshold_denominator,
                        args.min_separation,
                    )
                    counterexample["edges"] = [
                        [left, right]
                        for left, neighbors in enumerate(adj)
                        for right in neighbors if left < right
                    ]
                    break
            if counterexample is not None:
                break
        print(json.dumps({
            "event": "row",
            "row": row_index,
            "rows": len(rows),
            "best_ratio": best["objective"]["ratio"] if best else None,
            "counterexample_found": counterexample is not None,
            "elapsed_seconds": time.time() - started,
        }), flush=True)
        if counterexample is not None:
            break

    payload = {
        "claim_scope": "exact optimization finalists; floating point used only for ranking",
        "identity": "I(T;x)=(1+x)^ell product I(C_i;x)+x product I(C_i-r_i;x)",
        "seed": args.seed,
        "source_rows": len(rows),
        "distinct_components_loaded": len(cache),
        "restarts": args.restarts,
        "max_hub_leaves": args.max_hub_leaves,
        "threshold": args.threshold,
        "min_separation": args.min_separation,
        "evaluated_exact_finalists": evaluated_finalists,
        "counterexample_found": counterexample is not None,
        "counterexample": counterexample,
        "best": best,
        "elapsed_seconds": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "event": "done",
        "counterexample_found": counterexample is not None,
        "best_ratio": best["objective"]["ratio"] if best else None,
        "out": str(args.out),
        "elapsed_seconds": time.time() - started,
    }), flush=True)
    raise SystemExit(1 if counterexample is not None else 0)


if __name__ == "__main__":
    main()
