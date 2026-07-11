#!/usr/bin/env python3
"""Audit variation-diminishing candidates in the maximum-set expansion.

For a maximum independent set M and H=T-M, aggregate

    I_T(x) = sum_{X in I(H)} x^q (1+x)^(alpha-b),
    q=|X|, b=|N_M(X)|, delta=2q-b.

The script checks exact reconstruction and falsifies/tests three increasingly
coarse order claims: monotonicity of delta under inclusion, TP2 for individual
binomial columns ordered only by delta, and TP2 after aggregating all columns
with the same delta.  It is an audit driver, not a theorem certificate.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter
from itertools import combinations

from indpoly import independence_poly
from trees import trees


def independent(mask: int, adj: list[list[int]]) -> bool:
    return all(
        not (mask >> u & 1 and mask >> v & 1)
        for u, neighbors in enumerate(adj)
        for v in neighbors
        if u < v
    )


def maximum_sets(adj: list[list[int]]) -> list[int]:
    candidates = [mask for mask in range(1 << len(adj)) if independent(mask, adj)]
    alpha = max(mask.bit_count() for mask in candidates)
    return [mask for mask in candidates if mask.bit_count() == alpha]


def hall_counts(adj: list[list[int]], maximum_mask: int):
    m_vertices = [v for v in range(len(adj)) if maximum_mask >> v & 1]
    h_vertices = [v for v in range(len(adj)) if not (maximum_mask >> v & 1)]
    m_index = {vertex: index for index, vertex in enumerate(m_vertices)}
    neighbor_masks = []
    h_conflicts = []
    for vertex in h_vertices:
        neighbor_masks.append(sum(
            1 << m_index[neighbor]
            for neighbor in adj[vertex]
            if neighbor in m_index
        ))
        h_conflicts.append(sum(
            1 << h_vertices.index(neighbor)
            for neighbor in adj[vertex]
            if neighbor in h_vertices
        ))
    counts: Counter[tuple[int, int]] = Counter()
    independent_subsets = []
    for mask in range(1 << len(h_vertices)):
        if any(mask >> index & 1 and mask & h_conflicts[index]
               for index in range(len(h_vertices))):
            continue
        neighborhood = 0
        for index, neighbor_mask in enumerate(neighbor_masks):
            if mask >> index & 1:
                neighborhood |= neighbor_mask
        q, b = mask.bit_count(), neighborhood.bit_count()
        counts[(q, b)] += 1
        independent_subsets.append((mask, q, b))
    return len(m_vertices), counts, independent_subsets


def reconstruct(alpha: int, counts: Counter[tuple[int, int]]) -> list[int]:
    out = [0] * (alpha + 1)
    for (q, b), count in counts.items():
        for relative in range(alpha - b + 1):
            out[q + relative] += count * math.comb(alpha - b, relative)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def first_tp2_failure(columns: dict[int, list[int]]):
    labels = sorted(columns)
    if len(labels) < 2:
        return None
    rows = range(len(next(iter(columns.values()))))
    for left_index, left_label in enumerate(labels):
        for right_label in labels[left_index + 1 :]:
            left, right = columns[left_label], columns[right_label]
            for low in rows:
                for high in range(low + 1, len(left)):
                    determinant = left[low] * right[high] - right[low] * left[high]
                    if determinant < 0:
                        return {
                            "labels": (left_label, right_label),
                            "rows": (low, high),
                            "minor": determinant,
                            "entries": (left[low], right[low], left[high], right[high]),
                        }
    return None


def audit(adj: list[list[int]], maximum_mask: int):
    alpha, counts, subsets = hall_counts(adj, maximum_mask)
    assert reconstruct(alpha, counts) == independence_poly(len(adj), adj)
    assert all(b >= q for q, b in counts)

    # The center parameter delta=2q-b need not increase under inclusion.
    inclusion_failure = None
    for mask, q, b in subsets:
        delta = 2 * q - b
        for larger, q2, b2 in subsets:
            if mask != larger and mask & larger == mask and 2 * q2 - b2 < delta:
                inclusion_failure = (mask, (q, b, delta), larger, (q2, b2, 2 * q2 - b2))
                break
        if inclusion_failure:
            break

    # Individual columns: keep (q,b) separate, order only by delta.  Equal
    # deltas are deliberately omitted because delta alone cannot order them.
    parameter_columns = {}
    parameters_by_delta: dict[int, list[tuple[int, int]]] = {}
    for q, b in counts:
        parameters_by_delta.setdefault(2 * q - b, []).append((q, b))
    individual_failure = None
    parameters = list(counts)
    for left in parameters:
        for right in parameters:
            delta_left = 2 * left[0] - left[1]
            delta_right = 2 * right[0] - right[1]
            if delta_left >= delta_right:
                continue
            columns = {}
            for label, (q, b) in ((0, left), (1, right)):
                columns[label] = [
                    math.comb(alpha - b, rank - q)
                    if 0 <= rank - q <= alpha - b else 0
                    for rank in range(alpha + 1)
                ]
            failure = first_tp2_failure(columns)
            if failure:
                individual_failure = {
                    "parameters": (left, delta_left, right, delta_right),
                    **failure,
                }
                break
        if individual_failure:
            break

    # Aggregate all X with the same delta before testing TP2.
    aggregated: dict[int, list[int]] = {}
    for (q, b), count in counts.items():
        delta = 2 * q - b
        column = aggregated.setdefault(delta, [0] * (alpha + 1))
        for rank in range(alpha + 1):
            if 0 <= rank - q <= alpha - b:
                column[rank] += count * math.comb(alpha - b, rank - q)
    aggregate_failure = first_tp2_failure(aggregated)
    return {
        "alpha": alpha,
        "counts": dict(counts),
        "inclusion_failure": inclusion_failure,
        "individual_tp2_failure": individual_failure,
        "aggregate_tp2_failure": aggregate_failure,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=11)
    parser.add_argument("--all-maximum-sets", action="store_true")
    args = parser.parse_args()
    first = {"inclusion": None, "individual": None, "aggregate": None}
    audited = 0
    for n in range(1, args.max_n + 1):
        tree_count = 0
        for graph6, adj in trees(n):
            tree_count += 1
            masks = maximum_sets(adj)
            if not args.all_maximum_sets:
                masks = masks[:1]
            for maximum_mask in masks:
                row = audit(adj, maximum_mask)
                audited += 1
                for key, field in (
                    ("inclusion", "inclusion_failure"),
                    ("individual", "individual_tp2_failure"),
                    ("aggregate", "aggregate_tp2_failure"),
                ):
                    if first[key] is None and row[field] is not None:
                        first[key] = {
                            "n": n,
                            "graph6": graph6.decode() if isinstance(graph6, bytes) else graph6,
                            "maximum_mask": maximum_mask,
                            "degrees": sorted(map(len, adj), reverse=True),
                            "failure": row[field],
                            "counts": row["counts"],
                            "adj": adj,
                        }
        print({"n": n, "trees": tree_count, "audited": audited, "first": first}, flush=True)
    print({"complete": True, "audited": audited, "first": first}, flush=True)


if __name__ == "__main__":
    main()
