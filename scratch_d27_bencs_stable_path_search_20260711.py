#!/usr/bin/env python3
"""Search Bencs stable-path trees built from nonunimodal join graphs.

For a globally ordered graph ``G`` and root ``u``, Bencs' recursive stable-path
tree has one child for each neighbour ``u_i``; the corresponding recursive
state deletes ``u`` and the earlier neighbours ``u_1,...,u_(i-1)``.  The
identity

    I(T_(G,u)) / I(T_(G,u)-u) = I(G) / I(G-u)

lets us collect the tree polynomial as a product of irreducible factors of
independence polynomials of induced subgraphs, without materialising the often
huge tree.  Floating convolution is only a filter.  Any apparent valley must
be replayed exactly before it is reported as a counterexample.
"""

from __future__ import annotations

import argparse
import heapq
from collections import Counter
from functools import lru_cache
from math import comb

import numpy as np
from scipy.signal import fftconvolve
from sympy import Poly, symbols

from indpoly import _polyadd


def joined_cluster_graph(
    a: int, m: int, b: int, n: int, second_first: bool
) -> tuple[list[int], list[str]]:
    """Return bit-adjacency for ``m K_a join n K_b``.

    The source polynomial is ``(1+a*x)^m+(1+b*x)^n-1``.  Entire clusters
    from the second summand are ordered first when ``second_first`` is true.
    """
    block_specs = [*(('A', a) for _ in range(m)), *(('B', b) for _ in range(n))]
    if second_first:
        block_specs.sort(key=lambda item: item[0], reverse=True)
    blocks: list[tuple[str, list[int]]] = []
    types: list[str] = []
    for kind, size in block_specs:
        vertices = list(range(len(types), len(types) + size))
        types.extend([kind] * size)
        blocks.append((kind, vertices))

    adjacency = [0] * len(types)
    for _, vertices in blocks:
        for index, u in enumerate(vertices):
            for v in vertices[index + 1 :]:
                adjacency[u] |= 1 << v
                adjacency[v] |= 1 << u
    left = [v for v, kind in enumerate(types) if kind == 'A']
    right = [v for v, kind in enumerate(types) if kind == 'B']
    for u in left:
        for v in right:
            adjacency[u] |= 1 << v
            adjacency[v] |= 1 << u
    return adjacency, types


class StablePathProduct:
    def __init__(self, adjacency: list[int]):
        self.adjacency = adjacency
        self.order = len(adjacency)

    @lru_cache(maxsize=None)
    def component(self, mask: int, root: int) -> int:
        seen = 0
        frontier = 1 << root
        while frontier:
            seen |= frontier
            neighbors = 0
            work = frontier
            while work:
                bit = work & -work
                vertex = bit.bit_length() - 1
                work -= bit
                neighbors |= self.adjacency[vertex]
            frontier = neighbors & mask & ~seen
        return seen

    @lru_cache(maxsize=None)
    def independence_polynomial(self, mask: int) -> tuple[int, ...]:
        if not mask:
            return (1,)
        bit = mask & -mask
        vertex = bit.bit_length() - 1
        excluded = self.independence_polynomial(mask & ~bit)
        selected = self.independence_polynomial(
            mask & ~bit & ~self.adjacency[vertex]
        )
        return tuple(_polyadd(list(excluded), [0, *selected]))

    @lru_cache(maxsize=None)
    def rational_factors(
        self, mask: int, root: int
    ) -> tuple[tuple[tuple[tuple[int, ...], int], ...], int]:
        """Collect numerator/denominator graph-polynomial factors."""
        mask = self.component(mask, root)
        factors: Counter[tuple[int, ...]] = Counter()
        tree_order = 1
        removed = 1 << root
        for neighbor in range(self.order):
            if not (mask & (1 << neighbor) and self.adjacency[root] & (1 << neighbor)):
                continue
            child_data, child_order = self.rational_factors(
                mask & ~removed, neighbor
            )
            factors.update(dict(child_data))
            tree_order += child_order
            removed |= 1 << neighbor

        whole = self.independence_polynomial(mask)
        deletion = self.independence_polynomial(mask & ~(1 << root))
        if len(whole) > 1:
            factors[whole] += 1
        if len(deletion) > 1:
            factors[deletion] -= 1
        return tuple(sorted((poly, power) for poly, power in factors.items() if power)), tree_order

    def irreducible_factors(
        self, root: int
    ) -> tuple[Counter[tuple[int, ...]], int]:
        x = symbols('x')
        rational, tree_order = self.rational_factors((1 << self.order) - 1, root)
        factors: Counter[tuple[int, ...]] = Counter()
        for coefficients, outer_power in rational:
            expression = sum(value * x**index for index, value in enumerate(coefficients))
            content, rows = Poly(expression, x, domain='ZZ').factor_list()
            if content not in (1, -1):
                raise AssertionError(('nonunit content', content, coefficients))
            for factor, inner_power in rows:
                high_to_low = [int(value) for value in factor.all_coeffs()]
                low_to_high = tuple(reversed(high_to_low))
                if low_to_high[0] == -1:
                    low_to_high = tuple(-value for value in low_to_high)
                factors[low_to_high] += outer_power * inner_power
        factors += Counter()  # delete zero and negative entries only after audit
        negatives = {factor: power for factor, power in factors.items() if power < 0}
        if negatives:
            raise AssertionError(('uncancelled denominator', negatives))
        return factors, tree_order

    def materialize(self, root: int) -> list[list[int]]:
        """Materialize the ordered stable-path tree without DAG sharing."""
        tree: list[list[int]] = []

        def build(mask: int, vertex: int) -> int:
            mask = self.component(mask, vertex)
            here = len(tree)
            tree.append([])
            removed = 1 << vertex
            for neighbor in range(self.order):
                if not (
                    mask & (1 << neighbor)
                    and self.adjacency[vertex] & (1 << neighbor)
                ):
                    continue
                child = build(mask & ~removed, neighbor)
                tree[here].append(child)
                tree[child].append(here)
                removed |= 1 << neighbor
            return here

        assert build((1 << self.order) - 1, root) == 0
        return tree


def normalized_convolution(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    if min(len(left), len(right)) < 64:
        out = np.convolve(left, right)
    else:
        out = fftconvolve(left, right)
    # FFT noise in extreme tails is irrelevant to the filter and must not be
    # mistaken for a coefficient sign.
    out[np.abs(out) < 1e-15 * max(float(out.max()), 1e-300)] = 0.0
    out[out < 0] = 0.0
    total = float(out.sum())
    if not total:
        raise ArithmeticError('floating convolution underflow')
    return out / total


def power_distribution(
    coefficients: tuple[int, ...], exponent: int, activity: float
) -> np.ndarray:
    base = np.asarray(coefficients, dtype=np.float64)
    if np.any(base < 0):
        raise ValueError(('signed irreducible factor', coefficients))
    base *= np.power(activity, np.arange(len(base), dtype=np.float64))
    base /= base.sum()
    out = np.asarray([1.0])
    while exponent:
        if exponent & 1:
            out = normalized_convolution(out, base)
        exponent //= 2
        if exponent:
            base = normalized_convolution(base, base)
    return out


def product_distribution(
    factors: Counter[tuple[int, ...]], activity: float = 1.0
) -> np.ndarray:
    heap: list[tuple[int, int, np.ndarray]] = []
    serial = 0
    for coefficients, exponent in factors.items():
        if exponent <= 0:
            continue
        row = power_distribution(coefficients, exponent, activity)
        heapq.heappush(heap, (len(row), serial, row))
        serial += 1
    while len(heap) > 1:
        _, _, left = heapq.heappop(heap)
        _, _, right = heapq.heappop(heap)
        row = normalized_convolution(left, right)
        heapq.heappush(heap, (len(row), serial, row))
        serial += 1
    return heap[0][2] if heap else np.asarray([1.0])


def robust_valleys(
    distribution: np.ndarray, activity: float
) -> list[tuple[int, float, float, float]]:
    peak = float(distribution.max())
    rows = []
    for rank in range(1, len(distribution) - 1):
        left, middle, right = map(float, distribution[rank - 1 : rank + 2])
        if middle > peak * 1e-10:
            # If p_k is the distribution after tilting coefficients by
            # activity^k, undo the tilt in the two adjacent ratios.
            left_ratio = activity * left / middle
            right_ratio = right / (activity * middle)
            margin = min(left_ratio, right_ratio)
            if margin > 1.000001:
                rows.append((rank, left_ratio, 1.0, right_ratio))
    return rows


def source_polynomial(a: int, m: int, b: int, n: int) -> list[int]:
    return [
        (comb(m, rank) * a**rank if rank <= m else 0)
        + (comb(n, rank) * b**rank if rank <= n else 0)
        - (1 if rank == 0 else 0)
        for rank in range(max(m, n) + 1)
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-source-order', type=int, default=26)
    parser.add_argument('--max-tree-order', type=int, default=500_000)
    args = parser.parse_args()
    tested = 0
    for a in range(1, 11):
        for m in range(1, 11):
            for b in range(1, 11):
                for n in range(1, 11):
                    if (a, m) >= (b, n) or a * m + b * n > args.max_source_order:
                        continue
                    source = source_polynomial(a, m, b, n)
                    source_valleys = [
                        rank
                        for rank in range(1, len(source) - 1)
                        if source[rank - 1] > source[rank] < source[rank + 1]
                    ]
                    if not source_valleys:
                        continue
                    for second_first in (False, True):
                        adjacency, types = joined_cluster_graph(
                            a, m, b, n, second_first
                        )
                        stable = StablePathProduct(adjacency)
                        representatives = {
                            next(index for index, kind in enumerate(types) if kind == 'A'),
                            next(index for index, kind in enumerate(types) if kind == 'B'),
                        }
                        for root in representatives:
                            factors, tree_order = stable.irreducible_factors(root)
                            if tree_order > args.max_tree_order:
                                continue
                            valleys = []
                            tree_alpha = sum(
                                exponent * (len(factor) - 1)
                                for factor, exponent in factors.items()
                            )
                            # Tilting moves every coefficient window through
                            # the numerically stable bulk.  The corrected
                            # adjacent ratios are those of the original
                            # unweighted polynomial.
                            for activity in np.logspace(-6, 6, 73):
                                distribution = product_distribution(
                                    factors, float(activity)
                                )
                                rows = robust_valleys(distribution, float(activity))
                                if rows:
                                    valleys.extend(
                                        (float(activity), *row) for row in rows
                                    )
                                    break
                            tested += 1
                            print(
                                {
                                    'source': (a, m, b, n),
                                    'source_order': len(adjacency),
                                    'source_valleys': source_valleys,
                                    'second_first': second_first,
                                    'root_type': types[root],
                                    'tree_order': tree_order,
                                    'tree_alpha': tree_alpha,
                                    'robust_float_valleys': valleys,
                                },
                                flush=True,
                            )
                            if valleys:
                                raise SystemExit('candidate requires exact replay')
    print({'tested': tested, 'robust_float_valleys': 0}, flush=True)


if __name__ == '__main__':
    main()
