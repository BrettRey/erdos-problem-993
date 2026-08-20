#!/usr/bin/env python3
"""Exact checks for the erasure-shadow / forest-poset reduction."""

from __future__ import annotations

from math import comb


def projection_profile(code: tuple[int, ...], m: int) -> list[int]:
    """Return p_k = sum over k-coordinate projections of their sizes."""
    out = [0] * (m + 1)
    for mask in range(1 << m):
        out[mask.bit_count()] += len({word & mask for word in code})
    return out


def is_log_concave(a: list[int]) -> bool:
    return all(a[k] * a[k] >= a[k - 1] * a[k + 1]
               for k in range(1, len(a) - 1))


def is_binomially_normalized_lc(a: list[int]) -> bool:
    m = len(a) - 1
    return all(
        a[k] * a[k] * comb(m, k - 1) * comb(m, k + 1)
        >= a[k - 1] * a[k + 1] * comb(m, k) ** 2
        for k in range(1, m)
    )


def exhaustive_small_codes() -> None:
    for m in range(1, 5):
        words = tuple(range(1 << m))
        for code_mask in range(1, 1 << (1 << m)):
            code = tuple(w for w in words if (code_mask >> w) & 1)
            profile = projection_profile(code, m)
            assert is_log_concave(profile)
            assert is_binomially_normalized_lc(profile)


def transitive_closure(n: int, covers: tuple[tuple[int, int], ...]) -> list[list[bool]]:
    leq = [[False] * n for _ in range(n)]
    for i in range(n):
        leq[i][i] = True
    for i, j in covers:
        leq[i][j] = True
    for k in range(n):
        for i in range(n):
            if leq[i][k]:
                for j in range(n):
                    leq[i][j] = leq[i][j] or leq[k][j]
    return leq


def order_ideal_code(n: int, covers: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    leq = transitive_closure(n, covers)
    code = []
    for word in range(1 << n):
        # Values are isotone: i <= j forbids x_i=1, x_j=0.
        if all(not ((word >> i) & 1) or ((word >> j) & 1)
               for i in range(n) for j in range(n) if leq[i][j]):
            code.append(word)
    return tuple(code)


def poset_bipartite_independence_profile(
    n: int, covers: tuple[tuple[int, int], ...]
) -> list[int]:
    """Enumerate independent sets of B(P) by their size via ternary words."""
    leq = transitive_closure(n, covers)
    out = [0] * (n + 1)
    for ternary in range(3 ** n):
        q = ternary
        zeros: list[int] = []
        ones: list[int] = []
        size = 0
        for i in range(n):
            digit = q % 3
            q //= 3
            if digit:
                size += 1
                (zeros if digit == 1 else ones).append(i)
        if all(not leq[i][j] for i in ones for j in zeros):
            out[size] += 1
    return out


def antichain_profile(
    n: int, covers: tuple[tuple[int, int], ...]
) -> list[int]:
    leq = transitive_closure(n, covers)
    out = [0] * (n + 1)
    for mask in range(1 << n):
        chosen = [i for i in range(n) if (mask >> i) & 1]
        if all(not leq[i][j] and not leq[j][i]
               for at, i in enumerate(chosen) for j in chosen[at + 1:]):
            out[len(chosen)] += 1
    return out


def whisker_transform(a: list[int], n: int) -> list[int]:
    """Coefficients of (1+t)^n A(t/(1+t))."""
    out = [0] * (n + 1)
    for j, count in enumerate(a):
        for k in range(j, n + 1):
            out[k] += count * comb(n - j, k - j)
    return out


def check_poset_identity() -> None:
    examples = (
        (4, ((0, 1), (0, 2), (0, 3))),
        (5, ((0, 1), (2, 1), (2, 3), (4, 3))),
        (6, ((0, 1), (1, 2), (3, 2), (3, 4), (5, 4))),
    )
    for n, covers in examples:
        via_projection = projection_profile(order_ideal_code(n, covers), n)
        via_graph = poset_bipartite_independence_profile(n, covers)
        via_antichains = whisker_transform(antichain_profile(n, covers), n)
        assert via_projection == via_graph == via_antichains
        assert is_log_concave(via_graph)


def check_ulc_counterexample() -> None:
    code = (24, 47, 65, 67, 86, 93, 97, 99)
    profile = projection_profile(code, 7)
    assert profile == [1, 14, 80, 187, 220, 145, 52, 8]
    assert is_log_concave(profile)
    assert not is_binomially_normalized_lc(profile)
    left = profile[6] ** 2 * comb(7, 5) * comb(7, 7)
    right = profile[5] * profile[7] * comb(7, 6) ** 2
    assert (left, right) == (56_784, 56_840)


if __name__ == "__main__":
    exhaustive_small_codes()
    check_poset_identity()
    check_ulc_counterexample()
    print("PASS: erasure shadows, poset graph identity, and exact ULC counterexample")
