#!/usr/bin/env python3
"""Exact audit for notes/law_v_one_hub_proof_2026-08-20.md."""

from itertools import combinations_with_replacement
from math import comb


def add(a: list[int], b: list[int]) -> list[int]:
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0)
            + (b[i] if i < len(b) else 0) for i in range(n)]


def mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def sub(a: list[int], b: list[int]) -> list[int]:
    return add(a, [-x for x in b])


def cf(a: list[int], k: int) -> int:
    return a[k] if 0 <= k < len(a) else 0


def path_poly(n: int) -> list[int]:
    if n <= 0:
        return [1]
    return [comb(n + 1 - k, k) for k in range((n + 1) // 2 + 1)]


def prod_path(arms: tuple[int, ...], shift: int = 0) -> list[int]:
    out = [1]
    for length in arms:
        out = mul(out, path_poly(length - shift))
    return out


def lr_leq(a: list[int], b: list[int]) -> bool:
    """a <=_lr b, including support-boundary minors."""
    n = max(len(a), len(b))
    return all(cf(b, k) * cf(a, k - 1) >= cf(b, k - 1) * cf(a, k)
               for k in range(n + 1))


def one_hub(arms: tuple[int, ...]) -> tuple[list[int], list[int], list[int], list[int]]:
    q = prod_path(arms)
    r = prod_path(arms, 1)
    c = add(q, [0] + r)
    b = add(c, [0] + q)
    return q, r, c, b


def general_pair(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[list[int], list[int], list[int]]:
    qa, ra = prod_path(a), prod_path(a, 1)
    qb, rb = prod_path(b), prod_path(b, 1)
    qq = mul(qa, qb)
    c = add(qq, [0] + mul(ra, rb))
    h = sub(add(mul(ra, qb), mul(qa, rb)), mul(ra, rb))
    full_b = add(c, [0] + h)
    return c, h, full_b


def main() -> None:
    # The binomial formula agrees with the path recurrence.
    prev2, prev1 = [1], [1, 1]
    assert path_poly(0) == prev2 and path_poly(1) == prev1
    for n in range(2, 41):
        cur = add(prev1, [0] + prev2)
        assert cur == path_poly(n)
        prev2, prev1 = prev1, cur

    # Lemma 1 for individual paths.
    for n in range(1, 41):
        assert lr_leq(path_poly(n - 1), path_poly(n))

    arms_pool: list[tuple[int, ...]] = []
    for count in range(1, 5):
        for arms in combinations_with_replacement(range(1, 11), count):
            if sum(arms) <= 24:
                arms_pool.append(arms)

    checked_entries = 0
    for arms in arms_pool:
        q, r, c, b = one_hub(arms)
        assert lr_leq(r, q)
        top = max(len(q), len(r), len(c), len(b)) + 2
        for k in range(top):
            # Q log-concavity and the two summands in equation (5).
            assert cf(q, k) * cf(q, k - 1) >= cf(q, k + 1) * cf(q, k - 2)
            assert cf(r, k - 1) * cf(q, k - 1) >= cf(r, k) * cf(q, k - 2)
            cross = cf(c, k) * cf(q, k - 1) - cf(c, k + 1) * cf(q, k - 2)
            assert cross >= 0
            v = cf(c, k) * cf(b, k) - cf(c, k + 1) * cf(b, k - 1)
            assert v >= 0
            checked_entries += 1

    # Regression: the direct two-hub lift of (5) is false, although LAW-V holds.
    a, b = (1, 1), (1, 1, 5)
    c, h, full_b = general_pair(a, b)
    k = 6
    naive = cf(c, k) * cf(h, k - 1) - cf(c, k + 1) * cf(h, k - 2)
    v = cf(c, k) * cf(full_b, k) - cf(c, k + 1) * cf(full_b, k - 1)
    assert naive == -1
    assert v >= 0

    print(f"PASS: {len(arms_pool)} one-hub arm multisets, "
          f"{checked_entries} exact index checks")
    print("PASS: naive two-hub lift regression has cross margin -1 while V >= 0")


if __name__ == "__main__":
    main()
