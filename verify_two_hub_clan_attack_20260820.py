#!/usr/bin/env python3
"""Exact arithmetic checks for the two-hub clan cancellation attack."""

from math import comb


def single_surplus(p: int, q: int) -> int:
    assert 0 <= q <= 2 and p >= q + 2
    if q == 0:
        return comb(p, 1) - comb(p, 0)
    if q == 1:
        return comb(p - 1, 1) - comb(p - 1, 0)
    return comb(p, 2) - comb(p, 1)


def double_surplus(p: int, q: int) -> int:
    assert q >= 3 and p >= q + 2
    r = p + q - 2
    return comb(r, q) - comb(r, q - 1)


def transform_unit_arms(alpha: tuple[int, int, tuple[int, ...], tuple[int, ...]],
                         sides: str) -> tuple[int, int, tuple[int, ...], tuple[int, ...]]:
    """Length-one special case: zero named hubs and clone first unit arm."""
    u, v, ua, va = alpha
    ua, va = list(ua), list(va)
    if "u" in sides:
        u = 0
        first = ua.index(1)
        ua[first] = 2
    if "v" in sides:
        v = 0
        first = va.index(1)
        va[first] = 2
    return u, v, tuple(ua), tuple(va)


def main() -> None:
    for q in range(3):
        for p in range(q + 2, 101):
            assert single_surplus(p, q) > 0
    for q in range(3, 101):
        for p in range(q + 2, 151):
            got = double_surplus(p, q)
            assert got > 0
            r = p + q - 2
            assert comb(r, q) * q == comb(r, q - 1) * (p - 1)

    # Exact collision from the note: five v-leaves, three u-leaves.
    double_source = (1, 1, (1, 1, 1), (1, 1, 1, 1, 1))
    single_source = (0, 1, (2, 1, 1), (1, 1, 1, 1, 1))
    image_double = transform_unit_arms(double_source, "uv")
    image_single = transform_unit_arms(single_source, "v")
    assert image_double == image_single
    assert sum(double_source[:2]) + sum(double_source[2]) + sum(double_source[3]) == 10
    assert sum(single_source[:2]) + sum(single_source[2]) + sum(single_source[3]) == 10
    assert sum(image_double[:2]) + sum(image_double[2]) + sum(image_double[3]) == 10

    print("PASS: all local binomial surpluses positive for audited ranges")
    print("PASS: exact single/double map collision reproduced at total weight 10")


if __name__ == "__main__":
    main()
