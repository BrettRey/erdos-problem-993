#!/usr/bin/env python3
"""Exact checks for the universal Pascal-smoothing shadow lemma."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb


def is_downset(family_mask: int, n: int) -> bool:
    if not (family_mask & 1):
        return False
    for face in range(1 << n):
        if not ((family_mask >> face) & 1):
            continue
        bit = face
        while bit:
            atom = bit & -bit
            if not ((family_mask >> (face ^ atom)) & 1):
                return False
            bit ^= atom
    return True


def face_counts(family_mask: int, n: int) -> list[int]:
    out = [0] * (n + 1)
    for face in range(1 << n):
        if (family_mask >> face) & 1:
            out[face.bit_count()] += 1
    return out


def shadow_profile(a: list[int], ambient: int) -> list[int]:
    return [
        sum(count * comb(ambient - j, d) for j, count in enumerate(a)
            if ambient - j >= d)
        for d in range(ambient + 1)
    ]


def pascal_profile(a: list[int], ambient: int) -> list[Fraction]:
    b = [Fraction(0)] * (ambient + 1)
    for j, count in enumerate(a):
        b[j] = Fraction(count, comb(ambient, j))
    return [
        sum(Fraction(comb(m, j)) * b[j] for j in range(m + 1))
        for m in range(ambient + 1)
    ]


def check_family(a: list[int], ambient: int) -> None:
    b = [Fraction(0)] * (ambient + 1)
    for j, count in enumerate(a):
        b[j] = Fraction(count, comb(ambient, j))
    assert all(b[j] >= b[j + 1] for j in range(ambient))

    e = shadow_profile(a, ambient)
    s = pascal_profile(a, ambient)
    for d in range(ambient + 1):
        assert Fraction(e[d], comb(ambient, d)) == s[ambient - d]

    for m in range(1, ambient):
        assert 9 * s[m] * s[m] >= 8 * s[m - 1] * s[m + 1]

    for d in range(1, ambient):
        m = ambient - d
        if d * m <= 8 * (ambient + 1):
            assert e[d] * e[d] >= e[d - 1] * e[d + 1]


def exhaustive_complexes_through_four_vertices() -> int:
    checked = 0
    for n in range(1, 5):
        for family_mask in range(1 << (1 << n)):
            if not is_downset(family_mask, n):
                continue
            a = face_counts(family_mask, n)
            # Also check up to three constant/ghost coordinates.
            for extra in range(4):
                check_family(a, n + extra)
            checked += 1
    return checked


def exhaustive_decreasing_sequences() -> int:
    checked = 0
    # The Pascal lemma only uses monotonicity, so also check a finite grid
    # independently of the simplicial-complex enumeration.
    for length in range(3, 8):
        for seq in product(range(5), repeat=length):
            if any(seq[j] < seq[j + 1] for j in range(length - 1)):
                continue
            s = [
                sum(comb(m, j) * seq[j] for j in range(m + 1))
                for m in range(length)
            ]
            for m in range(1, length - 1):
                assert 9 * s[m] * s[m] >= 8 * s[m - 1] * s[m + 1]
            checked += 1
    return checked


def check_depth_three_constant() -> None:
    for ambient in range(4, 200):
        lhs = Fraction(8, 9) * Fraction(
            4 * (ambient - 2), 3 * (ambient - 3)
        )
        assert lhs == Fraction(32 * (ambient - 2), 27 * (ambient - 3))
        assert lhs > 1


if __name__ == "__main__":
    complexes = exhaustive_complexes_through_four_vertices()
    sequences = exhaustive_decreasing_sequences()
    check_depth_three_constant()
    print(
        "PASS: Pascal smoothing; "
        f"{complexes} complexes and {sequences} decreasing sequences"
    )
