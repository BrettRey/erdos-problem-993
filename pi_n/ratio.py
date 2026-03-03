from __future__ import annotations

import math

from .types import Ratio


def canon_ratio(r: Ratio) -> Ratio:
    g = math.gcd(abs(r.num), r.den)
    return Ratio(num=r.num // g, den=r.den // g)


def cmp_ratio(r1: Ratio, r2: Ratio) -> int:
    lhs = r1.num * r2.den
    rhs = r2.num * r1.den
    if lhs < rhs:
        return -1
    if lhs > rhs:
        return 1
    return 0
