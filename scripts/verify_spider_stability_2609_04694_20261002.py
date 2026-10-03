"""Certified check of Zhang & Tu, arXiv:2609.04694v1, Theorem 4, on small spiders.

Claim: every independence root of every spider S(l_1, ..., l_d) has strictly
negative real part.

For every spider on n <= N vertices (multisets of leg lengths summing to
n - 1), build I(S, z) = prod F_{l_j} + z prod F_{l_j - 1} exactly, where F_m
is the path polynomial (F_0 = 1, F_1 = 1 + z, F_m = F_{m-1} + z F_{m-2}),
isolate its roots with Arb (python-flint), and require each root's real-part
ball to lie strictly left of 0. A ball touching 0 counts as undecided, a ball
strictly right of 0 as a counterexample. Never float64 (project rule).

Run: venv/bin/python scripts/verify_spider_stability_2609_04694_20261002.py [N]
"""

from __future__ import annotations

import sys

import flint


def path_polys(m_max: int) -> list[list[int]]:
    F = [[1], [1, 1]]
    for m in range(2, m_max + 1):
        a, b = F[m - 1], [0] + F[m - 2]
        n = max(len(a), len(b))
        F.append([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
    return F


def mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def partitions(total: int, max_part: int | None = None):
    if max_part is None:
        max_part = total
    if total == 0:
        yield ()
        return
    for p in range(min(total, max_part), 0, -1):
        for rest in partitions(total - p, p):
            yield (p,) + rest


def spider_poly(legs: tuple[int, ...], F: list[list[int]]) -> list[int]:
    excl, incl = [1], [1]
    for l in legs:
        excl = mul(excl, F[l])
        incl = mul(incl, F[l - 1])
    incl = [0] + incl
    n = max(len(excl), len(incl))
    return [(excl[i] if i < len(excl) else 0) + (incl[i] if i < len(incl) else 0) for i in range(n)]


def main() -> None:
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 35
    F = path_polys(N)
    checked = undecided = bad = 0
    worst = None  # largest certified upper bound on Re(root)
    for n in range(2, N + 1):
        for legs in partitions(n - 1):
            p = spider_poly(legs, F)
            roots = flint.fmpz_poly(p).complex_roots()
            checked += 1
            for r, _m in roots:
                re = r.real
                upper = re.mid() + re.rad()
                if worst is None or upper > worst[0]:
                    worst = (upper, legs)
                if re < 0:
                    continue
                if re > 0:
                    bad += 1
                    print(f"COUNTEREXAMPLE: legs={legs} root={r}")
                else:
                    undecided += 1
                    print(f"undecided: legs={legs} root={r}")
    print(
        f"spiders on 2..{N} vertices: {checked}; counterexamples {bad}; undecided {undecided}; "
        f"max certified upper bound on Re(root) = {float(worst[0]):.6f} at legs={worst[1]}"
    )


if __name__ == "__main__":
    main()
