#!/usr/bin/env python3
"""Check the 32 rational inequalities of Section 4.1 (range 3 < H <= 16).

For fixed K, A_K(delta) = (1/2) sum_{i,j=-K}^{K} w_i w_j (i-j)^2 with w_0 = 1,
w_r = R_r = a^r prod_{j=1}^{r-1} (1 - j delta), w_{-r} = L_r =
(1-delta)^{-r} prod_{j=2}^{r+1} (1 - j delta), a = (1 - 2 delta)/(1 - delta),
is nonincreasing in delta on 0 < delta <= 1/(K+1), and so is the target
(3 + delta)/(4 delta^2). So on [lo, hi] with (K+1) hi <= 1 it suffices that
(3 + lo)/(4 lo^2) <= A_K(hi).

Cells: K = m on [1/(m+2), 1/(m+1)] for m = 4..15; K = 3 on 20 equal pieces
of [1/5, 1/4]. Also checks the remark that the symmetrized bound fails at
H = 7/2 (S T < Q(H)), and that the ratios printed in the paper are the
rounded-down values. Standard library only; exit status 0 means every check
passed.
"""
from __future__ import annotations

from fractions import Fraction as F
import sys

PRINTED_M = ["1.072", "1.439", "1.763", "2.050", "2.304", "2.531",
             "2.735", "2.920", "3.089", "3.243", "3.385", "3.516"]
PRINTED_K3_MIN = "1.0078"


def weights(d: F, K: int) -> dict[int, F]:
    a = (1 - 2 * d) / (1 - d)
    w = {0: F(1)}
    for r in range(1, K + 1):
        R = a ** r
        for j in range(1, r):
            R *= 1 - j * d
        L = 1 / (1 - d) ** r
        for j in range(2, r + 2):
            L *= 1 - j * d
        w[r], w[-r] = R, L
    return w


def A(d: F, K: int) -> F:
    w = weights(d, K)
    return F(1, 2) * sum(w[i] * w[j] * (i - j) ** 2 for i in w for j in w)


def target(d: F) -> F:
    return (3 + d) / (4 * d * d)


def floor_str(x: F, places: int) -> str:
    q = (x.numerator * 10 ** places) // x.denominator
    return f"{q // 10 ** places}.{q % 10 ** places:0{places}d}"


def main() -> int:
    ok = True
    ratios = []
    for m in range(4, 16):
        lo, hi = F(1, m + 2), F(1, m + 1)
        assert (m + 1) * hi <= 1
        r = A(hi, m) / target(lo)
        ok &= r >= 1
        ratios.append(floor_str(r, 3))
    print("m = 4..15 ratios (rounded down):", ", ".join(ratios))
    ok &= ratios == PRINTED_M
    lo0, hi0 = F(1, 5), F(1, 4)
    k3 = []
    for i in range(20):
        lo = lo0 + (hi0 - lo0) * i / 20
        hi = lo0 + (hi0 - lo0) * (i + 1) / 20
        assert 4 * hi <= 1
        k3.append(A(hi, 3) / target(lo))
    ok &= all(r >= 1 for r in k3)
    kmin = min(k3)
    print(f"K = 3, 20 pieces: min ratio {floor_str(kmin, 4)} on piece {k3.index(kmin)}")
    ok &= floor_str(kmin, 4) == PRINTED_K3_MIN and k3.index(kmin) == 19
    # remark: symmetrized bound fails at H = 7/2 (K = 3)
    H = F(7, 2)
    b = [F(1)]
    for s in range(1, 4):
        b.append(b[-1] * (1 - F(s) / H))
    S = 1 + 2 * sum(b[1:])
    T = 2 * sum(r * r * b[r] for r in range(1, 4))
    Q = (3 * H + 4) * (H + 1) / 4
    print("H = 7/2: S T - Q =", S * T - Q)
    ok &= S * T < Q
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
