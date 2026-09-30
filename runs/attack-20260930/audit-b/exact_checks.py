#!/usr/bin/env python3
"""Exact-integer side checks for the audit (Python ints / Fractions only).

1. A000055 (unlabelled trees) via Otter's formula from A000081 (rooted trees),
   to ground the statement about how far exhaustive enumeration can reach.
2. Patching (Section 8 / proof of Thm 1.1): the smallest alpha from which every
   possible valley index j (ceil(n/4) < j < ceil((2 alpha - 1)/3)) lies inside the
   log-concave window ceil(n/5) <= j <= floor(c alpha), for c = 17/25 (paper) and
   c = 64/95 (Lean), given alpha >= n/2.
3. Section 7 margin: log(13)/Q*(12) vs 17/25 is not rational; recorded as float diagnostic
   in explore_lammax.py.  Here only the Lean rational chain 64/95 = (64/25)/(19/5) is checked.
"""
from fractions import Fraction as Fr
import json


def rooted_trees(N):
    # A000081: a(n+1) = (1/n) sum_{k=1}^{n} (sum_{d|k} d a(d)) a(n-k+1)
    a = [0] * (N + 1)
    a[1] = 1
    s = [0] * (N + 1)  # s[k] = sum_{d|k} d a(d)
    for n in range(1, N):
        s[n] = sum(d * a[d] for d in range(1, n + 1) if n % d == 0)
        tot = sum(s[k] * a[n - k + 1] for k in range(1, n + 1))
        assert tot % n == 0
        a[n + 1] = tot // n
    return a


def free_trees(N):
    a = rooted_trees(N)
    t = [0] * (N + 1)
    for n in range(1, N + 1):
        # Otter: t(n) = a(n) - (1/2) [ sum_{i+j=n} a(i)a(j) - (a(n/2) if n even) ]
        conv = sum(a[i] * a[n - i] for i in range(1, n))
        if n % 2 == 0:
            conv -= a[n // 2]
        assert conv % 2 == 0
        t[n] = a[n] - conv // 2
    return t


def ceil_div(x, y):
    return -((-x) // y)


def patch_ok(n, alpha, c):
    lo_valley = ceil_div(n, 4) + 1
    hi_valley = ceil_div(2 * alpha - 1, 3) - 1
    if lo_valley > hi_valley:
        return True
    lc_lo = ceil_div(n, 5)
    lc_hi = (c * alpha).__floor__()
    return lc_lo <= lo_valley and hi_valley <= lc_hi


if __name__ == "__main__":
    t = free_trees(60)
    assert t[1:11] == [1, 1, 1, 2, 3, 6, 11, 23, 47, 106]
    assert t[15] == 7741 and t[20] == 823065 and t[25] == 104636890 and t[26] == 279793450
    out = {"A000055": {n: t[n] for n in range(1, 61)}}
    for name, c in [("17/25", Fr(17, 25)), ("64/95", Fr(64, 95))]:
        bad = [(n, al) for n in range(1, 400) for al in range((n + 1) // 2, n + 1) if not patch_ok(n, al, c)]
        out["patch_failures_" + name] = bad[:20]
        out["patch_max_failing_n_" + name] = max((b[0] for b in bad), default=None)
    json.dump(out, open("data/exact_checks.json", "w"), indent=1)
    for n in (32, 33, 36, 40, 45, 50, 60):
        print(n, t[n], f"{t[n]:.3e}")
    print("patch failures 17/25:", out["patch_failures_17/25"], "max n", out["patch_max_failing_n_17/25"])
    print("patch failures 64/95:", out["patch_failures_64/95"], "max n", out["patch_max_failing_n_64/95"])
