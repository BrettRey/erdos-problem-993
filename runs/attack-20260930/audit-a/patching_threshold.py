"""Exact integer check of the Section 8 patching step of arXiv:2609.20961 (p. 16) and of the
Lean `index_interleave` (Main.lean:151).  For each n and every alpha with n/2 <= alpha <= n
(a forest is bipartite, so alpha >= n/2), test the Lean criterion L <= P < m <= U with
L = ceil(n/5), P = ceil(n/4), m = ceil((2 alpha - 1)/3), U = floor(c_hi alpha),
c_hi = 17/25 (paper) or 64/95 (Lean).  Report the least n0 such that every n >= n0 passes
(checked up to 5000; beyond that the slack c_hi - 2/3 >= 2/285 grows linearly)."""
from fractions import Fraction as Fr
import math, json
def ceil_div(a, b): return -((-a) // b)
def ok(n, chi):
    L, P = ceil_div(n, 5), ceil_div(n, 4)
    for al in range(ceil_div(n, 2), n + 1):
        m = ceil_div(2 * al - 1, 3); U = math.floor(chi * al)
        if not (L <= P < m <= U): return False
    return True
out = {}
for name, chi in [('paper_17/25', Fr(17, 25)), ('lean_64/95', Fr(64, 95))]:
    bad = [n for n in range(1, 5001) if not ok(n, chi)]
    out[name] = dict(last_failing_n=max(bad) if bad else None, n0=(max(bad) + 1) if bad else 1,
                     failing_count=len(bad))
print(json.dumps(out, indent=1))
