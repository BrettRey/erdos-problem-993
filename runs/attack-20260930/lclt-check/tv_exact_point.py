#!/usr/bin/env python3
"""EXACT: TV(h_B)^2 v^3 at the float-scan maximiser (B=22, p=12/13) and neighbours,
compared with rational bounds 1.579^2 and 1.581^2 (brackets TV v^{3/2})."""
import json
from fractions import Fraction as Fr
from kernel_check_exact import tv_exact
res = []
for B in (20, 21, 22, 23, 24):
    TV, v, ch = tv_exact(B, Fr(12, 13)); val = TV * TV * v ** 3
    res.append(dict(B=B, gt_1579=bool(val > Fr(1579, 1000) ** 2), lt_1581=bool(val < Fr(1581, 1000) ** 2), float_diag=float(val) ** 0.5))
print(json.dumps(res)); json.dump(res, open("rerun/tv_exact_point.json", "w"), indent=1)
