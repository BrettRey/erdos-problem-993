"""orbit.analyse_spec vs slib.analyse on expanded trees: exact interval ends must agree."""
from fractions import Fraction as F
from orbit import analyse_spec, expand, H, MSH
from slib import analyse
specs = [H(10, 2), H(66, 5), H(20, 3), MSH(8, 10, 2), MSH(12, 10, 2), MSH(3, 5, 1),
         [(3, [(4, [(2, [])]), (2, [])]), (5, [])], [(2, [(3, [(2, [(1, [])])])]), (4, [(3, [])]), (1, [(1, [])])],
         [(6, [(3, [(5, [(2, [])])])])]]
for sp in specs:
    a = analyse_spec(sp); b = analyse(expand(sp))
    ok = a['n'] == b['n'] and F(int(a['low'][0]), int(a['low'][1])) == F(int(b['low'][0]), int(b['low'][1])) \
        and F(int(a['up'][0]), int(a['up'][1])) == F(int(b['up'][0]), int(b['up'][1])) and (a['n_viol23'] > 0) == (b['n_viol23'] > 0) \
        and abs(a['maxL23'] - b['maxL23']) < 1e-12
    print(ok, a['n'], round(a['low_f'], 6), round(a['up_f'], 6), a['maxL23'], b['maxL23'])
    assert ok
