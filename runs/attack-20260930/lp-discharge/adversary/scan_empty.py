"""Fine scan for the smallest single tree with an empty constant-sigma interval
(max lower end > min upper end, exact) among MSH(h;m,s), plus the lower-end trend."""
import sys, json
from fractions import Fraction as F
from orbit import analyse_spec, MSH
out = open(sys.argv[1], 'w')
for s in (1, 2, 3, 4):
    for m in range(6, 21):
        for h in range(8, 60):
            n = 1 + h * (1 + m * (s + 1))
            if n > 1700: break
            r = analyse_spec(MSH(h, m, s)); r.update(family='MSH', params=dict(h=h, m=m, s=s))
            r['empty'] = 'low' in r and 'up' in r and F(int(r['low'][0]), int(r['low'][1])) > F(int(r['up'][0]), int(r['up'][1]))
            out.write(json.dumps(r) + '\n'); out.flush()
