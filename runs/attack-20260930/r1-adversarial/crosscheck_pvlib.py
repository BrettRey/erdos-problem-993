"""Cross-check r1lib (rerooting, flint) against route-freecount/pv_lib (fresh DP per forest)."""
import sys, random, time
sys.path.insert(0, '../route-freecount')
import pv_lib
from r1lib import all_seqs, adj_from_edges, evaluate, B
rng = random.Random(1)
ok = 0
for trial in range(300):
    n = rng.randrange(2, 40)
    edges = [(i, rng.randrange(i)) for i in range(1, n)]
    adj = adj_from_edges(n, edges)
    I, J = all_seqs(adj)
    I2, J2, _, _ = pv_lib.tree_data(adj)
    assert I == I2, (I, I2)
    for v in range(n):
        a = J[v][:]; b = J2[v][:]
        while len(a) > 1 and a[-1] == 0: a.pop()
        while len(b) > 1 and b[-1] == 0: b.pop()
        assert a == b, (v, a, b)
    ok += 1
print("random trees agree:", ok)
def hub(ts):
    b = B()
    for t in ts:
        c = b.add(0)
        for _ in range(t): b.add(c)
    return b.adj()
r = evaluate(hub([5]*18)); print("hubstars_18_5 normB", r['B'], "(wave-1: -0.09919366569759212)", r['n'])
t=time.time(); r = evaluate(hub([2]*49)); print("n", r['n'], "time", time.time()-t, r['A'], r['B'], r['C'])
t=time.time(); rr = [(i, rng.randrange(i)) for i in range(1,150)]; r=evaluate(adj_from_edges(150, rr)); print("random n150 time", time.time()-t, r['A'], r['B'], r['C'])
