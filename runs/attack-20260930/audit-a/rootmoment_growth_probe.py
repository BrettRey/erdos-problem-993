"""DIAGNOSTIC (floats, heuristic search, not exhaustive): how large can the Prop. 4.1 quantity
q_v(1-q_v) delta_v^2 get on a rooted tree with n vertices at lambda = 12 (and 1/4)?
Bottom-up states (r, delta, n) with r = root occupation odds, delta = E(X|root in) - E(X|root out).
A new root with children multiset {s_i}: r = lambda / prod(1 + r_i), delta = 1 - sum q_i delta_i,
n = 1 + sum n_i (paper Eqs. 4.1-4.2).  Beam search: children drawn as d1 copies of state A plus
d2 copies of state B from the current beam; states bucketed by (n-scale, log r, sign delta) and
kept by |delta|.  Reports, per n-scale, the best q(1-q)delta^2 found and its log/log n slope.
A polynomial slope bounded away from 0 would be a lower bound on the true exponent a."""
import math, random, json, sys
random.seed(1)
def run(lam, rounds=14, beam=400, dmax=6, dlist=None):
    leaf = (lam, 1.0, 1)
    states = [leaf]
    best = {}
    for _ in range(rounds):
        new = []
        pool = states[:]
        for A in pool:
            for B in random.sample(pool, min(len(pool), 12)):
                for d1 in (dlist or range(1, dmax + 1)):
                    for d2 in range(0, 3):
                        n = 1 + d1 * A[2] + d2 * B[2]
                        if n > 10 ** 7: continue
                        qa = A[0] / (1 + A[0]); qb = B[0] / (1 + B[0])
                        lr = math.log(lam) - d1 * math.log1p(A[0]) - d2 * math.log1p(B[0])
                        if lr < -700: continue
                        r = math.exp(lr)
                        dl = 1 - d1 * qa * A[1] - d2 * qb * B[1]
                        new.append((r, dl, n))
        allst = states + new
        buckets = {}
        for s in allst:
            key = (int(math.log2(s[2])), round(math.log(s[0] + 1e-300), 0), s[1] > 0)
            if key not in buckets or abs(s[1]) > abs(buckets[key][1]): buckets[key] = s
        states = sorted(buckets.values(), key=lambda s: -abs(s[1]) * math.sqrt(s[0] / (1 + s[0]) ** 2))[:beam]
        for s in allst:
            q = s[0] / (1 + s[0]); val = q * (1 - q) * s[1] ** 2
            k = int(math.log2(s[2]))
            if k not in best or val > best[k][0]: best[k] = (val, s)
    return {k: dict(n=v[1][2], qqdelta2=v[0], slope=(math.log(v[0]) / math.log(v[1][2]) if v[1][2] > 1 and v[0] > 0 else None))
            for k, v in sorted(best.items())}
dl = [1, 2, 3, 4, 5, 6, 8, 11, 15, 20, 30]
out = {'12.0': run(12.0), '0.25': run(0.25), '0.25_wide_degrees': run(0.25, dlist=dl), '4.0_wide_degrees': run(4.0, dlist=dl), '12.0_wide_degrees': run(12.0, dlist=dl)}
for l, d in out.items():
    print('lambda', l, file=sys.stderr)
    for k, v in d.items(): print('  n~2^%d n=%d  max q(1-q)delta^2=%.4g  slope=%s' % (k, v['n'], v['qqdelta2'], None if v['slope'] is None else round(v['slope'], 3)), file=sys.stderr)
print(json.dumps(out, indent=1))
