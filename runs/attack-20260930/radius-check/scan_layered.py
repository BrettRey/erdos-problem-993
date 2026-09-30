"""For layered trees of depth 2..6: at each window level k, the sign of the
defect D at every depth, the widest run of consecutive positive depths, and
R = max over positive depths of the distance to the nearest negative depth
(a radius-r certificate, even tree-specific, needs R <= 2r). Exact integers."""
import itertools, json, time
from math import ceil
from layered import layered, coeffs
co = lambda p, k: p[k] if 0 <= k < len(p) else 0
def analyse(b):
    I, J, N, deg = layered(b); I = coeffs(I); Jc = [coeffs(x) for x in J]
    n = sum(N); a = len(I) - 1; q = ceil((2 * a - 1) / 3); lo = ceil(n / 4)
    worstR, worstW, pats = 0, 0, set()
    for k in range(lo, min(q, a - 1) + 1):
        D = [k * co(I, k - 1) * co(Jc[i], k) - (k + 1) * co(I, k) * co(Jc[i], k - 1) for i in range(len(b) + 1)]
        assert sum(N[i] * D[i] for i in range(len(D))) == k * (k + 1) * (co(I, k - 1) * co(I, k + 1) - co(I, k) ** 2)
        sg = ''.join('+' if d > 0 else '-' if d < 0 else '0' for d in D); pats.add(sg)
        neg = [i for i, d in enumerate(D) if d < 0]; pos = [i for i, d in enumerate(D) if d > 0]
        if pos:
            R = max(min(abs(i - j) for j in neg) for i in pos) if neg else 99
            w = max(len(r) for r in ''.join('+' if d > 0 else '.' for d in D).split('.'))
            worstR, worstW = max(worstR, R), max(worstW, w)
    return n, worstR, worstW, sorted(pats)
t0 = time.time(); res = {}
grids = {
    2: [(m, s) for m in list(range(2, 41)) + list(range(45, 201, 5)) for s in range(1, 13)],
    3: [(h, m, s) for h in list(range(2, 21)) + [24, 28, 32, 38, 44, 50, 60] for m in range(1, 21) for s in range(1, 5)],
    4: [(c, h, m, s) for c in range(2, 11) for h in range(1, 13) for m in range(1, 13) for s in range(1, 4)],
    5: [(a, c, h, m, s) for a in range(2, 7) for c in range(1, 7) for h in range(1, 7) for m in range(1, 9) for s in range(1, 4)],
    6: [(z, a, c, h, m, s) for z in range(2, 5) for a in range(1, 5) for c in range(1, 5) for h in range(1, 5) for m in range(1, 7) for s in range(1, 3)],
}
for depth, grid in grids.items():
    best = []; counts = {}
    for b in grid:
        n, R, W, pats = analyse(list(b))
        counts[R] = counts.get(R, 0) + 1
        best.append((R, W, n, b, pats))
    best.sort(key=lambda x: (-x[0], -x[1], x[2]))
    res[depth] = dict(trees=len(grid), R_histogram=counts, top=[(R, W, n, list(b), pats[:6]) for R, W, n, b, pats in best[:5]])
    print(f'depth {depth}: {len(grid)} trees; histogram of R (0 = no positive defect in window): {dict(sorted(counts.items()))}; worst: {[(R, W, n, b) for R, W, n, b, _ in best[:3]]} ({time.time()-t0:.0f}s)', flush=True)
json.dump(res, open('scan_layered.json', 'w'), indent=1, default=str)
