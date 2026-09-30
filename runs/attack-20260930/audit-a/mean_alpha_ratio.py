"""Diagnostic for improvement lever 3 (lowering lambda_max): over ALL trees with n <= NMAX
(nauty gentreeg; forests reduce to trees because E_lambda X and alpha are additive over
components), compute min E_lambda X / alpha exactly (fractions) for several lambda and test
E_lambda X > 2 alpha / 3.  Stars K_{1,m} give E/alpha -> lambda/(1+lambda), so lambda > 2 is
necessary asymptotically.  This is finite evidence, not a proof."""
import subprocess, sys, json
from fractions import Fraction as Fr
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 18
LAMS = [Fr(2), Fr(5, 2), Fr(3), Fr(4)]
def parse_g6(line):
    # gentreeg -p parent array: entry i (1-indexed vertex i) is its parent, 0 for the root
    par = [int(x) for x in line.split()]; n = len(par)
    adj = [[] for _ in range(n)]
    for i, pp in enumerate(par):
        if pp: adj[i].append(pp - 1); adj[pp - 1].append(i)
    return n, adj
def indep_poly(n, adj):
    order, parent, seen = [], [-1] * n, [False] * n
    st = [0]; seen[0] = True
    while st:
        v = st.pop(); order.append(v)
        for w in adj[v]:
            if not seen[w]: seen[w] = True; parent[w] = v; st.append(w)
    A = [[1] for _ in range(n)]; B = [[0, 1] for _ in range(n)]   # excluded / included
    def mul(p, q):
        r = [0] * (len(p) + len(q) - 1)
        for i, a in enumerate(p):
            if a:
                for j, b in enumerate(q): r[i + j] += a * b
        return r
    for v in reversed(order):
        for w in adj[v]:
            if w != parent[v]:
                s = [x + y for x, y in zip(A[w] + [0] * (len(B[w]) - len(A[w])), B[w])]
                A[v] = mul(A[v], s); B[v] = mul(B[v], A[w])
    Z = [x + y for x, y in zip(A[0] + [0] * (len(B[0]) - len(A[0])), B[0])]
    while Z and Z[-1] == 0: Z.pop()
    return Z
out = {}
for n in range(2, NMAX + 1):
    res = {str(l): [None, None] for l in LAMS}
    p = subprocess.run(['/opt/homebrew/bin/gentreeg', '-q', '-p', str(n)], capture_output=True, text=True)
    for line in p.stdout.splitlines():
        nn, adj = parse_g6(line); Z = indep_poly(nn, adj); al = len(Z) - 1
        for l in LAMS:
            num = sum(k * c * l ** k for k, c in enumerate(Z)); den = sum(c * l ** k for k, c in enumerate(Z))
            ratio = num / den / al
            cur = res[str(l)]
            if cur[0] is None or ratio < cur[0]: res[str(l)] = [ratio, line.strip()]
    out[n] = {l: dict(min_ratio=float(v[0]), exceeds_2_3=bool(v[0] > Fr(2, 3)), argmin_parent_array=v[1]) for l, v in res.items()}
    print(n, {l: round(v['min_ratio'], 5) for l, v in out[n].items()}, file=sys.stderr)
print(json.dumps(out, indent=1))
