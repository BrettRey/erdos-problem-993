"""Check the symrows closed forms against the generic forest DP on small cases."""
import sys
sys.path.insert(0, '.')
from symrows import H, MSH, coeffs
from forestdp import indpoly


def build_H(m, s):
    adj = [[] for _ in range(1 + m + m * s)]
    c = 1
    rep = {'hub': 0}
    for _ in range(m):
        sp = c; c += 1
        adj[0].append(sp); adj[sp].append(0); rep.setdefault('sup', sp)
        for _ in range(s):
            adj[sp].append(c); adj[c].append(sp); rep.setdefault('leaf', c); c += 1
    return adj, rep


def build_MSH(h, m, s):
    adj = [[] for _ in range(1 + h * (1 + m + m * s))]
    c = 1
    rep = {'cen': 0}
    for _ in range(h):
        hb = c; c += 1
        adj[0].append(hb); adj[hb].append(0); rep.setdefault('hub', hb)
        for _ in range(m):
            sp = c; c += 1
            adj[hb].append(sp); adj[sp].append(hb); rep.setdefault('sup', sp)
            for _ in range(s):
                adj[sp].append(c); adj[c].append(sp); rep.setdefault('leaf', c); c += 1
    return adj, rep


ok = True
for fam, args in [('H', (4, 3)), ('H', (7, 2)), ('MSH', (3, 4, 2)), ('MSH', (2, 5, 3)), ('MSH', (4, 3, 1))]:
    I, J, n = (H if fam == 'H' else MSH)(*args)
    A, rep = (build_H if fam == 'H' else build_MSH)(*args)
    this = coeffs(I) == indpoly(A, set(range(n)))
    for t, (Jp, d, nb) in J.items():
        v = rep[t]
        this &= coeffs(Jp) == indpoly(A, set(range(n)) - {v} - set(A[v]))
        this &= len(A[v]) == d
        this &= sorted((len(A[w]) for w in A[v])) == sorted(dd for (_, cnt, dd) in nb for _ in range(cnt))
    ok &= this
    print(fam, args, 'match' if this else 'MISMATCH')
print('ALL OK' if ok else 'FAIL')
