"""Independent recheck: closed-form class values vs generic tree DP (pv_lib) on explicit trees."""
import sys
sys.path.insert(0, "../route-freecount")
from pv_lib import tree_data, co as pco, window
import fam

def generic_D(adj, k):
    I, J, _E, alpha = tree_data(adj)
    return I, [k * pco(I, k - 1) * pco(Jv, k) - (k + 1) * pco(I, k) * pco(Jv, k - 1) for Jv in J]

def cmp(F, adj, rep):
    """rep: class name -> a representative vertex index"""
    I, J, _E, alpha = tree_data(adj)
    assert [co for co in I] == [fam.co(F.I, i) for i in range(F.I.degree() + 1)], F.label
    assert alpha == F.alpha
    lo, q = window(len(adj), alpha)
    for k in range(1, alpha + 1):
        D, L, *_ = F.level(k)
        Dg = [k * pco(I, k - 1) * pco(Jv, k) - (k + 1) * pco(I, k) * pco(Jv, k - 1) for Jv in J]
        for nm, v in rep.items():
            assert Dg[v] == D[nm], (F.label, k, nm)
            from fractions import Fraction
            Lg = Fraction(Dg[v], len(adj[v]) + 1) + sum(Fraction(Dg[w], len(adj[w]) + 1) for w in adj[v])
            assert Lg == L[nm], (F.label, k, nm, "L")
    return True

for m, s in [(3, 1), (4, 2), (9, 2), (5, 3), (4, 4)]:
    F = fam.hubstar(m, s); adj = fam.build_hubstar(m, s)
    print(F.label, cmp(F, adj, {"hub": 0, "mid": 1, "leaf": 2}))
for m, s in [(3, 2), (5, 1), (4, 3)]:
    F = fam.twohub(m, s); adj = fam.build_twohub(m, s)
    print(F.label, cmp(F, adj, {"hub": 0, "mid": 2, "leaf": 3}))
for h, m, s in [(2, 3, 2), (3, 2, 1), (4, 3, 2), (2, 2, 3)]:
    F = fam.starofhubs(h, m, s); adj = fam.build_starofhubs(h, m, s)
    print(F.label, cmp(F, adj, {"z": 0, "hub": 1, "mid": 2, "leaf": 3}))
