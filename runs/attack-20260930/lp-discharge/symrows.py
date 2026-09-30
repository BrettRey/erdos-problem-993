"""Closed-form defects for two symmetric families, via python-flint fmpz_poly.
H(m,s): hub -> m supports -> s leaves each.  Vertex types: hub, sup, leaf.
MSH(h;m,s): centre -> h hubs, each hub -> m supports -> s leaves. Types: cen, hub, sup, leaf.
J_type = I(T - N[v]) for a vertex v of that type. Returns I and {type: (J, degree, [(nbr type, count, nbr degree)])}."""
from flint import fmpz_poly
X = fmpz_poly([0, 1]); ONE = fmpz_poly([1]); P1 = fmpz_poly([1, 1])
def H(m, s):
    st = P1 ** s + X                              # rooted star K_{1,s} (support + s leaves)
    I = st ** m + X * P1 ** (s * m)
    J = {'hub': (P1 ** (s * m), m, [('sup', m, s + 1)]),
         'sup': (st ** (m - 1), s + 1, [('hub', 1, m), ('leaf', s, 1)]),
         'leaf': ((st ** (m - 1) + X * P1 ** (s * (m - 1))) * P1 ** (s - 1), 1, [('sup', 1, s + 1)])}
    return I, J, 1 + m + m * s
def MSH(h, m, s):
    st = P1 ** s + X
    EH, IH = st ** m, X * P1 ** (s * m)           # hub-star rooted at hub: hub excluded / included
    EHp, IHp = st ** (m - 1), X * P1 ** (s * (m - 1))
    I = (EH + IH) ** h + X * EH ** h
    J = {'cen': (st ** (h * m), h, [('hub', h, m + 1)]),
         'hub': ((EH + IH) ** (h - 1) * P1 ** (m * s), m + 1, [('cen', 1, h), ('sup', m, s + 1)]),
         'sup': (((EH + IH) ** (h - 1) + X * EH ** (h - 1)) * st ** (m - 1), s + 1, [('hub', 1, m + 1), ('leaf', s, 1)]),
         'leaf': (((EH + IH) ** (h - 1) * (EHp + IHp) + X * EH ** (h - 1) * EHp) * P1 ** (s - 1), 1, [('sup', 1, s + 1)])}
    return I, J, 1 + h * (1 + m + m * s)
def coeffs(p): return [int(c) for c in p.coeffs()]
