"""Closed-form (vertex-class) exact evaluation of D_v(k) and R1 loads
L_u(k) = sum_{v in N[u]} D_v(k)/(deg v + 1) for symmetric tree families.

Exact: all polynomial arithmetic is python-flint fmpz_poly (exact integers);
L_u compared via Fraction. Floats only in labelled diagnostic fields.
Internal consistency (asserted at every k): sum_v j^v_{k-1} = k i_k and
sum_v D_v(k) = k(k+1)(i_{k-1} i_{k+1} - i_k^2).
"""
from fractions import Fraction
from flint import fmpz_poly

X = fmpz_poly([0, 1])
ONE = fmpz_poly([1])


def co(p, k):
    return int(p[k]) if 0 <= k <= p.degree() else 0


def cdiv(a, b):
    return -((-a) // b)


class Fam:
    """classes: dict name -> (count, degree, J poly); nbhd: name -> {name: mult}"""
    def __init__(self, label, n, I, classes, nbhd):
        self.label, self.n, self.I, self.classes, self.nbhd = label, n, I, classes, nbhd
        assert sum(c for c, _, _ in classes.values()) == n, (label, "count")
        self.alpha = I.degree()
        self.lo = cdiv(n, 4)
        self.q = cdiv(2 * self.alpha - 1, 3)

    def level(self, k):
        I = self.I
        ik, ikm, ikp = co(I, k), co(I, k - 1), co(I, k + 1)
        D = {}
        s1 = 0
        for name, (cnt, deg, J) in self.classes.items():
            D[name] = k * ikm * co(J, k) - (k + 1) * ik * co(J, k - 1)
            s1 += cnt * co(J, k - 1)
        assert s1 == k * ik, (self.label, k, "double count")
        tot = sum(cnt * D[nm] for nm, (cnt, _, _) in self.classes.items())
        assert tot == k * (k + 1) * (ikm * ikp - ik * ik), (self.label, k, "sum D")
        L = {}
        Nabs = {}
        for u, nb in self.nbhd.items():
            L[u] = sum(Fraction(mult * D[v], self.classes[v][1] + 1) for v, mult in nb.items())
            Nabs[u] = sum(Fraction(mult * abs(D[v]), self.classes[v][1] + 1) for v, mult in nb.items())
        meanabs = Fraction(sum(cnt * abs(D[nm]) for nm, (cnt, _, _) in self.classes.items()), self.n)
        return D, L, Nabs, meanabs, tot

    def scan(self, ks=None):
        out = []
        ks = range(max(self.lo, 1), self.q + 1) if ks is None else ks
        for k in ks:
            D, L, Nabs, meanabs, tot = self.level(k)
            row = {"k": k, "lc_ok": tot <= 0,
                   "pv_fail": sorted(nm for nm in D if D[nm] > 0),
                   "r1_fail": sorted(u for u in L if L[u] > 0),
                   # diagnostics (floats): load normalised by the neighbourhood's own |D| mass, and by mean |D_v|
                   "r1_nbhd_norm": {u: (float(L[u] / Nabs[u]) if Nabs[u] else 0.0) for u in L},
                   "r1_mean_norm": {u: (float(L[u] / meanabs) if meanabs else 0.0) for u in L},
                   "pv_ratio": {nm: (float(Fraction(k * co(self.I, k - 1) * co(J, k), (k + 1) * co(self.I, k) * co(J, k - 1))) if co(J, k - 1) else None)
                                for nm, (_, _, J) in self.classes.items()}}
            out.append(row)
        return out


def AB(s):
    B = (1 + X) ** s
    return B + X, B


def hubstar(m, s):
    A, B = AB(s)
    I = A ** m + X * B ** m
    cl = {"hub": (1, m, B ** m), "mid": (m, s + 1, A ** (m - 1)),
          "leaf": (m * s, 1, (1 + X) ** (s - 1) * (A ** (m - 1) + X * B ** (m - 1)))}
    nb = {"hub": {"hub": 1, "mid": m}, "mid": {"mid": 1, "hub": 1, "leaf": s}, "leaf": {"leaf": 1, "mid": 1}}
    return Fam(f"H({m},{s})", 1 + m * (s + 1), I, cl, nb)


def twohub(m, s):
    A, B = AB(s)
    I = A ** (2 * m) + 2 * X * A ** m * B ** m
    cl = {"hub": (2, m + 1, B ** m * A ** m), "mid": (2 * m, s + 1, A ** (m - 1) * (A ** m + X * B ** m)),
          "leaf": (2 * m * s, 1, (1 + X) ** (s - 1) * (A ** (2 * m - 1) + X * A ** (m - 1) * B ** m + X * A ** m * B ** (m - 1)))}
    nb = {"hub": {"hub": 2, "mid": m}, "mid": {"mid": 1, "hub": 1, "leaf": s}, "leaf": {"leaf": 1, "mid": 1}}
    return Fam(f"TH({m},{s})", 2 * (1 + m * (s + 1)), I, cl, nb)


def starofhubs(h, m, s):
    """centre z adjacent to h hubs; each hub has m middles; each middle has s leaves."""
    A, B = AB(s)
    U = A ** m + X * B ** m
    Up = A ** (m - 1) + X * B ** (m - 1)
    I = U ** h + X * A ** (m * h)
    cl = {"z": (1, h, A ** (m * h)),
          "hub": (h, m + 1, B ** m * U ** (h - 1)),
          "mid": (h * m, s + 1, A ** (m - 1) * (U ** (h - 1) + X * A ** (m * (h - 1)))),
          "leaf": (h * m * s, 1, (1 + X) ** (s - 1) * (Up * U ** (h - 1) + X * A ** (m - 1) * A ** (m * (h - 1))))}
    nb = {"z": {"z": 1, "hub": h}, "hub": {"hub": 1, "z": 1, "mid": m},
          "mid": {"mid": 1, "hub": 1, "leaf": s}, "leaf": {"leaf": 1, "mid": 1}}
    return Fam(f"SH({h},{m},{s})", 1 + h * (1 + m * (s + 1)), I, cl, nb)


# explicit trees (adjacency lists) for independent recheck with pv_lib
def adj_hubstar(m, s):
    return _build([("hub", None)] + [])


def build_starofhubs(h, m, s):
    par = [-1]  # vertex 0 = centre
    for _ in range(h):
        hub = len(par); par.append(0)
        for _ in range(m):
            c = len(par); par.append(hub)
            for _ in range(s):
                par.append(c)
    return _adj(par)


def build_hubstar(m, s):
    par = [-1]
    for _ in range(m):
        c = len(par); par.append(0)
        for _ in range(s):
            par.append(c)
    return _adj(par)


def build_twohub(m, s):
    par = [-1, 0]
    for hub in (0, 1):
        for _ in range(m):
            c = len(par); par.append(hub)
            for _ in range(s):
                par.append(c)
    return _adj(par)


def _adj(par):
    n = len(par)
    adj = [[] for _ in range(n)]
    for i, p in enumerate(par):
        if p >= 0:
            adj[i].append(p); adj[p].append(i)
    return adj
