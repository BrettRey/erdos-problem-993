"""Exact tools for the Vatter-order certificate on trees (follow-up to arXiv:2608.22147).

Fix a vertex order. Classify independent sets by their least vertex:
    i_{k+1}(T) = sum_v i_k(T[L(v)]),   L(v) = later vertices not adjacent to v.
Hence rho_{k+1}(T) is the i_{k-1}(L(v))-weighted average of rho_k(L(v)), with
rho_k = i_k / i_{k-1}, and T is log-concave at k whenever every tail with
i_{k-1}(L(v)) > 0 has rho_k(L(v)) <= rho_k(T). An order with that property at k
is a Vatter order at k.

Everything here is exact (python-flint fmpz_poly). Forest polynomials are built
component by component, with component polynomials cached by an AHU canonical
form and repeated components combined by powering, so symmetric constructions
with hundreds of vertices are cheap.
"""

from __future__ import annotations

from collections import Counter
from math import ceil

from flint import fmpz_poly

X = fmpz_poly([0, 1])
ONE = fmpz_poly([1])


# ---------------------------------------------------------------- families

class Tree:
    def __init__(self, n: int):
        self.adj: list[list[int]] = [[] for _ in range(n)]

    def add(self, a: int, b: int) -> None:
        self.adj[a].append(b)
        self.adj[b].append(a)

    @property
    def n(self) -> int:
        return len(self.adj)


def _grow(adj: list[list[int]], parent: int) -> int:
    adj.append([])
    v = len(adj) - 1
    if parent >= 0:
        adj[v].append(parent)
        adj[parent].append(v)
    return v


def hubstar(m: int, s: int) -> list[list[int]]:
    """H(m,s): hub joined to m supports, each carrying s leaves (n = 1+m+ms)."""
    adj: list[list[int]] = []
    hub = _grow(adj, -1)
    for _ in range(m):
        sup = _grow(adj, hub)
        for _ in range(s):
            _grow(adj, sup)
    return adj


def msh(ms: list[int], s: int) -> list[list[int]]:
    """Centre joined to hubs; hub i carries ms[i] supports, each with s leaves.
    MSH(h;m,s) = msh([m]*h, s); S(h,m,s) of the R1 record is the same tree."""
    adj: list[list[int]] = []
    c = _grow(adj, -1)
    for m in ms:
        hub = _grow(adj, c)
        for _ in range(m):
            sup = _grow(adj, hub)
            for _ in range(s):
                _grow(adj, sup)
    return adj


def spider(legs: list[int]) -> list[list[int]]:
    adj: list[list[int]] = []
    c = _grow(adj, -1)
    for l in legs:
        prev = c
        for _ in range(l):
            prev = _grow(adj, prev)
    return adj


def path(n: int) -> list[list[int]]:
    return spider([n - 1]) if n > 1 else [[]]


def caterpillar(spine: int, legs: int) -> list[list[int]]:
    adj: list[list[int]] = []
    prev = -1
    for _ in range(spine):
        v = _grow(adj, prev)
        for _ in range(legs):
            _grow(adj, v)
        prev = v
    return adj


# ---------------------------------------------------------- polynomials

_comp_cache: dict[str, fmpz_poly] = {}
_pow_cache: dict[tuple[str, int], fmpz_poly] = {}


def _rooted_dp(adj, alive, root) -> tuple[fmpz_poly, str]:
    """Return (I(component), canonical form) for the component of root."""
    order, parent = [root], {root: -1}
    i = 0
    while i < len(order):
        v = order[i]
        i += 1
        for w in adj[v]:
            if w in alive and w not in parent:
                parent[w] = v
                order.append(w)
    comp = order
    # centre(s) by leaf peeling
    deg = {v: sum(1 for w in adj[v] if w in alive) for v in comp}
    layer = [v for v in comp if deg[v] <= 1]
    remaining = len(comp)
    removed = set()
    while remaining > 2:
        nxt = []
        for v in layer:
            removed.add(v)
            remaining -= 1
            for w in adj[v]:
                if w in alive and w not in removed:
                    deg[w] -= 1
                    if deg[w] == 1:
                        nxt.append(w)
        layer = nxt
    centres = [v for v in comp if v not in removed]

    def encode(r: int, block: int) -> str:
        # iterative AHU encoding of the tree rooted at r, not crossing to block
        par = {r: -1}
        seq = [r]
        j = 0
        while j < len(seq):
            v = seq[j]
            j += 1
            for w in adj[v]:
                if w in alive and w != block and w not in par:
                    par[w] = v
                    seq.append(w)
        code: dict[int, str] = {}
        for v in reversed(seq):
            kids = sorted(code[w] for w in adj[v] if w in alive and par.get(w) == v and w != block)
            code[v] = "(" + "".join(kids) + ")"
        return code[r]

    if len(centres) == 1:
        canon = encode(centres[0], -1)
    else:
        a, b = centres
        ea, eb = encode(a, b), encode(b, a)
        canon = "E" + min(ea + eb, eb + ea)
    if canon in _comp_cache:
        return _comp_cache[canon], canon
    out: dict[int, fmpz_poly] = {}
    inn: dict[int, fmpz_poly] = {}
    for v in reversed(comp):
        ex, ic = ONE, X
        for w in adj[v]:
            if w in alive and parent.get(w) == v:
                ex = ex * (out[w] + inn[w])
                ic = ic * out[w]
        out[v], inn[v] = ex, ic
    poly = out[root] + inn[root]
    _comp_cache[canon] = poly
    return poly, canon


def forest_poly(adj, alive) -> fmpz_poly:
    alive = set(alive)
    seen: set[int] = set()
    counts: Counter = Counter()
    for v in alive:
        if v in seen:
            continue
        poly, canon = _rooted_dp(adj, alive, v)
        # mark the component
        stack = [v]
        seen.add(v)
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w in alive and w not in seen:
                    seen.add(w)
                    stack.append(w)
        counts[canon] += 1
    total = ONE
    for canon, c in counts.items():
        key = (canon, c)
        if key not in _pow_cache:
            _pow_cache[key] = _comp_cache[canon] ** c
        total = total * _pow_cache[key]
    return total


def coeff(p: fmpz_poly, k: int) -> int:
    return int(p[k]) if 0 <= k <= p.degree() else 0


def ratio_le(small: fmpz_poly, big: fmpz_poly, k: int) -> bool | None:
    """rho_k(small) <= rho_k(big); None if small has i_{k-1} = 0 (zero weight)."""
    a1, a0 = coeff(small, k), coeff(small, k - 1)
    if a0 == 0:
        return None
    return a1 * coeff(big, k - 1) <= coeff(big, k) * a0


def window(n: int, alpha: int) -> range:
    return range(ceil(n / 4), ceil((2 * alpha - 1) / 3) + 1)


# ------------------------------------------------------------- checks

def first_vertex_filter(adj, ks) -> dict[int, list[int]]:
    """For each k, the vertices v with rho_k(T - N[v]) <= rho_k(T) (possible first vertices)."""
    n = len(adj)
    I = forest_poly(adj, range(n))
    tails = {}
    for v in range(n):
        tails[v] = forest_poly(adj, set(range(n)) - {v} - set(adj[v]))
    return {k: [v for v in range(n) if ratio_le(tails[v], I, k) is not False] for k in ks}


def order_tails(adj, order) -> list[fmpz_poly]:
    pos = {v: i for i, v in enumerate(order)}
    out = []
    for v in order:
        nb = set(adj[v])
        out.append(forest_poly(adj, [u for u in order if pos[u] > pos[v] and u not in nb]))
    return out


def check_order(adj, order, ks, assert_identity: bool = True) -> dict[int, list[int]]:
    """Return {k: [positions whose tail violates rho_k(tail) <= rho_k(T)]}."""
    n = len(adj)
    I = forest_poly(adj, range(n))
    tails = order_tails(adj, order)
    if assert_identity:
        for k in range(I.degree()):
            assert coeff(I, k + 1) == sum(coeff(t, k) for t in tails), ("identity", k)
    return {k: [i for i, t in enumerate(tails) if ratio_le(t, I, k) is False] for k in ks}


def bfs_order(adj, root: int) -> list[int]:
    order, seen = [root], {root}
    i = 0
    while i < len(order):
        v = order[i]
        i += 1
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                order.append(w)
    return order


def alpha(adj) -> int:
    return forest_poly(adj, range(len(adj))).degree()


def centre(adj) -> int:
    """A centre of the tree (first of the two if bicentral), by leaf peeling."""
    n = len(adj)
    deg = [len(a) for a in adj]
    layer = [v for v in range(n) if deg[v] <= 1]
    rem, gone = n, set()
    while rem > 2:
        nxt = []
        for v in layer:
            gone.add(v)
            rem -= 1
            for w in adj[v]:
                if w not in gone:
                    deg[w] -= 1
                    if deg[w] == 1:
                        nxt.append(w)
        layer = nxt
    return [v for v in range(n) if v not in gone][0]


def order_revbfs(adj) -> list[int]:
    """Leaves first: reverse of the BFS order from a centre."""
    return bfs_order(adj, centre(adj))[::-1]


def order_degasc(adj) -> list[int]:
    """Ascending degree, ties by vertex index."""
    return sorted(range(len(adj)), key=lambda v: (len(adj[v]), v))
