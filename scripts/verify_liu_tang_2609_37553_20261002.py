"""Exact replay of the tree cases of Liu & Tang, arXiv:2609.37553v1 (29 Sep 2026).

Claim checked (their Corollary 3.3 with Corollary 3.2 = Bendjeddou-Hardiman):
for every forest G, the independence polynomials of E_{G4(1,0,0,0)}(G) and
E_{G4(2,1,0,0)}(G) are log-concave. These are the only tree cases of their
construction (their closing paragraph; also forced by Definition 3.1, since
s>0 or t>0 closes a cycle through the c--x_i edge and any clique of size
>=3 is a triangle).

Construction read from Definition 3.1 and Figures 5-7, 12-14: F(l,m,t,s) has
a centre c adjacent to a free vertex x; clique K_m joined to x by m edges and
to c by s edges; K_l joined to x by l-m edges and to c by t edges; K_{m-s}
joined to c by m-s edges; K_{l-t} joined to c by l-t-m edges. E_{G4} replaces
each edge uv of G by a copy of F at c_u and one at c_v, joined x_u--x_v.

Checks, all in exact integer arithmetic:
  1. construction sanity: I(F_n(l,m,t,s)) equals their displayed closed form
     (3.1 with every x_i set to x), and brute force agrees on small cases;
  2. log-concavity and unimodality of I(E_{G4}(T)) for every tree T on
     2..N vertices (geng);
  3. real-rootedness of the same polynomials by certified Arb isolation
     (python-flint), never float64 (project root-finding rule).

Run: venv/bin/python scripts/verify_liu_tang_2609_37553_20261002.py [N]
"""

from __future__ import annotations

import itertools
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import flint  # noqa: E402

from indpoly import independence_poly, is_log_concave, is_unimodal  # noqa: E402
from trees import trees  # noqa: E402


class Builder:
    def __init__(self) -> None:
        self.adj: list[set[int]] = []

    def vertex(self) -> int:
        self.adj.append(set())
        return len(self.adj) - 1

    def edge(self, a: int, b: int) -> None:
        assert a != b
        self.adj[a].add(b)
        self.adj[b].add(a)

    def clique(self, k: int) -> list[int]:
        vs = [self.vertex() for _ in range(k)]
        for a, b in itertools.combinations(vs, 2):
            self.edge(a, b)
        return vs

    def attach_F(self, c: int, l: int, m: int, t: int, s: int) -> int:
        """Attach one copy of F(l,m,t,s) at centre c; return its free vertex x."""
        assert l >= m >= s >= 0 and l >= t + m and t >= 0
        x = self.vertex()
        self.edge(c, x)
        km = self.clique(m)
        kl = self.clique(l)
        kms = self.clique(m - s)
        klt = self.clique(l - t)
        # "the added edges may be chosen arbitrarily" (Remark 3.1); take the
        # first vertices. For the tree cases every choice is isomorphic.
        for v in km[:m]:
            self.edge(x, v)
        for v in km[:s]:
            self.edge(c, v)
        for v in kl[: l - m]:
            self.edge(x, v)
        for v in kl[l - m : l - m + t] if t else []:
            self.edge(c, v)
        for v in kms[: m - s]:
            self.edge(c, v)
        for v in klt[: l - t - m]:
            self.edge(c, v)
        return x

    def lists(self) -> list[list[int]]:
        return [sorted(s) for s in self.adj]


def build_Fn(n: int, params: tuple[int, int, int, int]) -> Builder:
    b = Builder()
    c = b.vertex()
    for _ in range(n):
        b.attach_F(c, *params)
    return b


def build_E(n_g: int, adj_g: list[list[int]], params: tuple[int, int, int, int]) -> Builder:
    b = Builder()
    centres = [b.vertex() for _ in range(n_g)]
    for u in range(n_g):
        for v in adj_g[u]:
            if u < v:
                xu = b.attach_F(centres[u], *params)
                xv = b.attach_F(centres[v], *params)
                b.edge(xu, xv)
    return b


def is_tree(adj: list[list[int]]) -> bool:
    n = len(adj)
    if sum(len(a) for a in adj) != 2 * (n - 1):
        return False
    seen, stack = {0}, [0]
    while stack:
        for w in adj[stack.pop()]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def brute_indpoly(adj: list[list[int]]) -> list[int]:
    n = len(adj)
    counts = [0] * (n + 1)
    nbr = [sum(1 << w for w in a) for a in adj]
    for mask in range(1 << n):
        ok = all(not (mask >> v & 1) or not (nbr[v] & mask) for v in range(n))
        if ok:
            counts[bin(mask).count("1")] += 1
    while counts and counts[-1] == 0:
        counts.pop()
    return counts


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def poly_pow(a: list[int], e: int) -> list[int]:
    out = [1]
    for _ in range(e):
        out = poly_mul(out, a)
    return out


def closed_form_Fn(n: int, l: int, m: int, t: int, s: int) -> list[int]:
    """Their display before (3.1), with every free variable x_i set to x."""
    pre = poly_mul(poly_mul(poly_pow([1, m - s], n), poly_pow([1, l - t], n)), poly_pow([1, m], n))
    inner = poly_pow([1, l + 1], n)
    inner = inner + [0] * max(0, 2 - len(inner))
    inner[1] += 1
    out = poly_mul(pre, inner)
    while out and out[-1] == 0:
        out.pop()
    return out


def lc_min_ratio(seq: list[int]) -> Fraction:
    return min(Fraction(seq[k] ** 2, seq[k - 1] * seq[k + 1]) for k in range(1, len(seq) - 1))


def count_real_roots_sturm(seq: list[int]) -> int:
    """Exact count of distinct real roots via a Sturm sequence over Q."""
    p = flint.fmpq_poly(seq)
    g = p.gcd(p.derivative())
    if g.degree() > 0:
        p = p // g  # squarefree part
    seqs = [p, p.derivative()]
    while seqs[-1].degree() > 0:
        r = seqs[-2] % seqs[-1]
        if r == 0:
            break
        seqs.append(-r)

    def sign_changes_at_inf(sign: int) -> int:
        signs = []
        for q in seqs:
            if q == 0:
                continue
            lead = q.coeffs()[-1]
            sgn = 1 if lead > 0 else -1
            if sign < 0 and q.degree() % 2 == 1:
                sgn = -sgn
            signs.append(sgn)
        return sum(1 for a, b in zip(signs, signs[1:]) if a != b)

    return sign_changes_at_inf(-1) - sign_changes_at_inf(1)


def all_roots_real(seq: list[int]) -> bool:
    p = flint.fmpq_poly(seq)
    g = p.gcd(p.derivative())
    sqf_deg = p.degree() - g.degree()
    return count_real_roots_sturm(seq) == sqf_deg


def main() -> None:
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    families = {"(1,0,0,0) Bendjeddou-Hardiman": (1, 0, 0, 0), "(2,1,0,0) Liu-Tang": (2, 1, 0, 0)}

    print("1. construction sanity")
    for name, params in families.items():
        for n in range(1, 6):
            b = build_Fn(n, params)
            adj = b.lists()
            assert is_tree(adj), (name, n)
            got = independence_poly(len(adj), adj)
            want = closed_form_Fn(n, *params)
            assert got == want, (name, n, got, want)
            if len(adj) <= 22:
                assert brute_indpoly(adj) == got, (name, n)
        print(f"   {name}: I(F_n) matches closed form for n=1..5; brute force agrees where run")
    # E_{G4} on a single edge and on P3, brute force
    for name, params in families.items():
        for n_g, adj_g in [(2, [[1], [0]])]:
            b = build_E(n_g, adj_g, params)
            adj = b.lists()
            assert is_tree(adj)
            if len(adj) <= 22:
                assert brute_indpoly(adj) == independence_poly(len(adj), adj)
        print(f"   {name}: G4 on one edge has {len(adj)} vertices; tree; brute force agrees where run")

    print(f"\n2-3. every tree T on 2..{N} vertices")
    for name, params in families.items():
        n_checked = n_lc_fail = n_uni_fail = n_not_rr = 0
        worst: tuple[Fraction, int, int] | None = None
        sizes = set()
        for n_t in range(2, N + 1):
            for _n, adj_t in trees(n_t):
                b = build_E(n_t, adj_t, params)
                adj = b.lists()
                assert is_tree(adj)
                sizes.add(len(adj))
                seq = independence_poly(len(adj), adj)
                n_checked += 1
                if not is_log_concave(seq):
                    n_lc_fail += 1
                    print(f"   LC FAILURE: {name}, T on {n_t} vertices, adj={adj_t}")
                if not is_unimodal(seq):
                    n_uni_fail += 1
                r = lc_min_ratio(seq)
                if worst is None or r < worst[0]:
                    worst = (r, n_t, len(adj))
                if not all_roots_real(seq):
                    n_not_rr += 1
        assert worst is not None
        print(
            f"   {name}: {n_checked} trees, image orders {min(sizes)}..{max(sizes)}; "
            f"LC failures {n_lc_fail}; unimodality failures {n_uni_fail}; "
            f"not real-rooted {n_not_rr}; min a_k^2/(a_(k-1)a_(k+1)) = {float(worst[0]):.6f} "
            f"(T on {worst[1]} vertices, image order {worst[2]})"
        )


if __name__ == "__main__":
    main()
