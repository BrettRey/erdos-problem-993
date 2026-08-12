#!/usr/bin/env python3
"""Atlas-truncation kill-test: bottom-level (Hyp) for tree independence
complexes.

Chan-Pak's combinatorial atlas (arXiv:2203.01533 expository /
arXiv:2110.10740 technical; full texts in the research vault) certifies
matroid log-concavity through a local-global principle whose entire
direct verification burden sits at the sink vertices (alpha, 0, 1),
carrying the matrix A(alpha, 1): the bordered "compatible-extensions"
matrix with entries

    A[x][y]    = 1 if alpha+x+y stays in the language (x != y)
    A[x][null] = 1 for every extendable x,   A[null][null] = 1.

Their Lemma 4.6 (expository paper) proves every sink satisfies (Hyp)
(= at most ONE strictly positive eigenvalue) via an equivalence
relation: x ~ y iff alpha+x+y leaves the language, whose TRANSITIVITY
is supplied by the matroid exchange property. The support-restricted
matrix is then a bordered complete multipartite matrix, which has one
positive eigenvalue.

For a GRAPH independence complex the same construction makes sense
verbatim, but the conflict relation at the empty word is graph
ADJACENCY, which is transitive iff the graph is a disjoint union of
cliques. A tree containing any path on 3 vertices is not, so the
bottom-level certificate fails at (empty, 0, 1) - a vertex contained in
the atlas of EVERY certifying index k. Index truncation therefore
cannot rescue the atlas for trees; the binding obstruction is the
missing exchange/transitivity structure, not the all-index
quantification.

This script verifies, in exact arithmetic:
  1. harness validity on matroid controls (free matroid; two parallel
     classes -> bordered complete multipartite): exactly 1 strictly
     positive eigenvalue;
  2. failure for paths, the claw, a small spider;
  3. an exhaustive sweep of all trees on 3..9 vertices: A(empty,1)
     satisfies (Hyp) for NO tree except those with <= 2 vertices'
     worth of conflict structure, and (Hyp)-holding coincides exactly
     with the conflict graph being a cluster graph;
  4. the count subtlety: strictly positive eigenvalue counts strip
     zero roots exactly (sympy count_roots includes interval
     endpoints, and these bordered matrices are rank-deficient).

Verdict recorded in notes/atlas_truncation_killtest_2026-08-12.md.
"""
import subprocess
import sympy as sp

GENG = "/opt/homebrew/bin/geng"


def strict_pos_eigs(A):
    M = sp.Matrix(A)
    cp = M.charpoly()
    p = sp.Poly(cp.as_expr(), cp.gen)
    coeffs = p.all_coeffs()[::-1]
    z = 0
    while z < len(coeffs) and coeffs[z] == 0:
        z += 1
    q = sp.Poly(coeffs[z:][::-1], p.gen)
    n = q.count_roots(0, sp.oo)
    if q.eval(0) == 0:
        n -= 1
    return n


def bordered_from_pairs(n, indep_pairs):
    A = [[0] * (n + 1) for _ in range(n + 1)]
    for x, y in indep_pairs:
        A[x][y] = A[y][x] = 1
    for x in range(n):
        A[x][n] = A[n][x] = 1
    A[n][n] = 1
    return A


def A_empty_1(n, edges):
    """A(empty,1) for the independence complex of graph (n, edges)."""
    es = {(min(a, b), max(a, b)) for a, b in edges}
    pairs = [(x, y) for x in range(n) for y in range(x + 1, n)
             if (x, y) not in es]
    return bordered_from_pairs(n, pairs)


def is_cluster_graph(n, edges):
    """Adjacency transitive <=> disjoint union of cliques."""
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    for x in range(n):
        for y in adj[x]:
            for z in adj[y]:
                if z != x and z not in adj[x]:
                    return False
    return True


def parse_graph6(line):
    bs = [ord(ch) - 63 for ch in line.strip()]
    n = bs[0]
    bits = []
    for x in bs[1:]:
        for i in range(5, -1, -1):
            bits.append((x >> i) & 1)
    edges = []
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                edges.append((i, j))
            idx += 1
    return n, edges


def main():
    ok = True

    def check(label, cond):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + label)
        ok = ok and cond

    # 1. matroid controls
    free = bordered_from_pairs(
        5, [(x, y) for x in range(5) for y in range(x + 1, 5)])
    check("free matroid control: exactly 1 positive eigenvalue",
          strict_pos_eigs(free) == 1)
    k33 = bordered_from_pairs(
        6, [(x, y) for x in range(3) for y in range(3, 6)])
    check("two-parallel-class control (bordered K_{3,3}): exactly 1",
          strict_pos_eigs(k33) == 1)

    # 2. named tree failures
    p6 = A_empty_1(6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)])
    check("P6 fails (Hyp): 3 positive eigenvalues",
          strict_pos_eigs(p6) == 3)
    claw = A_empty_1(4, [(0, 1), (0, 2), (0, 3)])
    check("claw fails (Hyp): 2 positive eigenvalues",
          strict_pos_eigs(claw) == 2)
    spider = A_empty_1(
        7, [(0, 1), (1, 2), (0, 3), (3, 4), (0, 5), (5, 6)])
    check("spider(2,2,2) fails (Hyp): 2 positive eigenvalues",
          strict_pos_eigs(spider) == 2)

    # 3. exhaustive sweep n = 3..9: (Hyp) <=> cluster conflict graph
    total = fails = mismatches = 0
    for n in range(3, 10):
        out = subprocess.run(
            [GENG, "-q", "-c", str(n), f"{n-1}:{n-1}"],
            capture_output=True, text=True).stdout
        for line in out.splitlines():
            nn, edges = parse_graph6(line)
            total += 1
            hyp = strict_pos_eigs(A_empty_1(nn, edges)) <= 1
            cluster = is_cluster_graph(nn, edges)
            if not hyp:
                fails += 1
            if hyp != cluster:
                mismatches += 1
    check(f"all {total} trees n=3..9: (Hyp) at (empty,0,1) <=> cluster "
          f"conflict graph (mismatches {mismatches})", mismatches == 0)
    check(f"every tree n=3..9 FAILS bottom-level (Hyp) "
          f"({fails}/{total} fail; connected non-clique trees cannot "
          f"be cluster graphs)", fails == total)

    print()
    print("KILL-TEST " + ("CONFIRMED" if ok else "INCONCLUSIVE"))
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
