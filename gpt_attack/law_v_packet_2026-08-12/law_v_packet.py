#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LAW-V / LAW-W PROOF PACKET  (Erdos Problem 993 campaign, 2026-08-12)
=====================================================================

You are being handed a frozen, graded proof target. Read this entire
docstring before doing anything. Run `python3 law_v_packet.py` first:
the self-test must pass before you trust the engine, and your final
deliverable must pass it too.

THE SETTING
-----------
P_n denotes the independence polynomial of the n-vertex path, as an
ascending integer coefficient list: P_{-1} = P_0 = [1], P_1 = [1,1],
P_n = P_{n-1} + x*P_{n-2}. For finite multisets a, b of positive
integers (pendant "arm" lengths on two hubs) define

    Q_a = prod_i P_{a_i}      R_a = prod_i P_{a_i - 1}
    B   = Q_a*Q_b + x*(R_a*Q_b + Q_a*R_b)
    C   = Q_a*Q_b + x*R_a*R_b

B is the independence polynomial of the tree with two adjacent hubs
carrying pendant paths a and b; C is that of the spider obtained by
contracting the hub edge. All coefficients are positive integers on
their supports.

THE TARGETS (single-index forms; out-of-range coefficients read 0)
------------------------------------------------------------------
    U_k = c_k b_k - c_{k-1} b_{k+1}
    V_k = c_k b_k - c_{k+1} b_{k-1}
    W_k = 2 c_{k-1} b_k - c_k b_{k-1} - c_{k-2} b_{k+1}

LAW-V:  V_k >= 0 for every arm pair (a, b) and every k.
LAW-W:  W_k >= 0 for every arm pair (a, b) and every k.

WHY THEY MATTER: a proved assembly lemma (see TOOLKIT) shows that
C's log-concavity (a published theorem: all spiders are strongly
log-concave, Li-Li-Yang-Zhang arXiv:2501.04245) plus LAW-V plus
U_m >= 0 at the relevant index plus LAW-W implies both partial-
synchronization base relations C ~p B and xC ~p B of the two-hub
program. Those base relations imply log-concavity of the two-hub tree
for EVERY connector length (notes/c2_bounded_pendant_core_2026-08-08).
Proving LAW-V and LAW-W (plus the known depth<=3 head analysis of U)
would deliver log-concavity, hence unimodality, for every tree with at
most two branch vertices: the C2 theorem.

EVIDENCE STATE (all exact arithmetic, replayable)
-------------------------------------------------
- Zero violations of LAW-V or LAW-W over: all arm pairs of total
  weight <= 24 (163,523 pairs); ~8,000 structured/random pairs to
  weight 200 with hill-climbing (margins plateau ~3.5%, never near 0);
  16,525-pair adversarial sweep to weight 6,002.
- U is NOT universally nonnegative: it fails only at reflected depth
  N - k <= 3 (N = deg B), never deeper, in every instance checked.
- The global tightest margin in the whole space is the top-corner
  xC ~p B inequality for arms (1,1) x (2h+1,2h+1): margin exactly
  h + 2 (proved for all h via polynomial certificates).
- WARNING from the wider campaign: the analogous forms for general
  tree edge contraction FAIL at unbounded depth when the tree pieces
  are non-log-concave. The C2 setting avoids this because all pieces
  here are products of path polynomials (log-concave, real-rooted).
  Do not "generalize" the target beyond path pieces; that version is
  refuted (notes/contraction_band_law_2026-08-11.md).

GRADED SUB-TARGETS (climb in order; each is valuable alone)
-----------------------------------------------------------
G0  Verify the engine: run this file; the self-test must pass.
G1  Single-hub core: prove V_k >= 0 for all k for the pair
    (C, B) = (Qa + x*Ra, Qa + x*Qa') for a single arm a = (m), i.e.
    arms a=(m), b=(). Closed forms via P-recurrences are available.
G2  Bottleneck family: arms (1,1) x (2h+1,2h+1). Prove U_k, V_k,
    W_k >= 0 for ALL (h, k). Head window (reflected depth <= 4 for
    sync margins, <= 7 for B's Turan gaps) is already proved by
    polynomial-in-h Sturm certificates; the body is open. Clean to
    h = 80 empirically, including U.
G3  One-hub-trivial family: a = (), b arbitrary. Then B = Qb + x(Qb +
    Rb), C = Qb + x*Rb. Prove LAW-V here.
G4  LAW-V in full (all arm pairs, all k).
G5  LAW-W in full.
G6  (Stretch) The U head strip: characterize exactly which arm pairs
    have U-failures and prove they are confined to reflected depth
    <= 3.

REFUTATION IS A SUCCESS MODE. A verified counterexample to LAW-V or
LAW-W at ANY weight — confirmed by this engine's exact arithmetic and
by rebuilding the two actual trees with an independent DP — is worth
as much as a proof. The adversarial searches above are evidence, not
proof, and the searches were budget-capped.

TOOLKIT (all statements sourced; do not re-derive from memory)
--------------------------------------------------------------
T1  HWZZ (Math. Inequal. Appl. 20 (2017) 91-103, on disk at
    notes/literature/hu_wang_zhao_zhao_2017_mia_partial_sync.txt):
    partial synchronicity (their Def 2.4) is closed under nonnegative
    linear combinations on both sides (Thm 3.5) and under common
    convolution by a log-concave no-internal-zero sequence (Thm 3.6);
    A ~p B implies uA + vB log-concave (Lemma 3.4).
T2  Assembly lemma (proved, notes/c2_single_index_laws_2026-08-11.md):
    with C log-concave, V_n >= 0 and U_m >= 0 give the (C,B) partial-
    sync inequality at every (m, n); with m >= n+1 the same gives the
    (xC, B) inequality, whose diagonal is exactly W.
T3  TP2 reformulation + K2 closure (notes/karlin_lr_order_toolkit_
    2026-08-12.md): LAW-V for (C,B) is consecutive-TP2 of the 2-row
    array [C; xB]; the property is preserved by common convolution
    with any log-concave no-internal-zero factor (Cauchy-Binet).
    Consequence: the failure geometry of the bilinear "dirty" term
    V(QaQb, xRaQb) equals that of the single-hub pair (Qa, xRa),
    transferred by the common factor Qb. G1 is therefore load-bearing.
T4  Bilinear localization (results/c2_single_index_laws_20260811.json):
    with C = g1 + g2, B = g1 + h2 + h3 (g1 = QaQb, g2 = xRaRb,
    h2 = xRaQb, h3 = xQaRb), termwise positivity FAILS: for V the
    dirty terms are (g1,h2), (g1,h3); for W everything except (g1,g1)
    is dirty somewhere. Any proof must dominate the dirty terms with
    the clean ones; pairwise-only arguments are refuted.
T5  Head certificates for G2 (scripts/bottleneck_head_certificates_
    20260811.py): reflected coefficients of B, C at fixed depth are
    explicit polynomials in h; all head margins positive for all h by
    exact Sturm isolation, reproducing the c2 note's h+2 margin.
T6  Path-polynomial identities: P_n = P_{n-1} + x P_{n-2};
    [x^k] P_n = binomial(n+1-k, k); Cassini-type identity
    P_{n-1} P_{n+1} - P_n^2 = (-1)^n x^{n+1} (verified in self-test).
    Products of P's are real-rooted (paths are claw-free), hence
    log-concave with no internal zeros: every K2/T1 hypothesis holds
    for every polynomial in this packet.

BINDING VERIFICATION RULES
--------------------------
1. Exact integer or exact rational arithmetic only. No floating-point
   root-finding, no numpy.roots, anywhere, for any claim.
2. Every claimed identity or inequality family must come with either
   a complete proof or a finite verification script using THIS
   engine, and the scope of finite evidence must be stated as such.
3. Label every statement THEOREM (proved here), CITED (with source),
   or EVIDENCE (finite computation). Do not blur the labels.
4. Self-grade at the end: FULL / PARTIAL / FAILED, with one sentence
   per graded target. An honest PARTIAL with one real proved step
   beats a grandiose FAILED-in-disguise. Overclaiming is the only
   disqualifying failure.

DELIVERABLE
-----------
A single markdown or text document containing: your self-grade table;
proofs (or refutations) for the targets you reached; and a runnable
verification script for every finite claim. Keep everything replayable
against this file's engine.
"""

# ======================================================================
# EXACT ENGINE (standard library only)
# ======================================================================

from fractions import Fraction
from itertools import combinations_with_replacement


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return out


def padd(a, b):
    m = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(m)]


_P = {}


def P(n):
    """Independence polynomial of the n-vertex path, ascending coeffs."""
    if n in _P:
        return _P[n]
    if n <= 0:
        return [1]
    a, b = [1], [1, 1]
    for _ in range(n - 1):
        a, b = b, padd(b, [0] + a)
    _P[n] = b
    return b


def cf(p, i):
    return p[i] if 0 <= i < len(p) else 0


def prodP(arms, shift=0):
    out = [1]
    for a in arms:
        out = pmul(out, P(a - shift))
    return out


def BC(a, b):
    """The pair (B, C) for arm multisets a, b."""
    Qa, Ra = prodP(a), prodP(a, 1)
    Qb, Rb = prodP(b), prodP(b, 1)
    QQ = pmul(Qa, Qb)
    B = padd(QQ, [0] + padd(pmul(Ra, Qb), pmul(Qa, Rb)))
    C = padd(QQ, [0] + pmul(Ra, Rb))
    return B, C


def UVW(a, b):
    """Lists of (k, value) for U, V, W over the full index range."""
    B, C = BC(a, b)
    N = len(B) - 1
    U = [(k, cf(C, k) * cf(B, k) - cf(C, k - 1) * cf(B, k + 1))
         for k in range(N + 2)]
    V = [(k, cf(C, k) * cf(B, k) - cf(C, k + 1) * cf(B, k - 1))
         for k in range(N + 2)]
    W = [(k, 2 * cf(C, k - 1) * cf(B, k) - cf(C, k) * cf(B, k - 1)
          - cf(C, k - 2) * cf(B, k + 1)) for k in range(N + 2)]
    return N, U, V, W


def violations(a, b):
    """Reflected depths of negative U, V, W entries."""
    N, U, V, W = UVW(a, b)
    return {"U": [N - k for k, v in U if v < 0],
            "V": [N - k for k, v in V if v < 0],
            "W": [N - k for k, v in W if v < 0]}


def psync_ok(A, B):
    """Full partial-synchronization check (HWZZ Def 2.4), exact."""
    D = max(len(A), len(B))
    for m in range(D + 1):
        for n in range(m + 1):
            lhs = cf(A, m) * cf(B, n) + cf(A, n) * cf(B, m)
            rhs = (cf(A, m + 1) * cf(B, n - 1)
                   + cf(A, n - 1) * cf(B, m + 1))
            if lhs < rhs:
                return False
    return True


def tree_ipoly(n, adj, root=0):
    """Independent DP for verification of any counterexample claim."""
    parent = [-1] * n
    order, seen, st = [], [False] * n, [root]
    seen[root] = True
    while st:
        v = st.pop()
        order.append(v)
        for u in adj[v]:
            if not seen[u]:
                seen[u] = True
                parent[u] = v
                st.append(u)
    f, g = [None] * n, [None] * n
    for v in reversed(order):
        fv, gv = [1], [0, 1]
        for u in adj[v]:
            if parent[u] == v:
                fv = pmul(fv, padd(f[u], g[u]))
                gv = pmul(gv, f[u])
        f[v], g[v] = fv, gv
    return padd(f[root], g[root])


def two_hub_tree(a, b):
    """Adjacency list of the two-hub tree for arms a, b (B's tree)."""
    adj = [[], []]
    adj[0].append(1)
    adj[1].append(0)

    def arm(hub, length):
        prev = hub
        for _ in range(length):
            adj.append([])
            v = len(adj) - 1
            adj[prev].append(v)
            adj[v].append(prev)
            prev = v
    for L in a:
        arm(0, L)
    for L in b:
        arm(1, L)
    return adj


def spider_tree(arms):
    """Single hub with the given arms (C's tree)."""
    adj = [[]]

    def arm(length):
        prev = 0
        for _ in range(length):
            adj.append([])
            v = len(adj) - 1
            adj[prev].append(v)
            adj[v].append(prev)
            prev = v
    for L in arms:
        arm(L)
    return adj


# ======================================================================
# SELF-TEST: grounded anchors. All must pass before and after your work.
# ======================================================================

def self_test():
    ok = True

    def check(label, cond):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + label)
        ok = ok and cond

    # A1: B and C really are the tree polynomials (independent DP).
    for a, b in [([1, 1], [3, 3]), ([2], [1, 4]), ([], [1, 3]),
                 ([1, 2, 3], [2, 2])]:
        Bp, Cp = BC(a, b)
        adj = two_hub_tree(a, b)
        check(f"A1 B==I(two-hub tree) for {a}x{b}",
              tree_ipoly(len(adj), adj) == Bp)
        sadj = spider_tree(list(a) + list(b))
        check(f"A1 C==I(spider) for {a}x{b}",
              tree_ipoly(len(sadj), sadj) == Cp)

    # A2: Cassini identity for path polynomials (T6).
    for n in range(1, 12):
        lhs = padd(pmul(P(n - 1), P(n + 1)),
                   [-x for x in pmul(P(n), P(n))])
        rhs = [0] * (n + 1) + [(-1) ** n]
        m = max(len(lhs), len(rhs))
        check(f"A2 Cassini n={n}",
              [cf(lhs, i) for i in range(m)] == [cf(rhs, i) for i in range(m)])
        if n >= 4:
            break  # spot-check low n plus one larger
    lhs = padd(pmul(P(9), P(11)), [-x for x in pmul(P(10), P(10))])
    check("A2 Cassini n=10",
          [cf(lhs, i) for i in range(12)] == [0] * 11 + [1])

    # A3: the (1,1)x(9,9) bottleneck margin is exactly 6 (h=4, h+2).
    # Top-corner diagonal of the (xC, B) partial-sync form sits at
    # absolute index m = N = deg B, where it equals
    # 2 C_{N-1} - B_{N-1} = h + 2 (since B_N = C_N = 1).
    Bp, Cp = BC([1, 1], [9, 9])
    N = len(Bp) - 1
    A = [0] + Cp
    m_ = N
    margin = (2 * cf(A, m_) * cf(Bp, m_)
              - cf(A, m_ + 1) * cf(Bp, m_ - 1)
              - cf(A, m_ - 1) * cf(Bp, m_ + 1))
    check("A3 (1,1)x(9,9) top-corner margin == 6 (== h+2 at h=4)",
          margin == 6)

    # A4: weight<=14 exhaustive re-check: no U/V/W violation deeper
    # than reflected depth 3, and V, W clean everywhere.
    deepest, vbad, wbad, npairs = -1, 0, 0, 0
    arms_pool = []
    for r in range(0, 4):
        for arms in combinations_with_replacement(range(1, 13), r):
            if sum(arms) <= 12:
                arms_pool.append(list(arms))
    for i, a in enumerate(arms_pool):
        for b in arms_pool:
            if sum(a) + sum(b) == 0 or sum(a) + sum(b) > 14:
                continue
            npairs += 1
            v = violations(a, b)
            if v["V"]:
                vbad += 1
            if v["W"]:
                wbad += 1
            for f in "UVW":
                if v[f]:
                    deepest = max(deepest, max(v[f]))
    check(f"A4 weight<=14 ({npairs} pairs): V clean", vbad == 0)
    check(f"A4 weight<=14: W clean", wbad == 0)
    check(f"A4 weight<=14: deepest U depth <= 3 (got {deepest})",
          deepest <= 3)

    # A5: the assembly implication holds end-to-end on samples:
    # where U, V, W are clean, both partial-sync relations hold.
    for a, b in [([1, 1], [5, 5]), ([2, 3], [4]), ([1], [2, 2, 2])]:
        v = violations(a, b)
        Bp, Cp = BC(a, b)
        if not v["U"] and not v["V"] and not v["W"]:
            check(f"A5 psync C~B and xC~B for {a}x{b}",
                  psync_ok(Cp, Bp) and psync_ok([0] + Cp, Bp))

    print()
    print("SELF-TEST " + ("PASSED" if ok else "FAILED"))
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if self_test() else 1)
