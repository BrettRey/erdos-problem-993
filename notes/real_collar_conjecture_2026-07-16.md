# The real-collar conjecture for tree independence polynomials
<!-- SUMMARY: Precise spectral conjectures distilled from the July 2026 campaign: a uniform real collar around the dominant zero for trees, an angle envelope, and an asymptotic bridge program to unimodality · status: conjecture note, literature-placed · updated: 2026-07-16 -->

Status: private research note. Statements are precise; support is
certified computation plus literature placement. No proof claims.

## 1. Known results (verified against sources 2026-07-16)

- **Csikvari (CPC 22/2013, "Note on the smallest root of the
  independence polynomial"):** for every graph G, the root of I(G,x) of
  smallest modulus, beta(G), is unique and real.
- **Prakash--Sharma (FSTTCS 2025, arXiv:2510.09197):** the gap from
  beta(G) to the next smallest root modulus is quantified for connected
  graphs; the general bound is exponentially small, (beta/n)^{O(n)},
  and their examples show paths have all-roots gap ratio only
  1 + Omega(1/n^2). Hence a class-uniform collar over ALL roots is
  false even for trees ~-- the crowding roots on paths are real.
- **Bencs--Csikvari (SIAM 2023, arXiv:2204.04868), Peters--Regts, de
  Boer et al. (arXiv:2111.06451):** zero-locus technology for
  bounded-degree graphs and tree recurrences; zeros are dense in the
  complement of the zero-free region.
- **Hibi--Kara--Vien (arXiv:2604.18824, April 2026):** constructs tree
  families (whisker/caterpillar-related) with symmetric AND unimodal
  independence polynomials, engaging the AMSE conjecture from the
  gamma-positivity side. Complementary to the spectral frame here; no
  collar-type statement.

## 2. Conjectures, after same-day kill-testing

**Refuted (2026-07-16, same session): the constant real collar.** The
first-draft conjecture ~-- a uniform delta > 0 with all zeros in
|z| <= (1 + delta) beta(T) real ~-- was killed within the hour by its
own kill-test (`scratch_collar_stress_20260716.py`). A random tree at
n = 72 already has a certified non-real pair at ratio 1.054, and
adversarial minimization drove the ratio to 1.01839 (n = 160, still
slowly descending, champion size drifting upward: the collar shrinks
like 1 + c/n, not a constant). Mechanism: two near-dominant REAL zeros
collide and split into a conjugate pair just off the axis; the
offending pair's angle is phi = 0.000021.

**Conjecture B' (cusp envelope, surviving).** For every tree, the
angular deviation phi of a non-real zero from the negative real axis is
bounded by a function of its dominance ratio that vanishes at 1:
phi <= F(|z|/beta) with F(1+) = 0. Measured points:
phi = 2.1e-5 at ratio 1.018; 0.054 at 1.2; 0.103 at 1.35; 0.181 at
1.5; 0.324 at 1.75; 0.449 at 2.0. Non-real zeros can approach beta
only along the real axis (collision cusps), so near-dominant zeros are
angularly trivial: their oscillation wavelength (2 pi / phi) vastly
exceeds any polynomial's degree.

**Conjecture C (positive-axis sector, robust at sampled scale).**
Non-real zeros of tree independence polynomials keep a sector around
the positive axis: measured minimum |arg z| = 0.90284 over all 28,998
certified trees, and a dedicated adversarial run (17,526 mutations,
n <= 160) failed to move that minimum at all ~-- the most rigid
invariant found in the entire program. Caveat from known theory: the
bounded-degree zero-locus results (Peters--Regts; de Boer et al.,
arXiv:2111.06451) suggest zeros of deep bounded-degree trees approach
the positive axis near the critical activity, so a FULLY uniform
sector over all trees is likely false asymptotically. The operative
conjecture is quantified: sigma = sigma(max degree, modulus window),
and the unimodality bridge needs the sector only along the saddle
contour for window-k coefficients. Reading task: extract from
arXiv:2111.06451 exactly where the tree zero-locus limit meets the
positive axis, and at what moduli.

None of these follow from Csikvari (first zero only) or Prakash--Sharma
(all-roots gap, no reality or angle information). The tree-specific
content is the coupling of angle to dominance (B') and the uniform
sector (C).

## 3. The discriminating control, and a corrected bridge

The split graph K_20 v E_6 has independent-set sequence
(1, 26, 15, 20, 15, 6, 1), a valley at k = 2, while its non-real zeros
sit at modulus ratios 30.0 and 63.0 ~-- far outside any collar. So
"collar implies unimodal" is FALSE as a pointwise claim: low-degree
valleys are controlled by the whole zero set, not the dominant zero.

The corrected program is asymptotic:

- Near-dominant zeros govern the TAIL third (alpha - k small), where
  coefficients are elementary symmetric functions of few reciprocal
  roots; the real collar should yield effective tail regularity.
- Mid-band coefficients (the open window k < 2 alpha / 3) are governed
  by saddle-point analysis: the oscillation suppressor there is the
  distance of non-real zeros from the positive real axis (the saddle
  contour), not from beta.
- Bridge target: collar (A) + envelope (B) + a positive-axis sector
  bound imply no valley for alpha >= alpha_0 effective; exhaustive
  verification (all trees n <= 29 unimodal) covers small cases. This
  is a program, not a theorem; the analytic step is Darboux/saddle
  asymptotics for polynomials with zero sets obeying A/B.

## 4. The tree-DP invariant, concretely (repo issue #1)

Work with the occupation ratio R_T(x) = I_root/E_root. Every rooted
tree is generated from the single vertex (R = x) by two maps:

    extend:  R  ->  x / (1 + R)          (new root above)
    merge:   (R_1, R_2)  ->  R_1 R_2 / x (identify roots)

and x0 is a zero of the full polynomial I(T; x) iff the reachable set
at x = x0 hits R = -1 (or E vanishes compatibly). So all three
spectral statements are reachability statements for this semigroup:

- Csikvari's theorem: the first x on the negative real ray where -1
  becomes reachable is a simple, isolated hit.
- Cusp envelope B': for x off the real axis but close to that first
  hit, -1 is NOT reachable unless arg(x) is tiny ~-- reachability
  leaks off the axis only tangentially.
- Sector C: for x in a sector around the positive axis (quantified by
  degree via the number of merges), -1 is never reachable.

Literature alignment (verified from the abstract of de Boer--Buys--
Guerini--Peters--Regts, arXiv:2111.06451, 2026-07-16): the zero-free
region for bounded degree "coincides with the normality region" of
exactly these occupation-ratio maps, the rescaled limit domain is a
cardioid, and they give "an exact formula for the boundary ... near
the positive real boundary point." So the invariant sought in issue #1
is, in this frame, a forward-invariant domain for extend/merge
avoiding -1 ~-- Asano/Lee--Yang contraction style. The genuinely new
work for #993 is (i) uniformity-in-degree bookkeeping for the merge
map (this is what breaks for general graphs), and (ii) the bridge:
converting the exact positive-boundary formula into a window
saddle-point bound at occupancy ~0.28.

Reading task (open): extract from the full paper the boundary
geometry near the negative real point (cusp confirmation for B') and
the exact positive-boundary formula with its modulus scaling in
Delta (input for the bridge inequality).

## 5. Kill-test results and next tests

Executed 2026-07-16 (`scratch_collar_stress_20260716.py`, 28,998
certified trees): the mass scan and adversarial minimization KILLED
the constant collar (see section 2) and left B' and C standing:
every near-dominant non-real pair found is angularly trivial
(cusp behavior), and the positive-axis sector minimum stayed at
0.90284 under both random and adversarial pressure.

Kill-tests 1 and 2 executed 2026-07-16, both survived rigidly:
1. B' attack (maximize angle within modulus band <= 1.2): 25,163
   adversarial mutations, zero improvement over the frontier value
   0.05404 (`results/cusp_attack_20260716.json`).
2. C attack (minimize |arg z|): 17,526 adversarial mutations, zero
   improvement over 0.90284 (`results/sector_attack_20260716.json`).
Contrast: the refuted collar conjecture moved within seconds under
identical adversarial pressure. Rigidity under adversarial search is
separating true invariants from artifacts.

Remaining kill-test:
3. Joint attack: minimize |arg z| subject to moderate modulus (the
   saddle-relevant zeros for window-k coefficients), ideally on
   high-degree deep trees where the Peters--Regts caveat predicts the
   sector must eventually narrow.

## 6. Provenance

Campaign data: notes/why_trees_resist_2026-07-16.md; DECISIONS.md
2026-07-15/16; results/pareto_root_frontier_20260716.json. Certified
arithmetic: python-flint / Arb; float64 root-finding is forbidden for
this purpose (claw-free control failure, DECISIONS.md 2026-07-16).
