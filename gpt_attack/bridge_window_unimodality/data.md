# Certified empirical inputs (2026-07-15/16 campaign)

All values from exact big-integer computation and certified Arb root
isolation; replay scripts in the repo root and scripts/. Full context:
notes/why_trees_resist_2026-07-16.md,
notes/real_collar_conjecture_2026-07-16.md.

## Coefficient-space barrier (why T3 refutation-by-construction is hard)

- Across ~1.4M exact evaluations (bouquets, hybrids, dumbbells,
  depth-3/4 stacking, subtractive joins, 540k joins of the 21 exact
  n=26/28 LC-failure trees, free mutation, independent Codex search to
  n=3,000): ZERO window valleys; under a genuine prior descent
  (threshold 0.001) the rebound deficit never beat ~11/n, constant
  GROWING with n; every rebound has rise distance c-b = 1.
- The T_{3,M,N} late shoulder slides out of the window as it scales
  (rise index converges to L from below, then crosses).

## Root-space frontier (certified)

Angle phi = angular deviation from the negative real axis; ratio
r = |z|/beta(T). Measured Pareto frontier (85,114 certified trees):

  r     : 1.05  1.1   1.2    1.35   1.5    1.75   2.0
  phi*  : 0     0     0.054  0.103  0.181  0.324  0.449

Cardioid cusp prediction (2(r-1))^{3/2}/6 (from expanding
C_infty = {-u e^{-u}} at u=1; OUR derivation, verify):

  r     : 1.2    1.35   1.5    2.0
  pred  : 0.042  0.098  0.167  0.471

- Constant real collar REFUTED: adversarial minimization reaches
  non-real pairs at r = 1.018 (angle 2.1e-5): real-root collisions,
  angularly trivial (cusp behavior).
- Cusp envelope adversarially rigid: 25,163 mutations, band r <= 1.2,
  zero improvement over 0.054.
- Positive-axis sector: min |arg z| = 0.90284 over 28,998 trees;
  17,526 adversarial mutations, zero improvement. Finite-size only:
  BBP Remark 5.5 implies deep spherically regular trees break any
  uniform sector.

## Control

Split graph K_20 v E_6: sequence (1, 26, 15, 20, 15, 6, 1), valley at
k = 2, non-real zeros at r = 30.0 and 63.0. Pointwise spectral bridges
are false; only asymptotic ones are candidates.

## Literature files (pdftotext conversions in notes/literature/)

- arxiv_2111_06451.txt  (Bencs--Buys--Peters: limit zero locus,
  cardioid, Gamma curve Cor 5.4, Remark 5.5)
- arxiv_2204_04868.txt  (Bencs--Csikvari: zero-free criterions,
  angular profiles Thm 8.1)
- arxiv_2510_09197.txt  (Prakash--Sharma: root gap (beta/n)^{O(n)};
  paths 1+Omega(1/n^2))
- scott_sokal_0309352.txt (Shearer radius, LLL)
- hibi_kara_vien_2604_18824.txt (symmetric+unimodal tree families)
- li_2026_unimodality_two_families.txt (T_{3,m,n} unimodal)
- li_yang_zhang_zhang_2025_symmetric_function_log_concavity.txt
  (spiders unimodal via 2-s-positivity)
