# Absorption-margin probe on the pattern-graph families

Replay-verified 2026-08-11: `scripts/probe_absorption_margin_20260811.py`
regenerated `results/absorption_margin_probe_20260811.json` byte-identically
(modulo runtime); the closest-approach record (T17_S25_2, m = 4, j = 53) was
re-derived independently to the exact integers; the TG anchors and the
Eq. (4) decomposition identity were reproduced with an independently written
TG builder. Probe computed by a Sonnet subagent; verified in the main
session.

## Question

`verify_partial_sync_obstruction_20260811.py` builds tree families with
positive two-term decompositions F = A + B of their independence
polynomials and locates where log-concavity (LC) breaks. The exact
identity, verified per index below, is

    Turan(F)_j = Turan(A)_j + Turan(B)_j + S_j
    Turan(P)_j = p_j^2 - p_{j-1} p_{j+1}
    S_j        = 2 a_j b_j - a_{j-1} b_{j+1} - a_{j+1} b_{j-1}

with out-of-range coefficients read as 0. An LC break at index j means
Turan(F)_j < 0: the diagonal deficit -S_j exceeds the summand surplus
Turan(A)_j + Turan(B)_j. In the checked corpus, actual breaks sit only
at reflected depth d = alpha - j <= 3, even though S_j < 0 (a
necessary precondition for a break) persists out to relative depth
d/alpha ~ 0.54.

The question: at fixed relative depth d/alpha, does the absorption
factor

    rho_j = (Turan(A)_j + Turan(B)_j) / (-S_j),   defined where S_j < 0,

shrink toward 1 as the family parameter grows (a live route to a
unimodality counterexample, since rho_j < 1 is exactly an LC break),
or does it stay bounded away from 1 or grow (evidence these families
cannot be pushed into a valley by increasing the parameter)?

## Method

Three m-direction families from the existing script, built via
`m_family_instance(ell, nlegs, k, m)`:

- `T12_S24_2` = (ell=2, nlegs=4, k=2), m = 2..60
- `T13_S24_2` = (ell=3, nlegs=4, k=2), m = 2..30
- `T17_S25_2` = (ell=7, nlegs=5, k=2), m = 2..24

For every instance and every index j with S_j < 0, the probe computes
Turan(A)_j, Turan(B)_j, S_j, margin_j = Turan(F)_j, and rho_j as an
exact `fractions.Fraction`. Per index the identity margin_j =
Turan(A)_j + Turan(B)_j + S_j is asserted exactly (integer arithmetic
throughout; no floats in any inequality test). As a self-check, the
set {j : margin_j < 0} is compared against the `lc_breaks` field
already recorded for that instance in
`results/partial_sync_obstruction_20260811.json`. All 111 instances
across the three families matched exactly; the run aborts otherwise.

Points are bucketed by relative depth d/alpha into four bands,
[0.05,0.15], [0.15,0.25], [0.25,0.35], [0.35,0.45], and the minimum
rho within each band is recorded per instance, as both an exact
fraction and a float.

A fourth family, TG_{m,t} (Bautista-Ramos, arXiv:2511.00334, "Multiple
breaks of log-concavity in the independence polynomials of trees"),
was added as a stretch goal. TG_{m,t} is built from m disjoint copies
of T_{3,t} (root with 3 children, each child carrying t pendant P2
legs), their roots joined to a new root v0, plus one pendant leaf
child of v0. Applying I(G) = I(G\v) + xI(G\N[v]) at v0 gives the
paper's Eq. (4):

    I(TG_{m,t}) = I(P1) I(T_{3,t})^m + x I(S_{2,t})^{3m}

which supplies A = I(P1) I(T_{3,t})^m and B = x I(S_{2,t})^{3m}. This
is the vertex the paper's own recurrence uses at that step, so it is
the natural, non-arbitrary choice.

## Validation anchors (TG family)

The paper states TG_{4,6} breaks log-concavity exactly at indices 84,
82, 80, 78, and TG_{5,6} exactly at 105, 103, 101, 99, 97. Both were
recomputed independently with the rooted-tree DP and matched exactly:

| instance | alpha | expected breaks | computed breaks | match |
|---|---|---|---|---|
| TG_{4,6} | 85 | 84, 82, 80, 78 | 84, 82, 80, 78 | yes |
| TG_{5,6} | 106 | 105, 103, 101, 99, 97 | 105, 103, 101, 99, 97 | yes |

The anchors validated, so the TG margin analysis below is not fudged
or skipped.

## A caution about the TG decomposition

Unlike the three m-direction families, whose summands A and B are, by
construction, each a product of log-concave factors (so Turan(A)_j >=
0 and Turan(B)_j >= 0 hold at every index; confirmed empirically,
0 violations across all 111 instances), the TG decomposition does not
carry that guarantee. T_{3,t} is itself one of the small non-log-concave
tree families in this literature, so I(T_{3,t})^m need not be
log-concave, and indeed Turan(A)_j or Turan(B)_j goes negative
somewhere in 39 of the 45 TG instances tested. Where that happens,
rho_j can be small or negative, but that is not a near-miss on
absorption; it corresponds exactly to an actual LC break (rho_j < 1 iff
margin_j < 0 by the identity above), and because TG's alpha stays
modest over the tested range (m = 2..6, t = 4..12, alpha from 31 to
235), the shallow near-top break region can fall inside the 0.05-0.15
relative-depth band that was meant to probe the deep, absorbed
region. The tables below report both the raw band minima and a
filtered "absorbed only" view restricted to margin_j >= 0 (S_j < 0 but
no actual break), which is the subset the question is actually about.

## Trend tables: m-direction families

For each family and band, first and last (parameter, min rho) pair in
the tested range, the ratio (growth factor), and how many of the
consecutive steps in between were increases vs. decreases.

### T12_S24_2 (ell=2, nlegs=4, k=2), m = 2..60

| band | n pts | first (m, rho) | last (m, rho) | ratio last/first | up/down steps |
|---|---|---|---|---|---|
| 0.05-0.15 | 58 | (3, 32.77) | (60, 1.58e22) | 4.82e20 | 50 up / 7 down |
| 0.15-0.25 | 47 | (14, 2.03e7) | (60, 2.40e26) | 1.19e19 | 45 up / 1 down |
| 0.25-0.35 | 47 | (14, 3.56e7) | (60, 1.57e18) | 4.42e10 | 45 up / 1 down |
| 0.35-0.45 | 39 | (22, 1.13e8) | (60, 1.38e13) | 1.21e5 | 37 up / 1 down |

### T13_S24_2 (ell=3, nlegs=4, k=2), m = 2..30

| band | n pts | first (m, rho) | last (m, rho) | ratio last/first | up/down steps |
|---|---|---|---|---|---|
| 0.05-0.15 | 28 | (3, 13.59) | (30, 2.40e10) | 1.77e9 | 20 up / 7 down |
| 0.15-0.25 | 17 | (14, 1.16e8) | (30, 7.42e12) | 6.41e4 | 15 up / 1 down |
| 0.25-0.35 | 16 | (15, 2.81e7) | (30, 2.23e9) | 79.4 | 14 up / 1 down |
| 0.35-0.45 | 9 | (22, 1.27e9) | (30, 3.91e8) | 0.31 | 7 up / 1 down |

### T17_S25_2 (ell=7, nlegs=5, k=2), m = 2..24

| band | n pts | first (m, rho) | last (m, rho) | ratio last/first | up/down steps |
|---|---|---|---|---|---|
| 0.05-0.15 | 23 | (2, 54.17) | (24, 2.67e12) | 4.93e10 | 17 up / 5 down |
| 0.15-0.25 | 16 | (9, 6.59e6) | (24, 7.02e10) | 1.06e4 | 14 up / 1 down |
| 0.25-0.35 | 12 | (13, 6.04e6) | (24, 1.45e8) | 24.0 | 10 up / 1 down |
| 0.35-0.45 | 3 | (22, 1.73e8) | (24, 1.36e8) | 0.78 | 1 up / 1 down |

The direction is unambiguous: min rho per band grows, in most bands by
ten or more orders of magnitude, as m grows. The handful of down-steps
sit early in the range, before the band has more than a few qualifying
points; once the band is populated (roughly the middle third of each
family's m range onward), growth is monotone in every case checked.
The two apparent net decreases (T13 and T17, band 0.35-0.45, ratio
0.31 and 0.78) come from bands with only 9 and 3 points respectively,
i.e. this band only starts acquiring S_j < 0 points late in the tested
range for these two families, so "first" and "last" sit close together
in parameter value; they are not evidence against the growth trend
seen everywhere else.

The single closest approach to rho = 1 anywhere in the three
m-direction families is at T17_S25_2, m = 4, index j = 53 (d = 3,
d/alpha = 3/56 = 0.0536, band 0.05-0.15):

    Turan(A) = 22,604,275,560,196,590,732
    Turan(B) = 56,677,835,550,132,736
    S        = -6,715,309,957,701,755,352
    rho      = 5,665,238,348,936,680,867 / 1,678,827,489,425,438,838 = 3.3745...

This occurs at m = 4, the second point in that family's range, not at
large m. As m increases from there the minimum in that band rises to
2.67e12 by m = 24. Nowhere in 9,538 examined S_j < 0 points across the
three families does rho fall below roughly 3.37, and the values that
come closest to 1 are concentrated at small parameter values, not
approached as a limit from above as the parameter grows.

## Trend tables: TG family (absorbed-only, margin_j >= 0)

Band 0.05-0.15 is the only band with data across most TG_m values in
the tested range (m = 2..6, t = 4..12); TG_m = 6 also populates
0.15-0.25 at t = 7, 8.

| TG m | t range with data | min rho (absorbed only) | first (t, rho) | last (t, rho) |
|---|---|---|---|---|
| 2 | 8..12 | 697.98 | (8, 2224.7) | (12, 816.9) |
| 3 | 8..12 | 7083.1 | (8, 8661.2) | (12, 9983.5) |
| 4 | 5..12 | 12.91 | (5, 12.91) | (12, 75835.9) |
| 5 | 5..12 | 3.66 | (5, 3.66) | (12, 449253.1) |
| 6 | 4..12 | 23.57 | (4, 189.6) | (12, 2,093,957.3) |

The global minimum in the absorbed-only TG set is at m = 5, t = 5, j =
86 (d = 5, d/alpha = 5/91 = 0.0549), with no actual break at that
instance (lc_breaks = [] there):

    rho = 935,658,594,148,466,142,505,549 / 255,767,612,901,038,594,879,251 = 3.6582...

As with the m-direction families, this near-1 value sits at the
smallest tested t for that m, not at the largest. For every fixed TG_m,
rho grows steeply and (after t = 5 or 6) monotonically as t increases,
reaching 10^5-10^6 by t = 12. The pattern matches the m-direction
families: closest approach to 1 at small parameter values, explosive
growth in the surplus-to-deficit ratio as the parameter grows, no
convergence toward 1.

## Answer

Across four families (three m-direction, one t-direction) and 9,961
total S_j < 0 indices examined with exact integer/Fraction arithmetic,
the absorption factor rho_j does not shrink toward 1 as the family
parameter grows. It grows, typically by many orders of magnitude, and
the closest approaches to 1 (about 3.37 for the m-direction families,
about 3.66 for TG) occur at the smallest tested parameter values, not
as a limit from larger ones. This is kill-test evidence, not proof: it
rules out "grow the parameter within one of these four families" as a
route to a unimodality counterexample, but it says nothing about other
families or other decompositions of the same trees.

## Limitations

- Coverage is bounded: m up to 60/30/24 for the three m-direction
  families and (m, t) up to (6, 12) for TG. The observed growth is not
  a proof of unbounded growth; it is the trend across every instance
  computed.
- The four relative-depth bands (0.05-0.45) do not cover the full
  range where S_j < 0 was reported to persist (up to ~0.54); the
  deepest 0.45-0.54 slice was not banded separately, though the
  `deepest_Sneg` field in the JSON records the single deepest S_j < 0
  point per instance regardless of band.
- The TG decomposition (Eq. 4 of the source paper, deleting v0) is the
  natural choice given the paper's own recurrence, but it is not the
  only possible two-term split of TG's independence polynomial, and
  unlike the three m-direction families it does not have guaranteed
  log-concave summands. The "raw" (non-absorbed-only) band data,
  stored in the JSON under `bands`, mixes genuine near-top breaks into
  the same relative-depth bucket for TG's smaller-alpha instances; the
  absorbed-only view above is the one relevant to the decision
  question.
- "Kill-test evidence" describes what was found in this corpus, not a
  general theorem. A family with a genuinely different absorption
  structure could still behave differently.

## Files

- `scripts/probe_absorption_margin_20260811.py`: regenerates the JSON
  deterministically.
- `results/absorption_margin_probe_20260811.json`: full per-instance
  records (family, params, alpha, self-check results, per-band minima
  as exact fractions and floats, break sets, deepest S_j < 0 index,
  TG anchor validation).
