# Break-depth lane: how far below the head can a tree LC break sit?
<!-- SUMMARY: New search lane (objective: minimize k_break - LK-threshold) opened after the Rota disproof read; first campaign finds a rigid dist=4 wall across 58.5M exact evaluations plus complete small mutation balls, and reduces dist=3 to a finite window (depth-2 break, alpha in {15,16}, n in [29,32]) · status: wall standing, finite target named · updated: 2026-08-14 -->

**Motivation:** `notes/matroid_nonunimodality_2026-08-14.md`. The
Divoux-Larson-Lowen-Wang disproof of Rota's conjecture (arXiv:2608.07342)
runs on three resources: severe LC failure, replication, and an in-category
tilt. Trees have the first two; the tilt is absent; and all known tree LC
damage sits in the tail that Levit-Kadrawi (arXiv:2603.17114, Lemma 2.19 /
Cor 2.20) prove harmless. So the binding resource for any #993 counterexample
is **break depth**: an LC break at index k below the threshold
`thr = ceil((2*alpha - 1)/3)` would breach the prefix that the reduction
needs. No prior lane had depth as its explicit objective (rebound deficit,
root angle, near-miss ratio, and absorption margin have all been optimized;
break position never).

**Metric.** For a break at k: `dist = k - thr`. Penetration of the prefix is
`dist < 0`. Secondary coordinate: reflected depth `d = alpha - k`.

**Artifacts.** `scripts/search_break_depth_20260814.py` (subcommands
`selftest` / `sweep` / `balls` / `evolve`; exact integer/Fraction arithmetic
throughout; every evaluated tree is also checked for unimodality, with any
valley dumped to an ALARM file as a #993 counterexample candidate).
Results: `results/break_depth_sweep_20260814.json`,
`results/break_depth_balls_20260814.json`,
`results/break_depth_search_20260814.json`,
island checkpoints in `results/break_depth_lane_20260814/`.

## The incumbent frontier: dist = 4, everywhere

- **Exhaustive small-n data.** All 21 LC-failure trees known from exhaustive
  enumeration (2 at n=26, 0 at n=27, 19 at n=28) have their break at
  k = alpha - 1 with alpha in {14, 15}: dist = 4 in every case
  (recomputed from `results/analysis_n26.json`,
  `results/analysis_n28_modal_lc_nm.json`; polys reconstructed from graph6
  where not stored).
- **Structured families (sweep: 982 trees, 128 with breaks).** TG_{m,t} grid
  (m <= 20, t <= 8), the KL pattern grammar (T_{1,l} : P_2)_ks over uniform
  and all mixed ks-tuples (length 2-3, entries 2-6, l in {2,3,4}), spider-
  unit pattern families, double-broom controls: **minimum dist = 4**,
  achieved only at the small corner (n = 26-46). Growing any family
  parameter moves breaks away from the threshold: in TG_{m,t} the deepest
  break sits at reflected depth 2m-1 while the threshold sits at reflected
  depth ~ (t+1)m, so the band can never reach the prefix within the family
  (2m-1 < (t+1)m for every t >= 1).
- **Deterministic mutation balls (`balls`).** The complete 1-mutation ball
  around the frontier trees (KL grammar instances + known failures; 339
  distinct polynomials) and the complete 2-mutation ball around all 21 known
  LC-failure trees (10,284 distinct polynomials), under pendant-leaf,
  pendant-P2, subdivision, leaf-deletion, and smoothing moves: minimum
  dist = 4, and **not one depth >= 2 break at alpha <= 18 anywhere**.
- **Evolutionary campaign (`evolve`, 30 min, 6 islands): 58,476,064 exact
  evaluations** (`results/break_depth_search_20260814.json`). Four islands
  optimized (-dist, then exact Turan ratio at indices deeper than the
  incumbent break) at n-caps 48/96/200 with full and descending-zone
  variants; two islands maximized the Turan ratio in the band
  [thr-2, thr+2] directly. Seeds included every known LC-failure tree, the
  KL grammar frontier, TG instances, and random trees; mutation set as in
  the ball audits plus subtree regraft, arm/broom growth, and T_{1,l}-unit
  grafts. No island beat dist = 4. Zero unimodality alarms. Coverage
  caveat: all four deepen islands converged onto the same n = 35,
  alpha = 19 basin (the klpat m = 3 shape, proxy plateau 0.783), so most of
  those evaluations probed one neighborhood intensively rather than the
  space broadly; a restart/diversity-forcing variant is the obvious second
  campaign if one is wanted.

## Why dist = 3 is a finite question right now

`dist = (alpha - d) - ceil((2*alpha - 1)/3)` depends only on (alpha, d).
Exact arithmetic (verified for d <= 7): **dist = 3 forces
alpha in {3d+8, 3d+9, 3d+10}**, and dist <= 2 forces alpha < 3d+8. So:

- d = 1: alpha in {11, 12, 13}, hence n <= 26: exhausted, no failures exist
  there at all (first failures are alpha = 14 at n = 26). Impossible.
- d = 2: alpha in {14, 15, 16}. alpha = 14 forces n <= 28: exhausted (all
  n <= 28 breaks are depth-1). **Remaining window: alpha in {15, 16} with
  n in [29, 32]** (n <= 2*alpha and n >= 29 since n <= 28 is exhausted).
- d = 3: alpha in {17, 18, 19}, n in [29, 38]. And so on: the required
  depth grows one per three alpha, while the observed depth-per-alpha line
  of the KL grammar (depth m-1 at alpha ~ 5m+9) loses two alpha-units per
  depth unit against it, i.e. the deficit widens with depth.

So the first place a dist-3 tree can exist is a **depth-2 break at
alpha in {15, 16}, order 29-32** ~-- just past the exhaustive boundary, in
trees that are near-perfect-matching (alpha barely above n/2). Everything
this campaign threw at that window (mutation balls of the alpha = 15
failures, four evolutionary islands whose n-48 population lives exactly at
these orders) failed to produce one, and the whole n <= 28 universe is
exhausted against it.

## Reading

Three coherent readings of the wall, in decreasing strength:

1. **Provable law?** "Every tree LC break has dist >= 4" is now an
   adversarially-tested statement in the campaign's rigid-conjecture sense
   (it did not move under 58.5M targeted evaluations plus complete small
   balls). If true, combined with Levit-Kadrawi it would come close to the
   prefix side of #993 for the head mechanism specifically ~-- worth checking
   whether the c2/contraction machinery says anything about depth-1-only
   breaks at small alpha.
2. **Finite check.** The d = 2 window is a bounded computation: enumerate
   trees with n in [29, 32], filter alpha in {15, 16} (cheap greedy alpha
   before the DP), check LC at reflected depth 2. The n = 29 unimodality
   sweep (5.47e9 trees) ran on Modal in 2026-03 with 1024 workers; this is
   the same shape of job with a strong prefilter (near-perfect-matching
   trees are a thin slice). If it comes back empty, dist >= 4 holds through
   n = 32 unconditionally, and the d >= 3 windows inherit the widening
   deficit.
3. **Null result.** At minimum: the Rota-recipe transfer question from
   2608.07342 is answered quantitatively for the only resource trees were
   missing a measurement of. Break depth is rigid; the damage cannot be
   steered toward the prefix by any move this campaign knows.

**Methodological note.** The band-objective islands confirmed that the
Turan ratio at the threshold band approaches 1 as 1 - Theta(1/n) (0.9480 at
n = 96, 0.9732 at n = 200, 0.9902 at n = 1042 in TG_{20,8}): pure mode
smoothness, no anomalous tightness. Any future depth lane should treat
near-mode ratio maximization as an attractor to design around, not a
signal.

## Provenance

Session 2026-08-14, Brett's direction ("run the depth-maximization lane")
following the 2608.07342 read. Self-test validates against stored
exhaustive polys, the TG_{4,8} break set {102,104,106,108}, and
500-mutation tree-invariant closure. No manuscript or public surface
touched.
