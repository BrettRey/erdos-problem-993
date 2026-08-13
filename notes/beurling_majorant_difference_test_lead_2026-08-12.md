# Beurling-tradition majorants for the adjacent-difference test: a lead
<!-- SUMMARY: Candidate lever for the C' bridge lane, killed same day ~-- band-limited one-sided majorants cannot resolve a width-1 valley on any tree (Selberg slack 2pi/w vs the torus cap w <= pi); kill-test in section 6 · status: KILLED (pure lane), hybrid residual not pursued · updated: 2026-08-12 -->

Status: speculative research lead, not a result. The mechanism below is a
construction by Claude (2026-08-12 session), prompted by the toolbox of a
newly filed paper, not by anything that paper proves. Refutation is a
success mode: one honest estimate showing the trade has no slack kills
this note and that is a fine outcome.

## 1. Trigger and provenance

Alex Cohen, "Fractal uncertainty in higher dimensions," Annals of
Mathematics 202 (2025), 265--307, is now filed at
`literature/cohen_2025_fractal_uncertainty_higher_dimensions.{pdf,md}`.
Its theorem (a fractal uncertainty principle for line-porous sets in
R^d) does NOT apply to #993: a full-text sweep found zero contact with
independence polynomials, unimodality, trees, or any combinatorics
(search terms recorded in section 5). The lead comes from its engine:
Cohen's main ingredient is a higher-dimensional Beurling--Malliavin
multiplier theorem, i.e. the tradition of constructing band-limited
functions with prescribed one-sided decay on a prescribed bad set. His
reference [25] (Mashregi--Nazarov--Khavin, "The Beurling-Malliavin
multiplier theorem: the seventh proof," Algebra i Analiz 17 (2005)) is
the classical 1-d form.

Two structural echoes, recorded for orientation but not actionable:

- **Conjecture B' is uncertainty-shaped.** The cusp envelope
  (`notes/real_collar_conjecture_2026-07-16.md`, section 2) says a
  non-real zero cannot be simultaneously near-dominant in modulus and
  non-tangential in angle: phi <= F(|z|/beta), F(1+) = 0. Two
  localization quantities that cannot both be extreme at once.
- **The C' enemy is line-concentration.** Cohen's porosity-on-lines
  notion exists because FUP fails when the set contains near-line
  structure; the sector profile C' fails (Jerrum--Patel) exactly when
  zeros accumulate on a line, the positive real axis, in C = R^2. Same
  degenerate geometry, no transferable lemma found.

## 2. The gap this aims at

The bridge lane's recorded obstruction
(`gpt_attack/bridge_window_unimodality/outcome_2026-07-16/bridge_lemma_report.md`,
line 616):

> The entropy argument proves the effective saddle bound in Lemma 8,
> but a standard relative-error local central limit theorem is too
> coarse to control the discrete second difference of adjacent saddle
> masses.

Same verdict in `notes/literature/literature_intake_2026-07-16.md`,
line 152: "central adjacent-mass cancellations occur at the same order
as the available approximation error, and local limits do not identify
the exact first strict descent."

So the known analytic route estimates each i_k and subtracts, and the
difference drowns in the per-coefficient error. The lead: change the
test, not the estimate.

## 3. The candidate mechanism

Write the first difference as a contour functional. With z = x e^{i theta}
on the saddle-radius circle,

    S_k := i_{k+1} - i_k
         = (x^{-k-1} / 2 pi) \int I(x e^{i theta}) e^{-i(k+1) theta}
           (1 - x e^{i theta}) d theta.

The mid-band goal (no valley in the window k < 2 alpha / 3) needs sign
control of S_k, or of staggered interval sums of the i_j, not accurate
values of individual i_k.

The Beurling--Selberg move: sharp interval sums against a sequence are
exactly what one-sided band-limited approximations were built for
(the large-sieve / Erdos--Turan application pattern). Replace the sharp
window indicator by a band-limited majorant M and minorant m with
m <= 1_{[a,b]} <= M and spectrum confined to width w. Then interval
sums of (i_j) are pinched between two contour integrals whose test
functions have Fourier support of width w, and the contribution of a
zero pair at angle phi is only uncontrolled when its oscillation
frequency lands inside that width. The Beurling--Malliavin theory is
precisely the quantified trade between spectral width w and achievable
one-sided pointwise control on the bad set (here: the near-axis zero
angles active on the saddle contour).

Why this could beat the local-limit route: the loss there is
per-coefficient additive error vs. an O(sigma^{-2}) relative
difference. Here the sharp object is the test (an indicator or sign
pattern), and the loss is an explicit spectral-width budget that can be
matched against the B' cusp envelope: near-dominant zeros have tiny
phi (slow oscillation, wavelength exceeding the degree), far zeros are
modulus-damped, so the dangerous band is intermediate, and a majorant
with width tuned to straddle that band is a designable object rather
than a fixed CLT error term.

Honest caveats, stated up front:

- Interval-sum control does not by itself forbid a valley; raw-sequence
  unimodality needs the pinch combined with a local step (e.g. bounded
  local variation, or log-concavity of the main term across the
  window). The assembly is undesigned.
- The trade may simply have no slack in the Jerrum--Patel regime, where
  C' is allowed to vanish. If the majorant width required to exclude
  the accumulating angles forces windows wider than the valley scale,
  the lead is dead. That is the first thing to test.

## 4. Next steps (either kills or promotes)

1. **Targeted literature probe** (outward, after this repo's inward
   check in section 5): queries of the shape "Beurling--Selberg
   majorant polynomial coefficients," "one-sided approximation
   unimodality coefficients," "Erdos--Turan discrepancy independence
   polynomial zeros," "band-limited majorant sign of Fourier
   coefficient." Also check whether the local-limit-from-zero-free-
   region line (possibly Michelen--Sahasrabudhe; attribution from
   memory, UNVERIFIED) has a difference-level refinement anywhere.
2. **Toy kill-test, certified arithmetic:** on the adversarial
   champions in `results/` and a phase-truncated Jerrum--Patel-style
   family, compute the actual angular spectrum of zeros on the saddle
   contour for a mid-band k, then check numerically whether ANY
   majorant width w separates the pinch from the valley scale. Exact
   coefficients; Arb intervals for any root data (float64 root-finding
   is forbidden, DECISIONS.md 2026-07-16). If no width works on the
   known hard family, record the kill and close the lane.

## 5. Search evidence

Claim: the extremal-function / Beurling-tradition toolkit is absent
from this repo's swept literature. Search run 2026-08-12:

    grep -rniE 'michelen|sahasrabudhe|local (central )?limit|local CLT|
    characteristic function|beurling|malliavin|selberg|
    extremal (function|minorant|majorant)'
    --include='*.md' --include='*.tex' --include='*.bib'  <repo root>

Zero hits for beurling, malliavin, selberg, extremal
function/minorant/majorant, michelen, sahasrabudhe. The local-limit
hits that do appear carry the "too coarse" verdicts quoted in
section 2. The 2026-07-16 intake's variant-search list
(`notes/literature/literature_intake_2026-07-16.md`, line 184) covers
local limits, ultra-log-concavity, adjacent ratios, and Turán deficits;
nothing from the extremal-function corner.

Claim: Cohen 2025 itself has no combinatorial contact. Two sweeps over
the full converted text (conversion verified complete through the
references and footnotes), nineteen terms: unimodal, independen, tree,
graph, log-concav, polynomial, erd, combinat, discrete, Newton,
interlac, real root, matching, Heilmann, coefficien, sequence, entrop,
Schwenk, Alavi. Zero hits. Its bibliography is entirely harmonic
analysis / quantum chaos.

## 6. Outcome (2026-08-12, same day): pure lane KILLED

The toy kill-test (section 4, step 2) ran the same day:
`scratch_beurling_band_killtest_20260812.py`, results in
`results/beurling_band_killtest_20260812.json`. Exact integer
coefficients (DP and the compressed word recurrence), exact rational
saddle bisection (integer sign tests only), Arb-certified root
isolation (max ball radius 1.9e-26 across all 16 cases; the P_101
control comes back real-rooted, 0 certified-non-real, as theory
requires).

**The mechanism dies structurally, on every tree, before zero geometry
even enters.** The accounting: a Selberg majorant/minorant pair for an
interval at bandwidth w carries mass slack ~ 2 pi / w in
coefficient-count units (optimal up to constants), so isolating a
single coefficient needs w > 2 pi. But a test kernel used against the
saddle-circle characteristic function lives on the frequency torus
[-pi, pi]: usable width is capped at pi, and any w > pi periodizes
over the whole torus, including the zero angles the method exists to
avoid. So even a completely zero-free spectrum resolves only
2-coefficient windows, and a width-1 valley is invisible. The
real-rooted control realizes exactly this bound (Delta_k = 2.00). The
uncertainty principle that inspired the lead is what kills it.

**Empirical margins, for the record.** On the July pareto champions
and the nine JP phase-truncated trees (n <= 189), the widest zero-free
angular band at modulus-activity eta = 0.5 gives interval resolution
Delta_k = 2 pi / w_avail between 2.19 and 4.74, worst at binary-h6,
k = 2 alpha / 3 (w = 1.326). So a hypothetical hybrid (majorant
interval pinches plus a local step covering <= 5 coefficients) is not
band-blocked at these finite sizes, but the resolution degrades
exactly toward the Jerrum--Patel accumulation regime, where w -> 0:
the hybrid inherits the same wall the C' profile already faces, on
top of its local step being undesigned. Not pursued.

**Scope of the kill.** What died is the mechanism this note proposed:
one-sided band-limited majorants tested against the saddle-circle
characteristic function. Untested and still conceivable: extremal
functions in other roles (multiplier damping applied to the polynomial
itself rather than the test; other dual pairs). No claim either way
about those; nobody should reopen this lane without an accounting that
beats the 2 pi / w slack or escapes the torus cap.
