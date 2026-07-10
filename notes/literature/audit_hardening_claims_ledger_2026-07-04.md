# Audit-Hardening Claims Ledger
Date: 2026-07-04

## Purpose

This ledger freezes the current claim status for the issue #5 signed reserve
work while the project is in audit-hardening mode. It is a reliability surface,
not a new proof note. While it is active, external-audit findings should be
patched before any extension of the proof route.

## Current Mode

Audit-hardening mode means:

1. Patch external-audit findings before extending the proof route.
2. Keep theorem, algebra, computation, conjecture, and disproved targets
   separated.
3. Add regression tests for each concrete mathematical error caught by audit.
4. Treat optimizer and smoke results as falsification evidence only.
5. Do not close issue #5 or claim the hub-included perturbation, the full
   hub-bouquet theorem, Case B, or Erdős 993.

### 2026-07-10 audit update

The one-sided theorem has now received an independent line-by-line audit. Its
Newton/localization argument is valid after adding an explicit support-domain
lemma: after zero parameters are removed, the first descent satisfies
`D <= m-1`, so `a_{D-1}`, `a_D`, and `a_{D+1}` are positive.

The audit also found that strict-first-descent reserve is not continuous under
a vanishing reflected perturbation when the endpoint law has a plateau.
Accordingly, any perturbative signed argument must use the first weak descent
as its one-sided endpoint and handle descent-index selection explicitly.

The bounded Skellam, one-reflected, and two-reflected extensions all survived
independent line-by-line audits. A subsequent Hillion--Johnson
cubic-curvature argument then proved the universal finite Poisson-binomial
quarter bound, closing the general signed bridge. Issue #5 nevertheless
remains open because the hub-included term has not been controlled.

## Theorem-Level Claims

### Local-Mode Mean Bound

Source:

```text
notes/literature/local_mode_mean_bound_proof_2026-07-04.md
```

Status: internally proved and externally audited earlier in this thread.

Claim:

```text
0 <= w_i <= 1, e_2 > 0, e_3 >= e_2
  => sum_i w_i/(1+w_i) >= 5/2.
```

Known caveats:

- This supplies the local-mode sublemma used in the current one-sided route.
- It does not prove the signed reserve theorem or issue #5.

### One-Sided Effective-Drop Bound

Sources:

```text
notes/literature/one_sided_localization_reduction_2026-07-04.md
notes/literature/one_sided_effective_drop_reduction_2026-07-03.md
```

Status: theorem-level within the current route after an independent
line-by-line audit and the support-domain repair. Do not promote it into the
manuscript or use it for issue closure without the remaining signed bridge.

Claim:

```text
For one-sided low-probability PB laws with V >= 1,
Delta_eff >= 1/(4V).
```

Known caveats:

- This is one-sided only.
- It is below the Poisson effective-drop ceiling `1/3`.
- It does not settle the signed case.
- Zero parameters must first be deleted, and the support-domain argument must
  be stated before the three adjacent ratios are defined.
- The same `1/(4V)` conclusion holds at the first weak descent. That is the
  plateau-safe endpoint for later signed perturbation.

### Fair-Binomial Signed Fallback

Source:

```text
notes/literature/fair_binomial_signed_fallback_2026-07-04.md
```

Status: theorem-level exact model calculation, externally audited and corrected.

Claim:

```text
X ~ Binomial(m, 1/2), Y ~ Binomial(n, 1/2), V=(m+n)/4 >= 1
  => V * Delta_eff >= 5/8
     V * reserve >= 3/4.
```

Sharpness:

```text
Both minima occur at total fair count m+n=4.
At m+n=5, V * Delta_eff = 3/4, not 1/2.
```

Known caveats:

- Exact fair-binomial only.
- Does not prove thresholded half-heavy stability.

### Skellam Effective-Drop Bound

Source:

```text
notes/literature/skellam_effective_drop_theorem_2026-07-10.md
```

Status: theorem-level as a limiting-family result after an independent
line-by-line audit.

Claim:

```text
Z = Pois(lambda)-Pois(eta), V=lambda+eta>=1, c_D>0
  => V * Delta_eff(D) > 1/4.
```

The proof combines a closed Darroch mode window with Baricz's strict Bessel
Turán lower bound and handles both zero-rate orientations separately. The
80-digit falsification harness refined 305 observed interior mode transitions
across 2,691 rows, with zero failures and observed minimum `1/3`.

Known caveats:

- This is the Poisson limit, not a finite signed PB comparison theorem.
- The only terminal first descent is the reflected `Pois(1)` endpoint; it is
  excluded from the displayed quotient.
- The high-precision grid is audit evidence, not part of the proof.

### One-Reflected-Bernoulli Effective-Drop Bound

Source:

```text
notes/literature/one_reflected_bernoulli_effective_drop_2026-07-10.md
```

Status: theorem-level within the signed-reserve route after an independent
line-by-line audit.

Claim:

```text
X low-probability PB, Var X>=1, Y=Bernoulli(q), 0<=q<=1/2
  => (Var X + q(1-q)) * Delta_eff(X-Y) >= 1/4.
```

The proof explicitly selects between the two possible strict-descents
`D in {S-1,S}`. It therefore remains valid across endpoint plateaus, where
strict-descent continuity fails. An exact-rational harness checked 47,850 rows
from 4,785 eligible profiles, with both branches and 47 threshold equalities
represented and no failures.

Known caveats:

- The proof relies on the already-audited one-sided localization
  `S+1<=4 Var X`.
- Its two-term coefficient transform is special to one reflected factor and
  cannot be naively iterated; the two-factor theorem below needs a distinct
  multi-term curvature argument.
- The exact finite grid is audit evidence, not part of the proof.

### Two-Reflected-Bernoulli Effective-Drop Bound

Source:

```text
notes/literature/two_reflected_bernoulli_effective_drop_2026-07-10.md
```

Status: theorem-level within the signed-reserve route after an independent
symbolic line audit.

Claim:

```text
X low-probability PB, Var X>=1,
Y=Bernoulli(q1)+Bernoulli(q2), 0<=q1,q2<=1/2
  => (Var X + Var Y) * Delta_eff(X-Y) >= 1/4.
```

The signed descent lies in `{S-2,S-1,S}`. The leftmost branch is an ordinary
shifted-Newton bound. The other branches follow from a six-term curvature
regrouping, monotonic reduction, and an exact nonnegative Bernstein-basis
certificate on `[0,1]^2`. A redundant audit reconstructed every derivative
identity and all three `5 x 5` Bernstein matrices exactly. The exact harness
then checked 97,488 rows from 2,708 eligible profiles, 584,928 component
bounds, and all 75 Bernstein coefficients with zero failures.

Known caveats:

- The proof again relies on `S+1<=4 Var X`.
- Naive iteration of the one-factor theorem is invalid because the required
  one-sided localization is not preserved after the first transform.
- This proves two reflected factors, not the general finite signed theorem.

### Universal Finite Poisson-Binomial Effective-Drop Bound

Source:

```text
notes/literature/universal_poisson_binomial_effective_drop_2026-07-10.md
```

Status: theorem-level after an independent proof audit and exact replay of all
certificates. The finite cells have both Sturm and independent positive-
Bernstein proofs.

Claim:

```text
For every finite Poisson-binomial law with V>=1 and a supported first strict
descent D,
  V * (1 - f[D-1]f[D+1]/f[D]^2) >= 1/4,
and therefore V * (1-f[D+1]/f[D]) >= 1/4.
```

The proof propagates Hillion--Johnson cubic curvature to force explicit
two-sided mass windows about the mode, combines a pairwise variance bound with
a sharp max-atom variance bound, and reduces the theorem to one scalar
inequality. The range `3<H<=16` is closed by exact Sturm certificates and,
independently, by full-degree Bernstein expansions with 275/275 positive
coefficients. The tail is closed by an analytic Bonferroni--Bernstein
argument. The main replay checks twelve compact cells and the complete
symbolic tail; the second replay checks all thirteen finite Bernstein
identities. A separate exact 11,320-vector PBD scan found no failure but is not
used by the proof.

Consequences:

- Shifting `X-Y` to an ordinary PB law proves the complete nonterminal finite
  signed effective-drop theorem; terminal raw reserve is immediate.
- The raw reserve follows because `1-R_+ >= 1-R_+/R_-` at a strict descent.
- For `A=(1+x)^sQ`, `s>=4`, and `V_A=s/4+V_Q`, every post-descent ratio is at
  most `1-1/(s+4V_Q)`.

Known caveats:

- This closes the product term, not perturbation by the hub-included `xR`.
- It does not prove a single-hub growing-arm theorem, the general Case-B mode
  bound, PNP, or Erdős 993.
- The effective quotient is undefined when the largest mode is the upper
  support endpoint; only the raw terminal reserve is asserted there.

## Algebraic Reductions

### Corrected Signed Conditional Reduction

Source:

```text
notes/literature/signed_conditional_reduction_2026-07-04.md
```

Status: finite algebra, externally audited after the missing X-side boundary
term was found; statement-level fixes have been applied. This route is now
superseded as a dependency by the universal theorem, but remains useful audit
history.

Key corrected X-side identity:

```text
c_{D-1}/c_D = E_pi[h_X] + beta_X,
beta_X = a_{n_X} b_{n_X+1-D}/c_D.
```

Key transfer statement:

```text
1 - R_+ >= Delta,
```

with strictness only when `R_+>0`.

Known caveats:

- Assumes `c_D>0` and contiguous side supports after deterministic shifts are
  stripped off.
- Terminal support descents are support-edge cases, not the interior
  obstruction targeted by the reduction.
- The algebra does not imply a signed reserve lower bound without additional
  control of dispersion and boundary terms, or a direct reserve fallback.

## Computational Evidence

### Corrected Side-Bound Audit

Source:

```text
notes/literature/corrected_side_bound_audit_2026-07-04.md
```

Status: computational falsification/evidence only.

Current conclusion:

```text
The side-bound-only target
  V * max(X_side_bound, Y_side_bound) >= 1/4
should be retired as a proof target.
```

The live shape is a two-clause target:

```text
Either the corrected conditional side-bound carries the descent,
or raw/effective reserve is already large by a direct fallback mechanism.
```

Known caveats:

- `0.75` is a search threshold, not a theorem constant.
- `0.8` is ruled out for the displayed half-heavy+dust fallback row by an
  exact rational regression test.
- The side-bound-only target now has a compact exact rational disproof:
  `X=Bernoulli(1/4)+Binomial(6,1/2)`, `Y=Binomial(4,1/2)` gives
  `V * max(X_side_bound,Y_side_bound) = 346967/1481760 < 1/4`.
- The large balanced `Binomial(500,1/2)` row has also been independently
  exact-certified; it is useful as an asymptotic explanation, but is no
  longer needed as the primary disproof.

## Disproved Or Superseded Targets

1. **Omitted-boundary X-side reduction.**
   The identity `c_{D-1}/c_D = E_pi[h_X]` is false without `beta_X`.

2. **Side-bound-only signed lemma.**
   The target `V * max(X_side_bound,Y_side_bound) >= 1/4` is too strong as a
   proof target; half-heavy balanced cases make the conditional bound tiny
   while actual reserve is large. The small half-heavy row above disproves it
   in exact rational arithmetic.

3. **Fair-binomial fallback constant `1/2`.**
   The statement `V * Delta_eff = 1/2` at `N=5` was an arithmetic error.
   Correct sharp constant is `5/8`, attained at `N=4`.

4. **Strict transfer `1-R_+ > Delta`.**
   Correct transfer is weak: `1-R_+ >= Delta`, strict iff `R_+>0`.

5. **Aggressive half-heavy fallback constants near `0.8`.**
   Half-heavy+dust rows show `V * Delta_eff` and `V * reserve` near `0.79`
   while the corrected side-bound is below `1/4`; the displayed rational row
   below `0.8` is now regression-covered exactly.

## Regression Tests Added

Current tests in `test_all.py` cover:

```text
X=Bernoulli(1/2), Y=Binomial(10,1/2):
  old omitted-boundary X expression = 4/11
  true Delta = 3/10
  beta_X = 42/55

Half-heavy+dust near misses:
  side-bound below 1/4
  actual V * Delta_eff and V * reserve above 0.75
  compact side-bound disproof certified exactly at 346967/1481760

Exact fair-binomial fallback:
  min_{N=4..40} V * Delta_eff = 5/8 at N=4
  min_{N=4..40} V * reserve = 3/4 at N=4

Prompt 4 hardening:
  first_descent uses a relative-drop tolerance so float-resolved plateaus are
  not treated as strict descents
  terminal descents are counted separately from malformed rows
  non-finite conditional terms and identity errors are fatal diagnostics
  a half-heavy+dust row below the 0.8 fallback threshold is checked exactly

Plateau safety:
  Binomial(5,1/2)-Bernoulli(1/1000) moves the strict descent from 4 to 3
  the perturbed drop is exact and converges to the one-sided weak-descent 1/2
```

## External Audit Status

Prompt 4, the code/formula and overclaim audit of

```text
notes/literature/corrected_side_bound_audit_2026-07-04.md
```

has been received and patched. It found exact code/formula alignment for the
ten audited quantities, plus one real numeric guard issue and five prose-scope
issues. The guard issue was the float-resolved plateau risk in descent
placement; the prose issues were missing `V>=1` scope in the fair-binomial
quote, overbroad topic sentences, missing finite-search caveats, missing
terminal-descent caveat, and too little separation between float certificates
and exact disproofs.

The 2026-07-10 follow-up audit of the one-sided route found no algebraic
defect. It supplied the missing support-domain lemma, a strict-descent
discontinuity example

```text
X=Binomial(5,1/2), Y=Bernoulli(q), q -> 0+,
```

and the first-weak-descent repair. It rules out any proof that treats the
strict first-descent reserve as a uniformly continuous function of the
reflected variance.

Independent follow-up audits of the Skellam and one-reflected-Bernoulli notes
found no substantive proof defect. Three off-by-one/presentation issues in the
one-Bernoulli note were corrected before this ledger update; they did not
change the theorem or its proof branches.

The two-reflected-Bernoulli proof then received a separate redundant audit.
It exactly reconstructed the six-term identity, four-step Newton factor, two
monotonicity derivatives, and all Bernstein matrices, with no proof defect.

The universal proof received two independent audits. One caught the need to
carry the endpoint convention `delta_0=delta_n=1`; without it, the abstract
recurrence admits truncated counterexamples. With the genuine PB endpoint
cubic inequalities included, both audits reproduced the mass windows, scalar
reduction, all compact Sturm root counts, and every analytic Bernstein
coefficient. The exact replay passed locally.

A subsequent independent exact conversion rewrote all thirteen finite-cell
numerators in the full Bernstein basis. All 275 coefficients are positive,
every polynomial identity reconstructs exactly, and the generated Lean 4.28
identity file compiles. This gives a second finite certificate and removes
Sturm theory from the preferred mechanical-formalization route.

## Next Permitted Work

Permitted in audit-hardening mode:

- patch any follow-up external-audit findings;
- add text or tests that prevent recurrence of known errors;
- reconcile issue comments and notes when a later correction supersedes an
  earlier statement;
- exact-certify or explicitly quarantine large float side-bound breakers;
- continue code/formula alignment audits after any script edits;
- use the certified product-term reserve to prove perturbation bounds for
  `I=A+xR`, first for fixed path arms and then for growing broom arms;
- certify the pre-descent and post-descent increment budgets separately;
- decide deliberately whether the universal PB theorem merits a standalone
  manuscript/revision; do not edit the submitted manuscript automatically.

Not permitted while audit-hardening mode is active:

- claim that issue #5, the hub-bouquet theorem, or Case B is closed;
- state a new half-heavy stability lemma as likely true without a breaker pass;
- close issue #5;
- import these notes into the manuscript as theorem-level claims.
