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
5. Do not close issue #5 or claim the hub-bouquet reserve, signed reserve, or
   Erdos 993.

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

Status: internal theorem-level claim within the current route, conditional on
the local-mode mean proof above. It has not yet received the same external
line-by-line audit as the signed-boundary and fair-binomial notes, so do not
promote it into the manuscript or issue closure until that audit is done.

Claim:

```text
For one-sided low-probability PB laws with V >= 1,
Delta_eff >= 1/(4V).
```

Known caveats:

- This is one-sided only.
- It is below the Poisson effective-drop ceiling `1/3`.
- It does not settle the signed case.

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

## Algebraic Reductions

### Corrected Signed Conditional Reduction

Source:

```text
notes/literature/signed_conditional_reduction_2026-07-04.md
```

Status: finite algebra, externally audited after the missing X-side boundary
term was found; statement-level fixes have been applied.

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
- The large `Binomial(500,1/2)` side-bound-only breaker is still a float
  certificate; exact certification would be needed before citing it as a
  theorem-level disproof.

## Disproved Or Superseded Targets

1. **Omitted-boundary X-side reduction.**
   The identity `c_{D-1}/c_D = E_pi[h_X]` is false without `beta_X`.

2. **Side-bound-only signed lemma.**
   The target `V * max(X_side_bound,Y_side_bound) >= 1/4` is too strong as a
   proof target; half-heavy balanced cases make the conditional bound tiny
   while actual reserve is large. Current certification is computational for
   the large balanced side-bound row.

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

Exact fair-binomial fallback:
  min_{N=4..40} V * Delta_eff = 5/8 at N=4
  min_{N=4..40} V * reserve = 3/4 at N=4

Prompt 4 hardening:
  first_descent uses a relative-drop tolerance so float-resolved plateaus are
  not treated as strict descents
  terminal descents are counted separately from malformed rows
  non-finite conditional terms and identity errors are fatal diagnostics
  a half-heavy+dust row below the 0.8 fallback threshold is checked exactly
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

## Next Permitted Work

Permitted in audit-hardening mode:

- patch any follow-up external-audit findings;
- add text or tests that prevent recurrence of known errors;
- reconcile issue comments and notes when a later correction supersedes an
  earlier statement;
- exact-certify or explicitly quarantine large float side-bound breakers;
- continue code/formula alignment audits after any script edits.

Not permitted while audit-hardening mode is active:

- claim progress toward proving issue #5 beyond the current ledger;
- state a new half-heavy stability lemma as likely true without a breaker pass;
- close issue #5;
- import these notes into the manuscript as theorem-level claims.
