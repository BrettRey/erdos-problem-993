# Fixed-`r` Proof Note Pressure Test

## Status

This is a new pressure-test document for:

```text
notes/fixed_r_proof_note_draft_2026-05-21.md
```

No existing proof notes or manuscript files were edited.

## Verdict

The draft is substantially better scoped than the earlier scratch: Theorem A
and Theorem B are separated, the `r=2,3` boundary issue is flagged, and the
corrected `a * perturbation` monotonicity condition is included.

It is not manuscript-ready yet.  The main remaining risks are:

```text
1. C1 in the certificate lemma still hides the global F-margin proof.
2. Theorem B says "computable threshold" but still assumes produced
   thresholds rather than proving the whole production pipeline in one place.
3. The small cases r=2,3 are left outside the arbitrary fixed-r criterion.
4. The sampled-lane proposition needs explicit certificate-file references
   if promoted to paper text.
```

## Quick Verification Run

Algebra/checks rerun during this pressure test:

```text
1. Path-moment identity for V+3K_3 checked through r=100.
2. First-order boundary cases checked through r=500: only r=2,3.
3. Corrected sampled certificate suite result read from JSON: all_ok=true.
```

These checks do not prove the lemmas by themselves, but they reduce the risk
of a transcription/sign error in the draft.

## Findings

### P1. Theorem A is safe, but needs artifact references

Draft location: lines 68-112.

The sampled-lane proposition is the strongest result in the document.  It is
properly framed as certificate-backed, not analytic.  The statement is safe if
the final version cites the exact artifacts:

```text
results/fixed_r_finite_route2_selected_a199.json
results/fixed_r_sampled_certificate_suite_corrected_monotonicity.json
```

and names the scripts used for the finite and asymptotic certificates.

Recommended change for a final proof note: add a short "Certificate artifacts"
paragraph immediately after Theorem A.

### P2. Certificate lemma has the right logic

Draft location: lines 114-195.

The C0-C3 implication is sound:

```text
C2 gives m - 4/3 + eta/a;
C3 costs at most eta/a;
therefore m - 4/3 > m - 3/2.
```

No issue with the algebraic implication.

The only caveat is that C1 is strong: it requires global mode verification,
not just adjacent mode verification.

### P3. Global `F`-margin bridge should become its own lemma

Draft location: lines 144-152.

The draft says the global margin follows from adjacent margins plus
real-rooted/log-concavity.  This is plausible and already part of the project
scaffold, but a polished proof note should state it explicitly:

```text
Lemma.  Let (c_k) be a positive log-concave sequence with mode m.  If
c_m-c_{m-1} >= eps c_m and c_m-c_{m+1} >= eps c_m, then
c_m-c_k >= eps c_m for every k != m.
```

Reason: log-concavity gives monotonicity of consecutive ratios on each side
of the mode.  Hence coefficients moving away from `m` cannot climb back
toward `c_m`.

This lemma would close the main hidden proof bridge in C1.

### P4. First-order shift rule is appropriately scoped

Draft location: lines 197-269.

This section now avoids overclaiming: it states the clean rule for `r>=4` and
marks `r=2,3` as boundary cases.

Pressure-test result: computation confirms that through `r=500`, only `r=2,3`
have non-unique first-order candidates.

Remaining issue: if Theorem B is stated for `r>=4`, this is fine.  If anyone
wants `r>=2`, the proof note must add separate `r=2,3` certificates or exact
formulas.

### P5. Hub-off reserve constant looks correctly stated

Draft location: lines 271-313.

The formula

```text
C_{r,q}=(9(V+3K_3)+24alpha+16)/12
```

is consistent with the earlier exact series extraction.  The pressure test
did not find a sign issue.

One wording improvement: "By the first-order shift interval, alpha >= -2/3"
should say:

```text
By the first-order shift interval, alpha lies in [-2/3,1/3].
```

Then the strict statement for `r>=4` follows because the lower endpoint is not
attained.

### P6. Path moment positivity is the strongest analytic piece

Draft location: lines 315-398.

The identity for `250F_{r+2}^3(V+3K_3)` was checked through `r=100`.  The
sign proof also looks structurally correct:

```text
odd r: bracket added;
even r>=6: use L_r<3F_r and F_r>=r;
small cases direct.
```

This section is close to manuscript-ready, except that the recurrence-solving
step is currently asserted.  A final version should say "Solving by induction"
or include a one-line verification that the displayed formulas satisfy the
recurrences and initial values.

### P7. Shifted-positive termination lemma is good but needs one caveat

Draft location: lines 400-419.

The derivative argument is sound for a polynomial with positive leading
coefficient.  The rational-function paragraph should explicitly say:

```text
discard common factors first, and ensure denominator has no zeros past the
chosen threshold by shifted positivity of the denominator.
```

Shifted positivity of the denominator already implies positivity and no zeros,
but stating it prevents a reader from objecting that rational functions may
have poles.

### P8. Hub-on perturbation section has the right correction

Draft location: lines 421-472.

The corrected monotonicity condition is the right one:

```text
a * perturbation decreases.
```

This should be retained prominently.  It fixes the exact weakness that would
otherwise make the asymptotic certificate logically incomplete.

Potential manuscript issue: the constant `C_r` in the mixture term is not
defined in the draft.  Add:

```text
C_r = 2 P_{r-1}(2)/P_r(1/2)
```

or whatever exact constant is used in the certificate.

### P9. Theorem B is a criterion, not yet a theorem with constructed A(r)

Draft location: lines 474-523.

The title correctly says "criterion."  Keep that word.  Do not retitle it as
"Fixed-r theorem" unless the proof note includes the full construction of all
five listed objects.

The current wording is mostly safe:

```text
Suppose the following objects are produced...
```

This is the right level.  It is a certificate theorem, not a standalone
analytic theorem.

### P10. Boundary and limitations sections are necessary and should stay

Draft location: lines 525-566.

These sections prevent the two most dangerous overclaims:

```text
1. arbitrary fixed-r for all r>=2 without small-case handling;
2. full d_leaf<=1 from spider-lane evidence.
```

Do not remove them in any tightened version.

## Recommended Next Edits in a New Draft

If creating a fifth/new version, make only these changes:

```text
1. Add artifact references under Theorem A.
2. Add a standalone global-F-margin lemma.
3. Define the mixture constant C_r in the hub-on perturbation lemma.
4. Add the denominator/no-poles sentence to shifted-positive termination.
5. Reword alpha range as alpha in [-2/3,1/3].
6. Add "verified by induction" language to the raw moment formulas.
```

Do not add more computations before these textual repairs.

## Current Call

The proof-note draft is viable as a controlled internal theorem note after one
more cleanup pass.  It should not yet be copied into `paper/main_v2.tex`.

The next work item should be a new "cleaned proof note" document that applies
the six recommended edits above, still leaving existing documents untouched.
