# Fixed-`r` Theorem Assembly Review

## Status

This is a new review document.  It audits:

```text
notes/fixed_r_theorem_assembly_scratch_2026-05-21.md
```

No existing manuscript or prior note was edited.

## Verdict

The scratch proof is a good assembly scaffold, but it should not yet be
promoted to manuscript text.  It is strongest as two separate claims:

```text
1. A certified sampled-lane proposition for
   r in {4,8,12,16,20,24,32,40,60,80}.

2. An arbitrary fixed-r computable-threshold schema, with several proof
   bridges still needing formal statements.
```

The main risk is overclaiming the arbitrary fixed-`r` theorem before the
certificate procedure is stated abstractly enough and before the small
boundary cases `r=2,3` are handled cleanly.

## Findings

### F1. The target theorem statement is too strong as written

Location: scratch lines 39-45.

The text says:

```text
For every fixed r>=2, there is a computable threshold A(r) ...
```

This is plausible and now well-supported, but the document has not yet proved
all pieces at manuscript standard.  In a theorem note, this should be split:

```text
Theorem A: sampled fixed-r lanes, fully certified.
Theorem B: fixed-r computable-threshold criterion/schema.
```

The arbitrary fixed-`r` result can be stated only after the certificate
procedure is formalized as a terminating construction.

### F2. The `r=2,3` boundary cases need their own lemma

Location: scratch lines 135-143 and 179-180.

The divisibility argument proves that first-order boundary cases occur only at
`r=2,3`.  The scratch says direct exact shifts select the lower candidate:

```text
r=2: {0:3, 1:1, 2:2}
r=3: {0:3, 1:4, 2:2}
```

But that is not yet a proof of eventual mode localization for those two
lanes.  A manuscript-ready version needs one of:

```text
1. exact closed formulas for the adjacent ratios in r=2 and r=3;
2. a shifted-coefficient certificate proposition for r=2 and r=3;
3. a reduction of r=2 to the pure spider case plus a direct r=3 check.
```

Without this, the arbitrary fixed-`r` theorem should probably start at
`r>=4`, with `r=2,3` listed as separate finite/certificate cases.

### F3. Adjacent `F` margins need to be connected to global mode margins

Location: scratch lines 309-330.

The hub-on mode step uses:

```text
2 max_k G_k < F_m * eta_mode/a.
```

This only suffices if there is a global lower bound

```text
F_m - F_k >= F_m * eta_mode/a
```

for every `k != m`, not just adjacent `k=m-1,m+1`.

The intended bridge is:

```text
P_r(x)(1+2x)^a is real-rooted with nonnegative coefficients,
so its coefficient sequence is log-concave/unimodal.
```

Then adjacent margins at the mode imply global margins.  This bridge is
already known in the project notes, but the scratch document needs to state it
as an explicit lemma if it becomes a theorem note.

### F4. The shifted-positive threshold procedure needs a termination lemma

Location: scratch lines 269-307.

The proof says eventual positivity is computable and uses shifted-coefficient
positivity.  To be theorem-ready, add the elementary lemma:

```text
If P(t) is a nonzero polynomial with positive leading coefficient, then
there exists T such that every coefficient of P(u+T) is positive.
```

Proof: the coefficient of `u^j` in `P(u+T)` is `P^{(j)}(T)/j!`, and each
nonzero derivative has positive leading coefficient after enough shifts.

For a rational function, apply this to numerator and denominator after fixing
positive leading signs.  This gives a genuine terminating certificate
procedure, not just an implementation pattern.

### F5. The hub-on perturbation monotonicity correction is important and should be foregrounded

Location: scratch lines 377-387.

This is a real repair.  The old condition, “perturbation decreases,” is weaker
than what the reserve comparison needs.  The correct condition is:

```text
a * perturbation decreases in each residue class.
```

This should be stated in any theorem note as a separate lemma or as part of
the hub-on perturbation certificate proposition, because it is the most likely
place for a reader to catch an error.

### F6. Hub-on remains the practical bottleneck

Location: scratch lines 389-405.

The new data support a useful diagnostic:

```text
hub-off/mode/lambda thresholds are small through r=120;
effective hub-on thresholds grow roughly linearly in the tested range.
```

This is not just engineering trivia.  It tells us what a stronger analytic
proof should target: reduce the crude mixture bound and the variance-based
fugacity-shift bound.

### F7. The sampled-lane proposition is the cleanest publishable result today

Location: scratch lines 408-437.

This part is ready to be turned into a rigorous proposition, because it has:

```text
finite exact checks;
asymptotic certificates;
corrected monotonicity rerun;
all_ok=True.
```

It should be stated separately from the arbitrary fixed-`r` schema.

## Proposed Next Document Shape

If we continue in a new document, use this structure:

```text
1. Definitions and polynomial split.
2. The sampled-lane proposition.
3. Certificate lemma: finite + mode + hub-off reserve + hub-on perturbation.
4. Analytic lemmas:
   a. first-order shift rule for r>=4;
   b. path moment positivity V+3K_3>=0;
   c. shifted-positive threshold termination;
   d. corrected a*perturbation monotonicity.
5. Fixed-r computable-threshold schema.
6. Remaining cases and limitations.
```

## Immediate Next Work

Do not run more large probes yet.  The next useful work is to create a
polished proof-note draft, still outside the manuscript, with:

```text
Theorem A. Sampled fixed-r certificate proposition.
Theorem B. Fixed-r computable-threshold criterion.
Lemma package. Shift rule, path moment positivity, shifted positivity,
hub-on perturbation.
```

The draft should explicitly mark `r=2,3` as separate boundary cases until
they are proved by exact formulas or certificates.
