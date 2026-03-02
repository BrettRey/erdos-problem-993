# Round 19, Prompt 3: Step-2 Closure Focus (All Deficits Concentrate Here)

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Empirical trigger

Through exhaustive scans (`n<=18`, boundary-correct prefix regime), every `X_k<0` case occurs at update step `t=2`.

This suggests a strategy:

- prove no-deficit (or easy dominance) for `t>=3`, and
- concentrate the hard inequality on step-2 structure.

Empirical seed from coarse full-data search (`n<=18`, all negatives):

`sum_err <= D + C10 + C01 + C11`

appears to hold with strong margin, where

- `C10 = sum_i Lambda_old(i) P(k-i-1)Q(k-i)`
- `C01 = sum_i Lambda_old(i) P(k-i)Q(k-i-1)`
- `C11 = sum_i Lambda_old(i) P(k-i-1)Q(k-i-1)`.

Calibrated scalar form also holds empirically through `n<=18`:

`sum_err <= D + lambda* (C10 + C01 + C11)` with `lambda* ≈ 0.0520138`.

---

## Step-2 normal form

At `t=2`, with leaf block size `ell` and first non-leaf child `(I_1,E_1)`:

- `A = (1+x)^ell * I_1`
- `B = E_1`
- second child gives factor pair `(P,Q)=(I_2,E_2)`.

Write `Lambda_k^{new}=D_k+X_k` and the diagonal-Abel decomposition in this normal form.

---

## Task 1: Prove / reduce “no new negatives after t=2”

Try to show one of:

1. `X_k >= 0` for all `t>=3` under chain STP2 + LC, or
2. stronger: odd residual term vanishes (or is auto-dominated) for `t>=3`.

If full proof unavailable, isolate the exact first inequality where the argument breaks.

---

## Task 2: Build a closed inequality for step-2

Using the step-2 normal form, derive an explicit bound of type

`Err_odd_step2 <= Reserve_step2`

where reserve is expressed via IH channels (`Lambda_old`-weighted terms) and/or a single extra testable hypothesis.

Prefer formulas that depend on small structural parameters (`ell`, child subtree sizes) so they can be checked parametrically.

If possible, connect `Reserve_step2` explicitly to the shifted-channel seed above.

---

## Task 3: Parity mechanism

Explain why odd diagonals are the only unresolved part in the Abel reduction and derive the sharp parity split formula (even part + odd residual).

Then show how your step-2 bound interacts with this split to yield `Lambda_k^{new}>=0`.

---

## Task 4: Output

Provide:

1. explicit theorem candidate for step-2,
2. proof sketch with all dependencies,
3. one finite computational checklist that would decisively validate/falsify the missing lemma.
