# Formalization request: the even-connector Laurent block

## Objective

Formalize in Lean 4 the remaining algebraic lemma in the two-branch-vertex clan proof.
This packet is independent of the still-active global clan-map partition audit. Prove
the uniform Laurent-polynomial theorem; do not replace it by bounded computation.

## Definitions

Work over Laurent polynomials with rational coefficients. Put

```text
X = z + z⁻¹
H n = z^n + z^(-n).
```

A Laurent polynomial has **centrally unimodal coefficients** when it is symmetric under
`m ↦ -m`, has nonnegative coefficients, and its coefficients weakly decrease as
`m ≥ 0` increases within each parity class. An equivalent precise definition is
acceptable, but it must imply

```text
coeff f m - coeff f (m+2) ≥ 0    for every integer m ≥ 0.
```

For natural numbers `r,s ≥ 1` and rationals `c,d ≥ 1`, define

```text
G(r,s,c,d)
  = c*d*X^(r+s+1)
      + d*X^s*H(r+1)
      + c*X^r*H(s+1)
      + H(r+s+1).
```

## Main theorem

Prove that `G(r,s,c,d)` has centrally unimodal coefficients for all
`r,s ≥ 1`, `c,d ≥ 1`. Consequently every difference

```text
coeff (G r s c d) m - coeff (G r s c d) (m+2)
```

is nonnegative for `m ≥ 0`.

## Suggested reduction to audit

Let

```text
A = X^(r+s+1)
B = X^s * H(r+1)
C = X^r * H(s+1)
D = H(r+s+1).
```

First verify the exact identity

```text
c*d*A + d*B + c*C + D
  = (c-1)*(d-1)*A
      + (c-1)*(A+C)
      + (d-1)*(A+B)
      + (A+B+C+D).
```

The first three summands are expected to be centrally unimodal. It therefore suffices
to prove the unit-constant case. Under the standard Laurent-to-ordinary coefficient
translation, that case has coefficient polynomial

```text
(1+y)^N
  + (1+y)^s * (1+y^(r+1))
  + (1+y)^r * (1+y^(s+1))
  + 1 + y^N,
N = r+s+1.
```

The source proof in `notes/c2_connector_clan_reduction_2026-08-20.md` uses a split at
`j<s`, `j=s`, and `s<j<N/2`, together with binomial-difference bounds and convolution
of symmetric unimodal sequences. Audit that proof; repair it if a boundary or parity
case is missing.

## Deliverable and soundness

- Standalone current-Mathlib Lean project.
- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, weakened theorem,
  or bounded replacement.
- A source-labelled main file and README traceability table.
- `lake build` plus `#print axioms` for the main theorem.
- If false, return the smallest exact rational counterexample and identify the failed
  proof step. If partial, state the strongest proved parameter range and the first exact
  missing theorem.
