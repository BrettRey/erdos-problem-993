# Formalization request: LAW-V when one hub is trivial

## Objective

Formalize the uniform theorem proved in `notes/law_v_one_hub_proof_2026-08-20.md`.
This is an independent fallback theorem in the two-hub program and does not depend on
the active global clan-partition job. No bounded replacement is acceptable.

## Definitions

Let `P n` be the independence polynomial of the path on `n` vertices, with coefficient

```text
[x^j] P n = choose (n + 1 - j) j,
```

and conventions `P 0 = P (-1) = 1` implemented with a natural-number formulation that
makes `P (a-1)` unambiguous. For a finite indexed family or list of positive arm
lengths `a`, put

```text
Q = ∏ᵢ P(aᵢ)
R = ∏ᵢ P(aᵢ-1)
C = Q + x*R
B = C + x*Q.
```

For coefficient sequences `q,r,c,b`, define

```text
V k = c k * b k - c (k+1) * b (k-1),
```

using zero outside natural support.

## Final theorem

For every finite family of positive arm lengths and every natural index `k`, prove

```text
0 ≤ V k.
```

Equivalently, prove `c_k b_k ≥ c_(k+1) b_(k-1)`. Include empty and one-arm boundary
cases explicitly rather than relying on division by possibly zero coefficients.

## Required proof layers

1. The path coefficient formula and adjacent likelihood-ratio inequality

   ```text
   p(n,j) * p(n-1,j-1) ≥ p(n,j-1) * p(n-1,j)
   ```

   with complete support-boundary handling.
2. Common-convolution preservation of likelihood-ratio order when the common factor is
   nonnegative, log-concave, and has no internal zeros. Prove the needed finite
   coefficient theorem directly or via a fully formal TP2/Cauchy–Binet argument.
3. Products of path polynomials are log-concave with no internal zeros, and replacing
   the factors one at a time gives `R ≤lr Q`.
4. Prove the spider polynomial `C = Q+xR` is log-concave for arbitrary arms. This is the
   known Li–Li–Yang–Zhang spider theorem used by the paper proof, but it must not enter
   as an axiom or an unproved hypothesis. Reuse a proved Mathlib result if one exactly
   matches; otherwise formalize enough of its proof for this coefficient statement.
5. Prove the shifted cross term

   ```text
   c_k*q_(k-1) - c_(k+1)*q_(k-2) ≥ 0
   ```

   and finish by decomposing `V k` into the Turán gap of `C` plus that cross term.

As an optional adversarial check, formalize the stated failure of the naive general
two-hub cross term for arms `(1,1)` and `(1,1,5)`; it must not affect the one-trivial-hub
theorem.

## Deliverable and grading

- Standalone current-Mathlib Lean project with source-labelled declarations and a
  README traceability table.
- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, hidden support
  assumptions, weakened theorem, or bounded computation.
- `lake build` and `#print axioms` for the final theorem and main closure lemmas.
- Grade `COMPLETE` only if the arbitrary-arm theorem is proved without assuming the
  spider log-concavity layer. If that layer is the only missing item, return a clearly
  labelled conditional core and the exact theorem signature still required. If false,
  return the smallest exact counterexample and failed step.
