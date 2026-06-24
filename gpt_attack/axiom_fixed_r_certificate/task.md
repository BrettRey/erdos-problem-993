# Task: Fixed-`r` Certificate Bridge

You are working on the Erdos Problem 993 repository.  Your task is not to
prove the full tree independence-polynomial unimodality conjecture.  Your task
is to prove or refute the bridge lemmas in this packet.

## Target

Start with `problem.lean`.  The first-pass arithmetic bridge lemmas there
already check in Lean.  Extend them toward the full certificate criterion, or
return a minimal counterexample showing that a stated hypothesis is
insufficient.

Current checked status:

- The adjacent-to-global margin lemma checks in Lean.
- The fugacity perturbation lemma checks in Lean.
- The finite-support Gibbs derivative identity checks in Lean.
- The concrete independence-polynomial mean-shift theorem checks in Lean,
  including positivity of the partition function from nonnegative coefficients
  and a positive constant term.
- The abstract composition theorem `fixed_r_certificate_composition` checks in
  Lean, with record wrapper `Route2Certificate`.
- The split `B` versus `F^-` composition theorem
  `fixed_r_certificate_composition_split` checks in Lean, with record wrapper
  `Route2SplitCertificate`.  This also consumes the hub-on mode-preservation
  hypothesis via `sum_mode_of_adjacent_margins_and_perturbation`.
- The family-level threshold split `route2_family_from_finite_and_tail`
  checks in Lean: exact finite checks handle `1 <= a < A`, and a function
  `a ↦ Route2SplitCertificateFor` handles `a >= A`.
- The bundled family record `Route2FamilyCertificate` checks in Lean, with
  unpacking theorem `route2_of_family_certificate`.

Preferred next order:

1. Formalize the spider-polynomial model: path polynomial recurrence,
   coefficient extraction for `P_r(x)(1+2x)^a` and
   `xP_{r-1}(x)(1+x)^a`, and the identity
   `I(T_{a,r}) = F + G`.
2. Connect the concrete script output to `Route2FamilyCertificate` instances
   once the script-side data format is frozen.

After those are settled, propose a Lean statement for the full fixed-`r`
certificate criterion in `fixed_r_certificate_target.tex`.

## Required Discipline

- Use Lean/mathlib v4.28.0.
- Keep all arithmetic exact.
- Do not cite floating-point scans as proof.
- Do not assume global log-concavity of tree independence polynomials.
- Do not assume ECMS, mode-mean, or Conjecture A.
- If a lemma needs an extra positivity or support hypothesis, state it
  explicitly and explain why the fixed-`r` certificate scripts provide it.

## Mathematical Context

The fixed-`r` lane is

```text
T_{a,r}=S(2^a,r)
I(T_{a,r};x)=P_r(x)(1+2x)^a + xP_{r-1}(x)(1+x)^a.
```

The certificate scripts separately check:

- finite exact Route-2 inequalities below a threshold;
- adjacent hub-off mode margins;
- hub-off reserve at the hub-off tie fugacity;
- hub-on mode domination;
- hub-on perturbation domination.

The proof gap is the bridge from these separate checks to a theorem.  Do not
try to replace the certificate machinery with a broad informal argument.

## Deliverable

Return:

```text
1. Lean code for the new proved lemmas, or a counterexample.
2. A short list of any extra hypotheses added.
3. A note explaining which existing certificate script supplies each
   hypothesis.
```
