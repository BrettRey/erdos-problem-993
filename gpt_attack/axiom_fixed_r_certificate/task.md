# Task: Fixed-`r` Certificate Bridge

You are working on the Erdos Problem 993 repository.  Your task is not to
prove the full tree independence-polynomial unimodality conjecture.  Your task
is to prove or refute the bridge lemmas in this packet.

## Target

Start with `problem.lean`.  The first-pass arithmetic bridge lemmas there
already check in Lean.  Extend them toward the full certificate criterion, or
return a minimal counterexample showing that a stated hypothesis is
insufficient.

Preferred order:

1. Replace the packaged Lipschitz assumption in
   `mean_shift_bound_from_lipschitz_certificate` with a finite-support Gibbs
   distribution proof.
2. State the full fixed-`r` certificate criterion in Lean.
3. Connect the exact certificate hypotheses from
   `fixed_r_certificate_target.tex` to the Route-2 conclusion.

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
