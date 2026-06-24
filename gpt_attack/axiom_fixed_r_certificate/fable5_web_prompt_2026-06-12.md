# Claude Fable 5 Web Prompt: Erdős 993 Fixed-`r` Lean Target

You are working on a narrow Lean 4 formalization target inside an Erdős Problem 993 project. The task is **not** to prove Erdős #993, not to revise the manuscript, and not to restart broad proof search.

## Files attached

Primary files:

- `gpt_attack/axiom_fixed_r_certificate/problem.lean`
- `gpt_attack/axiom_fixed_r_certificate/task.md`
- `gpt_attack/axiom_fixed_r_certificate/README.md`
- `gpt_attack/axiom_fixed_r_certificate/fixed_r_certificate_target.tex`
- `notes/leap_reopen_assessment_2026-06-05.md`

Lean environment files:

- `lean-toolchain`
- `lakefile.toml`
- `lake-manifest.json`

Optional context:

- `gpt_attack/conjecture_and_state.md`
- `notes/fixed_r_proof_note_clean_2026-05-21.md`

## Current state

The file `gpt_attack/axiom_fixed_r_certificate/problem.lean` currently checks with:

```bash
lake env lean gpt_attack/axiom_fixed_r_certificate/problem.lean
```

The intended next target is the finite-support Gibbs derivative identity used by the already-proved mean-shift bridge:

```lean
deriv mu t = finiteVariance N (p t) / t
```

At present this identity is assumed as `hderiv_id` in:

```lean
mean_shift_bound_from_finite_gibbs_distribution
```

The useful task is to replace that assumption with a Lean theorem for the concrete polynomial-weight Gibbs family, or to identify the exact additional hypotheses required.

## Target theorem shape

Prefer a theorem approximately of this form, adjusting the statement if Lean/mathlib needs more explicit hypotheses:

```lean
def gibbsZ (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) : ℝ :=
  ∑ k, w k * t ^ k.val

def gibbsProb (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) (k : Fin (N + 1)) : ℝ :=
  (w k * t ^ k.val) / gibbsZ N w t

def gibbsMean (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ) : ℝ :=
  finiteMean N (gibbsProb N w t)

theorem deriv_gibbsMean_eq_variance_div
    (N : Nat) (w : Fin (N + 1) -> ℝ) (t : ℝ)
    (ht : 0 < t)
    (hZ : gibbsZ N w t ≠ 0) :
    deriv (fun s => gibbsMean N w s) t =
      finiteVariance N (gibbsProb N w t) / t := by
  ...
```

If that statement is too weak or too strong, return the corrected theorem. In particular, make explicit whether Lean needs:

- `0 < gibbsZ N w t`;
- nonnegative weights `∀ k, 0 <= w k`;
- at least one positive weight;
- local nonvanishing of `gibbsZ`;
- an interval restriction rather than a pointwise theorem;
- a version stated with `HasDerivAt` before translating to `deriv`.

## Required discipline

- Use Lean 4 / mathlib compatible with the attached `lean-toolchain`.
- Keep all arithmetic exact.
- Do not introduce `axiom`, `admit`, `sorry`, or `by native_decide` placeholders.
- Do not change existing proved theorem statements unless necessary; if you must change one, explain why.
- Do not rely on floating-point computations.
- Do not assume global log-concavity, ECMS, the mode-mean conjecture, Conjecture A, or the full Erdős #993 conjecture.
- Prefer small helper lemmas that can be checked independently.

## Deliverable

Return one of the following:

1. A patch or replacement block for `problem.lean` that proves the Gibbs derivative identity and preserves the existing checks.
2. A smaller verified Lean lemma that is a genuine step toward the derivative identity, plus the exact remaining subgoal.
3. A counterexample or corrected theorem statement showing that the proposed hypotheses are insufficient.

Also include:

- the exact command I should run locally;
- any new imports required;
- a short explanation of the proof structure;
- a list of any extra hypotheses added and why they are mathematically needed.

Do not give a broad proof essay. I need Lean code or a precise Lean-facing obstruction.
