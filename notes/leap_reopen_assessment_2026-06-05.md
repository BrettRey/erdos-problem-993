# LEAP Reopen Assessment

Date: 2026-06-05

## Recommendation

Do **not** reopen the whole Erdős #993 project because of LEAP. Reopen a narrow formalization lane only if we can aim it at the fixed-`r` certificate bridge.

The LEAP paper changes the expected value of formalization workflows, not the mathematical state of the conjecture. It gives evidence that general-purpose models can become effective Lean formalizers when the task is decomposed into verified subgoals with a persistent proof DAG. That matches the packet already assembled in `gpt_attack/axiom_fixed_r_certificate/`.

## What LEAP Changes

LEAP makes it more plausible that a bounded Lean target can be pushed further than Aristotle did in March. The important design points are:

- parent reductions are Lean-checked before child lemmas are accepted;
- new `sorry`s are allowed only for explicit child subgoals;
- repeated subgoals are memoized in a DAG rather than rediscovered;
- an LLM reviewer rejects decompositions that merely restate the parent theorem.

This is especially relevant because the current packet already has a Lean-facing file that checks:

```bash
lake env lean gpt_attack/axiom_fixed_r_certificate/problem.lean
```

## What LEAP Does Not Change

LEAP does not supply a new #993 proof idea, and the released Google DeepMind repository appears to publish proof artifacts rather than a runnable LEAP agent. So the paper is not a reason to resume broad search, manuscript rewriting, ECMS, STP2 closure, or generic proof discovery.

It also does not make the old STP2 closure target attractive again. That theorem remains too close to the central open obstruction, and previous Aristotle runs already showed that formal tools can burn time there without making decisive progress.

## Best Reopen Target

The best target is the finite-support Gibbs mean-shift lemma in `gpt_attack/axiom_fixed_r_certificate/problem.lean`: replace the packaged Lipschitz assumption in `mean_shift_bound_from_lipschitz_certificate` with a direct finite-support Gibbs distribution proof.

This target is worth reopening because it is:

- narrow enough for an agentic Lean workflow;
- useful even if it only produces a sharper statement or a counterexample;
- aligned with the existing fixed-`r` certificate plan;
- less likely to collapse into an open-ended attack on Erdős #993.

## Proposed Reopen Scope

Reopen for one bounded session with this objective:

> Formalize or refute the finite-support Gibbs mean-shift lemma needed by the fixed-`r` certificate bridge, using exact rational/Lean statements and preserving the distinction between certificate facts and heuristic diagnostics.

Success should be one of:

- a Lean proof of the mean-shift lemma;
- a corrected Lean statement with explicit extra hypotheses;
- a small counterexample showing the current statement is insufficient;
- a LEAP-style subgoal DAG that isolates one genuinely smaller unresolved lemma.

## Stop Conditions

Stop if the work drifts into:

- full Erdős #993 proof search;
- STP2 closure;
- mode-mean or ECMS;
- global tree log-concavity;
- manuscript revision;
- numerical scans not tied to exact certificate replay.

## Decision

My recommendation is: **yes, reopen, but only as a bounded formalization experiment on the fixed-`r` certificate bridge.** Keep the main project parked otherwise.
