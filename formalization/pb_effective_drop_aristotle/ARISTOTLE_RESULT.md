# Aristotle Result: Endpoint-Aware Reserve Core

Date: 2026-07-10

- Aristotle project: `8d59a353-5c7b-4071-b837-9ab7bf561be3`
- Aristotle task: `eb147738-3627-4d2b-96a5-6922b28b92e6`
- Toolchain: Lean 4.28.0 with Mathlib 4.28.0

## Outcome

Aristotle filled all nine proof holes in `PBReserve/Core.lean` without changing
any theorem statement. Independent local checks confirmed:

- `lake build` succeeds;
- the file contains no `sorry`, `axiom`, `admit`, or `implemented_by`;
- all nine theorems depend only on Mathlib's standard
  `propext`, `Classical.choice`, and `Quot.sound` axioms;
- the curvature-propagation, endpoint-exclusion, crossing-ratio, and
  raw-from-effective arguments are mathematically sound.

Aristotle's generated summary said “eight” proof holes, but the submitted and
returned files each contain nine theorem statements, and the source diff shows
that all nine `sorry`s were replaced. The count in that summary was a reporting
typo, not a missing proof.

## Scope

This is a completed formalization of the first proof layer, not of the entire
universal Poisson-binomial theorem. In particular, it assumes the normalized
Hillion--Johnson recurrence as `hstep`. The remaining end-to-end layers include
the mass-window products, variance reductions, scalar certificate assembly,
bounded-density max-atom argument, and the Hillion--Johnson-to-PBD bridge.

The formal `curvature_propagation` hypothesis uses `(r+1)d < 1`, slightly
narrower than the informal standalone recurrence bound `rd < 1`. This is
sufficient for every mass window and endpoint exclusion used in the current
proof, whose radius is already selected by `(r+1)d < 1`.
