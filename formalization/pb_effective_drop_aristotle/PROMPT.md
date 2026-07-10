# Aristotle request: endpoint-aware Poisson-binomial reserve core

Work only in this minimal Lean 4 project. The toolchain is Lean `v4.28.0` with
Mathlib `v4.28.0`.

## Goal

Fill every `sorry` in `PBReserve/Core.lean` and make `lake build` pass with no
`sorry`, no new `axiom`, and no weakened theorem statement.

The file formalizes a bounded, high-value slice of a new proof:

- iteration of the normalized Hillion--Johnson curvature recurrence;

- exclusion of a support endpoint from the forced curvature window;

- the first-crossing ratio lower bound;

- conversion from effective Turán drop to raw drop at a strict descent.

Read `PROOF_CONTEXT.md` before editing.

## Scope boundary

The recurrence hypothesis `hstep` is intentionally assumed. Do not claim that
this packet formalizes Hillion--Johnson's cubic inequalities or the complete
Poisson-binomial theorem.

## Requirements

- Preserve all theorem statements unless one is false. If you believe a
  statement is false, do not weaken it silently: explain the issue and provide
  a concrete counterexample in the output summary.

- Prefer transparent algebraic/order proofs over heavy automation when the
  proof is short.

- Keep changes inside this project.

- Run `lake build` before returning the result.

## Deliverable

A compiling project in which `PBReserve/Core.lean` contains no `sorry` and no
new assumptions beyond the hypotheses already stated.
