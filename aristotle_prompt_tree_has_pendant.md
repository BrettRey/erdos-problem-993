# Formalization Request: `tree_has_pendant` in `Formal/P3.lean`

Please fill only the theorem `tree_has_pendant` in `Formal/P3.lean`.

## Project / dependency constraints

- Upload the repository root as the Lean project.
- The canonical Lean project is the repo root, not `lean/`.
- Dependencies should not exceed Lean `v4.28.0`.
- The root project already satisfies this:
  - `lean-toolchain`: `leanprover/lean4:v4.28.0`
  - `lakefile.toml`: only `mathlib` at `v4.28.0`

## Goal

Prove:

```lean
theorem tree_has_pendant (hT : G.IsTree) (hcard : 1 < Fintype.card V) :
    ∃ ℓ r : V, IsPendant G ℓ r := by
```

where

```lean
def IsPendant (ℓ r : V) : Prop :=
  G.Adj ℓ r ∧ ∀ w, G.Adj ℓ w → w = r
```

and this theorem lives in `Formal/P3.lean`.

## Mathematical proof

A finite tree with at least two vertices has a leaf `ℓ`.
Because the tree is connected and has more than one vertex, `ℓ` has a neighbor `r`.
Because `ℓ` is a leaf, every neighbor of `ℓ` equals `r`.
Therefore `IsPendant G ℓ r`.

## Requirements

- Work in the existing file/API; avoid unrelated refactors.
- Prefer existing Mathlib lemmas about finite trees, connected graphs, leaves, and degree-1 vertices.
- If Mathlib already has a theorem stating that a finite tree with at least two vertices has a leaf, use it.
- If a tiny helper lemma is needed, that is fine, but the target theorem itself should be completed.
- Do not touch `Formal/STP2Closure.lean` in this run.
- Keep the rest of the project compiling under Lean `v4.28.0`.

## Context

`Formal/P3.lean` is otherwise complete except for this standard tree-existence lemma.
This is intentionally a bounded formalization task, not an attempt to solve the open STP2 closure problem.
