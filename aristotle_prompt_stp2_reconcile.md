# Formalization Request: reconcile `Formal/STP2Closure.lean` with prior Aristotle progress

Please work only on `Formal/STP2Closure.lean` in the repo-root Lean project.

## Project / dependency constraints

- Upload the repository root as the Lean project.
- The canonical Lean project is the repo root.
- Dependencies must stay within Lean `v4.28.0`.
- Do not switch to the old `lean/` project as the main target.

## Goal

Bring `Formal/STP2Closure.lean` up to the level of the earlier Aristotle-assisted auxiliary progress, while **explicitly leaving the main open theorem open**.

The target outcome is:

1. `Formal/STP2Closure.lean` builds successfully.
2. The file contains the proved auxiliary/base-case lemmas that are actually valid.
3. The file **still leaves**
   - `stp2_conv_closure`
   - `stp2_multi_child_closure`
   as `sorry`.
4. False lemmas are removed or replaced with correct guarded/base-case statements.

## Important context

This is **not** a request to solve the open combinatorial problem.
The broad theorem

```lean
stp2_conv_closure
```

is still open in this project and should remain open.

An earlier Aristotle run already made useful progress on the auxiliary side:

- `conv_comm`
- `conv_nonneg`
- leaf-factor definitions / lemmas
- degree-0 base case
- structural observations about the current `isLC` definition using `ℕ` subtraction

That earlier run also identified that some intended lemmas in the old file are not valid as stated.

## What to do

Use the prior work as guidance, but make the canonical target be `Formal/STP2Closure.lean`.

In particular:

- Prove easy/valid lemmas like:
  - `conv_comm`
  - `conv_nonneg`
- If appropriate, add valid helper lemmas for:
  - zero second argument
  - leaf-factor STP2
  - degree-0 closure
- If `lc_conv` is false under the current `isLC` definition, do **not** force a bogus proof.
  Replace it with a correct weaker statement or remove it if unused.
- If `cb_expansion_form1` is false or incorrectly stated, do **not** fake it.
  Either correct the statement precisely or leave it out if it is not needed for the file to be a truthful record of progress.

## Existing local context

The repo may also contain earlier Aristotle scratch work under `lean/`.
You may inspect that as reference material, but the canonical file to edit is:

- `Formal/STP2Closure.lean`

## Non-goals

- Do not try to prove `stp2_conv_closure`.
- Do not try to prove `stp2_multi_child_closure`.
- Do not make broad unrelated refactors outside `Formal/STP2Closure.lean`.

## Deliverable

A truthful, compiling `Formal/STP2Closure.lean` that captures the real auxiliary progress and leaves only the genuinely open theorem/corollary as `sorry`.
