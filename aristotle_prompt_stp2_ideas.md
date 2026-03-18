# Formalization / ideas request: sharpen the open `STP2` closure step

Please work only in the repo-root Lean project (Lean `v4.28.0`) and focus on:

- `Formal/STP2Closure.lean`

Do **not** switch to any older scratch project.

## Goal

This is **not** a request to prove the full open theorem. The theorem

```lean
stp2_conv_closure
```

should remain open unless you genuinely finish it.

Instead, the goal is to find the strongest **correct, formalizable reduction** around the current gap. Please prioritize actual Lean progress over speculative prose.

## Current state

The file already builds and currently leaves exactly two `sorry`s:

- `stp2_conv_closure`
- `stp2_multi_child_closure`

The auxiliary cleanup is already done:

- `conv_comm`
- `conv_nonneg`
- `conv_zero`
- `conv_zero_left`
- `conv_zero_right`
- leaf-factor definitions and basic lemmas
- false lemmas like the old `lc_conv` have been removed/documented

## What would be high-value progress

Please try to produce as many of the following as are actually true.

### 1. Repair the multi-child wrapper

The current multi-child statement folds convolution from the zero sequence, which is not the right identity for nonempty products.

Please:

- define the correct convolution identity sequence `δ0` (delta at zero), if needed
- prove the identity lemmas
  - `conv f δ0 = f`
  - `conv δ0 f = f`
- restate / repair `stp2_multi_child_closure` so it is mathematically correct
- ideally make it a corollary that reduces the multi-child case to the 2-child case

If this is awkward in the current file, it is fine to replace the current theorem with a more honest theorem statement that clearly captures the intended induction reduction.

### 2. Derive an exact expansion or reformulation for the open binary step

Please try to derive a **correct** lemma expanding

```lean
ladder (conv I₁ I₂) (conv E₁ E₂) k
```

into a usable sum or equivalent reformulation.

The old `cb_expansion_form1` was not trustworthy as stated, so do **not** fake that theorem.

If you can prove a corrected expansion under finite-support hypotheses, that is valuable even if it does not close the inequality.

### 3. Prove nontrivial special cases of `stp2_conv_closure`

Examples of potentially useful special cases:

- one factor is the leaf pair `leafI, leafE`
- one child pair has support at most `1`
- one child pair has support at most `2`
- one factor is exactly `δ0`

If a natural special case is false, do not force it. Either skip it or replace it with the strongest true version.

### 4. Find the weakest extra assumptions that make closure provable

If you discover that closure becomes provable with an additional assumption, such as:

- a corrected boundary condition at `k = 0`
- a stronger support bound
- a monotonicity / ratio condition beyond current `isSTP2`
- a tree-specific realizability hypothesis that can be stated cleanly

then formalize that theorem if possible and leave a short comment explaining the role of the extra assumption.

## Non-goals

- Do not add bogus proofs.
- Do not reintroduce false lemmas just to fill the file.
- Do not make broad unrelated refactors outside `Formal/STP2Closure.lean` unless a tiny import/helper change is genuinely needed.

## Deliverable

A truthful improvement to `Formal/STP2Closure.lean` that leaves the real open step honest, but strengthens the surrounding formal structure as much as possible.

If the full theorem still remains open, that is fine. The goal is to improve the reduction and isolate the exact mathematical obstruction.
