# Perturbation Reduction for Degree-2 Leaf Steps

This note records the clean reduction and the current hard gap.

## Exact reduction

If `T` has a leaf `v` whose neighbor `u` has degree 2 (other neighbor `w`),
set:

  `T'' = T - {u, v}`.

Then:

  `I(T; x) = (1+x) I(T''; x) + x I(T'' - w; x)`.

Define:

  `f = I(T'')`, `g = I(T''-w)`, `h = (1+x)f + xg`.

So the degree-2 step reduces unimodality of `I(T)` to unimodality of a
structured perturbation `h`.

## What is true / false

### False route (too strong)

The tail ratio condition

  `g_{k+1}/f_{k+1} <= g_k/f_k` for `k >= d(f)`

is **not** generally true.

Computed with `scripts/perturb_ratio_check.py` through `n <= 20`
(`results/perturb_ratio_n20.json`):

  - `leaf_cases = 3,258,826`
  - `ratio_failures = 42,196`

Example witness:

  - `n = 9`, `graph6 = H?AE@`g`
  - `f = [1,7,15,11,2]`, `g = [1,6,10,5,1]`
  - `d(f)=3`, and `g_4/f_4 = 1/2 > 5/11 = g_3/f_3`.

So ratio monotonicity is not the right theorem target.

### Strong empirical invariants (holding so far)

From `results/perturb_ratio_n20.json`:

  - `mode_failures = 0` for `mode(g) <= d(f)` across all 3,258,826 degree-2 leaf cases.
  - `tail_h_failures = 0` for `h` being nonincreasing on `k >= d(f)+1`.
  - `pre_h_failures = 1,938,197`: early descent before `d(f)` is common.

Interpretation:

  - We should **not** require monotone increase up to `d(f)`.
  - The right target is “no tail rise after `d(f)+1`”.

## Mode-alignment evidence (full vertex scans)

Independent scan (all vertices, not just degree-2 reductions):

  - `results/mode_alignment_n18.json`: 3,553,678 vertex cases, 0 failures.
  - `results/mode_alignment_n20.json`: 22,502,445 vertex cases, 0 failures.

Extended complete scans:

  - `results/mode_alignment_n21_mod8_merged.json`: 45,034,605 vertex cases, 0 failures.
  - `results/mode_alignment_n22_mod8_merged.json`: 123,722,632 vertex cases, 0 failures.
  - `results/mode_alignment_n23_mod8_merged.json`: 341,045,702 vertex cases, 0 failures.

Combined evidence to `n <= 23`:

  - 535,859,062 vertex cases with no failure of `mode(I(T-w)) <= d(I(T))`.
  - `max |mode(I(T-w)) - mode(I(T))| = 1` in tested range.

Important negative fact:

  - Mode alignment between `I(T-w)` and `I(T-N[w])` is not usable.
    The gap can be much larger than 1 (up to 9 in earlier scans).

## Practical proof target

For the degree-2 reduction `h = (1+x)f + xg`, the most useful target is:

1) prove `mode(g) <= d(f)` (equivalently `Δg_k <= 0` for `k >= d(f)`), and
2) prove boundary control at `k = d(f)` for `h`.

Given (1), tail monotonicity of `h` from `d(f)+1` is straightforward.
The real unresolved step is the boundary index.

## Open tasks

1) Prove `mode(I(T-w)) <= d(I(T))` for all trees/vertices.
2) Derive a boundary inequality at `k=d(f)` that is compatible with observed
   early-descent behavior.
3) If (1) remains hard globally, prove it on the finite reduced kernel class
   from `notes/erdos_plan_onepage.md` and close by exhaustive verification.

## Leaf-step monotonicity reduction (new)

For a leaf attachment `T = H + leaf_at_u`, write

`I(T) = I(H) + x I(H-u)`.

In minimal-counterexample mode, `I(H)` is unimodal. Therefore, if one proves

`d(I(T)) >= d(I(H))`,

then automatically

`mode(I(H)) = d(I(H))-1 <= d(I(T))-1 < d(I(T))`,

which gives the strict leaf branch
`mode(I(T-w)) < d(I(T))` for leaf `w`.

Empirical support:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n16.json`
  reports zero failures of `d(I(T)) >= d(I(T-w))` across `244,692` leaf-cases
  (`n<=16`).
- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n18_combined.json`
  extends this to `1,723,516` leaf-cases through `n<=18` with zero failures.
