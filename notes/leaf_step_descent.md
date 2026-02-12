# Leaf-Step Descent Monotonicity (Local Subproblem)

## Setup

Let `T` be a tree, `w` a leaf, and `u` the unique neighbor of `w`.
Set `H := T-w`.

Write

- `f := I(T)`,
- `g := I(H) = I(T-w)`,
- `q := I(H-u) = I(T-N[w])`,
- `r := I(H-N_H[u])`.

Exact identities:

1. `f = g + x q`,
2. `g = q + x r`,
3. `f = (1+x)q + x r`.

(2)-(3) are the leaf normal form.

## Target local inequality

The distilled leaf-step target is:

`d(f) >= d(g)` where `d(P)` is first descent index (`deg+1` if none).

Why this matters:

- In minimal-counterexample mode, `H` is unimodal, so
  `mode(g) = d(g)-1`.
- Therefore `d(f) >= d(g)` implies
  `mode(g) < d(f)`, i.e. strict leaf branch
  `mode(I(T-w)) < d(I(T))`.

So proving this single local monotonicity would settle all leaf-vertex cases
in the mode-alignment step.

Conjecture L (sharpened local form):
for every leaf attachment `T = H + leaf_at_u`,

`d(I(T)) - d(I(H)) in {0,1}`.

This implies the required monotonicity `d(I(T)) >= d(I(H))`.

## Equivalent difference form

For all `k >= 0`:

`Delta f_k = Delta g_k + Delta q_{k-1}`

with `Delta p_k := p_{k+1}-p_k`.

Hence early-descent transfer `d(f) < d(g)` would require some
`k < d(g)` with

`Delta q_{k-1} < -Delta g_k`.

So any proof of `d(f) >= d(g)` can focus on ruling out this negative
`Delta q` overcompensation window before the first descent of `g`.

### Critical prewindow scan (first-descent indexing used by scripts)

In script indexing (`d(P)=min{i>=1: p_i<p_{i-1}}`), the only indices that can
force `d(f)<d(g)` are

`k = 0,1,...,d(g)-2`

since descent at index `i` corresponds to `Delta` at `k=i-1`.

Artifact:
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_critical_window_n18_combined.json`
(exhaustive leaf-cases through `n<=18`, `1,723,516` cases).

Observed:

- `0` failures of prewindow nonnegativity:
  `Delta g_k + Delta q_{k-1} >= 0` for all `k<=d(g)-2`.
- Boundary prewindow margin is always strictly positive:
  `Delta g_{d(g)-2} + Delta q_{d(g)-3} > 0` in all cases.
- Gap law confirmed again:
  `d(f)-d(g) in {0,1}` only.

So the leaf-step law is fully concentrated in this prewindow inequality family.
Empirically there is no need to control `k>=d(g)-1` to get `d(f)>=d(g)`.

### Conditional lemma (prefix-monotone criterion)

If `Delta q_t >= 0` for all `t <= d(g)-2`, then `d(f) >= d(g)`.

Proof:
for `k < d(g)`, we have `Delta g_k >= 0` by definition of first descent.
By the hypothesis, `Delta q_{k-1} >= 0` for all such `k`.
Therefore

`Delta f_k = Delta g_k + Delta q_{k-1} >= 0`

for every `k < d(g)`. So `f` has no descent before `d(g)`, i.e.
`d(f) >= d(g)`.

QED.

Corollary (descent-index form):
if `d(q) >= d(g)-1`, then `d(f) >= d(g)`.

Reason: by definition of first descent, `d(q) >= d(g)-1` implies
`Delta q_t >= 0` for all `t <= d(g)-2`, so the criterion applies.

Empirical frequency through `n<=16` (leaf cases), artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_prefix_monotone_n16.json`:
this criterion holds in `48,938 / 244,692` cases (about `20%`), so it is
useful but far from the full statement.

For the weaker descent-index condition (artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_dq_dg_condition_n16.json`):
`d(q) >= d(g)-1` holds in `155,619 / 244,692` cases (about `63.6%`).

New empirical pattern (artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_dq_vs_dg_n16.json`):
through `n<=16`, always

`d(q) <= d(g)`.

Observed distribution of `d(q)-d(g)`:
`{-6,-5,-4,-3,-2,-1,0}` with maximum `0`.

Extended check through `n<=18` (combined summary
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_dq_vs_dg_n18_combined.json`):

- `leaf_cases = 1,723,516`,
- `0` violations of `d(q) <= d(g)`.

Conjecture (auxiliary):
for leaf attachment normal form `g=q+xr`, one always has `d(q) <= d(g)`.
If proved, this would significantly constrain the only failure mechanism
`Delta q_{k-1} < -Delta g_k` before `d(g)`.

Finer sign split at the first post-descent `q`-difference from the new
critical-window artifact (`leaf_critical_window_n18_combined.json`):

- `Delta q_{d(g)-1} > 0` occurs in `35,039` cases, all with `d(f)-d(g)=1`.
- No case with `d(f)-d(g)=0` has `Delta q_{d(g)-1} > 0`.

So empirical implication:

`Delta q_{d(g)-1} > 0  =>  d(f)=d(g)+1`.

## Branch-product normal form at the attachment vertex

Let `c_1,...,c_m` be neighbors of `u` in `H`, and let `S_i` be the component
of `H-u` containing `c_i`.
Define

- `A_i := I(S_i)`,
- `B_i := I(S_i-c_i)`.

Then:

- `q = prod_i A_i`,
- `r = prod_i B_i`,
- `g = q + x r`,
- `f = (1+x)q + x r`.

So the leaf-step inequality `d(f) >= d(g)` is a pure coefficient statement for
these two branch products.

This removes graph topology from the local step: only the paired branch
polynomials `(A_i,B_i)` matter.

## Exhaustive evidence

Artifacts:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n16.json`
- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_descent_monotone_n18.json`
- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n18_combined.json`

Combined through `n <= 18`:

- `leaf_cases = 1,723,516`,
- `0` failures of `d(I(T)) >= d(I(T-w))`,
- `0` failures of strict `mode(I(T-w)) < d(I(T))`,
- observed `d(I(T)) - mode(I(T-w)) in {1,2}` (from the `n<=16` full
  leaf diagnostic).

Additional distribution detail through `n<=16`:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_descent_gap_n16.json`
  gives
  - `d(I(T)) - d(I(T-w)) in {0,1}`,
  - `d(I(T)) - mode(I(T-w)) in {1,2}`.

Extended distribution through `n<=18`:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_descent_gap_n18_combined.json`
  gives
  - `d(I(T)) - d(I(T-w)) = 0` in `933,089` leaf-cases,
  - `d(I(T)) - d(I(T-w)) = 1` in `790,427` leaf-cases,
  - no other values observed.

Degree split at the attachment site (`deg_H(u)`), artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_descent_split_by_deg_n18.json`:

- `deg_H(u)=1`: `452,658` cases, still only gaps `{0,1}`.
- `deg_H(u)>=2`: `1,270,858` cases, still only gaps `{0,1}`.

## Route that fails (explicit)

The naive condition `mode(q) >= d(g)-1` is false in about 80% of tested leaf
cases through `n<=16`, so mode-position arguments on `q` alone are not enough.

Also, abstract coefficientwise domination is insufficient:
`r <= q` does **not** imply `d((1+x)q + xr) >= d(q + xr)` in general.

Concrete non-graph counterexample:

- `q = [1,7,4,8,5]`,
- `r = [0,4,2,5,1]` (entrywise `r<=q`),
- `g=q+xr = [1,7,8,10,10,1]`, so `d(g)=5`,
- `f=(1+x)q+xr = [1,8,15,14,18,6]`, so `d(f)=3<5`.

Even adding log-concavity on both `q,r` is still insufficient in abstract:

- `q=[3,8,7,2]` (LC),
- `r=[0,1,6,2]` (LC, and `r<=q`),
- `g=[3,8,8,8,2]` with `d(g)=4`,
- `f=[3,11,16,15,4]` with `d(f)=3<4`.

So any proof must use tree/forest structure beyond generic nonnegative or
log-concave coefficient constraints.
