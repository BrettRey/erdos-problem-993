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

### Creative route H_R (ratio monotonic bridge)

Define

`R_k := g_k / f_k`  (equivalently check `g_{k+1} f_k <= g_k f_{k+1}`).

Hypothesis H_R:
for leaf cases, `R_k` is nonincreasing on the prewindow
`k=0,...,d(g)-2`.

If H_R holds, then `d(f)>=d(g)` follows immediately:

- for `k<=d(g)-2`, we have `Delta g_k >= 0`;
- `R_{k+1}<=R_k` means
  `g_{k+1} f_k <= g_k f_{k+1}`;
- since `g_{k+1}>=g_k>0` in this prewindow,
  this forces `f_{k+1}>=f_k`, i.e. `Delta f_k>=0`.

Therefore no descent of `f` can occur before `d(g)`.

So H_R is a single-point creative reduction:
prove prewindow monotonicity of `g_k/f_k`, and leaf-step monotonicity is done.

Exhaustive evidence for H_R through `n<=18`
(artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_ratio_hypothesis_n18.json`):

- `leaf_cases = 1,723,516`,
- checked prewindow pair inequalities: `9,189,052`,
- violations on prewindow: `0`.

Diagnostic (same artifact):
the same ratio monotonicity fails on the full common prefix in `44,481` cases,
so this appears genuinely **prewindow-local**, not globally true.

Random stress (artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_ratio_hypothesis_random.json`):

- `4,000` random Prüfer trees (`n in [25,220]`), one random leaf each: `0` failures;
- `800` random Prüfer trees (`n in [25,180]`), all leaves (`29,854` leaf-cases):
  `0` failures.

### Hard-core covariance reformulation (conceptual bridge)

Let `lambda>0`, and hard-core measure on independent sets of `T`:

`P_lambda(S) proportional to lambda^{|S|}`.

With leaf `w`, define:

- `f=I(T)`, `g=I(T-w)`,
- `Theta(lambda)=g(lambda)/f(lambda)=P_lambda(w notin S)`,
- `mu_f(lambda)=lambda f'(lambda)/f(lambda)=E_lambda[|S|]`,
- `mu_g(lambda)=lambda g'(lambda)/g(lambda)=E_lambda[|S| | w notin S]`.

Then

`Cov_lambda(1_{w notin S}, |S|) = Theta(lambda) * (mu_g(lambda)-mu_f(lambda))`
and equivalently
`lambda Theta'(lambda)=Cov_lambda(1_{w notin S}, |S|)`.

So covariance sign is exactly the sign of `mu_g-mu_f` (since `Theta>0`).

#### Empirical verdict

Global covariance negativity is false:
`mu_g <= mu_f` fails often at larger `lambda`.

Windowed version is strikingly true:

If `mu_f(lambda) <= d(g)-1`, then `mu_g(lambda) <= mu_f(lambda)`.

Exhaustive artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_covariance_window_n18.json`
(`n<=18`, `1,723,516` leaf-cases, 49-point `lambda` grid from `1e-6` to `1e6`):

- global checks: `84,452,284`, global violations: `862,446`,
- window checks (`mu_f<=d(g)-1`): `41,573,584`,
- window violations: `0`.

Extended one-level run at `n=19`:
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_covariance_window_n19.json`

- `leaf_cases = 2,902,599`,
- global checks: `142,227,351`, global violations: `1,524,805`,
- window checks: `70,286,945`,
- window violations: `0`.

Combined (`n<=18` + `n=19`) artifact:
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_covariance_window_n19_combined.json`

- `leaf_cases = 4,626,115`,
- global checks: `226,679,635`, global violations: `2,387,251`,
- window checks: `111,860,529`,
- window violations: `0`.

This is the current best conceptual synthesis:
the right covariance monotonicity is **windowed by the descent boundary**, not
global in `lambda`.

### Bridge problem (covariance -> coefficient ratio)

Current leaf bottleneck can now be stated as two candidate statements:

- `H_R` (coefficient side): `R_k=g_k/f_k` is nonincreasing for `k<=d(g)-2`.
- `H_CW` (ensemble side): for all `lambda` with `mu_f(lambda)<=d(g)-1`,
  `mu_g(lambda)<=mu_f(lambda)`.

Empirically both hold through exhaustive ranges above.

The conceptual gap is to prove a **windowed MLR bridge**:

`H_CW  =>  H_R` (on the same boundary window).

This is not true as a generic statement for arbitrary sequences, so any proof
must exploit tree/branch structure of `(f,g)` (leaf-realizable pairs).

Practical consequence:
if this bridge is proved for leaf-realizable pairs, then leaf-step monotonicity
`d(f)>=d(g)` follows (via the earlier `H_R` implication), closing the leaf
branch of the mode-alignment proof spine.

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
