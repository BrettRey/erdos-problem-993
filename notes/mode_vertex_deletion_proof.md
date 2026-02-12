# Mode-Vertex Deletion: Proved Lemmas Only

## Setup and Notation

Let `T` be a tree and `w in V(T)`.

Write

- `f(x) = I(T; x) = sum_{k>=0} f_k x^k`,
- `g(x) = I(T-w; x) = sum_{k>=0} g_k x^k`,
- `q(x) = I(T-N[w]; x) = sum_{k>=0} q_k x^k`.

The vertex-deletion recurrence is

  `f(x) = g(x) + x q(x)`,

hence coefficientwise

  `f_k = g_k + q_{k-1}` (with `q_{-1} := 0`).

We use:

- `alpha(G)`: independence number,
- `mode(P)`: last index of a maximum coefficient of polynomial `P`,
- `d(P)`: first descent index (`min{k>=1: p_k < p_{k-1}}`, or `deg(P)+1` if no descent).

Target invariant:

  `M(T,w): mode(I(T-w)) <= d(I(T))`.

---

## Proof Spine (Exact Dependency)

Let `m = mode(f)` and `d = d(f)`.

### Lemma 0 (Last-mode upper bound)

For every nonnegative coefficient sequence with last mode index `m`,

  `d <= m+1`.

#### Proof

Because `m` is the **last** mode index, `f_{m+1} < f_m`, so a descent occurs at
index `m+1`. Therefore first descent satisfies `d <= m+1`.

QED.

### Lemma 0.1 (Mode-shift implication)

If `mode(g) = m+1` and `M(T,w)` holds, then `d = m+1`.

#### Proof

`M(T,w)` gives `mode(g) <= d`, so `m+1 <= d`. Lemma 0 gives `d <= m+1`.
Hence `d = m+1`.

QED.

This isolates the key unresolved statement: prove `M(T,w)` globally.

---

## Theorem 1 (Alpha Characterization at a Deleted Vertex)

Let `N(w) = {u_1, ..., u_r}`. In `T-w`, let `C_i` be the component containing `u_i`.
Define

- `a_i = alpha(C_i)`,
- `b_i = alpha(C_i - u_i)`,
- `delta_i = a_i - b_i`.

Then:

1. `alpha(T-w) = sum_i a_i`.
2. The maximum size of an independent set of `T` containing `w` is

   `alpha_w = 1 + sum_i b_i`.

3. Consequently,

   `alpha(T) = max(sum_i a_i, 1 + sum_i b_i)`.

4. Therefore,

   `alpha(T) = alpha(T-w)` iff `sum_i delta_i >= 1`.

### Proof

(1) In `T-w`, components `C_i` are disjoint, so a maximum independent set is the union of componentwise maxima.
Hence `alpha(T-w) = sum_i alpha(C_i) = sum_i a_i`.

(2) If an independent set contains `w`, it cannot contain any `u_i`. On each component `C_i`, allowed vertices are exactly in `C_i-u_i`, giving contribution `b_i`. Including `w` adds `1`.
Hence `alpha_w = 1 + sum_i b_i`.

(3) Every independent set of `T` either excludes `w` (size at most `alpha(T-w)`) or includes `w` (size at most `alpha_w`). Taking the better of the two gives the formula.

(4) Substitute (1) and (2) into (3):

`alpha(T) = alpha(T-w)` iff `sum_i a_i >= 1 + sum_i b_i` iff `sum_i (a_i-b_i) >= 1` iff `sum_i delta_i >= 1`.

QED.

### Lemma 1.1 (Branch drop is binary)

For every `i`, `delta_i in {0,1}`.

#### Proof

By definition, `delta_i = alpha(C_i) - alpha(C_i-u_i) >= 0`.
Also, for any graph `H` and vertex `u`, every independent set of `H` either
avoids `u` (size at most `alpha(H-u)`) or contains `u` (then deleting `u`
leaves an independent set in `H-u`), so `alpha(H) <= alpha(H-u)+1`.
Applying this with `H=C_i`, `u=u_i` gives `delta_i <= 1`.
Hence `delta_i in {0,1}`.

QED.

### Lemma 1.2 (Degree-gap and support window)

Set

- `a := deg(g) = alpha(T-w)`,
- `b := deg(q) = alpha(T-N[w])`,
- `s := a-(b+1) = g_w`.

Then:

1. `f_k = g_k` for all `k >= b+2`.
2. `Delta f_k = Delta g_k` for all `k >= b+2`.
3. If `s >= 1`, then for each integer `j` with `0 <= j <= s-1`,
   `f_{a-j} = g_{a-j}` (top-`s` coefficients agree).

#### Proof

From `f_k = g_k + q_{k-1}` and `deg(q)=b`, we have `q_{k-1}=0` for all
`k-1 > b`, i.e. `k >= b+2`; this gives (1).
For (2), use Lemma 2's difference identity and note `q_k=q_{k-1}=0` when
`k >= b+2`.
For (3), if `0 <= j <= s-1` then
`a-j >= a-(s-1) = b+2`, so (1) applies at `k=a-j`.

QED.

Corollary (explicit top-degree split):

- `s >= 1` implies `f_a = g_a` (top coefficient comes from `g` only).
- `s >= 2` implies `f_a = g_a` and `f_{a-1} = g_{a-1}`.

This isolates where `xq` can affect mode/descent: only at indices `k <= b+1`.

### Lemma 1.3 (Leaf-vertex consequence for the gap)

If `deg_T(w)=1`, then `s=g_w <= 0`.

#### Proof

When `deg_T(w)=1`, `T-w` has exactly one branch component (`r=1` in Theorem 1).
So `s = g_w = delta_1 - 1` with `delta_1 in {0,1}` by Lemma 1.1.
Hence `s in {-1,0}` and in particular `s <= 0`.

QED.

Corollary:

`g_w >= 1` implies `deg_T(w) >= 2`.

### Lemma 1.4 (Leaf normal form)

Assume `deg_T(w)=1`, and let `u` be the unique neighbor of `w`.
Set `H := T-w` and `r(x) := I(H-N_H[u];x)`.
Then:

1. `g(x) = q(x) + x r(x)`,
2. `f(x) = (1+x) q(x) + x r(x)`.

#### Proof

In `H`, apply vertex deletion at `u`:

`I(H) = I(H-u) + x I(H-N_H[u])`.

But `H-u = T-N_T[w]`, so `I(H-u)=q`.
Also `I(H)=g` by definition, and `I(H-N_H[u])=r`.
Thus `g = q + x r`, proving (1).
Now use `f = g + xq`:

`f = (q + xr) + xq = (1+x)q + xr`,

proving (2).

QED.

This gives an exact decomposition for all leaf-deletion cases; any proof that
handles `(1+x)q + xr` uniformly in this branch would settle `M(T,w)` for
all leaf vertices.

### Proposition 1.5 (Leaf branch reduction to descent monotonicity)

Let `w` be a leaf of `T`, and set `H := T-w` (`g = I(H)`).
Assume:

1. `H` is unimodal,
2. `d(I(T)) >= d(I(H))`.

Then:

`mode(I(T-w)) < d(I(T))`.

#### Proof

By unimodality of `H`,

`d(I(H)) = mode(I(H)) + 1 = mode(g)+1`.

From (2),

`d(I(T)) >= d(I(H)) = mode(g)+1`,

hence `mode(g) < d(I(T))`.
Since `g=I(T-w)`, this is exactly the claim.

QED.

This reduces the leaf branch in a minimal-counterexample proof to one local
monotonicity statement:

`d(I(H + leaf_at_u)) >= d(I(H))` for all trees `H` and vertices `u`.

---

## Lemma 2 (Exact Boundary Identity / Agent 1 Implication Lemma)

For every index `k >= 0`,

  `Delta g_k = Delta f_k - (q_k - q_{k-1})`,

where `Delta p_k := p_{k+1} - p_k`.

Equivalently,

  `g_k >= g_{k+1}` iff `q_k - q_{k-1} >= f_{k+1} - f_k`.

### Proof

From `f_k = g_k + q_{k-1}` and `f_{k+1} = g_{k+1} + q_k`, subtract:

`f_{k+1} - f_k = (g_{k+1} - g_k) + (q_k - q_{k-1})`.

Rearrange to get the first identity, and move terms to get the equivalence.

QED.

---

## Corollary 3 (Sufficient Boundary Check at a Fixed Index)

Fix `m >= 0`. If

  `q_m - q_{m-1} >= f_{m+1} - f_m`,

then `g_m >= g_{m+1}`.

This condition is **sufficient**, not necessary, for the local one-step descent of `g` at `m`.

### Empirical note

In the checked subset with `alpha(T) = alpha(T-w)` (through `n <= 16`), this boundary inequality fails in `2,644` cases while `mode(I(T-w)) <= d(I(T))` still holds; see `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/mode_alignment_status.md`.

So the boundary inequality is not a necessary condition for mode alignment.

---

## Lemma 3 (Tail-Dominance Sufficient Condition for Mode Alignment)

If for all `k >= d(f)` we have

  `q_k - q_{k-1} >= f_{k+1} - f_k`,

then `mode(g) <= d(f)`.

### Proof

By Lemma 2, the displayed inequality is equivalent to `g_k >= g_{k+1}` for all
`k >= d(f)`. So `g` is nonincreasing from index `d(f)` onward. Therefore no
new maximum can occur at index `> d(f)`, and `mode(g) <= d(f)`.

QED.

Remarks:

- This is a clean sufficient route to `M(T,w)`.
- It is stronger than necessary (the boundary-only condition already fails as a
  necessary criterion by the `2,644` alpha-preserved counterexamples noted above).

---

## Verified Computational Status (for the target invariant)

Target invariant:

  `mode(I(T-w)) <= d(I(T))`.

Verified artifacts:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n18.json`:
  - `vertex_cases = 3,553,678`, failures `0`.
- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n20.json`:
  - `vertex_cases = 22,502,445`, failures `0`.

Combined through `n <= 20`: `26,056,123` vertex-cases, failures `0`.

Additional complete scan at `n = 21` (partition-merged):

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n21_mod8_merged.json`:
  - `trees = 2,144,505`,
  - `vertex_cases = 45,034,605`,
  - `total_failures = 0`.

Combined through `n <= 21`: `71,090,728` vertex-cases, failures `0`.

Additional complete scan at `n = 22` (partition-merged):

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n22_mod8_merged.json`:
  - `trees = 5,623,756`,
  - `vertex_cases = 123,722,632`,
  - `total_failures = 0`.

Combined through `n <= 22`: `194,813,360` vertex-cases, failures `0`.

Additional complete scan at `n = 23` (partition-merged):

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n23_mod8_merged.json`:
  - `trees = 14,828,074`,
  - `vertex_cases = 341,045,702`,
  - `total_failures = 0`.

Combined through `n <= 23`: `535,859,062` vertex-cases, failures `0`.

Leaf-focused exhaustive diagnostic through `n <= 16`:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n16.json`
  - `leaf_cases = 244,692`,
  - `fail_mode_le_d = 0`,
  - `fail_mode_lt_d = 0`,
  - `fail_d_monotone = 0` for `d(I(T)) >= d(I(T-w))`,
  - observed `d(I(T)) - mode(I(T-w)) in {1,2}`.

Leaf extension through `n <= 18` (combined artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n18_combined.json`):

- `leaf_cases = 1,723,516`,
- `0` failures of `mode(I(T-w)) <= d(I(T))`,
- `0` failures of strict `mode(I(T-w)) < d(I(T))`,
- `0` failures of `d(I(T)) >= d(I(T-w))`.
- exact gap split
  (`results/leaf_descent_gap_n18_combined.json`):
  `d(I(T)) - d(I(T-w)) in {0,1}` only.

## One-sided empirical laws (n<=20)

From `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/theorem_mining_one_sided_n20.json`,
with `E := [mode(I(T-w)) = d(I(T))]`:

- No-counterexample necessary implications:
  - `E => [d(I(T)) = mode(I(T)) + 1]`,
  - `E => [d(I(T)) = mode(I(T)) + 1 and deg(w) >= 2]`,
  - `E => [g_d > g_{d+1}]`.

- No-counterexample sufficient implications (in tested range):
  - `[d=mode+1 and deg(w)>=2 and g_d>=g_{d-1} and g_d>g_{d+1}] => E`,
  - `[d=mode+1 and deg(w)>=2 and g_d>g_{d-1} and g_d>g_{d+1}] => E`.

These are empirical guides, not proved statements.

## Current Gap (Formal)

What is already formal in this note:

1. Vertex-deletion alpha decomposition (`Theorem 1`).
2. Exact local boundary identity (`Lemma 2`).
3. Exact implication in the `mode(g)=mode(f)+1` branch (`Lemma 0.1`).

What remains unproved:

- The global invariant `M(T,w): mode(I(T-w)) <= d(I(T))` for all trees and vertices.

At this stage, all reduction paths still bottleneck on that statement.

Leaf-subproblem distilled:

- By Proposition 1.5, in minimal-counterexample mode it is enough to prove
  descent-index monotonicity under leaf attachment:
  `d(I(H + leaf_at_u)) >= d(I(H))`.
- Empirical support for this local statement:
  - no failures through `n<=16` in
    `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n16.json`,
  - no failures in an additional exhaustive `n=17,18` run
    (summary artifact
    `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n18_combined.json`).

## g_w-Sign Stratification (Empirical)

Using `g_w = alpha(T-w) - alpha_w` from
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/alpha_vertex_characterization.md`,
artifact:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_gw_n16_summary.json`

through `n <= 16`:

- total vertex cases: `497,379`,
- failures of `mode(I(T-w)) <= d(I(T))`: `0`,
- tight equalities `mode(I(T-w)) = d(I(T))`: `3,009`.

Sign split:

- `g_w <= 0`: `231,797` cases, `0` tight equalities (in this `n<=16` sample).
- `g_w >= 1`: `265,582` cases, all `3,009` tight equalities occur here.

Extended check through `n<=18` (artifact
`/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_gw_n18_summary.json`)
finds `2` tight equalities with `g_w<=0` versus `19,009` with `g_w>=1`.

So a sharper theorem target is:

1. prove strict inequality `mode(I(T-w)) < d(I(T))` when `g_w <= 0`,
2. prove non-strict inequality `mode(I(T-w)) <= d(I(T))` when `g_w >= 1`.

That split may be easier than a single uniform argument.
