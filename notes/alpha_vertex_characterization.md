# Vertex-Deletion Alpha Characterization

Let `T` be a tree and `w` a vertex. Let `N(w) = {u_1, ..., u_r}`.
In `T-w`, let `C_i` be the component containing `u_i`.

Define:

- `a_i = alpha(C_i)`
- `b_i = alpha(C_i - u_i)`
- `delta_i = a_i - b_i` (binary: `delta_i in {0,1}`)

Then:

  `alpha(T-w) = sum_i a_i`,

and the best independent set in `T` that includes `w` has size

  `alpha_w = 1 + sum_i b_i`.

So:

  `alpha(T) = max(sum_i a_i, 1 + sum_i b_i)`.

Why `delta_i in {0,1}`:

- `delta_i >= 0` is immediate from `C_i-u_i` being an induced subgraph of `C_i`.
- `delta_i <= 1` because every independent set of `C_i` either avoids `u_i`
  (size at most `alpha(C_i-u_i)`) or contains `u_i` (then removing `u_i`
  leaves an independent set in `C_i-u_i`).
  Hence `alpha(C_i) <= alpha(C_i-u_i)+1`.

## Consequence

`alpha(T) = alpha(T-w)` iff

  `sum_i a_i >= 1 + sum_i b_i`

iff

  `sum_i delta_i >= 1`.

This is equivalent to saying `w` is not in every maximum independent set
(`w` not in the core of `T`).

Equivalent branch-local form:

`alpha(T) = alpha(T-w)` iff at least one branch root `u_i` is forced in
maximum independent sets of `C_i` (`delta_i = 1` for some `i`).

## Gap parameter

Define

  `g_w = alpha(T-w) - alpha_w = sum_i delta_i - 1`.

Since each `delta_i in {0,1}`, `sum_i delta_i` counts exactly how many branch
roots `u_i` are forced in maximum independent sets of their branch components.

Interpretation:

- `g_w >= 0` iff `alpha(T) = alpha(T-w)`.
- Larger `g_w` means a wider gap between maximum sets avoiding vs including `w`.
- `g_w >= 1` iff at least two branch contributions force exclusion of `w` in
  any maximum-set comparison.

This parameter is structurally useful for mode/descent analysis because the
second term in

  `I(T) = I(T-w) + x I(T-N[w])`

has degree exactly `alpha_w`.

## Use in the mode-alignment proof

Write

- `f = I(T)`,
- `g = I(T-w)`,
- `q = I(T-N[w])`,
- `f = g + x q`.

Then `deg(g) = alpha(T-w)` and `deg(xq) = alpha_w`.
So `g_w = alpha(T-w) - alpha_w` is exactly the top-degree separation between
the two summands in the deletion recurrence.

This gives a clean structural split:

- `g_w < 0`: top degree is driven by the `xq` term.
- `g_w = 0`: both summands meet at top degree.
- `g_w > 0`: top degree is driven by `g`.

The unresolved mode statement `mode(g) <= d(f)` appears empirically stable
across all three regimes; proving it likely needs case analysis by `g_w`.

## Local computational sanity check

Verified by exhaustive checks through `n <= 12`:

- identity `alpha(T) = max(sum a_i, 1+sum b_i)` held in all `11,005` vertex cases.
- equivalence `alpha(T)=alpha(T-w)` iff `sum delta_i >= 1` held in all cases.

Additional stratified scan with the mode invariant (`d`-convention: first descent,
or `deg+1` if no descent), artifact:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_gw_n12_dconv.json`

Observed through `n <= 12`:

- No failures of `mode(I(T-w)) <= d(I(T))` in any `g_w` class.
- Equality `mode(I(T-w)) = d(I(T))` appears only for `g_w >= 1`
  (none observed for `g_w <= 0` in this range).
