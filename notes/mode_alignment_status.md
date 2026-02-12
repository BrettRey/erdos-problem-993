# Mode Alignment Status (Current)

Target invariant:

  `mode(I(T-w)) <= d(I(T))` for all trees `T` and vertices `w`.

## Exhaustive verified artifacts

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n18.json`
  - `trees = 205,003`
  - `vertex_cases = 3,553,678`
  - `failures = 0`
  - `max_abs_mode_diff = 1`

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n20.json`
  - `trees = 1,141,020` (n = 19..20)
  - `vertex_cases = 22,502,445`
  - `failures = 0`
  - `max_abs_mode_diff = 1`

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n21_mod8_merged.json`
  - `trees = 2,144,505`
  - `vertex_cases = 45,034,605`
  - `total_failures = 0`
  - `max_abs_mode_diff = 1`
  - `is_complete = true`

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n22_mod8_merged.json`
  - `trees = 5,623,756`
  - `vertex_cases = 123,722,632`
  - `total_failures = 0`
  - `max_abs_mode_diff = 1`
  - `is_complete = true`

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n23_mod8_merged.json`
  - `trees = 14,828,074`
  - `vertex_cases = 341,045,702`
  - `total_failures = 0`
  - `max_mode_gap = 0`
  - `max_abs_mode_diff = 1`
  - `is_complete = true`

Combined through `n <= 23`:

- `trees = 23,942,358`
- `vertex_cases = 535,859,062`
- `failures = 0`
- `max_abs_mode_diff = 1`

Leaf-only diagnostic (new):

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n16.json`
  - exhaustive through `n <= 16`,
  - `leaf_cases = 244,692`,
  - `mode(I(T-w)) <= d(I(T))`: `0` failures,
  - strict `mode(I(T-w)) < d(I(T))`: `0` failures,
  - leaf-step descent monotonicity `d(I(T)) >= d(I(T-w))`: `0` failures,
  - observed gap `d(I(T)) - mode(I(T-w))` is always `1` or `2`.
  - finer split (`results/leaf_descent_gap_n16.json`):
    `d(I(T)) - d(I(T-w)) in {0,1}`.
  - auxiliary split (`results/leaf_dq_vs_dg_n16.json`):
    `d(I(T-N[w])) - d(I(T-w)) <= 0` in all leaf-cases.

Leaf extension to `n <= 18`:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/leaf_mode_descent_n18_combined.json`
  - `leaf_cases = 1,723,516`,
  - `0` failures of `mode(I(T-w)) <= d(I(T))`,
  - `0` failures of strict `mode(I(T-w)) < d(I(T))`,
  - `0` failures of `d(I(T)) >= d(I(T-w))`.
  - with exact split
    (`results/leaf_descent_gap_n18_combined.json`):
    `d(I(T)) - d(I(T-w))` is always `0` or `1`.
  - and auxiliary relation
    (`results/leaf_dq_vs_dg_n18_combined.json`):
    `d(I(T-N[w])) <= d(I(T-w))` in all checked leaf-cases.

## In-progress run (n = 24, mod 16 partitions)

Current partial merge artifact:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_n24_partial_mod16_merged.json`
  - `is_complete = false`
  - `partitions_found = [0,1,2,3,4,5,6,7,9,10,12,13]` (latest snapshot)
  - `total_trees = 28,156,954`
  - `total_vertex_cases = 675,766,896`
  - `total_failures = 0`
  - `max_mode_gap = 0`
  - `max_abs_mode_diff = 1`

## Local diagnostic checks

One-off exhaustive checks through `n <= 16`:

- Cases with `mode(I(T-w)) = mode(I(T)) + 1`:
  - `3,009` cases
  - `0` violations of `d(I(T)) = mode(I(T)) + 1`.

- Subset with `alpha(T) = alpha(T-w)`:
  - `265,582` vertex cases
  - `0` failures of `mode(I(T-w)) <= d(I(T))`
  - boundary inequality
    `q_m - q_{m-1} >= f_{m+1} - f_m`
    fails in `2,644` cases (sufficient, not necessary).

Alpha-gap (`g_w`) stratified check through `n <= 12`:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_gw_n12_dconv.json`
- `11,005` vertex cases, `0` failures of `mode(I(T-w)) <= d(I(T))`.
- In this range, equality `mode(I(T-w)) = d(I(T))` is observed only when `g_w >= 1`.

Extended `g_w`-sign summary through `n <= 16`:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_gw_n16_summary.json`
- `497,379` vertex cases, `0` failures of `mode(I(T-w)) <= d(I(T))`.
- `3,009` tight equalities `mode(I(T-w)) = d(I(T))`, all with `g_w >= 1`.
- no tight equalities observed with `g_w <= 0` in this range.

Extended `g_w`-sign summary through `n <= 18`:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/mode_alignment_gw_n18_summary.json`
- `3,553,678` vertex cases, `0` failures of `mode(I(T-w)) <= d(I(T))`.
- `19,011` tight equalities total:
  - `19,009` with `g_w >= 1`,
  - `2` with `g_w <= 0`.

So `g_w >= 1` captures essentially all tight cases in this range, but not all.

## Equality-condition tests (`E := mode(I(T-w)) = d(I(T))`)

Artifact:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/tight_mode_equivalence_n20.json`

Tested biconditional:

  `E` iff `(alpha(T)=alpha(T-w) and d(I(T))=mode(I(T))+1 and deg(w)>=2)`.

Result:

- `equivalence_confirmed = false`
- `lhs_true_rhs_false = 6`
- `lhs_false_rhs_true = 10,929,098`

## One-sided mined laws (empirical, n<=20)

Artifact:

- `/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/results/theorem_mining_one_sided_n20.json`

No-counterexample necessary implications:

- `E => (d = mode(I(T)) + 1)` (`N_D`)
- `E => (d = mode(I(T)) + 1 and deg(w) >= 2)` (`N_DG2`)
- `E => (g_d > g_{d+1})` (`N_DEC_STR`)

No-counterexample sufficient implications:

- `(d=mode+1 and deg(w)>=2 and g_d>=g_{d-1} and g_d>g_{d+1}) => E`
  (`S_DG2_INC_DECSTR`)
- `(d=mode+1 and deg(w)>=2 and g_d>g_{d-1} and g_d>g_{d+1}) => E`
  (`S_DG2_INCSTR_DECSTR`)

These are theorem-mining candidates only, not proofs.
