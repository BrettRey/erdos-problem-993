# Mode-Tie Leaf Bridge Progress (2026-02-18)

We tested the decomposition identity route for the focused tie inequality
`mu(lambda_m) >= m-1` (`m = mode(I_T)` at `lambda=1`, `lambda_m=i_{m-1}/i_m`).

## Bridge functional

For a leaf `l` with support `s`, define:

- `A = T-l`
- `B = T-{l,s}`
- `Phi_q(F;lambda) = lambda I'_F(lambda) - (q-1)I_F(lambda)`.

At `lambda=lambda_m(T)`:

`Phi_m(T) = Phi_m(A) + lambda * Phi_{m-1}(B)`.

So a sufficient local condition is:

`Phi_m(A) >= 0` and `Phi_{m-1}(B) >= 0`.

## Staged frontier results (d_leaf<=1)

Using `conjecture_a_mode_tie_leaf_bridge_scan.py` in staged passes:

- `results/whnc_mode_tie_leaf_bridge_stage2_dleaf_n20.json`
- `results/whnc_mode_tie_leaf_bridge_stage3_dleaf_n21.json`
- `results/whnc_mode_tie_leaf_bridge_stage3_dleaf_n22.json`
- `results/whnc_mode_tie_leaf_bridge_stage3_dleaf_n23.json`
- aggregate: `results/whnc_mode_tie_leaf_bridge_dleaf_n23_staged_summary.json`

Aggregate (`n<=23`, `931,596` trees):

- **existence failures** (`exists leaf with both nonnegative`): `0`
- `all_leaves` failures (not all leaves good): `4,754`

Heuristic leaf choices:

- `min_parent_deg`: `0` failures
- `parent_deg2_else_max`: `0` failures
- `first_leaf`: `342` failures
- `max_parent_deg`: `4,726` failures

So the robust choice is: **pick a leaf whose support has minimum degree**.

## General-tree contrast

On all trees (`n<=16`), existence already fails (`2` failures):

- `results/whnc_mode_tie_leaf_bridge_stage2_all_n16.json`

Hence this bridge appears specific to `d_leaf<=1`, matching the conjecture class.

## Degree-2 support local test

New scanner: `conjecture_a_mode_tie_deg2_support_scan.py`.

Hypothesis tested:

For every leaf `l` whose support `s` has `deg(s)=2`,

`Phi_m(T-l;lambda_m(T)) >= 0` and `Phi_{m-1}(T-{l,s};lambda_m(T)) >= 0`.

Checks performed (all `d_leaf<=1`):

- `n<=20`: `334,140` degree-2-support leaf checks, `0` failures
- `n=21`: `458,915` checks, `0` failures
- `n=22`: `1,100,831` checks, `0` failures
- `n=23`: `2,649,484` checks, `0` failures

Total: `4,543,370` checks, `0` failures.

## Structural mode profile for chosen leaf (min support degree)

One-leaf-per-tree scan through full `n<=23` (`931,596` trees):

- `Phi` failures: `0`
- `mode(A)-m in {-1,0}` always
- `mode(B)-(m-1) in {0,1}` for tested range (`n<=20` exact profile).

## High-ROI proof target

A clean two-lemma route now looks plausible:

1. **Existence lemma (combinatorial):** Every `d_leaf<=1` tree has a support vertex
   of degree `2`.

   Quick counting proof sketch:
   let `S` be supports, `L` leaves (`|L|=|S|`), `C=V\(L\cup S)`.
   If all supports had degree `>=3`, then
   `2(n-1) = sum deg >= |L| + 3|S| + 2|C| = 2n`, contradiction.

2. **Local bridge lemma (open):** If `l` is a leaf with `deg(s)=2`, then at
   `lambda=lambda_m(T)`, both bridge terms are nonnegative.

If (2) is proved, (1) gives a canonical good leaf, and the focused tie route
becomes a concrete induction candidate.
