# Mode Rigidity for `m <= 3` in the Min-`u` Canonical Setup

## Setup
Let `T` be a tree in the gated class `d_leaf <= 1`, with independence polynomial
`I(T;x) = sum_{k>=0} i_k x^k`.
Let `m` be the leftmost mode index and define `lambda := i_{m-1}/i_m` (for `m>=1`).

In the accepted min-`u` canonical decomposition, `I(T;x) = (1+2x)P(x) + (1+x)Q(x)` and
`i1 = [x]I = [x]P + 3`, so `N := [x]P = i1 - 3`.

## Lemma 1 (`lambda < 1`)
If `m` is the leftmost mode index, then `0 <= lambda < 1`.

Reason: leftmost mode implies `i_{m-1} < i_m`, and coefficients are nonnegative.

## Lemma 2 (No Unique Admissible Triplet for `d_leaf <= 1`)
Let `(ell,s,u)` be admissible if `ell` is a leaf, `deg(s)=2`, and `N(s)={ell,u}`.
Then admissible triplets are in bijection with degree-1 vertices of `Prune(T)` (delete all leaves).
Hence the admissible triplet count is never exactly 1.

Reason: forests have either 0 leaves (all isolated components) or at least 2 leaves.

## Proposition 3 (Tree Identity for `i2`)
For any tree:
`i2 = (i1-1)(i1-2)/2`.

Reason: independent 2-sets are non-edges:
`i2 = C(n,2) - e`, with `n=i1`, `e=n-1`.

## Corollary 4 (`m=2` Rigidity)
If `m=2`, then `lambda = i1/i2` uniquely determines `i1`, hence determines `N=i1-3`.

Reason: `lambda = 2n/((n-1)(n-2))` with `n=i1`, and this map is strictly decreasing for integers `n>=3`.

## Proposition 5 (Tree Formula for `i3`)
Let `W(T) := sum_v C(deg(v),2)` (wedge count). Then
`i3 = C(n,3) - (n-1)(n-2) + W(T)`, where `n=i1`.

## Lemma 6 (Sharp Bounds on `W` under `d_leaf <= 1`)
For a tree on `n>=2`:
1. `W(T) >= n-2` (equality for path `P_n`).
2. If `d_leaf <= 1`, leaf count `L <= floor(n/2)`.
3. If `d_leaf <= 1`,
   `W(T) <= W_max(n) := C(floor(n/2),2) + (n-floor(n/2)-1)`,
   i.e.
   - odd `n`: `W_max = (n^2-1)/8`
   - even `n`: `W_max = (n^2+2n-8)/8`
   and this bound is attained.

## Theorem 7 (`m=3` Rigidity)
Assume `m=3`. Then `lambda = i2/i3` uniquely determines `n=i1`, hence determines `N=i1-3`.

Proof sketch:
1. `i2 = (n-1)(n-2)/2`.
2. `i3 = (n-1)(n-2)(n-6)/6 + W(T)`.
3. For `n<5`, `i3=0`, so `m=3` is impossible.
4. For `n>=5`, use `n-2 <= W <= W_max(n)` to obtain:
   - `L(n) <= lambda <= U(n)`,
   - `U(n) = 3(n-1)/((n-3)(n-4))`,
   - `L(n) = 3 / ( (n-6) + 6W_max(n)/((n-1)(n-2)) )`.
5. Show `U(n+1) < L(n)` for all `n>=5`, so intervals are pairwise disjoint.
6. Therefore each admissible `lambda` corresponds to a unique integer `n`.

## Corollary 8 (`(m,lambda) -> N` for `m <= 3`)
In the min-`u` canonical setup:
- `m=1`: `lambda = i0/i1 = 1/i1`, so `i1` is fixed.
- `m=2`: by Corollary 4.
- `m=3`: by Theorem 7.

Hence `(m,lambda)` determines `N` for `m<=3`.

## Current Frontier (`m >= 4`)
`(m,lambda) -> N` is false (explicit split exists at `(m,lambda)=(5,7/9)`).
Open direction:
- decide whether `(m,lambda,rho) -> N` or `(m,lambda,sigma) -> N`
  on the gated class with min-`u` canonicalization.

Empirical status:
- no split through `n<=25` for either `(m,lambda,rho)` or `(m,lambda,sigma)`
  (artifact: `results/canonical_projection_battery_minu_mlambda_rho_sigma_mge4_n25_exact.json`).
