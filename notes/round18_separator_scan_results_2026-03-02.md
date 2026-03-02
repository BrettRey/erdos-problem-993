# Round 18 Local Separator Scan Results (2026-03-02)

Script: `scan_round18_separator.py`  
Primary output: `results/round18_separator_scan_n18.json`

## Scope and conventions

- Trees exhaustive via `geng` for `n=3..18`.
- Support-step processing:
  - absorb leaf children first (`(1+x)^ell`),
  - then process non-leaf children by increasing rooted-subtree size.
- Prefix regime per step: `k < mode(I_new)`, smallest-index mode.
- Boundary-correct CB indexing includes index `-1`.

## Core totals (n <= 18)

- Trees: 205,002
- Support rootings: 1,103,584
- Steps with `t>=2`: 510,433
- Prefix checks: 2,657,523
- `X_k < 0`: 135,976 (rate 5.1166%)

This matches the Round 17 external report at the level of key aggregate counts.

## Extremals

- Minimum `X_k`:
  - `X = -613,470`
  - `D = 7,913,763`
  - witness: `n=18`, `root=0`, `step=2`, `k=6`,
    graph6 `Q?????????????????C??o?A~~_`

- `D/|X|` on `X<0` cases:
  - minimum: `85/22 = 3.8636`
  - median: `163/8 = 20.375`

## Candidate invariant family checks

### 1) Negative gap-mass budget

Define:
- `neg_gap = sum_g max(0, -gap_sum(g))`

Observed:
- `neg_gap <= D` in 2,656,819 / 2,657,523 checks (99.9735%).
- So `neg_gap <= D` is **almost** true but has rare failures.
- Worst observed ratio:
  - `neg_gap / D = 24017/18600 = 1.2912`
  - witness: `n=18`, `root=16`, `step=2`, `k=7`.

Conclusion: `neg_gap <= D` is not a universal invariant.

### 2) Diagonal-Abel difference sign-change structure

For each diagonal `s`, define
`D_i^(s) = W_i^(s) - W_{i+1}^(s)` with
`W_i^(s) = P(k-i)Q(k-s+i)`.

On all `X<0` cases through `n<=18`:
- diagonal sequences checked: 1,311,324
- diagonals with `>1` sign change: 0
- maximum sign changes seen: 1
- `X<0` cases with any diagonal `>1` sign changes: 0

Conclusion: the one-sign-change property is a strong surviving empirical invariant.

### 3) `w_m(I,E) >= 0` (adjacent minor) [separate local gate]

This candidate was tested separately and is **false** on trees, including in prefix-only checks.
Earliest witness:
- `n=7`, graph6 `F??Fw`, root `r=6`, `m=2`
- `I=[1,7,15,20,15,6,1]`, `E=[1,6,15,20,15,6,1]`
- `w_2 = -15`

So this route is dead.

## Practical implication for Round 18 theory prompts

The strongest new empirical lever is:

1. `X_k` can be negative, but
2. `D/|X|` stays comfortably > 1 on all scanned negatives, and
3. diagonal-Abel difference sequences exhibit at most one sign change.

This suggests focusing proofs on diagonal-Abel aggregation with a unimodal/single-crossing
kernel argument, not on pairwise term positivity or adjacent-minor positivity.

