# `(m,lambda)` Overlap Profile (min-u, gated, `m>=4`, `n<=20`)

## Command
- `nice -n 10 python3 scripts/mlambda_overlap_profile_minu.py --min-n 3 --max-n 20 --m-min 4 --max-examples 30 --out results/mlambda_overlap_profile_minu_mge4_n20.json`

## Result summary
From `results/mlambda_overlap_profile_minu_mge4_n20.json`:
- `checked_total = 77059`
- `keys_total = 17316` for key `(m,lambda)`
- `keys_with_multiple_N = 12`
- `adjacent_pairs_total = 12` (i.e., same `(m,lambda)` with `N` and `N+1`)
- `adjacent_pairs_even_d_eq_halfN = 7`
- ratio `7/12 = 0.5833`

## Interpretation
- Adjacent odd/even overlaps are already present for `(m,lambda)` alone.
- The condition `d_even = N_even/2` holds for many adjacent pairs, but not all (`58.3%` here).
- So this condition is not a universal separator for `(m,lambda)` overlaps by itself.
- This does **not** contradict the current `(m,lambda,rho)` agenda:
  - current exact scans still report no `(m,lambda,rho)` split through `n<=26`;
  - the pending `n<=27` adjacent scan remains the decisive computation.

## Representative adjacent `(m,lambda)` pair (known split)
- `m=5`, `lambda=7/9`
- odd: `N=13`, `d=7`, `g6=O??????_A?C?E?@_WG@j?`
- even: `N=14`, `d=7`, `g6=P????A?OD?E?E?B??o?E?OO?`
- here `d_even = N_even/2 = 7`.
