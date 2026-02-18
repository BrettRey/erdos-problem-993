# Mode-Mean Tiepoint Scan (2026-02-18)

This note records a direct probabilistic scan aimed at proving
`mode <= ceil(mu)` without requiring full log-concavity.

## Setup

Let

`I_T(x) = sum_k i_k x^k`.

For each adjacent pair (`k>=1`, `i_{k-1},i_k>0`), define the tie fugacity

`lambda_k = i_{k-1}/i_k`,

so

`i_{k-1} lambda_k^{k-1} = i_k lambda_k^k`.

Let

`mu(lambda) = E_lambda[|S|]`

for the hard-core tilted distribution `P_lambda(|S|=j) propto i_j lambda^j`.

## Packaged script

- `conjecture_a_tie_mean_scan.py`

Run:

```bash
python3 conjecture_a_tie_mean_scan.py \
  --min-n 3 --max-n 23 \
  --out results/whnc_tie_mean_scan_n23_modefocus.json
```

## Full `n<=23` results (`931,596` d_leaf<=1 trees)

From `results/whnc_tie_mean_scan_n23_modefocus.json`:

- adjacent tie points checked: `11,172,104`,
- all were jump-like (`mode_set={k-1,k}` at `lambda_k`): `11,172,104`,
- global lower tie bound failures:
  - `mu(lambda_k) >= k-1` failed in `2` cases,
  - worst margin `mu(lambda_k)-(k-1) = -0.04470192405278617`
    at `n=22`, `k=12`, `lambda=33`.
- global upper tie bound:
  - `mu(lambda_k) <= k` had `0` failures.
- at `lambda=1`:
  - `mode <= ceil(mu)` had `0` failures,
  - Darroch set-membership (`mode in {floor(mu), ceil(mu)}`) had `0` failures.

## Focused mode-tie condition (promising)

For each tree, let `m = mode(I_T)` at `lambda=1` (leftmost mode in script),
and test only the tie `k=m`:

`mu(lambda_m) >= m-1`, where `lambda_m = i_{m-1}/i_m`.

Result:

- checked: `931,596`,
- failures: `0`,
- minimum margin: `0.36225903819956606`.

Witness for minimum margin:

- `n=21`,
- `g6=T???????C?G?G?C?@??G??_?@??@???_B~o?`,
- `m=7`,
- `lambda_m=0.879383429672447`,
- `mu(lambda_m)=6.362259038199566`.

## Dedicated focused scan (packaged)

To isolate the proof target and track extremals, we added:

- `conjecture_a_mode_tie_focused_scan.py`

Runs:

```bash
# d_leaf<=1 frontier
python3 conjecture_a_mode_tie_focused_scan.py \
  --min-n 3 --max-n 23 \
  --out results/whnc_mode_tie_focused_dleaf_n23.json

# all-tree extension
python3 conjecture_a_mode_tie_focused_scan.py --all-trees \
  --min-n 3 --max-n 22 \
  --out results/whnc_mode_tie_focused_alltrees_n22.json
```

Results:

- `d_leaf<=1` (`n<=23`): checked `931,596`, failures `0`, minimum margin
  `0.36225903819956606`.
- all trees (`n<=22`): checked `9,114,283`, failures `0`, minimum margin
  `0.36225903819956606`.

Both minima are attained by the same witness at `n=21`:

- `g6=T???????C?G?G?C?@??G??_?@??@???_B~o?`,
- degree signature `{1:10, 2:10, 10:1}` (balanced length-2 spider `S(2^10)`),
- `m=7`, `lambda_m=0.879383429672447`,
- `mu(lambda_m)=6.362259038199566`.

Additional per-`n` scan (d_leaf<=1) shows `n=23` minimum is again a balanced
length-2 spider (`{1:11,2:11,11:1}`), while even `n` minima can come from
near-spider mixed signatures.

## Interpretation

The broad tiepoint claim (`mu(lambda_k) >= k-1` for all `k`) is false, so that
route is not viable globally.

But the **single-tie-at-the-actual-mode** condition is robust on the full
frontier and gives a sharper target:

1. prove `mu(lambda_m) >= m-1` for `m = mode(I_T)` at `lambda=1`,
2. use monotonicity of `mu(lambda)` in `lambda` to lift to `mu(1) >= m-1`,
3. conclude `m <= ceil(mu(1))`.

This would bypass full LC while preserving the mode-vs-mean bridge needed
after the Steiner mean bound.
