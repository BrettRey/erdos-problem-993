# Sub-claim A structure scan (2026-02-18)

Target: prove for mixed spiders `S(2^k,1^j)`:

`margin(k,j+2) >= margin(k,j)` for all `k>=6`, `j>=0`.

## New scan artifact

- Script: `analyze_subclaim_a_structure.py`
- Output: `results/subclaim_a_structure_k1200_j120.json`
- Run: `python3 analyze_subclaim_a_structure.py --k-max 1200 --j-max 120 --out results/subclaim_a_structure_k1200_j120.json`

Checked pairs: `144,595`.

## Empirical lemmas with 0 failures in this range

1. **Mode step is rigid**:
   `mode(k,j+2) - mode(k,j) = 1`.

2. **Tie fugacity is nondecreasing on parity subsequences**:
   `lambda(k,j+2) - lambda(k,j) >= 0`.

3. **Sub-claim A itself**:
   `margin(k,j+2) - margin(k,j) >= 0`.

Minima from scan:

- min `lambda` step: `1.699e-6` at `(k,j)=(1200,119)`.
- min margin step: `4.51668e-5` at `(k,j)=(8,120)`.

## Mean-step decomposition

Let `mu_j(lambda) := mu(S(2^k,1^j), lambda)` and `lambda_j := lambda(k,j)`.

Then

`mu_{j+2}(lambda_{j+2}) - mu_j(lambda_j)`
`= [mu_{j+2}(lambda_j) - mu_j(lambda_j)]`
`+ [mu_{j+2}(lambda_{j+2}) - mu_{j+2}(lambda_j)]`.

The scan records minima of these two pieces:

- min fixed-`lambda_j` step: `0.7709778` at `(k,j)=(8,0)`.
- min `lambda`-gain piece: `0.000224698` at `(k,j)=(6,119)`.
- min total mean step: `1.00004517` at `(k,j)=(8,120)`.

Since `mode(j+2)-mode(j)=1`, this is equivalent to a positive margin step.

## Useful exact algebra identity at fixed lambda

Write

- `a = 2k*lambda/(1+2lambda)`,
- `u = lambda/(1+lambda)`,
- `b = 1 + k*lambda/(1+lambda)`,
- `r_j = lambda*(1+lambda)^(k-j)/(1+2lambda)^k`,
- `s = (1+lambda)^2`,
- `Y_j = (a + j*u) - b`.

Then

`mu_j(lambda) = (a + j*u + r_j*b)/(1+r_j)` and

`mu_{j+2}(lambda)-mu_j(lambda)`
`= [Y_j*r_j*(s-1) + 2*s*u*(r_j+1)] / ((1+r_j)*(r_j+s))`.

This identity may be useful for an algebraic lower bound on the fixed-`lambda` piece.
