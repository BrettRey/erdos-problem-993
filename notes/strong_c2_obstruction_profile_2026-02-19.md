# STRONG C2 Obstruction Profile (2026-02-19)

This note profiles the only apparent obstruction term in the STRONG C2 split
for canonical degree-2 leaves.

## Setup

At the bridge leaf `(l,s)` with `deg(s)=2`, let:

- `B = T-{l,s}`
- `P = dp_B[u][0]`, `Q = dp_B[u][1]`, `I(B)=P+Q`
- `m = mode(I(T))`
- indices:
  - `p0=p_{m-2}`, `p1=p_{m-1}`
  - `q0=q_{m-2}`, `q1=q_{m-1}`
  - `b0=b_{m-2}`, `b1=b_{m-1}`, `b2=b_m`

Define:

- `mismatch = p0*b1 - p1*b0`
- `neg = -mismatch` (for mismatch-negative cases)
- `rise = b1*(b1-b0)`

## Key identity

For every checked case:

`rise - neg = p1*(b1-b0) + b1*(q1-q0).`

So `neg <= rise` is equivalent to:

`p1*(b1-b0) + b1*(q1-q0) >= 0.`

When `q1<q0` and `b1>b0`, this is the ratio form:

`(- (q1-q0) / (b1-b0)) <= p1/b1.`

Interpretation:

- LHS = normalized `Q`-drop rate near `m-1`,
- RHS = `P`-mass fraction at level `m-1` inside `B`.

## New scanner and artifacts

- script: `conjecture_a_strong_c2_obstruction_scan.py`
- staged artifacts:
  - `results/whnc_strong_c2_obstruction_n22.json`
  - `results/whnc_strong_c2_obstruction_n23.json`
  - `results/whnc_strong_c2_obstruction_staged_summary_n23.json`

## Full frontier summary (`d_leaf<=1`, `n<=23`, 931,596 trees)

From `results/whnc_strong_c2_obstruction_staged_summary_n23.json`:

- `checked = 931,596`
- `mismatch_neg = 129`
- `combined_neg = 0`
- `rise_fail = 0`

Mismatch-negative split by `mode(B)-(m-1)`:

- shift `0`: `91`
- shift `1`: `38`
- other: `0`

`Q`-drop (`q1<q0`) is extremely rare:

- `q_drop_all = 2` (across all checked trees)
- both are mismatch-negative and shift `0`:
  - `q_drop_mismatch_neg = 2`
  - `q_drop_mismatch_neg_shift0 = 2`
  - `q_drop_mismatch_neg_shift1 = 0`

Ratio-form slack on those two `Q`-drop witnesses:

- `max_transfer_ratio_neg_qdrop = max((-dq)/db) = 0.0267707752`
- `min_need_ratio_neg_qdrop = min(p1/b1) = 0.9283765825`
- `min_ratio_gap_neg_qdrop = min( (p1/b1) - ((-dq)/db) ) = 0.9016058073`

So even in the only observed `Q`-drop cases, the ratio condition has a huge margin.

Witness details (both canonical-leaf bridge):

1. `n=20`, `g6 = S???????C?G?G?C?@??G??_?@??@?F~_?`
   - degree signature: `{1:10, 2:9, 10:1}`
   - `dq = -14`, `db = 658`
   - `(-dq)/db = 0.0212766`, `p1/b1 = 0.9846154`

2. `n=23`, `g6 = V???????????_?O?C??_?A??C??C??A???_?{E??^g??`
   - degree signature: `{1:11, 2:10, 6:1, 7:1}`
   - `dq = -48`, `db = 1793`
   - `(-dq)/db = 0.0267708`, `p1/b1 = 0.9283766`

## Immediate reduction

To prove the empirical rise compensation algebraically, it is enough to prove
the ratio form:

`(- (q1-q0) / (b1-b0)) <= p1/b1`

whenever `mismatch<0` and `q1<q0` (the only nontrivial regime in this profile).

Given only 2 such witnesses on the full frontier and large ratio slack, this
looks substantially narrower than the original determinant inequality.

## Commands used

```bash
python3 conjecture_a_strong_c2_obstruction_scan.py \
  --min-n 4 --max-n 22 \
  --out results/whnc_strong_c2_obstruction_n22.json

python3 conjecture_a_strong_c2_obstruction_scan.py \
  --min-n 23 --max-n 23 --mod 8 \
  --out results/whnc_strong_c2_obstruction_n23.json
```
