# STRONG C2 Route A: Q-drop vs shift-0 (2026-02-19)

## Goal

Investigate whether, in the canonical degree-2 bridge setup for `d_leaf <= 1` trees,

`q_{m-1} < q_{m-2}` (Q-drop at `m-1`)

forces

`mode(B) = m-1` (shift-0).

Here:
- `m = mode(I(T))` (leftmost mode),
- `B = T - {leaf, support}`,
- `P = dp_B[u][0]`,
- `Q = dp_B[u][1] = x * R`,
- `R = product(dp0[child])` at hub `u` in `B`.

---

## New verification scripts

- `verify_strong_c2_qdrop_shift0_2026_02_19.py`
  - Exhaustive scan with witness dump.
  - Verifies Q-index identity `Q = xR` at witness indices.
  - Tracks whether any Q-drop has nonzero shift.

- `analyze_strong_c2_qdrop_witnesses_2026_02_19.py`
  - Deep analysis of Q-drop witnesses.
  - Decomposes
    `P = product(F_i + G_i) = sum_{t=0}^d H_t`
    by number `t` of `G_i=dp1[child_i]` factors.
  - Quantifies which `t`-layers drive growth near `m-1`.

Artifacts:
- `results/verify_strong_c2_qdrop_shift0_2026_02_19.json`
- `results/analyze_strong_c2_qdrop_witnesses_2026_02_19.json`

---

## Main frontier result (`d_leaf <= 1`, canonical leaf, `n <= 23`)

From `results/verify_strong_c2_qdrop_shift0_2026_02_19.json`:

- `seen = 23,942,356`
- `considered = 931,596`
- `checked = 931,596`
- `qdrop_total = 2`
- `qdrop_shift0 = 2`
- `qdrop_shift1 = 0`
- `qdrop_shift_other = 0`
- `qdrop_shift_nonzero = 0`

So on the full checked frontier:

`Q-drop => shift-0` holds with `0` exceptions.

Also verified at witness indices:
- `qdrop_q_index_identity_fail = 0` for `q_{m-2}=R_{m-3}`, `q_{m-1}=R_{m-2}`.

---

## Witnesses and local data

The only Q-drop witnesses are:

1. `n=20`, `g6 = S???????C?G?G?C?@??G??_?@??@?F~_?`
2. `n=23`, `g6 = V???????????_?O?C??_?A??C??C??A???_?{E??^g??`

For both:
- `mode(B)=m-1` (shift-0),
- `mode(R)=m-3`,
- therefore `mode(Q)=m-2`,
- `mode(P)=m-1`.

Numerically:

- `n=20`: `dq=-14`, `db=658`, `dp=672`, `p1/b1=0.984615`, `(-dq/db)=0.021277`.
- `n=23`: `dq=-48`, `db=1793`, `dp=1841`, `p1/b1=0.928377`, `(-dq/db)=0.026771`.

So even in Q-drop:
- the negative `Q` slope is tiny,
- `P` mass fraction is very large,
- and the ratio gap is huge (`~0.96`, `~0.90`).

---

## What drives Q-drop in these two trees?

From `results/analyze_strong_c2_qdrop_witnesses_2026_02_19.json`:

- `n=20` witness: hub child subtree sizes are `[2,2,2,2,2,2,2,2,1]`.
  - This gives `R = (1+x)^8`, `Q = x(1+x)^8`.
  - Q-drop is exactly binomial right-tail descent at the checked index.

- `n=23` witness: child sizes are `[2,2,2,13,1]`.
  - One larger branch plus small arms; still `R` has mode at `m-3`.

In both witnesses, at the same window (`k=m-2 -> m-1`):
- `delta_R < 0` (the excluded product decreases),
- `delta_P > 0` (the full product increases),
- hence `delta_S = delta(P-R)` is strongly positive.

Values:
- `n=20`: `delta_R=-28`, `delta_P=+672`, `delta_S=+700`.
- `n=23`: `delta_R=-404`, `delta_P=+1841`, `delta_S=+2245`.

So the `dp1`-enabled part `S=P-R` more than compensates:
- compensation factor `delta_S/|delta_R| = 25.0` (n=20),
- `= 5.56` (n=23).

---

## dp1-contribution structure (key point for Question 3)

Using `P = sum_t H_t` where `H_t` uses exactly `t` child `dp1` factors:

At the same window (`m-2 -> m-1`):

- `n=20`:
  - `H_0 + H_1 = -154` (non-positive block),
  - `sum_{t>=2} H_t = +826`,
  - total `delta_P = +672`.

- `n=23`:
  - `H_0 + H_1 = -936`,
  - `sum_{t>=2} H_t = +2777`,
  - total `delta_P = +1841`.

Interpretation:
- The rise of `P` is not coming from `R=H_0` (it is decreasing).
- It is not coming from only single-`dp1` terms either (`H_1` is net negative here).
- The decisive effect is **multi-child (`t>=2`) dp1 interaction**, which shifts mass to higher indices strongly enough to overwhelm the `R` drop.

This is the concrete mechanism behind:
- `R` decreasing while `P` increases,
- and why `B=P+Q` still peaks at `m-1` in these witnesses.

---

## Is this structural or coincidence?

### Within the target regime (`d_leaf <= 1`, canonical bridge)

Empirically on the full checked frontier `n<=23`, it is structural:

`Q-drop => shift-0` with `0` exceptions.

### Outside the target regime (sanity check)

Running the same script on **all trees** up to `n<=20` gives:
- `qdrop_total = 20005`
- `qdrop_shift1 = 2442`

So `Q-drop => shift-0` is **not** universally true for arbitrary trees.
Any proof must exploit the `d_leaf<=1` bridge structure, not just generic tree DP algebra.

---

## Commands used

```bash
python3 verify_strong_c2_qdrop_shift0_2026_02_19.py \
  --min-n 4 --max-n 23 \
  --out results/verify_strong_c2_qdrop_shift0_2026_02_19.json

python3 analyze_strong_c2_qdrop_witnesses_2026_02_19.py \
  --out results/analyze_strong_c2_qdrop_witnesses_2026_02_19.json

# sanity check outside d_leaf<=1 (not part of main claim)
python3 verify_strong_c2_qdrop_shift0_2026_02_19.py \
  --all-trees --min-n 4 --max-n 20 \
  --out /tmp/verify_qdrop_shift0_alltrees_n20.json
```

---

## Practical conclusion for Route A

For the intended class, Q-drop appears to be confined to configurations where:
- `Q` is already one step left (`mode(Q)=m-2`),
- `P` is at `m-1`,
- and multi-child `dp1` terms strongly dominate the local slope.

This is exactly the pattern needed for shift-0 closure in the hard lane.
