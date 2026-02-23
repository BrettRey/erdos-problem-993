# STRONG C2 Route C: Structural Analysis of `p_{m-1} >= q_{m-1}`

Date: 2026-02-19

## 1. Setup

For a `d_leaf <= 1` tree `T`, choose a leaf `l` whose support `s` has degree `2`, and let `u` be the other neighbor of `s`.
Set `B = T - {l,s}` and root `B` at `u`.

Let the children of `u` in `B` be `c_1, ..., c_d`.
Define

- `f_i := I(T_{c_i}) = dp0[c_i] + dp1[c_i]`
- `g_i := dp0[c_i]`
- `h_i := dp1[c_i]`

Then

- `P = prod_i f_i`
- `P' = prod_i g_i`
- `Q = x P'`
- `I(B) = P + Q`

and at coefficients:

- `p_k = [x^k] P`
- `q_k = [x^k] Q = [x^{k-1}] P'`

Let `m = mode(I(T))` (leftmost mode).
Target: prove `p_{m-1} >= q_{m-1}`.

---

## 2. Product-expansion decomposition (exact)

Because `P = prod_i (g_i + h_i)`:

`P = P' + E`, where

`E = sum_{empty != S subseteq [d]} (prod_{i in S} h_i)(prod_{j notin S} g_j)`.

Hence `E_k >= 0` for all `k`.
At index `m-1`:

`p_{m-1} - q_{m-1}`
`= P_{m-1} - P'_{m-2}`
`= (P'_{m-1} - P'_{m-2}) + E_{m-1}`.

This is the standard Route C split:

- ascent/descent term of `P'` at `m-1`, plus
- nonnegative excess from `dp1`-terms.

---

## 3. Check of the proposed shortcut `mode(P') >= m-1`

The shortcut would be:

1. `P >= P'` coefficient-wise, so `P_{m-1} >= P'_{m-1}`.
2. If `P'_{m-1} >= P'_{m-2}`, then `P_{m-1} >= P'_{m-2} = q_{m-1}`.

So the key subclaim is `mode(P') >= m-1`.

This subclaim is **false**.

Counterexample (first occurs at `n=8`):

- graph6: ``G?`@F_``
- `m=3`
- `P' = [1,2,1]`, so `P'_{m-1}=P'_2=1 < P'_1=2=P'_{m-2}` (descending)
- yet `P=[1,5,8,4]`, `Q=[0,1,2,1]`, so `p_{m-1}=8 >= 2=q_{m-1}`.

So Route C cannot be closed by monotonicity of `P'` alone.

---

## 4. Correct structural reduction

A better reduction is through `P` itself.

### Lemma (exact reduction)
If `P_{m-1} >= P_{m-2}` (equivalently `mode(P) >= m-1`), then `p_{m-1} >= q_{m-1}`.

### Proof

`p_{m-1} - q_{m-1}`
`= P_{m-1} - P'_{m-2}`
`= (P_{m-1} - P_{m-2}) + (P_{m-2} - P'_{m-2})`.

The second term is nonnegative because `P = P' + E` with `E >= 0` coefficient-wise, so
`P_{m-2} >= P'_{m-2}`.

If additionally `P_{m-1} >= P_{m-2}`, then both terms are nonnegative, hence
`p_{m-1} - q_{m-1} >= 0`.

QED.

So Route C reduces to one structural lemma:

`mode(P) >= m-1`.

---

## 5. Full-frontier computation (`n<=23`) supporting the reduction

Command run:

```bash
python3 verify_route_c_p_dominance_2026_02_19.py --min-n 4 --max-n 23
```

Observed:

- Total checked: `931,596`
- `p_{m-1} >= q_{m-1}`: `931,596 / 931,596`
- `mode(P) >= m-1`: `931,596 / 931,596`
- `mode(P') >= m-1`: `835,396 / 931,596` (`89.7%`)
- `P'` descending at `m-1`: `95,362` (`10.2%`)
- All descending cases compensated by excess `E`: `95,362 / 95,362`
- First-order excess `E1 = sum_i h_i prod_{j!=i} g_j` also compensates all descents:
  `95,362 / 95,362`
- Minimum ratios on descending cases:
  - `E_{m-1} / descent = 5.0`
  - `E1_{m-1} / descent = 4.0`

Thus empirically the correct invariant is the `P`-mode condition, not the `P'`-mode condition.

---

## 6. Route C status

What is now structurally proved:

1. Exact expansion formula for `p_{m-1} - q_{m-1}`.
2. The proposed `mode(P') >= m-1` route is invalid (explicit counterexample).
3. A corrected reduction: `mode(P) >= m-1` is sufficient to prove `p_{m-1} >= q_{m-1}`.

What remains for a fully symbolic closure of Route C:

- A structural proof that `mode(P) >= m-1` for the canonical degree-2 bridge decomposition in `d_leaf <= 1` trees.

Current computational evidence for this remaining lemma is complete through `n<=23` on the checked frontier.
