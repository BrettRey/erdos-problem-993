# Round 19, Prompt 2: Corrected Diagonal-Abel Theorem (with Odd Residual)

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context and constraints

We now have a refined empirical picture:

- Even-diagonal channel bound appears exact/robust in data (`n<=18`): no failures.
- Odd-diagonal residual is nonzero in many cases.
- Global `sum_s Err_s <= D_k` is close but not always true (small finite failure set).
- Empirical shifted-channel reserve seed:
  `R_shift = C10 + C01 + C11` (with `C10,C01,C11` defined from shifted
  `Lambda_old` channels) appears to dominate full tail error with wide margin
  on exhaustive `n<=18` negatives.

So the prior heuristic “odd diagonals free” is false.

Do not attempt dead routes (`X_k>=0`, pairwise `S>=0`, `F>=0`, abstract `w_m(I,E)>=0`).

---

## Setup

For one update step:

- `A=E_acc`, `B=J_acc`, `P=I_c`, `Q=E_c`
- `Lambda_k^{new} = sum_{i,j} Delta_{i,j}(A,B) P(k-i)Q(k-j)`
- `D_k = sum_i Lambda_old(i) P(k-i)Q(k-i)`
- `X_k = Lambda_k^{new} - D_k`

Diagonal decomposition (`s=i+j`) and Abel form:

- `u_i^{(s)} = A(i)B(s-i)`
- `W_i^{(s)} = P(k-i)Q(k-s+i)`
- `D_i^{(s)} = W_i^{(s)}-W_{i+1}^{(s)} = Delta_{k-i,k-s+i}(P,Q)`
- `S_k^{(s)} = sum_i u_i^{(s)} D_i^{(s)}`
- `Lambda_k^{new} = sum_s S_k^{(s)}`

---

## Task 1: Prove the structural pieces cleanly

1. Exact diagonal-Abel identity with finite-support boundary handling.
2. Complementary half-diagonal sign facts from STP2:
   - `u_i-u_{i-1} >= 0` on left half,
   - `D_i >= 0` on right half.
3. If possible, derive one-sign-change for `D_i` from LC of `P,Q` (or state exact needed assumptions).

---

## Task 2: Correct the reduced obligation

The previous reduced claim “odd Err=0” is false.

Derive a corrected theorem template:

- even diagonals: prove/justify a bound of type
  `Err_{2t} <= Channel_t` (where `Channel_t` is an explicit nonnegative IH channel);
- odd diagonals: isolate residual term `R_{2t+1}` and prove
  `Err_{2t+1} <= R_{2t+1}` with explicit formula.

Then reduce global closure to

`sum_t Channel_t + sum_t R_{2t+1} <= D_k + Extra_k`

with `Extra_k` explicit and ideally already nonnegative from IH-level invariants.

---

## Task 3: Propose one minimal extra hypothesis H_odd

If full closure does not follow from current hypotheses, state one minimal, testable tree-specific hypothesis `H_odd` that controls odd residuals.

Requirements:

- `H_odd` must be local/diagonal (not a global black-box claim).
- It must be computationally testable on finite scans.
- Show exactly how `H_odd` plugs into the corrected theorem above.

No vague “assume better alignment”; write the inequality explicitly.

---

## Task 4: Deliverable format

Provide either:

1. full proof skeleton with clear theorem/lemma/conjecture separation, or
2. reduced obligation list with exactly one missing lemma (`H_odd`), plus a precise computational test recipe.
