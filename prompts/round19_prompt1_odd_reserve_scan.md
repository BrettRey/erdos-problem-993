# Round 19, Prompt 1: Odd-Reserve Scan and Separator Search

**Target model:** Codex 5.3 (exhaustive computation)

---

## Hard facts carried forward

All scans use boundary-correct CB indexing and support-step prefix regime (`k < mode(I_new)`, smallest mode index).

Established through `n <= 18`:

- `X_k < 0` occurs in 135,976 / 2,656,341 instances (5.119%).
- All `X_k < 0` instances occur at `step t = 2`.
- Pairwise termwise routes are dead: `S_{i,j}(k) >= 0` false; `F(i,j) >= 0` false.
- `STP2(I,E)` all roots still has 0 failures through `n<=18`.
- Diagonal-Abel sign changes: each `D_i^{(s)} = W_i^{(s)}-W_{i+1}^{(s)}` has at most one sign change on all tested `X_k<0` cases.

New local gate (`scan_round19_tail_obligations.py`, `n<=18`):

- Even-channel bound held with 0 failures:
  - for even `s=2t`, `Err_s <= Lambda_old(t)*P(k-t)*Q(k-t)`.
- But global closure `sum_s Err_s <= D_k` fails in 85 cases (max ratio `1.071973...`).
- Odd-diagonal error is real:
  - odd `Err_s > 0` in 38,219 / 135,976 negative-X cases.

So the missing piece is an **odd-diagonal reserve**.

Additional empirical seed (full negatives, coarse grid search):

- with shifted channels
  - `C10 = sum_i Lambda_old(i) P(k-i-1)Q(k-i)`
  - `C01 = sum_i Lambda_old(i) P(k-i)Q(k-i-1)`
  - `C11 = sum_i Lambda_old(i) P(k-i-1)Q(k-i-1)`
- the choice `R_shift = C10 + C01 + C11` gave worst-case
  `sum_err / (D + R_shift) = 0.5696` on all 135,976 negative cases (`n<=18`),
  versus baseline worst `sum_err / D = 1.07197`.

Treat this as a strong candidate reserve to validate/refine.

Additional calibrated scalar result:

- minimal global scalar for this reserve through `n<=18`:
  `sum_err <= D + lambda*(C10+C01+C11)` with
  `lambda* = 0.05201381704686925`.

So Round 19 should prioritize proving a small-constant shifted-channel odd reserve.

---

## Task 1: Reproduce the Round 19 obligation diagnostics

Recompute through `n<=18`:

1. Even-channel bound fail count (expect 0).
2. Global fail count for `sum_s Err_s <= D_k` (expect 85).
3. Distribution of odd-error mass:
   - per-case `odd_sum_err / D_k`,
   - count with `odd_sum_err > 0`.
4. Confirm `X_k<0` only at `t=2`.

If mismatched, report exact convention drift.

---

## Task 2: Search candidate odd reserves `R_odd` that close globally

For each negative case, define

- `Err_odd := sum_{s odd} Err_s`
- `Err_even := sum_{s even} Err_s`
- `D_k = sum_i Lambda_old(i) P(k-i)Q(k-i)`.

Test candidate families for `Err_odd <= R_odd` and/or `Err_even + Err_odd <= D_k + R_odd`.

### Family A: shifted diagonal channels

Use channels

`C_{a,b}(k) := sum_i Lambda_old(i) P(k-i-a) Q(k-i-b)`

for small shifts `(a,b) in {0,1}^2` (valid coefficient convention outside range = 0).

Search nonnegative combinations:

`R_odd = alpha*C_{1,0} + beta*C_{0,1} + gamma*C_{1,1}`

on coarse then refined grids for `(alpha,beta,gamma)` minimizing worst-case ratio

`(Err_even+Err_odd)/(D_k + R_odd)`.

Start by confirming the seed point `(alpha,beta,gamma)=(1,1,1)` and then refine around it.
Also test the scalarized form directly to recover/beat `lambda*` above.

### Family B: odd-diagonal local channels

For each odd `s=2t+1`, define local channel candidates from neighboring `t,t+1`, e.g.

- `R_s^min = min(Lambda_old(t)P(k-t)Q(k-t), Lambda_old(t+1)P(k-t-1)Q(k-t-1))`
- `R_s^avg = 0.5*(...)`
- `R_s^sum = (...) + (...)`.

Test whether `Err_s <= R_s^*` holds and aggregate tightness.

---

## Task 3: Step-2 specialization profiling (critical regime)

Since all negatives are at `t=2`, profile by child-size pair `(a,b)` with `a<=b`:

1. counts of negative cases by `(a,b)`,
2. worst `X`, worst `Err_odd/D`, worst `(Err_even+Err_odd)/D`,
3. whether a simple reserve constant by `(a,b)` closes all cases.

Goal: identify whether proof can be reduced to finite parametric family in step-2.

---

## Output requirements

Provide:

1. summary table (core counts + reproduced diagnostics),
2. best odd-reserve candidate(s) with coefficients and worst-case ratios,
3. top-10 witness cases for remaining failures,
4. machine-readable JSON artifact with all aggregates and candidate-fit stats.
