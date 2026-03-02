# Round 18, Prompt 2: Diagonal-Abel Route to `D_k + X_k >= 0`

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Hard constraints from Round 17

Please do **not** attempt to prove any of the following (they are false):

- `X_k >= 0` termwise.
- Pairwise symmetric nonnegativity `S_{i,j}(k) >= 0`.
- Symmetrized bracket nonnegativity `F(i,j) >= 0`.
- Deriving `STP2(I,E)` from only abstract `STP2(E,J) + LC + J<=E`.

The target is only:

`D_k + X_k = Lambda_k^{new} >= 0`.

---

## Setup

At one support-step update:

- `A = E_acc`, `B = J_acc`
- `P = I_c`, `Q = E_c`
- `Lambda_k^{new} = sum_{i,j} Delta_{i,j}(A,B) P(k-i)Q(k-j)`
- `D_k = sum_i Lambda_i^{old} P(k-i)Q(k-i)`
- `X_k = sum_{i!=j} Delta_{i,j}(A,B) P(k-i)Q(k-j)`

Use boundary-inclusive indexing (`...,-1,0,1,...`) so identities are exact.

---

## Task 1: Write a clean diagonal decomposition and Abel transform

For each diagonal `s=i+j`, define:

- `u_i^{(s)} := A(i)B(s-i)`
- `W_i^{(s)} := P(k-i)Q(k-s+i)`
- `Delta_{i,s-i}(A,B) = u_i^{(s)} - u_{i-1}^{(s)}`
- `D_i^{(s)} := W_i^{(s)} - W_{i+1}^{(s)} = Delta_{k-i,k-s+i}(P,Q)`

Derive rigorously:

`S_k^{(s)} := sum_i Delta_{i,s-i}(A,B)W_i^{(s)} = sum_i u_i^{(s)}D_i^{(s)}`

and

`Lambda_k^{new} = sum_s S_k^{(s)}`

(with precise finite-support boundary handling).

---

## Task 2: State the exact complementary STP2 geometry on each diagonal

Prove the two half-diagonal sign facts:

1. From `STP2(A,B)`: `u_i^{(s)} - u_{i-1}^{(s)} >= 0` on one half (`i <= s/2`, with exact integer convention).
2. From `STP2(P,Q)`: `D_i^{(s)} >= 0` on the opposite half (`i >= s/2`).

Make index inequalities explicit; no informal quadrant language.

---

## Task 3: Attempt a *provable* lower bound for each diagonal sum

Goal: derive a usable lower bound of the form

`S_k^{(s)} >= -Err_k^{(s)}`

where `Err` is explicit and summable, and then show `sum_s Err_k^{(s)} <= D_k` (or another nonnegative reserve term from IH).

You may introduce one additional tree-specific hypothesis **only if clearly stated** and only if it is computationally testable (e.g. adjacent-minor positivity for `(I,E)`).

Do not assert unproved inequalities. If a step blocks, isolate the exact inequality and provide a minimal checkable conjecture.

---

## Task 4: Bridge to a proof strategy (or precise blocker)

Deliver one of:

1. A complete argument from accepted hypotheses to `Lambda_k^{new} >= 0`, or
2. A reduced proof obligation of the form:
   - finite family of inequalities on diagonal profiles, or
   - one explicit extra invariant `H` such that
     `STP2 + LC + H => D_k + X_k >= 0`.

The output must clearly separate theorem, lemma, conjecture, and heuristic.

