# Round 11 Triage (2026-03-01)

**Models:** Codex 5.3 (1 instance), GPT 5.2 Pro (2 instances)

## Codex P1: Diagonal/Cross Verification + (Agg) + STP2

**CRITICAL CORRECTIONS:**

### Diagonal/Cross decomposition REFUTED
- D_k < 0 at n=7 (4684 cases = 7.55% of checks in n≤14)
- X_k < 0 at n=11 (28 cases = 0.045%)
- w_k itself always ≥ 0 (the SUM is fine, individual parts are not)
- Signs are ORDER-DEPENDENT on child processing order
- Prior agent result (0 failures) was wrong: used k < mode(I_T) instead of k < mode(I_new)
- **Do NOT use** diagonal/cross as a proof target

### (Agg) condition DEAD
- r_{k-i}·J_{k-1} ≤ r_{k-1-i}·J_k fails ~47%

### STP2: PROMISING (0 failures)
- r_{m+1}·q_n ≤ r_m·q_{n+1} for m > n ≥ 0
- 0 failures across 50,917 factor pairs
- Most promising new algebraic condition

### Additional findings
- At tightest w_k: D ≤ X always (cross term dominates at critical k)
- Child ordering: D_k, X_k signs depend on processing order; w_k is invariant

## GPT P2: s=2 Binomial Smoothing Case

**KEY CORRECTIONS:**

### Polynomial factor cancellation NOT valid
- f ≽ g does NOT imply H·f ≽ H·g (only Karlin direction holds)
- The "E_1 cancels" simplification in Round 11 Prompt 2 was only sufficient, not equivalent

### w_k[f^(1)] CAN be negative for k < mode
- Cancellation between w_k[f^(0)] and w_k[f^(1)] is necessary
- Cannot prove each part ≥ 0 separately

### E≽J FAILS GLOBALLY at n=32 (CONFIRMED)
- Construction: support vertex r, ℓ=1, child 1 = P_2, child 2 = T_{3,4} broom
- n = 1 + 1 + 2 + 28 = 32
- d_15 = -3498 (our verification; GPT said -4482, minor construction difference)
- Mode = 10, so prefix E≽J (k < mode) still holds
- E is NOT LC at this tree; J IS LC

### Clean x-shift identity
- w_k[f^(1)] = H_k·C_k - H_{k-1}·C_{k+1} (shifted minor)

### Expectation-disintegration of correction
- d_{k-1}(B,C)/(C_k·C_{k-1}) = difference of conditional expectations of ratio profile ρ under weights μ_k(i)
- Potentially useful for bounding, not yet exploited

## GPT P3: Cross-Term Algebraic Proof

**KEY INSIGHT — THE CORRECT SPLIT:**

### Karlin split (CORRECT decomposition)
```
w_k = d_k(E_acc·g, J_acc·g) + d_k(x·E_acc·h, J_acc·g)
```
where g = E_t, h = J_t (so f = I_t = g + xh).

- **First term (Karlin): PROVABLY ≥ 0** by Karlin/TP2 convolution theorem
  - E_acc ≽ J_acc (by induction), g is LC, so g is PF2
  - Karlin: PF2 kernel preserves ≽ under convolution
  - Therefore d_k(E_acc·g, J_acc·g) ≥ 0 for all k

- **Second term (correction): sign-indefinite**
  - d_k(x·E_acc·h, J_acc·g) can be negative
  - Must be bounded by the Karlin term

### Explicit counterexamples for D/X
- D_k < 0 at n=10 (GPT's explicit example)
- X_k < 0 at n=12 (GPT's explicit example)
- Confirms Codex P1's findings independently

### "Pure-g convolution minor" is the right positive piece
- d_k(E_acc·g, J_acc·g) ≥ 0 ALWAYS (standard antisymmetric CB with PF2 kernel)
- This is NOT the same as D_k (diagonal of the old decomposition)
- The correct split: Karlin + correction, NOT diagonal + cross

## Strategic Assessment

### Confirmed alive
1. **STP2**: 0 failures across 50K+ factor pairs. Would resolve the mixed bracket obstruction if true.
2. **Karlin split**: First term provably ≥ 0. Correction needs bounding.
3. **Prefix E≽J**: Still holds at n=32 (failure only at k=15 > mode=10). Proof target should be prefix, not universal.
4. **W-form**: Still verified 0 failures. s=1 PROVED. s≥2 open.

### Confirmed dead
1. **Diagonal/Cross both ≥ 0**: FALSE. Do not use.
2. **(Agg) condition**: ~47% failure. Dead.
3. **Universal E≽J**: FALSE at n=32. Restrict to prefix.
4. **Polynomial factor cancellation**: f≽g ⇏ Hf≽Hg. Only Karlin direction valid.
5. **w_k[f^(0)] ≥ 0 and w_k[f^(1)] ≥ 0 separately**: FALSE for f^(1). Must use cancellation.

### Priorities for Round 12
1. **Test STP2 more thoroughly** (larger n, more factor pairs, look for failures)
2. **Bound the correction term** d_k(x·E_acc·h, J_acc·g) using STP2 or other conditions
3. **Verify prefix E≽J exhaustively** at n=23-27 (is it truly prefix-universal?)
4. **Exploit expectation-disintegration** for analytic bound on correction
5. **s=2 explicit**: balanced double-star ratio → 3/2, use this for concrete bound

### Open questions
- Does prefix E≽J hold at ALL support vertices of ALL trees? (verified n≤22, confirmed at n=32)
- Does STP2 hold for all tree-realizable factor pairs? (0 failures so far)
- Can the correction term be bounded using the x-divisibility of (q-r)?
