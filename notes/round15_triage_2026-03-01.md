# Round 15 Triage (2026-03-01)

**Models:** Codex 5.3 (1 instance), GPT 5.2 Pro (2 instances)

## Instance 1 (Codex): Step-Ratio Monotonicity & Extensions

### Step-ratio monotonicity: FALSE
- **22 finite-finite violations** at n=15-17 where step t+1 ratio < step t ratio.
- All violations are "soft" (ratios remain well above 1).
- Step 2 remains the **global minimum** across all n values.
- Step 2 min ratio always comes from **s=2** SVs (regardless of total s at the SV).
- **Conclusion:** "step t ratio monotonically nondecreasing" is false as a strict claim, but step 2 is still always hardest.

### Step-2 ratio independence of total s: CONFIRMED
- Min step-2 ratio is the same whether the SV has s=2, s=3, or s=4 non-leaf children.
- Step 2 depends only on (ℓ, child_1, child_2), not on children 3+.

### LC gap vs ratio: WEAK correlation
- Normalized LC gap of J_acc increases at later steps (products of more LC factors).
- Correlation with ratio is positive but weak. Not a clean explanatory variable.

### Ordering independence at step 2: CONFIRMED for global min
- Ordering affects individual step-2 ratios but NOT the global minimum.
- All orderings give the same global min at step 2.

### σ/M(F) matrix: VERIFIED computationally
- M(F·G) = M(F)·M(G) confirmed for all tree-realizable factors n≤12.
- STP2 diagonal condition corresponds to: σ(E)(n)·J(n+1) - E(n)·σ(J)(n+1) ≥ 0.
- Matches the 2×2 minor condition on the combined matrix [E-part | J-part].

### Chain STP2 extended to n=22: 0 failures
- Exhaustive n=18 (123K trees): 0 failures.
- Sampled n=19-22 (10K each): 0 failures.
- Diagonal form (Lemma A) agrees with full STP2 in all cases (as expected under LC).

## Instance 2 (GPT): σ/M(F) → TN2 Product Closure

### Full TN2 of concatenated Φ: FALSE
- **Counterexample**: rooted path of length 3 (root→v1→v2→v3).
  E=(1,3,3,1), J=(1,2,1). σ(E)=(3,3,1), σ(J)=(2,1).
  Φ = [3,3,1,0,...,2,1,0,...; 1,3,3,1,...,1,2,1,...]
  det[col_0(E), col_0(J)] = 3·1 - 1·2 = 1 ≥ 0 ✓
  det[col_1(E), col_0(J)] = 3·1 - 3·2 = -3 < 0 ✗
- Full TN2 requires ALL 2×2 minors ≥ 0; this fails.
- STP2 only needs **specific cross-minors** (col_n(E) × col_{n+1}(J)), not all.

### Correct target: "Ladder-minor" condition (partial TN2)
- The STP2(E,J) + LC(E) + LC(J) conditions use only:
  - E-block: consecutive minors (= LC gaps of E) ✓
  - J-block: consecutive minors (= LC gaps of J) ✓
  - Cross: det[col_n(E), col_{n+1}(J)] ≥ 0 (= STP2 diagonal) — SPECIFIC pairs only
- This is weaker than full TN2. The "specified minors" or "ladder-minor" framework from algebraic combinatorics may apply.

### 4×4 update matrix N_c (KEY NEW IDEA)
- Augmented state vector: (σ(E_acc), E_acc, σ(J_acc), J_acc) as polynomial entries.
- When processing child c with (I_c, E_c, J_c):
  ```
  σ(E_acc · I_c) = σ(E_acc)·I_c + E_acc·σ(I_c) = σ(E_acc)·I_c + E_acc·(σ(E_c)+J_c)
  E_acc · I_c = E_acc · I_c
  σ(J_acc · E_c) = σ(J_acc)·E_c + J_acc·σ(E_c)
  J_acc · E_c = J_acc · E_c
  ```
- This gives a 4×4 matrix N_c (entries are polynomials) such that state_new = N_c · state_old.
- If the "ladder-minor" condition on the state vector's coefficient matrix is preserved under left-multiplication by N_c, and N_c satisfies appropriate positivity, closure follows.
- **Karlin's TN product theorem**: if N_c's Toeplitz embedding is TN, and state is TN, product is TN. But we need only partial TN (ladder minors), not full TN.

### σ(I) = σ(E) + J: CONFIRMED as key structural identity
- This is the tree-specific relation linking E, J, I in the σ framework.
- M(I) = M(E) + [xJ, 0; J, 0] — rank-1 additive correction.

### Planar network interpretation (SPECULATIVE)
- TN matrices arise from planar networks (Lindström's lemma).
- If M(F) coefficients can be realized as path weights in a planar graph, Karlin's closure follows from concatenating networks.
- Not developed beyond the suggestion.

### Assessment: partial TN2 closure is the right target
- Full TN2 is too strong (FALSE).
- The correct condition is a "specified minors" or "ladder-minor" positivity.
- The 4×4 N_c matrix formulation is the most concrete path forward.
- Need: (1) verify ladder-minor preservation computationally, (2) find conditions on N_c that guarantee it.

## Instance 3 (GPT): Direct Step-2 Proof

### 4-term decomposition of w_k at step 2 (DERIVED)
Using I_1 = E_1 + xJ_1, the step-2 w_k decomposes as:
```
w_k = d_k(β·C, C) + d_k(x·β·J_1·g, C) + d_k(x·β·E_1·h, C) + d_k(x²·β·J_1·h, C)
```
where β = (1+x)^ℓ, C = E_1·E_2, g = E_2, h = J_2.
- **Term 1** ≥ 0: d_k((1+x)^ℓ · C, C) ≥ 0 since C is LC (binomial smoothing lemma).
- **Terms 2-4**: sign-indefinite individually.

### Star+star ALL-K PROOF (MAJOR RESULT)
**Method: derivative trick.**

At the star+star, J_1 = J_2 = [1], so Terms 3 and 4 simplify dramatically:
- B = E_acc (since h = [1])
- C = (1+x)^{a+b}

Define W(x) = Σ_k w_k x^k (generating function of w_k sequence).
Then:
```
W(x) = (1+x)·E'(x) - D·E(x)
```
where E = E_acc, D = deg(C) = a+b, and E'(x) = dE/dx.

**Key insight:** Write E = (1+x)^ℓ · ((1+x)^a + x) · (1+x)^b where the middle factor handles the extra x from I_1 = (1+x)^a + x. Then:
```
W(x) = β · G(x)
```
where β = (1+x)^{ℓ-1} · (something positive) and G(x) has all nonneg coefficients.

The proof works for all ℓ ≥ 1 (at least one leaf at the support vertex).

**This completes the star+star case for ALL k, not just k=2.**

### J-complexity monotonicity: FALSE
- **Counterexample provided** (GPT): a specific pair (E_1, J_1) where w_k(t) at the interpolation J_1(t) = (1-t)·[1] + t·J_1 is NOT monotone in t.
- The "reduce to star+star base case via interpolation" strategy is dead.

### Factor-level scalar rescue ratios: MIRAGE
- The rescue ratio Karlin/|correction| varies by k value within the same tree.
- A single scalar bound per factor would need to work for all k simultaneously.
- This is impossible because the extremal k shifts as factors change.
- "The rescue ratio is a diagnostic, not a proof mechanism."

### "Condition C" framework (PROMISING)
The curvature identity reformulation:
```
C_k · w_k = C_k · Karlin_k + C_{k+1} · d_{k-1}(B,C) + B_k · c_k(C)
```
- Term 1 ≥ 0 (Karlin)
- Term 3 ≥ 0 (curvature of LC polynomial C)
- Term 2 can be negative

**Condition C**: −d_{k-1}(B,C) ≤ θ_k · (B_k / C_{k+1}) · c_k(C) for some θ_k ∈ [0,1].

If Condition C holds, then:
```
C_k · w_k ≥ C_k · Karlin_k + (1 - θ_k) · B_k · c_k(C) ≥ 0
```

At the star+star extremal, θ_k < 1 always (verified). The question is whether Condition C is preserved under the tree DP.

### Best path assessment (GPT ranking):
1. **Curvature identity + STP2 sign control + extremal argument** — most concrete
2. **4×4 N_c matrix with ladder-minor preservation** — most algebraically clean
3. **Direct step-2 inequality via 4-term decomposition** — requires bounding 3 sign-indefinite terms

## Strategic Assessment

### What's now PROVED or CONFIRMED
1. Star+star ALL-K: w_k ≥ 0 for all k (derivative trick proof)
2. Chain STP2 extended to n=22 (0 failures)
3. σ/M(F) matrix identity verified computationally
4. Step-2 always hardest (global min at step 2, s=2)
5. Step-2 ratio independent of total s and child ordering (for global min)

### What's DEAD (new this round)
1. **Step-ratio monotonicity**: FALSE (22 violations n≤17)
2. **Full TN2 of concatenated Φ**: FALSE (counterexample: path length 3)
3. **J-complexity interpolation**: FALSE (counterexample from GPT)
4. **Factor-level scalar rescue ratios**: diagnostic only, not a proof mechanism

### THE SINGLE REMAINING GAP (unchanged)
**Chain STP2 multi-child closure**: Given chain STP2 for each child, show STP2(∏I_i, ∏E_i) holds.

### Most promising next directions (ranked, updated)

1. **4×4 N_c matrix + ladder-minor preservation** (GPT R15 P2):
   Formalize the update matrix, check computationally whether ladder-minor positivity
   is preserved, then prove it via structural properties of N_c.

2. **Condition C framework** (GPT R15 P3):
   Prove −d_{k-1}(B,C) ≤ θ · (B_k/C_{k+1}) · c_k(C) at each step.
   Relates to ratio control via STP2. Star+star is the extremal.

3. **Direct step-2 proof via 4-term decomposition + aggregation**:
   Prove w_k ≥ 0 at step 2, then show later steps are easier (not via monotonicity
   of ratios, but via growth of Karlin term relative to correction).
