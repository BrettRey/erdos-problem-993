# Round 15, Prompt 1: Step-by-Step Ratio Monotonicity

**Target model:** Codex 5.3 (exhaustive computation)

---

## Context

The W-form at each incremental step gives w_k = Karlin(k) + correction(k), where Karlin ≥ 0.

Round 14 confirmed: the global minimum W ratio (Karlin/|correction| when correction < 0) across ALL trees and ALL s values is always achieved at **step = 2** and **k = 2** (for n ≥ 13). This suggests step 2 is the hardest step.

If we can show:
1. w_k ≥ 0 at step 2 (direct proof)
2. The W ratio is nondecreasing in step t (later steps are easier)

then w_k ≥ 0 follows for all steps.

**Known results:**
- s=1 (step 1): PROVED
- Step 2 extremal: star+star with ratio → 1+√3 ≈ 2.732 (bounded away from 1)
- w_k ≥ 0: 0 failures exhaustive n≤17, sampled to n=24

## Setup

At a support vertex with ℓ leaves and s non-leaf children c_1, ..., c_s (ordered by increasing subtree size):

```
Step 0: E_acc = (1+x)^ℓ, J_acc = [1]
Step 1: E_acc *= I_1, J_acc *= E_1    (s=1 case, PROVED)
Step 2: E_acc *= I_2, J_acc *= E_2    (first genuine multi-factor interaction)
Step t: E_acc *= I_t, J_acc *= E_t
```

At each step t ≥ 2, define:
```
A = E_acc · E_t,  B = E_acc · J_t,  C = J_acc · E_t
Karlin_k = A[k+1]·C[k] - A[k]·C[k+1]
correction_k = B[k]·C[k] - B[k-1]·C[k+1]
w_k = Karlin_k + correction_k
ratio_k = Karlin_k / |correction_k|  (when correction_k < 0)
```

## Tasks

### Task 1: Ratio monotonicity in step t

For each tree at n=9-17 (exhaustive), at each support vertex with s ≥ 3, record the minimum W ratio at each step t = 2, 3, ..., s.

**Key question:** Is min_ratio(step t) ≤ min_ratio(step t+1) always? I.e., does the ratio never decrease as you process more children?

Report:
- For each n, the minimum ratio at step 2, step 3, step 4, ...
- Any violations of monotonicity (if step t+1 ratio < step t ratio)
- The tree/rooting where the worst step-3 ratio occurs

### Task 2: Profile the step-2 ratio across all s values

At a support vertex with s ≥ 3 non-leaf children, step 2 processes only the 2nd child. The accumulated pair at step 2 depends only on (ℓ, child_1), not on children 3+.

For each n=9-17, among ALL support vertices (regardless of total s), record the minimum step-2 ratio.

**Key question:** Is the step-2 minimum ratio the same regardless of whether the SV has s=2, s=3, or s=4 non-leaf children?

If yes, this confirms that step 2 is universally the hardest and the star+star analysis covers all cases.

### Task 3: What makes step 2 harder than step 3+?

At step 2: J_acc = E_1 (a single factor).
At step 3: J_acc = E_1 · E_2 (a product of two factors, hence "more LC", "smoother").

Hypothesis: the accumulated J_acc becomes more LC (larger LC gaps) at later steps, making the Karlin term relatively larger.

For each n=9-17 and s ≥ 3:
- Compute the normalized LC gap c_k(J_acc)/J_acc(k)² at each step
- Compute the normalized Karlin/correction ratio at each step
- Check correlation: does larger LC gap → larger ratio?

### Task 4: Minimum ratio at step 2 with general child ordering

Previous scans used increasing-subtree-size ordering (smallest child first). Check:
- At step 2, does the ordering matter for the minimum ratio?
- Try: largest child first, random orderings
- Report the minimum ratio across ALL orderings at step 2

If ordering doesn't affect the minimum ratio at step 2, this simplifies the proof target.

### Task 5: σ/M(F) matrix verification

The GPT Instance 3 from Round 14 proposed:
```
σ(F)(k) = F(k+1)  (shift operator, for F(0)=1)
M(F) = [[F, 0], [σ(F), 1]]
```
with M(F·G) = M(F)·M(G).

Verify computationally:
1. For each tree-realizable E_t with n ≤ 12, compute M(E_t) as a polynomial matrix.
2. For products E_1·E_2, check M(E_1·E_2) = M(E_1)·M(E_2) (matrix multiplication of polynomial matrices, where entries are polynomials multiplied by convolution).
3. The STP2 diagonal condition J(k+1)·E(k-1) ≤ J(k)·E(k) (under LC(E)) should correspond to some condition on the M matrices. Check: does it correspond to a 2×2 minor of a Toeplitz-like matrix built from M(E) and M(J)?

This is exploratory — the goal is to identify what positivity condition on M(F) matrices encodes STP2.

### Task 6: Extend chain STP2 to n=22

Round 14 checked chain STP2 through n=18 (exhaustive n≤17 + 50K at n=18). Extend:
- Exhaustive n=18 (123K trees)
- Sample 10K trees at n=19-22
- Report any failures

Also check: for each tree-realizable pair, compute both the diagonal STP2 condition (Lemma A: g(k+1)f(k-1) ≤ g(k)f(k)) and the full STP2 condition. Verify they agree (they should under LC).

## Deliverables

1. Step-ratio monotonicity table: min ratio at step 2, 3, 4, ... by n
2. Step-2 ratio independence of total s (confirmation or refutation)
3. LC gap vs ratio correlation analysis
4. Ordering independence at step 2
5. σ/M(F) matrix verification results
6. Chain STP2 extension to n=22
7. Assessment: is "step 2 is always hardest" provable from structural properties?
