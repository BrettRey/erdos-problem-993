# Round 16, Prompt 3: Condition C Closure Under the Tree DP

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving the independence set sequence of every tree is unimodal. The proof reduces to showing w_k ≥ 0 at each incremental step of the support-vertex DP, where:
```
w_k = Karlin_k + correction_k
```
with Karlin_k ≥ 0 (proved) and correction_k sign-indefinite.

The curvature identity (Round 12) gives:
```
C_k · w_k = C_k · d_k(A,C) + C_{k+1} · d_{k-1}(B,C) + B_k · c_k(C)
```
where:
- A = E_acc · g, B = E_acc · h, C = J_acc · g
- g = E_t (current child's exclude-root poly), h = J_t (include-root poly)
- Term 1 = C_k · d_k(A,C) ≥ 0 (Karlin, proved)
- Term 3 = B_k · c_k(C) ≥ 0 (C is LC, so c_k(C) ≥ 0; B_k ≥ 0)
- Term 2 = C_{k+1} · d_{k-1}(B,C) can be negative

**Condition C** asserts that Term 2 is controlled by Term 3:
```
-d_{k-1}(B,C) ≤ θ · (B_k / C_{k+1}) · c_k(C)
```
for some θ < 1. This gives:
```
C_k · w_k ≥ C_k · d_k(A,C) + (1 - θ) · B_k · c_k(C) ≥ 0
```

At the star+star extremal (J_1 = J_2 = [1]), θ < 1 always (verified). The rescue ratio at the extremal converges to 1+√3 ≈ 2.732, which corresponds to θ converging to some value strictly less than 1.

## Known results

1. **s=1 PROVED**: one non-leaf child, handles 63% of support vertices.
2. **Star+star ALL-K PROVED**: W(x) = β·G with G nonneg coefficients (derivative trick).
3. **w_k ≥ 0**: 0 failures exhaustive n≤17, sampled n≤24.
4. **Chain STP2**: STP2(E_t, J_t) and STP2(I_t, E_t) hold at all tree-realizable factors, 0 failures n≤22.
5. **STP2 diagonal**: under LC(E_t), STP2(E_t, J_t) ⟺ J_t(k+1)·E_t(k-1) ≤ J_t(k)·E_t(k).
6. **J ≤ E coefficientwise**: proved.
7. **J and E are LC**: proved (products of LC/PF2 factors).

## Your Tasks

### Part 1: Unpack d_{k-1}(B, C)

```
d_{k-1}(B,C) = B_k · C_{k-1} - B_{k-1} · C_k
```

Substitute B = E_acc · h and C = J_acc · g where h = J_t, g = E_t:
```
B_k = [E_acc · J_t](k) = Σ_i E_acc(i) · J_t(k-i)
C_{k-1} = [J_acc · E_t](k-1) = Σ_j J_acc(j) · E_t(k-1-j)
```

Expand d_{k-1}(B,C) as a double sum:
```
d_{k-1}(B,C) = Σ_{i,j} E_acc(i) · J_acc(j) · [J_t(k-i) · E_t(k-1-j) - J_t(k-1-i) · E_t(k-j)]
```

The bracket [J_t(k-i)·E_t(k-1-j) - J_t(k-1-i)·E_t(k-j)] involves 2×2 minors of the (J_t, E_t) factor at shifted indices.

When i = j: the bracket is d_{k-i-1}(J_t, E_t) — a minor of the factor pair.

When i ≠ j: the bracket involves off-diagonal entries.

**Key question:** Can you express d_{k-1}(B,C) in terms of:
- STP2 minors of (E_t, J_t)
- LC gaps of E_acc and J_acc
- Cross terms involving (E_acc, J_acc) and (E_t, J_t)

### Part 2: Sign analysis using STP2

STP2(E_t, J_t) gives: for m > n, J_t(m+1)/J_t(m) ≤ E_t(n+1)/E_t(n).

Diagonal form: J_t(k+1)·E_t(k-1) ≤ J_t(k)·E_t(k) for all k.

In the bracket [J_t(k-i)·E_t(k-1-j) - J_t(k-1-i)·E_t(k-j)]:
- Let a = k-i, b = k-1-j.
- If i < j (so a > b+1): bracket = J_t(a)·E_t(b) - J_t(a-1)·E_t(b+1).
  This is a STP2 minor with m-n = a-b-1 > 0. By STP2(E_t, J_t): this is ≤ 0.
- If i > j (so a < b+1): bracket = J_t(a)·E_t(b) - J_t(a-1)·E_t(b+1).
  Now a-1 < b, so this is a "reversed" minor. Sign?
- If i = j: bracket = J_t(a)·E_t(a-1) - J_t(a-1)·E_t(a) = -d_{a-1}(E_t, J_t).
  By E_t ≽ J_t (prefix): this is ≤ 0 in the prefix.

So the diagonal terms (i=j) contribute negatively. The off-diagonal terms with i < j also contribute negatively (STP2). The terms with i > j have uncertain sign.

**The sum** d_{k-1}(B,C) is a weighted sum (weights E_acc(i)·J_acc(j) ≥ 0) of these brackets. The negative contributions (i ≤ j) are balanced against the positive contributions (i > j).

**Can you bound** the negative sum by the positive sum? Or by c_k(C)?

### Part 3: The curvature rescue c_k(C)

```
c_k(C) = C_k² - C_{k-1}·C_{k+1} = [J_acc · E_t]_k² - [J_acc · E_t]_{k-1} · [J_acc · E_t]_{k+1}
```

This is the LC gap of C = J_acc · E_t. Since both J_acc and E_t are LC, their product C is LC, so c_k(C) ≥ 0.

By the Cauchy-Binet identity for LC gaps of products:
```
c_k(C) = Σ_j c_j(J_acc) · E_t(k-j)² + Σ_j J_acc(j)² · c_{k-j}(E_t)
        + Σ_{j<j'} [cross terms involving both J_acc and E_t LC gaps]
```

Actually, the correct Cauchy-Binet for LC gaps of convolutions is:
```
c_k(FG) = Σ_{i,j} [F(i)·F(j) - F(i-1)·F(j+1)] · [G(k-i)·G(k-j) - G(k-i-1)·G(k-j+1)] / ...
```

This is getting complicated. **Derive** the exact Cauchy-Binet expansion of c_k(C) in terms of c_j(J_acc) and c_m(E_t).

Then compare term-by-term with the expansion of d_{k-1}(B,C) from Part 1. If each term in d_{k-1}(B,C) can be bounded by a corresponding term in c_k(C) (with a ratio < 1), Condition C follows.

### Part 4: Condition C at step 2

At step 2:
- E_acc = (1+x)^ℓ · I_1, J_acc = E_1
- g = E_2, h = J_2
- B = (1+x)^ℓ · I_1 · J_2
- C = E_1 · E_2

Since J_acc = E_1 (a single factor), the LC gaps c_k(C) = c_k(E_1 · E_2) have a simpler structure.

**Specialize** the Condition C analysis to step 2:
1. d_{k-1}(B,C) with B = (1+x)^ℓ · I_1 · J_2 and C = E_1 · E_2
2. The (1+x)^ℓ factor adds binomial smoothing
3. The I_1 = E_1 + x·J_1 structure means B has an "extra J_1" component

At the star+star extremal (J_1 = J_2 = [1]):
- B = (1+x)^ℓ · ((1+x)^a + x)
- C = (1+x)^{a+b}
- d_{k-1}(B,C) involves only binomials and the "extra x" term

**Can you prove Condition C (θ < 1) at step 2 for general (E_1, J_1, E_2, J_2)?**

### Part 5: Inductive preservation of Condition C

Suppose Condition C holds at step t with parameter θ_t. At step t+1, the state updates:
```
E_acc → E_acc · I_{t+1}
J_acc → J_acc · E_{t+1}
```

The new Condition C parameter θ_{t+1} depends on:
- The old θ_t (how well-controlled the current state is)
- The factor (I_{t+1}, E_{t+1}, J_{t+1})
- The LC structure of the new C = J_acc · E_{t+1} · E_{t+2} · ...

**Key question:** Does θ decrease (improve) at later steps? The computational evidence says the W ratio (related to 1/θ) increases at later steps. This suggests θ decreases.

**If θ_2 < 1 is proved (step 2, the hardest), does θ_t ≤ θ_2 follow for t ≥ 3?**

One mechanism: at later steps, J_acc is a product of more E factors, hence "more LC" (larger LC gaps relative to coefficients). This makes c_k(C) relatively larger, giving more curvature rescue.

**Formalize:** if c_k(J_acc · g) is "monotonically larger relative to |d_{k-1}(E_acc · h, J_acc · g)|" as more factors are accumulated, then θ decreases.

### Part 6: STP2 as a sufficient condition for Condition C

STP2(E_t, J_t) controls the ratio structure of the factor. Condition C needs the negative part of d_{k-1}(B,C) to be controlled by c_k(C).

**Conjecture:** STP2(E_t, J_t) + LC(E_t) + LC(J_t) + J_t ≤ E_t ⟹ Condition C at each step.

This would close the gap: chain STP2 is the inductive hypothesis, and Condition C at each step gives w_k ≥ 0, which gives E_new ≽ J_new, which gives STP2(E_new, J_new) (since STP2 follows from the prefix E ≽ J structure at support vertices... does it?).

Actually, the logical flow needs to be:
1. Assume chain STP2 for all children (IH)
2. At each step, use STP2(E_t, J_t) + properties to prove Condition C
3. Condition C gives w_k ≥ 0
4. w_k ≥ 0 gives E_new ≽ J_new (prefix)
5. After all steps: E ≽ J at this vertex
6. From E ≽ J + LC + J≤E: derive STP2(E, J) at this vertex
7. From I = E + xJ + STP2(E,J): derive STP2(I, E) at this vertex
8. This completes the induction

**Check:** does step 6 work? I.e., does E ≽ J + LC(E) + J ≤ E imply STP2(E, J)?

We know STP2 does NOT follow from abstract axioms (Round 14 counterexample). But the counterexample doesn't have J ≤ E coefficientwise AND E ≽ J prefix AND LC(E) AND LC(J) all simultaneously. Check: does the counterexample E=(1,19,77,72,50,22), J=(1,17,21,24,25) satisfy ALL these conditions?

## Deliverables

1. Expansion of d_{k-1}(B,C) as double sum with sign analysis
2. Which terms are controlled by STP2, which are not
3. Cauchy-Binet expansion of c_k(C) for comparison
4. Step-2 specialization of Condition C
5. Inductive preservation mechanism (does θ decrease?)
6. STP2 as sufficient condition for Condition C — logical feasibility
7. Logical flow check: does the induction close (steps 1-8)?
8. Overall: can Condition C be proved at step 2 and then propagated?
