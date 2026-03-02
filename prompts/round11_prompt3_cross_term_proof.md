# Round 11, Prompt 3: Algebraic Proof of the Cross-Term Non-Negativity

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving E ≽ J (ratio dominance) at every support vertex of every tree, which implies unimodality of the IS sequence (Erdos Problem #993).

At each incremental step processing a non-leaf child t of the support vertex:
- E_acc, J_acc: accumulated polynomials from previous children (E_acc ≽ J_acc by induction)
- f = I_t = E_t + xJ_t, g = E_t: DP pair from child t (the IS poly and exclude-root poly of child subtree)
- E_new = E_acc · f (convolution), J_new = J_acc · g (convolution)

The P2 minor to verify:
```
w_k = d_k(E_new, J_new) = E_new_{k+1} · J_new_k - E_new_k · J_new_{k+1} ≥ 0
```

## The Diagonal/Cross Decomposition

Expanding E_new = E_acc * f and J_new = J_acc * g:

```
w_k = Σ_{i,j} E_acc_i · J_acc_j · φ(i,j)
```

where φ(a,b) = f_{k+1-a}·g_{k-b} - f_{k-a}·g_{k+1-b}.

Split: w_k = D_k + X_k where

**Diagonal (i=j):** D_k = Σ_i E_acc_i · J_acc_i · d_{k-i}(f, g)

Here d_p(f,g) = f_{p+1}·g_p - f_p·g_{p+1} is the factor-level LR minor. Since f = I_t and g = E_t, this is d_p(I_t, E_t), which can be negative (SCC/I≽E fails for some subtrees, e.g. T_{3,4} broom at n=28).

**Cross (i≠j):** X_k = Σ_{i≠j} E_acc_i · J_acc_j · φ(i,j)

### CRITICAL COMPUTATIONAL FINDING:

Exhaustive verification for n = 3..18 (410K support vertices with s ≥ 2, all incremental steps, all k):

- **X_k ≥ 0 in ALL cases (0 failures)**
- **D_k ≥ 0 in ALL cases (0 failures)**

Both parts are independently non-negative. Since w_k = D_k + X_k and both are ≥ 0, the proof of w_k ≥ 0 reduces to proving EITHER one.

## Available Structural Facts

At each incremental step, the following are known:

1. **E_acc ≽ J_acc** (ratio dominance of accumulated pair, by induction)
   - Equivalently: E_acc_i/J_acc_i is nondecreasing in i
   - For i < j: E_acc_i · J_acc_j ≤ E_acc_j · J_acc_i

2. **E_acc, J_acc have nonneg coefficients** (IS counting polynomials)

3. **g = E_t is LC** (log-concave, verified for all trees n ≤ 27)

4. **f = I_t = E_t + xJ_t** (IS polynomial of child subtree)
   - f is LC for most trees (but can fail)
   - f = g + xJ_t where J_t ≤ g coefficientwise (J ≤ E is proved)

5. **f ≽ g does NOT hold in general** (this is the SCC condition I ≽ E, which fails at n=28)

6. **g ≽ J_t** (E_t ≽ J_t at the child, by induction at support vertices)
   - But the child root may not be a support vertex. More carefully: this holds when the child subtree root happens to be a support vertex of the child subtree, which isn't always the case.

## Your Tasks

### Part 1: Prove D_k ≥ 0 (Diagonal Non-Negativity)

D_k = Σ_i E_acc_i · J_acc_i · d_{k-i}(f, g)

The weights E_acc_i · J_acc_i are all ≥ 0. The factor-level minors d_p(f,g) = d_p(I_t, E_t) can be negative.

But d_p(I_t, E_t) = J_{t,p}·E_{t,p} - J_{t,p-1}·E_{t,p+1} (by direct computation). This is positive when J has "enough curvature" relative to E.

Can you show that the convolution Σ_i w_i · d_{k-i}(f,g) is ≥ 0 for any LC weight sequence w_i ≥ 0? If not, what property of (E_acc, J_acc) beyond non-negativity is needed?

Note: the product E_acc_i · J_acc_i is the "pointwise product" of two sequences related by ratio dominance. Does this give the product sequence any special properties (e.g., stronger LC, TP2)?

### Part 2: Prove X_k ≥ 0 (Cross-Term Non-Negativity)

Group by pairs (i,j) with i < j:

```
X_k = Σ_{i<j} [E_acc_i · J_acc_j · φ(i,j) + E_acc_j · J_acc_i · φ(j,i)]
```

Let P(i,j) = E_acc_i · J_acc_j · φ(i,j) + E_acc_j · J_acc_i · φ(j,i) be the pair contribution.

From E_acc ≽ J_acc: for i < j, E_acc_j · J_acc_i ≥ E_acc_i · J_acc_j. So the "larger" weight multiplies φ(j,i).

Can you show P(i,j) ≥ 0 for each pair (i,j)?

Or if pointwise P(i,j) ≥ 0 fails (check!), can you show a coupling across different pairs that makes the total non-negative?

Explicitly:
```
φ(i,j) = f_{k+1-i}·g_{k-j} - f_{k-i}·g_{k+1-j}
φ(j,i) = f_{k+1-j}·g_{k-i} - f_{k-j}·g_{k+1-i}
```

Let a = k-i, b = k-j (so a > b since i < j):
```
φ(i,j) = f_{a+1}·g_b - f_a·g_{b+1}        (a "long-range" LR minor)
φ(j,i) = f_{b+1}·g_a - f_b·g_{a+1}        (another "long-range" LR minor)
```

The first bracket involves columns (a+1, b) in the (f,g) matrix.
The second involves columns (b+1, a).

If f and g are related by TP2 (f ≽ g), then both brackets have controlled signs. But f = I_t ≽ g = E_t is NOT known (SCC fails).

However, f = g + xJ_t gives f a specific structure. Can this be exploited?

### Part 3: Connection to Known Inequalities

The cross-term X_k ≥ 0 says something like: "the off-diagonal bilinear form Σ_{i≠j} a_i b_j M_{ij} ≥ 0" where the weight matrix (a_i b_j) comes from a ratio-dominant pair and M_{ij} comes from the factor's coefficient structure.

Does this resemble:
- **FKG inequality** (positively correlated random variables on a lattice)?
- **Rearrangement inequality** (similarly ordered sequences have maximal inner product)?
- **Schur convexity** (symmetric functions of the eigenvalues)?
- **Covariance inequality** (E[XY] ≥ E[X]E[Y] for positively associated X,Y)?

The "HPC projection" interpretation: individual factor-level LR minors can fail, but the aggregate product structure (via E_acc ≽ J_acc weighting) always restores non-negativity. This is a "positive association" phenomenon.

### Part 4: Which Target is More Tractable?

Given the algebra, which is easier to prove: D_k ≥ 0 or X_k ≥ 0?

D_k involves a 1D convolution (sum over single index i). X_k involves a 2D sum over pairs. Intuitively D_k is simpler, but it requires understanding the sign pattern of factor-level LR minors d_p(I_t, E_t).

X_k uses the stronger 2D structure of E_acc ≽ J_acc, which might give more "handles" for algebraic manipulation.

Your assessment of which target to prioritize, and what additional structural data would help.

## Deliverables

1. Proof of D_k ≥ 0, or identification of why it fails algebraically and what's missing
2. Proof of X_k ≥ 0, or analysis of pair contributions P(i,j)
3. Connection to known inequality families (FKG, rearrangement, etc.)
4. Assessment: which target (D or X) is more promising?
5. If neither can be proved: what additional computational data would help? What conditions on (f, g, E_acc, J_acc) would suffice?
