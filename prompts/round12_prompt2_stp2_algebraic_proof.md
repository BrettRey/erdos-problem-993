# Round 12, Prompt 2: STP2 Algebraic Analysis and Sign Control

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving prefix E ≽ J (ratio dominance for k < mode) at support vertices of trees. The W-form induction gives at each step:

```
w_k = d_k(E_new, J_new) where E_new = E_acc·f, J_new = J_acc·g
```

with f = I_t = E_t + xJ_t (child IS poly) and g = E_t (child exclude-root poly).

The Karlin split:
```
w_k = d_k(E_acc·g, J_acc·g) + d_k(x·E_acc·h, J_acc·g)
      \_________________/     \______________________/
       Karlin term (≥ 0)        correction (sign-indefinite)
```
where h = J_t. The Karlin term is non-negative by TP2 convolution (E_acc ≽ J_acc, g is LC).

The Cauchy-Binet expansion of W_k involves mixed brackets:
```
L_{i,j} = r_{k-i}·q_{k-1-j} - r_{k-1-i}·q_{k-j}
```
where q = E_t = g and r = J_t = h. The mixed bracket L_{i,j} is NOT sign-definite in general, and this is precisely the obstruction to proving w_k ≥ 0.

## The STP2 Condition

**STP2**: r_{m+1}·q_n ≤ r_m·q_{n+1} for all m > n ≥ 0.

Equivalently: the ratio r_m/q_m is nonincreasing in m (when both are positive). This says the Toeplitz-like matrix with rows indexed by r and columns by q has all 2×2 minors non-negative in a specific pattern.

**Computational status**: 0 failures across 50,917+ factor pairs at support vertices through n ≤ 20. This is the most promising algebraic condition found.

## Your Tasks

### Part 1: What Does STP2 Mean Structurally?

Here r = J_t and q = E_t are the DP polynomials of a child subtree rooted at vertex t.

STP2 says J_t(m+1)·E_t(n) ≤ J_t(m)·E_t(n+1) for m > n.

Interpret: the ratio J_t(k)/E_t(k) is nonincreasing. We already know:
- E_t ≽ J_t at support vertices (by induction): E_t(k+1)·J_t(k) ≥ E_t(k)·J_t(k+1)
- This says E_t(k+1)/E_t(k) ≥ J_t(k+1)/J_t(k), i.e., E_t's ratios dominate J_t's ratios

STP2 says something different: J_t(m+1)/J_t(m) ≤ E_t(n+1)/E_t(n) for m > n. This is a CROSS-INDEX comparison: the ratio of J at a LATER index is bounded by the ratio of E at an EARLIER index.

Since E_t is LC, E_t(n+1)/E_t(n) is nonincreasing. And if J_t is also LC, J_t(m+1)/J_t(m) is also nonincreasing. STP2 asks: do J_t's ratios at large indices stay below E_t's ratios at small indices?

This is automatically true when m = n+1 (that's just E_t ≽ J_t). The content is for m > n+1: the "gap" between the index positions makes J_t's ratio even smaller relative to E_t's ratio.

**Question**: Can STP2 be proved from:
1. E_t ≽ J_t (ratio dominance)
2. E_t and J_t both LC (log-concave)
3. J_t ≤ E_t coefficientwise

Or does it require additional tree-specific structure?

### Part 2: STP2 Implies Sign Control on L_{i,j}

Compute the sign of L_{i,j} under STP2.

L_{i,j} = r_{k-i}·q_{k-1-j} - r_{k-1-i}·q_{k-j}

Set a = k-i, b = k-j. For i < j (so a > b):
```
L_{i,j} = r_a·q_{b-1} - r_{a-1}·q_b
```

STP2 with m = a-1 > n = b-1 (since a > b) gives:
```
r_a·q_{b-1} ≤ r_{a-1}·q_b
```

So L_{i,j} ≤ 0 when i < j (under STP2).

This means the correction term J_{k+1}·L_{i,j} ≤ 0 for i < j. The correction is systematically NEGATIVE.

But the Karlin term involves K_{i,j} = q_{k+1-i}·q_{k-j} - q_{k-i}·q_{k+1-j}. For i < j, K_{i,j} ≤ 0 (from LC of q). And the coefficient (f_i·g_j - f_j·g_i) can have either sign.

**Key question**: Under STP2, can we show that the Karlin contributions (via K_{i,j}) dominate the correction contributions (via L_{i,j})?

Specifically, for each pair (i,j) with i < j, we need:
```
|J_k·K_{i,j}| ≥ |J_{k+1}·L_{i,j}|
```
after weighting by the appropriate (f,g) coefficients. Or if pointwise domination fails, can we show it in aggregate?

### Part 3: STP2 Proof Attempt

Try to prove STP2 for tree-realizable (E_t, J_t) pairs.

**Induction structure**: At a leaf, E = [1] and J = [1]. STP2 is vacuous.

At an internal vertex v with children c_1, ..., c_d:
```
E_v = product of I_{c_i} = product of (E_{c_i} + x·J_{c_i})
J_v = product of E_{c_i}
```

If STP2 holds for each child's (E_{c_i}, J_{c_i}), does it hold for (E_v, J_v)?

The product structure might preserve STP2 via a Karlin-type argument: if two pairs (q_1, r_1) and (q_2, r_2) each satisfy STP2, does the convolution pair (q_1·q_2·...·, r_1·r_2·...·) satisfy STP2?

Note: E_v = ∏I_{c_i} is NOT the product of E_{c_i}'s alone — it involves x·J_{c_i} terms. So this is not a simple product closure.

### Part 4: Alternative — Direct Correction Bound

If STP2 cannot be proved algebraically, consider bounding the correction directly.

The correction at step t is:
```
corr_k = d_k(x·E_acc·h, J_acc·g) = Σ_i E_acc_i · [h * J_acc]_{k-i} · g_{???} - ...
```

This is complicated. But empirically:
- The Karlin term always dominates the correction through n ≤ 20
- The min ratio (Karlin/|correction|) is 1.195 at n=20 (pendant-star, s=1, already PROVED)
- For s ≥ 2, the min ratio INCREASES: 1.26 (s=2), 1.32 (s=3), etc.

**Key structural fact**: h = J_t has degree one less than g = E_t. And x·h shifts h up by one degree. So x·E_acc·h has its mass shifted one position right relative to J_acc·g.

Can you bound d_k(x·E_acc·h, J_acc·g) using:
- The x-shift (h coefficients start one position later)
- J_t ≤ E_t coefficientwise (so h ≤ g)
- E_acc ≽ J_acc (by induction)
- LC of all polynomials involved

### Part 5: The Prefix Restriction

We actually only need w_k ≥ 0 for k < mode(I_new), not for all k. E≽J fails at n=32 only in the tail (k=15, mode=10).

Does restricting to k < mode help the algebra? For k < mode:
- I_k is increasing, so the coefficients are growing
- E_k and J_k are also near their growth phase
- The ratio J_k/E_k might be better behaved for small k

If STP2 gives L_{i,j} ≤ 0 (Part 2), the bound |J_{k+1}·L_{i,j}| involves J_{k+1}. For k < mode, J_{k+1}/J_k is larger (we're in the growth phase), which makes the correction relatively larger. So the prefix restriction might NOT help directly.

But the Karlin term also involves growth-phase coefficients. Analyze the relative scaling.

## Deliverables

1. Structural interpretation of STP2 (what property of trees ensures it?)
2. Sign of L_{i,j} under STP2 (confirm L_{i,j} ≤ 0 for i < j)
3. Bound on |correction|/|Karlin| under STP2 (pointwise or aggregate)
4. Proof of STP2 for tree-realizable pairs, or identification of what's missing
5. Assessment: does STP2 suffice to complete the proof? What else is needed?
