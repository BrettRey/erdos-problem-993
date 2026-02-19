# Algebraic Analysis of the Degree-2 Bridge Decomposition

Date: 2026-02-18

## Setup

Let T be a d_leaf ≤ 1 tree on n vertices. By the existence lemma (counting proof), T
has a leaf l whose support s has degree 2. Let u be the other neighbour of s.

Define:
- A = T - l (remove the leaf)
- B = T - {l, s} (remove the pendant edge)
- P = dp_B[u][0] (IS polynomial of B's subtree at u, with u excluded)
- Q = dp_B[u][1] (IS polynomial of B's subtree at u, with u included)

## Polynomial identities

Since s has degree 2 and is adjacent to l and u:

```
I(T) = (1+2x)P + (1+x)Q
I(A) = (1+x)P + Q
I(B) = P + Q
```

Equivalently: I(T) = (1+x)I(B) + xP.

**Verification**: exact match on all 77,141 d_leaf ≤ 1 trees through n=20.

## Three-term decomposition of Φ_m

The bridge identity gives Φ_m(T; λ) = Φ_m(A; λ) + λ·Φ_{m-1}(B; λ). Since
I(A) = (1+x)I(B) + xP, and Φ_m is linear in the coefficients:

```
Φ_m(A; λ) = Φ_m((1+x)I_B; λ) + Φ_m(xP; λ)
```

For the first term, if F = (1+x)G has coefficients f_k = g_k + g_{k-1}:
- Z_F(λ) = (1+λ)Z_G(λ)
- μ_F(λ) = μ_G(λ) + λ/(1+λ)
- Φ_m(F; λ) = (1+λ)Φ_m(G; λ) + λZ_G(λ)

For the second term, if H = xP has coefficients h_k = p_{k-1}:
- Z_H(λ) = λ Z_P(λ)
- μ_H(λ) = μ_P(λ) + 1
- Φ_m(H; λ) = λ Z_P(μ_P + 1 - (m-1)) = λ Z_P(μ_P - (m-2)) = λ Φ_{m-1}(P; λ)

Combining with the bridge identity:

```
Φ_m(T; λ) = [(1+λ)Φ_m(B; λ) + λZ_B(λ) + λΦ_{m-1}(P; λ)] + λΦ_{m-1}(B; λ)
           = (1+λ)Φ_m(B; λ) + λZ_B(λ) + λΦ_{m-1}(P; λ) + λΦ_{m-1}(B; λ)
```

Wait -- let me redo this carefully. We have:
- Φ_m(T) = Φ_m(A) + λΦ_{m-1}(B)       ... (bridge)
- Φ_m(A) = (1+λ)Φ_m(B) + λZ_B + λΦ_{m-1}(P)  ... (A decomposition)

So:
```
Φ_m(T; λ) = (1+λ)Φ_m(B; λ) + λZ_B(λ) + λΦ_{m-1}(P; λ) + λΦ_{m-1}(B; λ)
```

This is a **four**-term decomposition. Let me group differently:

```
Φ_m(T) = (1+λ)Φ_m(B) + λΦ_{m-1}(B) + λZ_B + λΦ_{m-1}(P)
        = Φ_m(B) + λ[Φ_m(B) + Φ_{m-1}(B)] + λZ_B + λΦ_{m-1}(P)
```

Note: Φ_m(B) + Φ_{m-1}(B) = Z_B(μ_B - (m-1)) + Z_B(μ_B - (m-2)) = Z_B(2μ_B - 2m + 3).

Actually, let me take a step back. The **original** two-term bridge identity is:

```
Φ_m(T; λ) = Φ_m(A; λ) + λ·Φ_{m-1}(B; λ)
```

The diagnostic confirmed: Φ_m(A) ≥ 0 and Φ_{m-1}(B) ≥ 0 for all d_leaf ≤ 1 trees
(with degree-2 support leaf choice). **Each of these two terms is independently ≥ 0.**

The three-term decomposition of Φ_m(A) = (1+λ)Φ_m(B) + λZ_B + λΦ_{m-1}(P) gives
additional structure.

## Diagnostic findings (n ≤ 20, 77,141 trees)

| Quantity | Failures | Notes |
|----------|----------|-------|
| Identity verification | 0 | P,Q decomposition exact |
| no_deg2_support | 0 | Every d_leaf≤1 tree has one |
| Φ_{m-1}(P; λ) < 0 | **0** | P always has μ_P(λ) ≥ m-2 |
| P not log-concave | **0** | P is always LC |
| Q not log-concave | **0** | Q is always LC |
| P not unimodal | **0** | P is always unimodal |
| (1+λ)Φ_m(B) < 0 | 63,617 (82%) | μ_B(λ) < m-1 is common |
| pendant bonus > 0 | 77,141 (100%) | λZ_B always positive (trivial) |
| pendant ratio ≥ 1.65 | 77,141 (100%) | Pendant always compensates |
| STRONG C2 fail | **0** | λ_m^T ≥ λ_{m-1}^B always |
| cross-tree fail | **0** | λ_m^A ≥ λ_{m-1}^B always |
| total Φ_m(T) < 0 | **0** | Main conjecture holds |

## Key structural observations

### 1. Φ_{m-1}(P; λ) ≥ 0 always (STRONG)

P is the IS polynomial of the "u-off" subtree in B. It is always log-concave and
its weighted mean μ_P(λ_m) ≥ m-2 always holds. The minimum value of Φ_{m-1}(P; λ)
across all 77,141 trees is **0.375** (at n=4, the path P_4).

This means Φ_m(A; λ) = (1+λ)Φ_m(B) + λZ_B + λΦ_{m-1}(P) has only one potentially
negative term: (1+λ)Φ_m(B).

### 2. The pendant bonus λZ_B always compensates

When (1+λ)Φ_m(B) < 0, the pendant bonus λZ_B combined with λΦ_{m-1}(P) always
compensates. The worst-case ratio is pendant_bonus/|negative_term| ≈ 1.65.

### 3. Sufficient condition: μ_B(λ) ≥ m - 1 - λ/(1+λ)

Since Φ_{m-1}(P) ≥ 0, the condition Φ_m(A) ≥ 0 reduces to:

```
(1+λ)Φ_m(B) + λZ_B ≥ 0
(1+λ)Z_B(μ_B - (m-1)) + λZ_B ≥ 0
(1+λ)(μ_B - m + 1) + λ ≥ 0
μ_B ≥ m - 1 - λ/(1+λ)
```

Since λ ∈ (0, 1) for mode-tied trees, λ/(1+λ) ∈ (0, 1/2), so this requires only
μ_B > m - 3/2, which is weaker than μ_B ≥ m-1.

### 4. The tightest witness

The tree with minimum STRONG C2 margin (0.114) at n=20 has degree signature
{1:9, 2:8, 3:1, 5:2} -- a near-caterpillar with two degree-5 hubs. The tree
with minimum pendant ratio (1.65) has degree signature {1:10, 2:9, 10:1} --
a star-like tree with one degree-10 hub.

## Proof strategy assessment

### Lane A (STRONG C2 + induction): VIABLE

STRONG C2 (λ_m^T ≥ λ_{m-1}^B) holds with 0 failures and decent margin. The algebraic
condition is:

```
a_{m-1}·b_{m-1} ≥ a_m·b_{m-2}
```

where a_k, b_k are coefficients of A = (1+x)B + xP and B. Substituting
a_k = b_k + b_{k-1} + p_{k-1}:

```
(b_{m-1} + b_{m-2} + p_{m-2})·b_{m-1} ≥ (b_m + b_{m-1} + p_{m-1})·b_{m-2}
```

Expanding:
```
b_{m-1}² - b_m·b_{m-2} + b_{m-2}·b_{m-1} - b_{m-1}·b_{m-2} + p_{m-2}·b_{m-1} - p_{m-1}·b_{m-2} ≥ 0
b_{m-1}² - b_m·b_{m-2} + p_{m-2}·b_{m-1} - p_{m-1}·b_{m-2} ≥ 0
```

This is: LC_defect(B, m-1) + p_{m-2}·b_{m-1} - p_{m-1}·b_{m-2} ≥ 0.

The first term is the log-concavity surplus of B at index m-1 (positive when B is LC).
The second pair compares the "shape" of P and B near m.

Since B is always LC (through n=26 for d_leaf≤1), this reduces to showing:
```
b_{m-1}² - b_m·b_{m-2} ≥ p_{m-1}·b_{m-2} - p_{m-2}·b_{m-1}
```

i.e., the LC surplus of B exceeds the P/B ratio mismatch.

### Lane B (pendant bonus): ALSO VIABLE

We need: (1+λ)Φ_m(B) + λZ_B ≥ 0, equivalently μ_B(λ) ≥ m - 1 - λ/(1+λ).

This is a weaker condition than μ_B ≥ m-1 and might follow from induction on the
tree structure: removing the pendant shifts the mean down by at most 1, but the
"slack" of λ/(1+λ) gives breathing room.

### Lane C (Φ_{m-1}(P) ≥ 0 + Φ_{m-1}(B) ≥ 0 separately): CLEANEST

The original bridge says Φ_m(T) = Φ_m(A) + λΦ_{m-1}(B). Both terms are independently
non-negative (0 failures). The Φ_m(A) ≥ 0 condition decomposes further but is itself
verified directly.

The key question: can we prove Φ_m(A) ≥ 0 and Φ_{m-1}(B) ≥ 0 by induction?

For Φ_{m-1}(B): B has n-2 vertices. If mode(B) ≥ m-1 (verified, 0 failures), and B
is also d_leaf ≤ 1 (it may not be!), then we could apply the induction hypothesis.
But B = T - {l,s} may introduce a new leaf at u, and u may now have d_leaf > 1 in B.

**Critical question**: Is B always d_leaf ≤ 1?

## Critical check: d_leaf ≤ 1 closure

**B = T-{l,s} is NOT always d_leaf ≤ 1**: 17,653 of 77,141 trees (23%) produce a B
that violates d_leaf ≤ 1. Similarly, A = T-l is not always d_leaf ≤ 1 (14,953 = 19%).

This means a naive induction *within the d_leaf ≤ 1 class* does not close. We would
need Φ non-negativity for a broader class of trees, or a different induction structure.

However, the Φ conditions (both bridge terms ≥ 0) hold with **0 failures** even for
those trees where B or A exit the d_leaf ≤ 1 class. This suggests the positivity is
actually a broader phenomenon.

### Induction obstacle

The issue: when we remove the pendant edge (l, s), vertex u may gain additional
leaf neighbours in B (vertices that were degree-2 in T but became degree-1 in B
after removing s from their neighbourhood). If u already had one leaf in T, it
now has two, violating d_leaf ≤ 1.

### Alternative strategies

**Strategy 1 (broader class):** Prove Φ_m ≥ 0 for all trees, not just d_leaf ≤ 1.
We already know this fails for general trees at n ≤ 16 (2 existence failures). But
the bridge Φ_m(A) + λΦ_{m-1}(B) decomposition with the degree-2 support leaf choice
might still hold for a broader class.

**Strategy 2 (STRONG C2 as the inductive invariant):** Instead of inducting on
Φ ≥ 0, induct on the STRONG C2 condition λ_m^T ≥ λ_{m-1}^B. This is 0 failures
through n=23 (4.5M degree-2 leaves) and has an algebraic form that interacts well
with the P,Q structure. The inductive step reduces to:

```
LC_surplus(B, m-1) ≥ p_{m-1}·b_{m-2} - p_{m-2}·b_{m-1}
```

Since P is always LC with bounded coefficients, this might be tractable.

**Strategy 3 (pendant bonus + mean bound):** Prove the weaker condition
μ_B(λ_m) ≥ m - 3/2 (instead of m-1). This has even more slack than the
exact condition. Might follow from the Steiner peeling apparatus already proved.

## Extended results: ALL trees (not just d_leaf ≤ 1)

Ran diagnostic on all trees through n=16 (29,355 trees with a degree-2 support leaf).
**All findings replicate exactly:**

| Quantity | Failures (all trees, n≤16) |
|----------|---------------------------|
| STRONG C2 | **0** (min margin 0.083) |
| Φ_{m-1}(P) < 0 | **0** |
| P not log-concave | **0** |
| total Φ < 0 | **0** |
| pendant ratio < 1 | **0** (min 1.78) |

Trees skipped (no deg-2 support) = 3,150, which are exactly the d_leaf > 1 trees.
The decomposition applies to every tree that has a degree-2 support leaf.

## Structural analysis of P

P = dp_B[u][0] where u is the hub vertex. When B is rooted at u:

```
P = ∏_{c ∈ children(u)} I(subtree_c)
```

So P is a **product of IS polynomials of subtrees**. Since:
1. IS polynomials of trees are log-concave (through n ≤ 26 computationally)
2. Products of LC sequences with positive terms are LC

P is always LC. This explains the observation P_not_lc = 0.

## Key open question: why μ_P(λ_m) ≥ m - 2

P = ∏ I(T_c). The weighted mean of a product at fugacity λ is:

```
μ_P(λ) = Σ_c μ_{T_c}(λ)
```

(because the independent sets from different subtrees are "independent" -- the
product structure means the total size is the sum of sizes from each component).

So μ_P(λ_m) ≥ m - 2 iff Σ_c μ_{T_c}(λ_m) ≥ m - 2.

This connects to the Steiner peeling framework: each subtree T_c contributes
its own weighted mean. The question is whether the sum of these means, evaluated
at the tree's mode fugacity, is at least m - 2.

Similarly, μ_B(λ) = μ_P(λ) + (probability u is in IS at fugacity λ) · 1. Wait, that's
not quite right. B = T - {l,s} and I(B) = P + Q. At fugacity λ:

```
μ_B(λ) = [μ_P(λ) · Z_P(λ) + μ_Q(λ) · Z_Q(λ)] / Z_B(λ)
```

But Q = dp[u][1] = x · ∏_c dp[c][0], so Q represents sets that include u. If we
denote the probability that u is in the IS as p_u = Z_Q / Z_B, then:

```
μ_B = (1-p_u) μ_P + p_u μ_Q = μ_P + p_u
```

(since μ_Q = μ_{xΠdp0} = 1 + Σ μ_{dp0[c]} by the same product-of-means argument).

Actually this is getting complicated. The simpler fact is:

## Sufficient condition summary

To prove Φ_m(T; λ_m) > 0 for all d_leaf ≤ 1 trees, it suffices to prove
**any one** of:

1. **Φ_m(A) ≥ 0 and Φ_{m-1}(B) ≥ 0** (both bridge terms non-negative)
   - 0 failures through n=23 (931K trees, 4.5M degree-2 leaves)

2. **Φ_{m-1}(P) ≥ 0** (the P-positivity condition, implies Φ_m(A) ≥ 0 via pendant bonus)
   - 0 failures through n=20 (77K trees)
   - Plus Φ_{m-1}(B) ≥ 0 (verified independently)

3. **STRONG C2 (λ_m^T ≥ λ_{m-1}^B)** (implies Φ_{m-1}(B) ≥ 0, then combine with Φ_m(A) ≥ 0)
   - 0 failures through n=23 (4.5M checks)

4. **μ_B(λ_m) ≥ m - 1 - λ/(1+λ)** (weaker mean bound on B, implies Φ_m(A) ≥ 0 directly)
   - Has ~0.5 units of slack beyond what's needed

## Next steps

1. Investigate the product-of-means structure for μ_P(λ_m) ≥ m-2
2. Try to connect STRONG C2 to the Steiner peeling cavity field framework
3. Consider whether the pendant bonus approach gives a clean one-page proof
