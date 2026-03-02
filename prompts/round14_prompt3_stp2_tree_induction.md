# Round 14, Prompt 3: STP2 for Tree-Realizable Sequences

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

STP2 is a key structural condition on tree-realizable polynomial pairs that controls the sign of mixed brackets in the Cauchy-Binet expansion of the W-form.

### Definition

**STP2(f, g)**: g(m+1)·f(n) ≤ g(m)·f(n+1) for all m > n ≥ 0.

Equivalently: the ratio g(k)/f(k) is nonincreasing in k (for positive entries).

### Computational status

- **STP2(E_t, J_t)**: 0 failures through n=20 (1.6M+ factor instances, 11M+ comparisons)
- **STP2(I_t, E_t)**: status unknown (needs verification — see Prompt 1)
- **Chain STP2**: both (E_t, J_t) and (I_t, E_t) simultaneously — unverified

### Key negative result (Round 13)

**STP2 is NOT product-closed for general sequences.** Counterexample:
```
E₁ = (1, 12, 1), J₁ = (1, 2, 1)
E₂ = (1, 16, 1), J₂ = (1, 2, 1)
```
Each pair satisfies STP2, but (E₁·E₂, J₁·J₂) VIOLATES STP2.

So a proof of STP2 for tree-realizable pairs CANNOT use generic product closure. It must exploit the specific structure of tree DP.

### Tree DP recurrence

At a vertex v with children c₁, ..., c_d:
```
E_v = ∏ᵢ I_{c_i} = ∏ᵢ (E_{c_i} + x·J_{c_i})
J_v = ∏ᵢ E_{c_i}
I_v = E_v + x·J_v
```

Base case (leaf): E = [1], J = [1], I = [1, 1] = 1+x.

### Why STP2 matters for the proof

In the CB expansion of w_k, the mixed brackets are:
```
L_{ij} = h(a)·g(b) - h(a-1)·g(b+1)    (for appropriate indices a > b)
```
where g = E_child, h = J_child.

Under STP2(g, h): h(a)·g(b) ≤ h(a-1)·g(b+1), so L_{ij} ≤ 0 with controlled magnitude.

Without STP2, mixed brackets have unpredictable sign, making the CB analysis intractable.

## What We Know About Tree-Realizable (E, J)

These are the properties we can use in the proof:

1. **E(0) = J(0) = 1** (empty independent set)
2. **J ≤ E coefficientwise** (PROVED by induction)
3. **E is LC** (product of LC polynomials I_{c_i} = E_{c_i} + x·J_{c_i})
4. **J is LC** (product of E_{c_i}, each LC)
5. **E - J has nonneg coefficients** (from property 2)
6. **deg(E) = α(T_v)** (independence number of subtree)
7. **deg(J) = α(T_v) - 1** (including v costs one)
8. **I = E + xJ** is the IS polynomial of the subtree (positive coefficients, unimodal for trees)
9. **E ≽ J in prefix** (verified 0 failures n ≤ 22, provably for s = 1)

The key constraint that general sequences lack: **E and J are built from the SAME factors** (the children's E_{c_i}), and I_{c_i} = E_{c_i} + x·J_{c_i} is an IS polynomial (not arbitrary).

## Your Tasks

### Part 1: Single-Child Step

When v has one child c (degree-2 vertex):
```
E_v = I_c = E_c + x·J_c
J_v = E_c
```

**Show STP2(E_v, J_v)**, i.e., J_v(m+1)·E_v(n) ≤ J_v(m)·E_v(n+1) for m > n.

Substituting: E_c(m+1)·(E_c(n) + J_c(n-1)) ≤ E_c(m)·(E_c(n+1) + J_c(n))

Rearranging:
```
[E_c(m+1)E_c(n) - E_c(m)E_c(n+1)] ≤ E_c(m)J_c(n) - E_c(m+1)J_c(n-1)
```

LHS = -d_n(E_c, E_c) with m,n shifted = -(LC minor of E_c) ≤ 0 (since E_c is LC and m > n).

Wait, more carefully: LHS = E_c(m+1)E_c(n) - E_c(m)E_c(n+1). For m > n with m ≥ n+2, this is a Toeplitz minor. For m = n+1, this is the LC gap c_n(E_c) with reversed sign? No: if m = n+1, LHS = E_c(n+2)E_c(n) - E_c(n+1)², which is -c_{n+1}(E_c) ≤ 0. Good.

For general m > n: by repeated application of the TP2 property of LC sequences (Efron 1965), E_c(m+1)E_c(n) ≤ E_c(m)E_c(n+1). So LHS ≤ 0.

RHS = E_c(m)J_c(n) - E_c(m+1)J_c(n-1). By STP2(E_c, J_c) with indices (m, n-1) where m > n-1: J_c(m+1)E_c(n-1) ≤ J_c(m)E_c(n). Hmm, that's not quite the same...

**Question:** Does STP2(E_c, J_c) imply RHS ≥ 0? If so, then LHS ≤ 0 ≤ RHS and we're done.

RHS ≥ 0 means E_c(m)/E_c(m+1) ≥ J_c(n-1)/J_c(n), i.e., the E ratio at position m dominates the J ratio at position n-1. Since m > n > n-1, and both ratios are "growing" (due to LC), this should hold if E's ratios dominate J's ratios at corresponding positions.

Complete this argument. Show that STP2(E_c, J_c) + LC(E_c) implies STP2(I_c, E_c).

### Part 2: Multi-Child Step (The Hard Part)

When v has children c₁, ..., c_d with d ≥ 2:
```
E_v = I₁ · I₂ · ... · I_d
J_v = E₁ · E₂ · ... · E_d
```

Assume inductively: STP2(E_{c_i}, J_{c_i}) and STP2(I_{c_i}, E_{c_i}) for all i.

**Show:** STP2(E_v, J_v), i.e., STP2(∏I_i, ∏E_i).

The generic product closure fails (counterexample above). But tree-realizable pairs have special structure: each (I_i, E_i) is an IS-polynomial/exclude-root pair from the SAME subtree.

**Key structural constraint:** I_i = E_i + x·J_i, so (I_i)(k) = (E_i)(k) + (J_i)(k-1). The "extra" contribution from I_i vs E_i at degree k is exactly J_i(k-1), the include-root term.

**Approach A (incremental):** Process children one at a time.

Define F_t = I₁·...·I_t and G_t = E₁·...·E_t. Then:
```
F_{t+1} = F_t · I_{t+1},  G_{t+1} = G_t · E_{t+1}
```

Base case: STP2(F_1, G_1) = STP2(I_1, E_1) (given inductively).
Step: Show STP2(F_t·I_{t+1}, G_t·E_{t+1}) from STP2(F_t, G_t) + STP2(I_{t+1}, E_{t+1}).

This requires a "conditional product closure" theorem: if (F,G) satisfies STP2 and (I,E) satisfies STP2, AND both pairs are tree-realizable, then (F·I, G·E) satisfies STP2.

The counterexample to generic closure uses pairs like (1,12,1) and (1,2,1) which are NOT tree-realizable (no tree has E = (1,12,1)). So the counterexample doesn't apply.

**Can you identify what property of tree-realizable pairs makes product closure work?**

**Approach B (Toeplitz matrix):** The STP2 condition says the matrix M with M_{m,n} = J(m)E(n) (for m > n) has a sign pattern controlled by the Toeplitz structure. For products, M_{m,n} = (J₁·J₂)(m)·(E₁·E₂)(n) = Σ_{i,j} J₁(i)J₂(m-i)E₁(j)E₂(n-j).

This is a sum of products of STP2 matrices. The Hadamard product of positive semidefinite matrices is PSD (Schur product theorem). Can you reformulate STP2 as a PSD condition and apply Schur?

**Approach C (ratio characterization):** STP2 says r(k) = J(k)/E(k) is nonincreasing. For products: r_product(k) = (J₁·J₂)(k)/(E₁·E₂)(k). This is a ratio of convolutions, NOT the product r₁(k)·r₂(k). But if r₁ and r₂ are both nonincreasing AND the underlying polynomials have LC structure, the convolution ratio might inherit monotonicity.

Known result (Efron 1965, Karlin 1968): convolution preserves PF2 (LC). Does it preserve ratio monotonicity?

### Part 3: Does STP2 Follow from Simpler Conditions?

Check whether STP2 is equivalent to, or follows from, any combination of:
(a) E ≽ J (prefix ratio dominance)
(b) J ≤ E coefficientwise
(c) LC of both E and J
(d) E(0) = J(0) = 1

Construct a counterexample (non-tree-realizable) satisfying (a)-(d) but violating STP2. If such examples exist, STP2 is genuinely tree-specific.

If no counterexamples exist for "reasonable" coefficient sizes (say degree ≤ 10, coefficients ≤ 1000), STP2 might follow from (a)-(d) after all.

### Part 4: Direct STP2 Proof via Tree DP

Instead of induction via products, try proving STP2 directly from the tree DP.

At vertex v with subtree T_v:
```
E_v(k) = #{independent sets of size k in T_v not containing v}
J_v(k) = #{independent sets of size k in T_v containing v}  (so J_v = I_v[1:] shifted)
```

Wait, J_v(k) counts IS of size k+1 containing v (after dividing out x). So:

STP2 says: for m > n,
```
J_v(m+1) · E_v(n) ≤ J_v(m) · E_v(n+1)
```

i.e., #{IS of size m+2 with v} · #{IS of size n without v} ≤ #{IS of size m+1 with v} · #{IS of size n+1 without v}

This is a log-concavity-like condition across the (include v, exclude v) partition.

**Combinatorial interpretation:** Consider the bipartite graph on IS's with v vs IS's without v, where we connect sets that differ by exactly one element. Does a flow argument or injection show the ratio J(k)/E(k) is nonincreasing?

This is related to the **hard Lefschetz theorem** for independence complexes (Adiprasito-Huh-Katz type results). The ratio J/E being nonincreasing might follow from the "flipping v" bijection interacting with the independence complex structure.

### Part 5: Prove Chain STP2

If you can prove STP2(E_v, J_v) for all v, then try:

**Chain STP2:** Both STP2(E_v, J_v) AND STP2(I_v, E_v) hold for all tree vertices v.

STP2(I_v, E_v) says E_v(m+1)·I_v(n) ≤ E_v(m)·I_v(n+1), i.e., the ratio E_v(k)/I_v(k) is nonincreasing.

Since I_v(k) = E_v(k) + J_v(k-1), this is:
```
E_v(m+1) · [E_v(n) + J_v(n-1)] ≤ E_v(m) · [E_v(n+1) + J_v(n)]
```

which is LC(E_v) + a correction involving J_v. The correction needs STP2(E_v, J_v) to be non-negative.

This looks like Part 1 but for general v (not just single-child).

## Deliverables

1. Single-child STP2 proof (complete the argument from Part 1)
2. Multi-child STP2: analysis of which approach (A, B, or C) is most viable
3. Counterexample search: does STP2 follow from (a)-(d)?
4. Assessment of the direct combinatorial approach (Part 4)
5. Chain STP2 analysis (Part 5)
6. Overall assessment: can STP2 be proved for tree-realizable sequences? What's the main obstacle?
