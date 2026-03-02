# Task: Prove product closure of Strong Condition C via incremental factor induction

## The problem

Prove that the independence polynomial of every tree is unimodal.

The proof reduces to a single algebraic claim: **Strong Condition C is preserved when multiplying by a single tree-realizable factor.** By iterating, this gives full product closure. This is verified computationally (0 failures, 2.8M products) but not proved.

Your task is to prove the **one-factor extension**: given that the accumulated product satisfies Strong Condition C, multiplying by one more tree-realized factor preserves it.

## Full setup (self-contained)

### Reduction to Condition C at support vertices

Every tree has a support vertex r (vertex adjacent to a leaf). At r with ℓ leaf children and non-leaf children c₁,...,c_s:

    E(x) = (1+x)^ℓ · A(x),    J(x) = B(x),
    I(T; x) = E(x) + x·J(x)

where A = ∏ I_{c_i} (product of subtree IS polys) and B = ∏ E_{c_i} (product of exclude-root polys).

- P3 (tail domination): **PROVED** at all support vertices.
- P2 (prefix ratio dominance): Reduces to Strong Condition C for the product (A, B).

### The key identity (PROVED, trivial algebra)

    b_{k-1}·Δ_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k

where Δ_k = P2 target, d_k = a_{k+1}b_k - a_kb_{k+1} (LR minors), c_k = b_k² - b_{k-1}b_{k+1} (LC gaps).

Strong Condition C = the RHS is ≥ 0, equivalently (1+x)·I ratio-dominates E.

### Strong Condition C (formal definition)

A pair (I, E) with I = E + x·J and nonneg integer coefficients satisfies **Strong Condition C** if for all k ≥ 1:

    b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0     (*)

where a_k = [x^k]I, b_k = [x^k]E, d_k = a_{k+1}b_k - a_kb_{k+1}, c_k = b_k² - b_{k-1}b_{k+1}.

**WARNING:** Replacing a_{k-1} with b_{k-1} (the WEAK version) FAILS at n=17. The amplification a_{k-1} ≥ b_{k-1} is essential.

### Tree-realizable triples

A triple (I, E, J) is **tree-realizable** if ∃ rooted tree T, root r such that:
- E = ∏_{c child of r} I(subtree(c); x)
- J = ∏_{c child of r} dp[c][0]  (exclude-root polys of subtrees)
- I = E + x·J

**Key structural property:** The SAME children c produce both E (via their IS polys I_c) and J (via their exclude-root polys E_c), with I_c ≥ E_c coefficientwise for each c.

**Verified properties (0 failures):**
1. J ≤ E coefficientwise (61.8M checks)
2. E is log-concave (11.9M checks)
3. Strong Condition C (11.9M factor + 59.9M root-level checks)

## The incremental product framework

### State after k factors

After processing children c₁,...,c_k, define:
- E^{(k)} = ∏_{i=1}^k I_{c_i}  (product of subtree IS polys)
- J^{(k)} = ∏_{i=1}^k E_{c_i}  (product of exclude-root polys)
- I^{(k)} = E^{(k)} + x·J^{(k)}

**Base case (k=0):** E^{(0)} = 1, J^{(0)} = 1, I^{(0)} = 1 + x. Condition C trivially holds.

**Base case (k=1):** E^{(1)} = I_{c₁}, J^{(1)} = E_{c₁}, I^{(1)} = I_{c₁} + x·E_{c₁}. Strong Condition C holds by the inductive hypothesis on subtree(c₁).

### The one-factor transition (k → k+1)

Adding factor (I_c, E_c) (abbreviating I_{c_{k+1}}, E_{c_{k+1}}):

    E^{(k+1)} = E^{(k)} · I_c
    J^{(k+1)} = J^{(k)} · E_c
    I^{(k+1)} = E^{(k+1)} + x · J^{(k+1)} = E^{(k)}·I_c + x·J^{(k)}·E_c

**The key asymmetry:** E is multiplied by I_c (the FULL subtree IS poly), while J is multiplied by E_c (only the EXCLUDE-root poly). Since I_c ≥ E_c, the E-product grows faster than the J-product at each step, maintaining J^{(k)} ≤ E^{(k)}.

### J^{(k)} ≤ E^{(k)} holds at all stages

J^{(k)} = ∏ E_{c_i} while E^{(k)} = ∏ I_{c_i}. Since E_{c_i} ≤ I_{c_i} coefficientwise for each i, and products of nonneg sequences preserve coefficientwise dominance, J^{(k)} ≤ E^{(k)} at every stage k. This is automatically induction-friendly.

### The one-factor extension claim

**Theorem (One-Factor Extension).** Suppose:
- (I^{(k)}, E^{(k)}) satisfies Strong Condition C
- E^{(k)} is LC
- I^{(k)} ≥ E^{(k)} coefficientwise
- J^{(k)} ≤ E^{(k)} coefficientwise
- (I_c, E_c) is a tree-realizable triple satisfying Strong Condition C, I_c ≥ E_c, E_c is LC, J_c ≤ E_c

Then (I^{(k+1)}, E^{(k+1)}) satisfies Strong Condition C.

By iterating from k=0 to k=s, this gives Strong Condition C for (A, B) = (∏I_{c_i}, ∏E_{c_i}).

### Why this might be easier than two-factor closure

Instead of proving closure for two arbitrary product triples (which involves double Cauchy products and four-index sums), you prove closure for one accumulated product times one tree-realized factor. The new factor has specific structure:
- deg(I_c) = α(subtree(c)) (bounded degree)
- I_c = (1+x)^{ℓ_c} · A_c + x · B_c (recursive structure)
- The simplest case: I_c = 1+x (leaf factor), E_c = 1

### Suggested attack: start with the leaf case

**Leaf factor case.** When I_c = 1+x, E_c = 1:

    E^{(k+1)} = E^{(k)} · (1+x)
    J^{(k+1)} = J^{(k)}
    I^{(k+1)} = (1+x)·E^{(k)} + x·J^{(k)}

New coefficients: b_k^{new} = b_k + b_{k-1} (convolution with 1+x). The J polynomial is unchanged. Strong Condition C for the new pair should follow from Condition C for the old pair plus the smoothing from (1+x). This is essentially what the key identity captures when ℓ = 1.

**P₃ factor case.** When I_c = 1+3x+x², E_c = 1+2x+x²:

Factor-level: d_0 = 1, d_1 = -1, c_0 = 1, c_1 = 3. Condition C at k=1: 1·(-1) + 2·1 + 1·3 = 4 ≥ 0. Even with d_1 < 0, the curvature bonus compensates.

If you can handle the leaf case and the P₃ case, the general case should follow a similar pattern.

## What I need from you

### Part 1: Prove the leaf case

When adding a leaf factor (I_c = 1+x, E_c = 1):

E^{new} = (1+x)·E^{old}, J^{new} = J^{old}, I^{new} = (1+x)·E^{old} + x·J^{old}.

Express the new d_k^{new}, c_k^{new}, and the Strong Condition C expression in terms of old quantities. Prove (*) for the new pair, using (*) for the old pair plus LC of E^{old}.

### Part 2: Prove the general one-factor extension

For arbitrary tree-realizable (I_c, E_c):

E^{new} = E^{old}·I_c, J^{new} = J^{old}·E_c, I^{new} = E^{old}·I_c + x·J^{old}·E_c.

Express d_k^{new} and c_k^{new} in terms of factor-level and accumulated quantities. The Cauchy product structure means:

    a_k^{new} = Σ_{i+j=k} (b_i + j_{i-1})·(b_j^c + j_{j-1}^c)  ... (messy)

where b^c, j^c are coefficients of E_c, J_c. The key structural constraint: the SAME E_c governs both J^{new} (as a multiplicand) and E^{new} (through I_c = E_c + x·J_c).

Show that the one-factor extension preserves (*), using:
- Factor-level Condition C for (I_c, E_c)
- J_c ≤ E_c (factor-level)
- E^{old} LC + E_c LC ⟹ E^{new} LC (if I_c is LC -- see Part 3)
- The asymmetric growth: E grows by I_c, J grows by E_c ≤ I_c

### Part 3: The LC closure subtlety

E^{(k)} = ∏I_{c_i}. For E^{(k)} to be LC, each I_{c_i} must be LC (product of nonneg LC sequences is LC, but the converse is needed: if a factor isn't LC, the product might not be).

**The issue:** I(T) fails LC at n = 26 (exactly 2 trees). For the induction on n: when proving unimodality for a tree on n vertices, each subtree has < n vertices. By induction, their IS polys are unimodal and satisfy Condition C, but might NOT be LC (for n ≥ 27 with subtrees at n ≥ 26).

**Possible resolutions:**
1. The induction also proves LC as a side result (it IS true for all trees n ≤ 25, and at n = 26 only 2 out of 279M trees fail). Perhaps Condition C implies LC? Unlikely, since LC is strictly stronger.
2. The induction only needs LC of E_c (which always holds), not LC of I_c. Reformulate the one-factor extension to use LC of E_c instead of LC of E^{(k)}.
3. The product ∏I_{c_i} might be LC even when individual I_{c_i} aren't. This requires a separate argument.

**Please analyze which resolution works and adjust the theorem statement accordingly.**

**Note:** E_c (the exclude-root polynomial dp[c][0]) is ALWAYS LC (0 failures, 11.9M checks across all trees n ≤ 20). This is a stronger fact than LC of I_c and might be the right induction hypothesis.

## Computational verification

| Check | Result |
|-------|--------|
| One-factor extension (tree factors) | 0 fails / 2.8M checks |
| Factor-level Condition C | 0 fails / 11.9M checks |
| E_c always LC (factor-level) | 0 fails / 11.9M checks |
| I_c always LC (factor-level) | FAILS at n = 26 (2 trees) |
| J_c ≤ E_c at factor level | 0 fails / 61.8M checks |
| J^{(k)} ≤ E^{(k)} at product level | 0 fails (follows from factor-level) |

## Dead ends (do not pursue)

- HWZZ partial synchronicity: FALSE for tree (I, E) at n = 12+.
- Factor-level ratio dominance (d_k ≥ 0): FALSE in ~31% of factors.
- Global LC of I(T): FALSE at n = 26.
- Weak Condition C: FALSE at n = 17.
- Generic product closure (J≤E + LC + Cond C without tree-realizability): FALSE (2/50K).

## Mechanism insight

When d_k < 0 at the product level, the curvature term T3 = a_{k-1}·c_k handles 100% of compensations (alone in 79.5% of cases). Memory (T2 = b_k·d_{k-1}) alone: 0%. The LC of E (producing c_k ≥ 0) and the amplification a_{k-1} ≥ b_{k-1} (from I ≥ E) are the primary rescue mechanisms.

## Notation summary

| Symbol | Definition |
|--------|-----------|
| I = E + xJ | IS polynomial decomposition at root |
| a_k, b_k, j_k | Coefficients of I, E, J |
| d_k = a_{k+1}b_k - a_kb_{k+1} | LR minors |
| c_k = b_k² - b_{k-1}b_{k+1} | LC gaps of E |
| Strong Cond C | b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0 |
| E^{(k)} = ∏_{i≤k} I_{c_i} | Accumulated exclude-root product |
| J^{(k)} = ∏_{i≤k} E_{c_i} | Accumulated include-root product |
