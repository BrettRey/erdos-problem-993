# Task: Prove product closure of Strong Condition C via Cauchy-Binet curvature domination

## The problem

Prove that the independence polynomial of every tree is unimodal.

The proof reduces to a single algebraic claim: **Strong Condition C is preserved under Cauchy products of tree-realizable factor pairs.** This is verified computationally (0 failures, 2.8M products) but not proved.

Your task is to prove product closure by directly analyzing the Cauchy product structure of the LR minors and LC gaps, and showing that the curvature bonus dominates the negative cross-terms.

## Reduction to Condition C at support vertices

Every tree has a support vertex r (vertex adjacent to at least one leaf). At r with ℓ leaf children and non-leaf children c₁,...,c_s:

    E(x) = (1+x)^ℓ · A(x),    J(x) = B(x),
    I(T; x) = E(x) + x·J(x)

where A = ∏ I_{c_i} and B = ∏ E_{c_i}, products over the non-leaf children.

- P3 (tail domination): **PROVED** at all support vertices (leaf-swap injection).
- P2 (prefix ratio dominance): Reduces to Strong Condition C for (A, B).

## The key identity (PROVED)

Define:
- a_k = [x^k]I,  b_k = [x^k]E  (coefficients of the IS poly and exclude-root poly)
- d_k = a_{k+1}·b_k - a_k·b_{k+1}  (LR minors)
- c_k = b_k² - b_{k-1}·b_{k+1}  (LC gaps of E)
- Δ_k = (a_{k+1}+a_k)·b_k - (a_k+a_{k-1})·b_{k+1}  (P2 target for ℓ=1)

Then:

    b_{k-1}·Δ_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k

The three terms: T1 = b_{k-1}·d_k (current minor), T2 = b_k·d_{k-1} (memory), T3 = a_{k-1}·c_k (curvature bonus).

## Strong Condition C

**Definition.** A pair (I, E) with I = E + x·J satisfies Strong Condition C if for all k ≥ 1:

    b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0

where a_k = [x^k]I, b_k = [x^k]E, d_k = a_{k+1}·b_k - a_k·b_{k+1}, c_k = b_k² - b_{k-1}·b_{k+1}.

By the identity, Strong Condition C implies Δ_k ≥ 0, hence P2.

**WARNING:** A WEAKER version (replacing a_{k-1} with b_{k-1}) FAILS at n=17. The amplification a_{k-1} ≥ b_{k-1} (from I ≥ E coefficientwise) is ESSENTIAL.

## Tree-realizable factor triples

A triple (I, E, J) is **tree-realizable** if there exists a rooted tree T with root r such that:
- E = ∏_{c child of r} I(subtree(c); x)
- J = ∏_{c child of r} E_{subtree(c)}
- I = E + x·J

The **same children** c govern both products (the "shared-factor" constraint).

**Properties (all verified computationally, 0 failures):**
1. I = E + x·J,  E[0] = 1
2. J ≤ E coefficientwise (61.8M checks)
3. E is log-concave (11.9M checks)
4. Strong Condition C holds (11.9M factor checks + 59.9M root-level checks)

## The product operation

Given two triples (I₁, E₁, J₁) and (I₂, E₂, J₂):
- E' = E₁ · E₂  (Cauchy product)
- J' = J₁·E₂ + E₁·J₂ + x·J₁·J₂
- I' = E' + x·J' = I₁·I₂

**Product closure claim:** If both factor triples satisfy properties 1-4, then (I', E', J') satisfies Strong Condition C.

Verified: 0 failures / 2.8M tree-derived products.

**Key negative result:** Generic polynomial triples with J≤E + LC + Cond C can FAIL product closure (2 failures / 50K synthetic pairs). Tree-realizability provides essential additional structure. Your proof must use it.

## Mechanism analysis (your proof guide)

We tested 2,001,000 pairwise products from tree factors (n ≤ 13), examining 25.99M (pair, k) checks.

When d_k^{prod} < 0 (1,939,383 events out of 25.99M):

| Compensator | Count | Fraction |
|-------------|-------|----------|
| T3 (curvature a_{k-1}·c_k) alone | 1,542,631 | 79.5% |
| Both T2 + T3 | 396,752 | 20.5% |
| T2 (memory b_k·d_{k-1}) alone | 0 | 0.0% |

**Curvature (T3) participates in 100% of rescues.** Memory alone never suffices.

Curvature fraction of positive load: min = 0.472, median = 1.000, p10 = 0.969.

The tightest margin: T1 = -7552, T2 = -224256, T3 = +319600, total = +87792 (margin 0.159). Margins shrink with subtree size but remain positive.

d_k < 0 is concentrated in the tail (k ≥ mode), driven by **off-diagonal** Cauchy-Binet cross-terms. The diagonal contribution is always nonneg.

**What this tells you:** The proof should focus on bounding the off-diagonal negativity and showing that the curvature bonus c_k^{prod} · a_{k-1}^{prod} absorbs it. The a_{k-1}/b_{k-1} amplification is essential.

## What I need from you

### Step 1: Decompose product-level LR minors

Write d_k^{prod} = a_{k+1}^{prod} · b_k^{prod} - a_k^{prod} · b_{k+1}^{prod} in terms of factor-level quantities.

The Cauchy product gives a_k^{prod} = Σ_{i+j=k} a_i^{(1)} a_j^{(2)} and similarly for b_k^{prod}. So:

    d_k^{prod} = Σ_{i₁+j₁=k+1, i₂+j₂=k} a_{i₁}^{(1)}a_{j₁}^{(2)} b_{i₂}^{(1)}b_{j₂}^{(2)}
               - Σ_{i₁+j₁=k, i₂+j₂=k+1} a_{i₁}^{(1)}a_{j₁}^{(2)} b_{i₂}^{(1)}b_{j₂}^{(2)}

Separate the **diagonal terms** (i₁ = i₂, j₁ = j₂ + 1 or vice versa) from the **off-diagonal cross-terms**. Show that:
- Diagonal terms are nonneg (they reduce to factor-level d_k terms)
- Off-diagonal terms can be negative but are bounded

### Step 2: Bound off-diagonal negativity

Use the factor-level properties to bound the negative off-diagonal terms:
- J_c ≤ E_c coefficientwise implies a_k^{(c)} - b_k^{(c)} = j_{k-1}^{(c)} ≤ b_{k-1}^{(c)}, so the "excess" is bounded
- Factor-level Strong Condition C gives: when d_k^{(c)} < 0, the curvature bonus at the factor level compensates
- Factor-level E is LC, so c_k^{(c)} ≥ 0

### Step 3: Show curvature dominates

The product-level curvature bonus is T3 = a_{k-1}^{prod} · c_k^{prod}.

- c_k^{prod} = (b_k^{prod})² - b_{k-1}^{prod} · b_{k+1}^{prod}. Since E₁ and E₂ are LC with nonneg coefficients, E₁·E₂ is also LC (Newton's inequality / Cauchy product preserves LC for nonneg sequences). So c_k^{prod} ≥ 0.

- a_{k-1}^{prod} ≥ b_{k-1}^{prod} (product of I ≥ E coefficientwise).

The question is whether a_{k-1}^{prod} · c_k^{prod} is large enough to absorb the total negative contribution from T1 + T2. The mechanism data says yes (margin ≥ 0.159), but you need to prove it.

### Step 4: Assemble the proof

Combine steps 1-3 to show:

    b_{k-1}^{prod} · d_k^{prod} + b_k^{prod} · d_{k-1}^{prod} + a_{k-1}^{prod} · c_k^{prod} ≥ 0

### Useful tools

**Binomial-minor expansion** (from a collaborator): Define shifted minors M_k^{(t)} = a_{k+1-t}·b_k - a_{k-t}·b_{k+1}. Then:
- M_k^{(0)} = d_k
- b_{k-1} · M_k^{(t)} = b_k · M_{k-1}^{(t-1)} + a_{k-t} · c_k
- Strong Condition C = Σ_{t=0}^{1} C(1,t) · M_k^{(t)}... wait, actually the binomial expansion is:
  Δ_k^{(ℓ)} = Σ_{t=0}^ℓ C(ℓ,t) · M_k^{(t)} for the (1+x)^ℓ case.

The shifted minor recursion b_{k-1} · M_k^{(t)} = b_k · M_{k-1}^{(t-1)} + a_{k-t} · c_k shows that each shift "harvests" one LC gap c_k. This might be useful for bounding the negative excursion.

## Constraints (hard boundaries)

- Do NOT assume d_k ≥ 0 at factor level. Fails in ~31% of factors.
- Do NOT assume global LC of I(T). Fails at n = 26.
- Do NOT assume J ≤ E is product-closed. It is NOT.
- Do NOT assume HWZZ partial synchronicity. FALSE for tree (I, E) pairs at n = 12+.
- The WEAK version of Condition C (with b_{k-1} instead of a_{k-1}) FAILS at n = 17. The amplification by a_{k-1} is essential.
- Generic pairs with J≤E + LC + Cond C can fail product closure. Tree-realizability is essential.

## Notation summary

| Symbol | Definition |
|--------|-----------|
| a_k = [x^k]I | Coefficients of I = E + xJ |
| b_k = [x^k]E | Coefficients of E (exclude-root poly) |
| j_k = [x^k]J | Coefficients of J (include-root poly / x) |
| d_k = a_{k+1}·b_k - a_k·b_{k+1} | LR minors of I vs E |
| c_k = b_k² - b_{k-1}·b_{k+1} | LC gaps of E |
| Δ_k = e_{k+1}·j_k - e_k·j_{k+1} | P2 target (e = coeffs of (1+x)^ℓ A) |
| Strong Cond C | b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0 |
