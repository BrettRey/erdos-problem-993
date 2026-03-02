# GPT 5.2 Round 6 Triage (2026-02-28)

## Instance 1: Bivariate Cancellation (MAJOR ADVANCE)

### Key identity (★)

Using the tree factor structure P = Q + xR (with R ≤ Q), the bivariate
identity F_{eP,bQ} simplifies via a cancellation: the -e(y)b(x)xR(x)Q(y)
term from F_{e,b}·xRQ cancels with +e(y)b(x)xR(x)Q(y) from the correction
term e(y)b(x)·F_{P,Q}.

Result: **exact two-term formula for the product SCC**:

    Δ⁺_k = Δ_k(e*Q, b*Q) + T_k

where:
- Δ_k(e*Q, b*Q) ≥ 0 by TP2 closure (SCC of (e,b) + PF2 of Q)
- T_k = (e*R)_k · (b*Q)_k − (e*R)_{k-1} · (b*Q)_{k+1}

### Significance

**Eliminates off-diagonal long minors entirely.** No K_{u,v} appears, no factor-level
curvature budget needed. Everything reduces to 1D convolutions with Q and R.

This explains WHY the factor-level curvature budget can fail while product SCC holds:
the "huge off-diagonal minors" are artifacts of an inferior decomposition.

### Reduction to Y_k monotonicity

Define Y_k = (e*R)_k / (b*Q)_{k+1}. Then T_k ≥ 0 iff Y_k is nondecreasing.

**Product-level bridge lemma** (the key remaining conjecture):
If e ≥_{lr} b (accumulated SCC), Q is PF2, P = Q + xR with R ≤ Q,
then Y_k = (e*R)_k / (b*Q)_{k+1} is nondecreasing in k.

If true → T_k ≥ 0 → Δ⁺_k ≥ 0 → SCC closed under one-factor extension.

### Weak form

Even without Y_k monotonicity, it suffices to show:
T_k ≥ −Δ_k(e*Q, b*Q)

This is a product-level comparison: both sides built from e,b,Q,R via 1D convolutions only.

## Instance 2: Weighted Span Telescoping

### Why factor-level curvature budget fails

The Toeplitz minor τ_{p,q} telescopes as:
τ_{p,q} = Σ_{t=p}^{q-1} c_t · (b_{p-1}·b_{q-1})/(b_{t-1}·b_t)

The weights b_{p-1}·b_{q-1}/(b_{t-1}·b_t) can exceed 1 in the tail regime.
This is exactly why |L_{8,10}| = 192 > Σc = 180 at n=19: the tail weights > 1.

So the "correct" span budget uses weighted curvatures, not unweighted sums.

### Product-level recursion telescopes span

Each application of the M_k^{(t)} recursion:
    b_{k-1}·M_k^{(t)} = b_k·M_{k-1}^{(t-1)} + a_{k-t}·c_k

reduces span by 1 and harvests one curvature term with multiplicative b-ratio prefactor.
For span-s objects, s reduction steps needed, accumulating a chain of b-ratio weights.

The Cauchy-Binet form of c_k^{prod} is organized in exactly these weighted "span currencies"
τ^{(i)}_{u,v}. This is the clean reconciliation: the recursion supplies weighted telescoping,
and Cauchy-Binet supplies the correct weighted span-resolved curvature budget.

### Term-by-term pairing

The (u,v) term in d_k^{prod} decomposes as:
    L^{(1)}_{u,v} · [τ^{(2)}_{k;u,v} + (b,j)-minor]

The extra (b,j)-minor is the discrepancy between using a^{(2)} vs b^{(2)} in the weight.
Domination by the single CB term requires controlling this discrepancy, but the n=19
violation shows this can't be done with naive unit weights.

The "memory + curvature" rescue pattern: most spans handled at depth 0 (curvature only),
remaining need depth 1 (memory + curvature).

## Instance 3: Sym² Dead, Bivariate Semigroup Alive

### Sym² is dead

In standard Sym² basis, any 3×3 minor using one row from each diagonal block
factors as a product of three Toeplitz entries (block upper-triangular structure).
**Cannot encode Δ_k.** No basis twist α helps for diagonal-block minors.

### Bivariate semigroup IS the right framework

The bivariate GF gives a genuine 2×2 matrix semigroup:

    [F_{AP,BQ}]   [P(y)Q(x)    F_{P,Q}(x,y)] [F_{A,B}]
    [G_{AP,BQ}] = [0           P(x)Q(y)     ] [G_{A,B}]

where G_{A,B}(x,y) = A(x)B(y). This M(P,Q) matrix is multiplicative:
M(P₂,Q₂)·M(P₁,Q₁) = M(P₁P₂, Q₁Q₂).

**SCC preservation = cone invariance:** show that M(P,Q) preserves the cone
of pairs (F,G) where superdiagonal coefficients of F are nonneg and G has nonneg
coefficients.

This is more faithful than Sym² because the obstruction is intrinsically
antisymmetric/wedge-like (F = B(x)A(y) − A(x)B(y)), not quadratic.

## Convergence between Instance 1 and Instance 3

Instance 1's clean cancellation + Instance 3's semigroup framework give the same picture:
- The "multiplicative" part (G evolving by P(x)Q(y)) is automatically nonneg
- The "antisymmetric" part (F evolving) has a correction F_{P,Q} that, for tree
  factors P = Q + xR, simplifies to the T_k obstruction
- The first-term Δ_k(e*Q, b*Q) ≥ 0 is the TP2 closure property of M(Q,Q)
- The obstruction T_k is the (P,Q)≠(Q,Q) asymmetry contribution

## Relation to my T1+T3 ≥ 0 finding

My computational finding (22M events, 0 failures) that T1+T3 = b_{k-1}·d_k + α_{k-1}·c_k ≥ 0
at ALL incremental product stages is CONSISTENT with Instance 1's formula.

The decomposition T1+T3 = c_k·(α_{k-1}+j_{k-1}) - b_k·r_k (verified algebraically,
0 identity failures) with ratio c_k·(α+j)/(b_k·r_k) → 2 from above.

The Instance 1 formula operates at the (e,b) = ((1+x)I, E) level, while my scan
checks (I_acc, E_acc) directly. The connection needs careful index tracking but both
capture the same phenomenon: product-level curvature dominates ratio-dominance deficits.

## Action items

1. **Verify Instance 1's ★ identity computationally** (highest priority)
2. **Check Y_k monotonicity** computationally at all intermediate stages
3. **Probe the bridge lemma** at larger n: is Y_k nondecreasing always?
4. **Write T_k ≥ 0 verification** into the scan infrastructure
5. Record Instance 2's weighted telescope insight (useful for paper exposition)

## Key dead ends confirmed

- Sym² Toeplitz minors: DEAD (block-triangular forces factoring)
- Factor-level curvature budget: DEAD (weighted tail ratios > 1)
- HWZZ partial synchronicity: already known DEAD
