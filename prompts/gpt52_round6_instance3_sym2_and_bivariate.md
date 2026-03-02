# Task: Sym² representation and bivariate generating function approach

## Context from your Round 5 output

Your semigroup encoding U(I,E,J) = [[I, xJ], [0, E]] with U₁U₂ = U(product) is clean and correct. Your diagnosis of why the naive TN₃ route fails (block-triangular 3×3 minors factor as b_{k-1} × 2×2 minor) is precise. Your suggestion to search for a 3-dimensional representation ρ of the 2×2 upper-triangular semigroup — specifically Sym² — is the most promising direction.

## New results since your Round 5 output

### Factor-level curvature budget FAILS

The "tree-only bound" from Instance 2 (|L_{p,q}| ≤ Σ c_t at factor level) is FALSE starting at n=15 (ratio exceeds 1 at n=19). So any TN₃ proof cannot rely on factor-level total positivity alone — it must use the product structure.

### Product-level SCC holds universally

89.9M intermediate product stages (n≤22), zero failures. Min nontrivial margin = 6. So the closure IS true; the question is purely about finding the right proof framework.

## Two directions to pursue (pick one or both)

### Direction A: Sym² computation

You proposed Sym² as the natural 3-dim representation. Here is the explicit setup.

For U = [[a, c], [0, b]] (upper triangular), Sym² gives the 3×3 matrix:

Sym²(U) = [[a², 2ac, c²], [0, ab, bc], [0, 0, b²]]

(using the standard basis {e₁², e₁e₂, e₂²} for Sym²(ℝ²).)

In our setting: a = I(x), b = E(x), c = xJ(x). So:

Sym²(U) = [[I², 2I·xJ, x²J²], [0, I·E, xJ·E], [0, 0, E²]]

The Toeplitz lift T(Sym²(U)) is a 3-block upper-triangular matrix. Its 3×3 minors involve adjacent Toeplitz entries from the blocks T(I²), T(2IxJ), T(x²J²), T(IE), T(xJE), T(E²).

**Specific question:** Take a 3×3 minor using one row from the T(I²) block, one from the T(IE) block, and one from the T(E²) block, at indices k-1, k, k+1. Compute this minor explicitly. Does it equal (or relate to) the Strong C quantity b_{k-1}·Δ_k or the SCC integer form S_k?

**Alternative basis:** Instead of {e₁², e₁e₂, e₂²}, try {e₁², e₁e₂ + αe₂², e₂²} for some constant α. The off-diagonal entries change, and the 3×3 minors change accordingly. Is there an α that makes a Toeplitz minor equal Δ_k?

**Or:** Try the basis {(1+x)I, something, E} instead of {I², IE, E²}. Since Strong C is about (1+x)I vs E ratio dominance, putting (1+x)I in the top row and E in the bottom might directly produce the right minors.

### Direction B: Bivariate generating function as a matrix

Strong Condition C is equivalent to: the bivariate antisymmetric generating function

F_{(1+x)I, E}(x, y) = E(x)·(1+x)I(y) - (1+x)I(x)·E(y)

has nonneg diagonal coefficients: [x^k y^{k+1}] F ≥ 0 for all k.

Now, under the product operation I → I·P, E → E·Q:

F_{(1+x)IP, EQ}(x,y) = EQ(x)·(1+x)IP(y) - (1+x)IP(x)·EQ(y)

Expand:
= E(x)Q(x)·(1+x)I(y)P(y) - (1+x)I(x)P(x)·E(y)Q(y)

Insert ± E(x)·(1+x)I(y)·P(x)·Q(y):

= [E(x)·(1+x)I(y) - (1+x)I(x)·E(y)] · P(x)Q(y)
  + (1+x)I(y)·E(x) · [Q(x)P(y) - P(x)Q(y)]

Wait, that's not quite right. Let me redo:

= E(x)Q(x)·(1+x)I(y)P(y) - (1+x)I(x)P(x)·E(y)Q(y)

Factor:
= [E(x)·(1+x)I(y)]·[Q(x)P(y)] - [(1+x)I(x)·E(y)]·[P(x)Q(y)]

= F_{(1+x)I, E}(x,y) · Q(x)P(y)  +  (1+x)I(x)·E(y) · [Q(x)P(y) - P(x)Q(y)]

Hmm, still not clean. Use Instance 1's identity (3.5) adapted for (e, b) = ((1+x)I, E):

F_{eP, bQ}(x,y) = F_{e,b}(x,y)·P(x)Q(y)  +  e(y)b(x)·F_{P,Q}(x,y)

where e = (1+x)I, b = E.

The diagonal coefficients of the first term are nonneg (old SCC times Cauchy product of P,Q with nonneg coefficients). The second term e(y)b(x)·F_{P,Q}(x,y) is the correction.

**Key observation:** F_{P,Q}(x,y) = Q(x)P(y) - P(x)Q(y) has diagonal coefficients = L_{p,q} (the long minors of the factor). The weighting by e(y)b(x) means the k-th diagonal coefficient of the correction is:

[x^k y^{k+1}] e(y)b(x)·F_{P,Q}(x,y) = Σ_{u,v} b_u · e_v · [x^{k-u} y^{k+1-v}] F_{P,Q}

This is a weighted sum of factor long minors, with weights from the accumulated product.

**Question:** Can you package this bivariate identity as a MATRIX operation? The bivariate generating function F_{A,B}(x,y) is an element of R[[x,y]]. The product operation acts linearly. Is there a matrix M over R[[x]] such that the diagonal coefficients of F are the entries of M·v for some vector v, and the product operation becomes matrix multiplication? This would connect Direction B to Direction A.

## What I'm hoping for

Either:
1. An explicit Sym² (or twisted Sym²) basis where a 3×3 Toeplitz minor equals Δ_k; or
2. A proof that no 3×3 minor of any polynomial-matrix representation equals Δ_k (with a clean reason why); or
3. A reformulation where the bivariate generating function approach gives a matrix factorization compatible with product closure.

## Dead ends

- Naive T(U) or T(Ũ) with (1+x)I: no 3×3 minor equals Δ_k (you verified this)
- Factor-level total positivity: FALSE at n≥15 (ratio > 1 at n=19)
- HWZZ partial synchronicity: FALSE at n=12+

## Notation

| Symbol | Definition |
|--------|-----------|
| U(I,E,J) | [[I, xJ], [0, E]] (2×2 polynomial matrix) |
| e | (1+x)I (shifted I) |
| Δ_k | e_{k+1}b_k - e_kb_{k+1} (Strong C target) |
| S_k | b_{k-1}·Δ_k = b_{k-1}d_k + b_kd_{k-1} + a_{k-1}c_k (integer form) |
| F_{A,B}(x,y) | B(x)A(y) - A(x)B(y) (bivariate antisymmetric GF) |
| K_{u,v} | p_{u+1}q_v - p_u q_{v+1} (factor cross-minor) |
