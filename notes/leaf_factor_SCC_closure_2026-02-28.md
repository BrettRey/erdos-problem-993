# Leaf-Factor SCC Closure Analysis (2026-02-28)

## Setup

Attaching a leaf to root r: E' = (1+x)E, J' = J. So:
- a'_k = a_k + b_{k-1}
- b'_k = b_k + b_{k-1}

## Key Identities

**Cross-determinant d'_k:**
```
d'_k = d_k + e_k    where e_k = b_{k-1}·j_k - b_k·j_{k-1}
```
This is a 2×2 minor of the (b, j) matrix: e_k = det(b_{k-1}, b_k; j_{k-1}, j_k).

**LC gap c'_k:**
```
c'_k = c_k + c_{k-1} + γ_k    where γ_k = b_k·b_{k-1} - b_{k-2}·b_{k+1}
```
γ_k ≥ 0 when E is LC (by ratio monotonicity of LC sequences).

**SCC decomposition S'_k:**
```
S'_k = S_k + (α_k + β'_k) + β_k + γ_{1,k} + γ_{2,k}
```
where:
- S_k = b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k  [original SCC, ≥ 0]
- α_k + β'_k = b_{k-2}·d'_k + b_{k-1}·d'_{k-1} + a_{k-1}·c_{k-1}  [shifted SCC-like, ≥ 0 in 99.99% of cases]
- β_k = b_{k-1}·e_k + b_k·e_{k-1}  [P2-dependent, ≤ 0 in 95.3% of cases]
- γ_{1,k} = b_{k-2}·(c_k + c_{k-1})  [≥ 0 always]
- γ_{2,k} = a'_{k-1}·γ_k  [≥ 0 when E is LC]

## Computational Verification (n ≤ 15)

| Metric | Value |
|--------|-------|
| Total nontrivial checks | 523,395 |
| β < 0 | 498,979 (95.3%) |
| α+β' < 0 | 60 (0.01%) |
| Tightest ratio (nonneg / |neg|) | **3.0** (at P_3, root=2, k=1) |
| S'_k < 0 | **0** |

## Significance

1. The leaf-factor case has ENORMOUS margin (ratio ≥ 3). The S_k + γ₁ + γ₂ terms
   always dominate the negative β term by at least 3×.

2. The β_k term is almost always negative (95.3%), meaning P2 holds broadly.
   This is consistent with the computational verification of P2 at support vertices.

3. The α+β' group is almost always nonneg (99.99%). The 60 exceptions are tiny
   and easily absorbed by other terms.

4. For an algebraic proof, the cleanest route would be:
   - Show S_k + γ_{2,k} ≥ |β_k| (since γ_{2,k} provides the curvature bonus)
   - Or show α+β'+γ₁ ≥ 0 always, reducing to S_k ≥ |β_k| - γ₂

5. The d'_k = d_k + e_k identity is beautiful and may extend to general factors
   via convolution (Cauchy products of e_k-like cross-determinants).
