# Degree-2 Chain Peeling (Path-Attachment Reduction)

Goal: record an exact, proved algebraic reduction for long degree-2 chains.
This is a path-attachment decomposition (not a leaf-star decomposition).

## Setup

Let `T` be a tree. A degree-2 chain is a maximal path

  `v_0 - v_1 - ... - v_s`

where `v_0` is a leaf, `v_1,...,v_{s-1}` have degree `2` in `T`, and `v_s`
has degree `!= 2` (or is the opposite leaf if `T` is a path).

Define `T_0 := T`, and for `i >= 1` define

  `T_i = T - {v_0, v_1, ..., v_{i-1}}`.

Then `T_s` is `T` with `v_0,...,v_{s-1}` removed, and `T_{s+1}=T_s-v_s`.

## One-step identity (s = 1)

If `v_0` is a leaf and `v_1` has degree `2` with other neighbor `w`, then:

  `I(T) = I(T - v_0) + x I(T - N[v_0])`
       `= I(T_1) + x I(T_2)`,

and since `v_1` is a leaf in `T_1`, we can rewrite:

  `I(T) = (1+x) I(T_2) + x I(T_2 - w)`.

This is the base path-attachment step (`s=1`).

## Iterated peeling (exact Fibonacci form)

For every `1 <= i <= s`, vertex `v_{i-1}` is a leaf in `T_{i-1}`, so:

  `I(T_{i-1}) = I(T_i) + x I(T_{i+1})`.      (R)

Define Fibonacci polynomials in `x` by
`F_0=1`, `F_1=1`, `F_t = F_{t-1} + x F_{t-2}` for `t>=2`.

### Theorem (Exact chain decomposition)

For every chain length `s>=1`,

  `I(T_0) = F_s(x) I(T_s) + x F_{s-1}(x) I(T_{s+1})`.

#### Proof

Claim by backward distance from the chain end:
for each `j=1,...,s`,

  `I(T_{s-j}) = F_j I(T_s) + x F_{j-1} I(T_{s+1})`.      (C_j)

Base `j=1`:
`I(T_{s-1}) = I(T_s) + x I(T_{s+1})` by (R), which is `(C_1)` since
`F_1=F_0=1`.

Step `j -> j+1` (`j+1 <= s`):
from (R),
`I(T_{s-(j+1)}) = I(T_{s-j}) + x I(T_{s-(j-1)})`.
Apply `(C_j)` and `(C_{j-1})`:

`I(T_{s-(j+1)})`
`= (F_j + xF_{j-1}) I(T_s) + x(F_{j-1}+xF_{j-2}) I(T_{s+1})`
`= F_{j+1} I(T_s) + x F_j I(T_{s+1})`.

So `(C_{j+1})` holds. By induction, `(C_s)` gives the theorem.

QED.

For `s>=2`, since `I(P_m;x)=F_{m+1}(x)`, this can be rewritten as

  `I(T) = I(P_{s-1};x) I(T_s) + x I(P_{s-2};x) I(T_{s+1})`.

So long degree-2 chains are exactly path-kernel attachments.

## Boundary implication (open)

The leaf-attachment boundary lemma does not apply directly. However:

  - Path polynomials `I(P_m; x)` are PF∞ and log-concave.
  - The decomposition above suggests a “path-attachment” analogue of the
    leaf-attachment lemma: if s is large, the ratio sequence of
    I(P_{s-1}) A + x I(P_{s-2}) B should be tail-monotone.

Two possible routes:

1) Prove tail ratio monotonicity for
   I(P_{s-1};x) A(x) + x I(P_{s-2};x) B(x) when s is large.
2) Treat long degree-2 chains via a finite-kernel check, if a bound on s
   is established by other means.

## Next step

Establish a path-attachment boundary lemma:

  I(T) = I(P_{s-1};x) A(x) + x I(P_{s-2};x) B(x)

with explicit s0(d,e) (depending on deg A, deg B) such that the boundary
indices are nonincreasing. This would yield an explicit bound on degree-2
chain length in any minimal counterexample.
