# D24 Bencs coefficient extraction and obstruction packet

## Verdict

Bencs' Christoffel--Darboux identity gives an exact induced-subtree
decomposition of the prefix-GSB deficit, but neither of the two natural
positive extractions survives the first unchecked rank.

1. Dropping the explicitly positive outer-pair correction from the diagonal
   energy identity already fails on `K_(1,9)` at prefix rank `r=4`.
2. Polarizing first and assigning the square term by the canonical identity
   `F'=sum_v I(T-N[v])` survives every possible rank-`>=4` prefix summand
   through order `10`, but fails at order `11` on an induced three-vertex
   path whose residual polynomial is `(1+x)^7`.

The second witness is decisive for this lane: the residual polynomial is
real-rooted, log-concave, unimodal, and already decreasing from its
Levit--Mandrescu tail threshold.  Residual tail control cannot make the
individual Bencs summands positive.  A viable continuation would have to
group complete color-switching orbits or otherwise compensate different
induced subtrees globally.

The exact replay is
`scratch_d24_bencs_deficit_20260711.py`.

## 1. Bivariate identity and coefficient extraction

Write

```text
F(x)=I(T;x)=sum_k a_k x^k.
```

For every nonempty connected induced subtree `H`, choose its bipartition
sizes `p>=q`, put `d=p-q`, and write

```text
G_H(x)=I(T-N[H];x)=sum_j b_j(H)x^j.
```

The bivariate form of Bencs' identity is

```text
x F'(x)F(y)-yF(x)F'(y)
 = sum_H d (x^p y^q-x^q y^p)G_H(x)G_H(y).       (CD2)
```

This is Theorem 1.7 of Bencs' primary paper,
[Christoffel--Darboux type identities for independence polynomial](https://arxiv.org/abs/1409.2527),
with the global sign fixed by the one-vertex graph.

Define the contribution of `H` to bidegree `(i,j)` by

```text
C_H(i,j)
 = d[b_(i-p)b_(j-q)-b_(i-q)b_(j-p)].
```

Coefficient extraction from `(CD2)` gives the exact identity

```text
(i-j)a_i a_j = sum_H C_H(i,j).                   (1)
```

Differentiating `(CD2)` in `y`, setting `y=x`, and multiplying by `x`
recovers the diagonal identity from the task:

```text
x^2F'(x)^2-x^2F''(x)F(x)-xF'(x)F(x)
 = -sum_H d^2 x^|H| G_H(x)^2.
```

Equivalently, if

```text
E_m(F)=1/2 sum_(i+j=m)(i-j)^2 a_i a_j,
```

then

```text
E_m(F)=sum_H d^2 [x^(m-|H|)]G_H(x)^2.            (2)
```

The replay checks every coefficient of `(1)` on every nonisomorphic tree
through order `9`: `3,616` exact coefficient equalities.

## 2. Exact GSB polarization

At rank `r`, put

```text
Delta_r
 = (r+1)a_(r+1)^2+a_r a_(r+1)-(r+2)a_r a_(r+2).
```

Taking `(i,j)=(r+1,r)` and `(r+2,r)` in `(1)` gives

```text
Delta_r
 = (r+1)a_(r+1)^2
   + sum_H [C_H(r+1,r)-(r+2)C_H(r+2,r)/2].       (3)
```

Thus the doubled local contribution is the integer

```text
Phi_H(r)=2C_H(r+1,r)-(r+2)C_H(r+2,r),

2Delta_r=2(r+1)a_(r+1)^2+sum_H Phi_H(r).         (4)
```

Formula `(3)` is the requested exact coefficient-level expression.  It is
strictly finer than the diagonal energy identity because it isolates the
adjacent product and the two-step product separately.

## 3. The strongest separated two-diagonal energy core fails

For compactness put

```text
A_t=a_(r-t)a_(r+2+t),
B_t=a_(r-t)a_(r+1+t),
C_m=[x^m]F(x)^2.
```

Then

```text
C_(2r+2)=a_(r+1)^2+2 sum_(t>=0) A_t,
E_(2r+2)=4 sum_(t>=0)(t+1)^2 A_t,

C_(2r+1)=2 sum_(t>=0) B_t,
E_(2r+1)=sum_(t>=0)(2t+1)^2 B_t.
```

Solving exactly for the central products while leaving every omitted outer
term nonnegative gives

```text
Delta_r = K_r
 + sum_(t>=1)[(3r+4)(t+1)^2-2(r+1)]A_t
 + sum_(t>=2)[((2t+1)^2-9)/8]B_t,                (5)
```

where

```text
K_r
 = (r+1)C_(2r+2)+(9/16)C_(2r+1)
   -((3r+4)/4)E_(2r+2)-(1/8)E_(2r+1).           (6)
```

The coefficients in both correction sums in `(5)` are nonnegative.  The
odd coefficient `-1/8` in `(6)` is forced at its largest possible value by
the first outer pair `B_1`; the even energy coefficient is forced by the
coefficient of `A_0`.  This is therefore the coefficientwise strongest core
of this two-diagonal form whose discarded correction is manifestly
nonnegative.

The convolution terms are themselves Bencs terms after adjoining one
isolated vertex: for `J=(1+x)F`, the singleton component gives
`xF(x)^2` in the diagonal identity.  Thus `(6)` is a genuine isolated-vertex
polarization of `(2)`, not an unrelated convolution estimate.

On `T=K_(1,9)`, whose polynomial is

```text
[1,10,36,84,126,126,84,36,9,1],
```

we have `alpha=9`, tail start `L=6`, and `r=4=L-2`.  Exact evaluation gives

```text
K_4                       = -180471,
positive outer correction = 212223,
Delta_4                    =  31752.
```

So the true GSB deficit is positive, but the separated energy core is
negative.  The outer correction is essential even on the star at the first
possible unchecked rank.

The same obstruction is visible directly in `(4)`.  For a leaf singleton
`H`, the residual forest is eight isolated vertices, with

```text
G_H=(1+x)^8=[1,8,28,56,70,56,28,8,1].
```

At `r=4`,

```text
C_H(5,4)=1764,
C_H(6,4)=2352,
Phi_H(4)=2*1764-6*2352=-10584.
```

There are nine such leaf terms.  Their total is `-95256` in the doubled
normalization, exactly compensated by the doubled square term `158760`.

## 4. Canonical square-term charging and its smallest-order failure

The derivative identity

```text
F'(x)=sum_(v in V(T)) I(T-N[v];x)
```

implies

```text
(r+1)a_(r+1)^2
 = a_(r+1) sum_v [x^r]I(T-N[v];x).               (7)
```

This gives a canonical attempt to repair the negative singleton terms in
`(4)`: assign the `v`-summand of `(7)` to `H={v}`.  In the doubled
normalization define

```text
Psi_H(r)=Phi_H(r)+2a_(r+1)[x^r]G_H
```

for singleton `H`, and `Psi_H=Phi_H` otherwise.  Then

```text
2Delta_r=sum_H Psi_H(r).                          (8)
```

This charging repairs every term of the star witness.  Exact enumeration
finds all `521` charged summands nonnegative at every still-open prefix rank
`r>=4` on every tree through order `10`.

It fails at order `11`.  Take the tree with edges

```text
01, 09, 9-10, 12, 13, 14, 15, 16, 17, 18.
```

Equivalently, start with a seven-leaf star and extend its remaining arm by
two vertices.  Its independence polynomial is

```text
[1,11,45,105,161,161,105,43,10,1].
```

Again `alpha=9`, `L=6`, and `r=4=L-2`.  Let `H` be the induced path on
vertices `0,9,10`.  Then `(p,q,d)=(2,1,1)`, while deleting `N[H]` leaves
seven isolated vertices:

```text
G_H=(1+x)^7=[1,7,21,35,35,21,7,1].
```

The exact polarized terms are

```text
C_H(5,4)=490,
C_H(6,4)=784,
Psi_H(4)=Phi_H(4)=2*490-6*784=-3724.             (9)
```

Because `H` is not a singleton, it receives no charge from `(7)`.  The whole
tree is nevertheless far inside GSB:

```text
Delta_4=54096>0.
```

The exhaustive order-`<=10` replay plus this order-`11` witness makes `(9)`
the smallest-order obstruction to the canonical charged-summand statement
at the still-open ranks `r>=4`.

## 5. Why the residual tail cannot repair the local statement

For the witness in `(9)`, the residual independence number is `7`, so its
Levit--Mandrescu tail starts at

```text
ceil((2*7-1)/3)=5.
```

The residual coefficients used in `(9)` satisfy all the strongest familiar
shape properties:

```text
G_H=(1+x)^7
```

is a Poisson-binomial generating polynomial, is real-rooted and log-concave,
and has `b_5=21<b_4=35` at the certified tail boundary.  Nevertheless,

```text
b_3^2-b_2b_4 = 490,
b_3b_4-b_2b_5 = 784,
490-3*784    = -1862.
```

Hence no hypothesis depending only on residual unimodality, log-concavity,
real-rootedness, or Levit--Mandrescu tail monotonicity can imply
`Phi_H(r)>=0` or the charged version `Psi_H(r)>=0`.

## 6. Surviving interpretation

Bencs' proof groups ordered pairs of independent sets by connected
components of their symmetric difference.  Switching the two color classes
of a component changes the rank difference by its imbalance `d=p-q`.
The diagonal identity is positive only after summing a complete switching
orbit; the witnesses above show that extracting one marked component destroys
that positivity.

Therefore this lane does not supply a branchwise positive proof of prefix
GSB.  The only materially different Bencs continuation is an orbit-level
inequality that retains all component switches simultaneously and couples
the even-total and odd-total orbits appearing in `Delta_r`.  Such a statement
would be new theorem-strength work; it is not furnished by the existing CD
identity or by the residual tail theorem.

## 7. Replay

Run

```bash
python3 scratch_d24_bencs_deficit_20260711.py
```

The script:

- verifies `3,616` exact bivariate coefficient equalities on all trees
  through order `9`;
- checks all `521` canonical charged summands at every possible open prefix
  rank through order `10`;
- reconstructs the star and broom trees from edge lists;
- checks both full independence polynomials, independence numbers, prefix
  boundaries, residual polynomials, local CD terms, energy-core values, and
  global GSB deficits using integer or rational arithmetic.

The final output ends with `certificate: passed`.
