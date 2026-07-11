# D14 aggregate-energy obstruction packet

## Frozen statement

The proposed D14 energy estimate

```text
Dir_r(e) <= E_r(e+h)/r                                      (ED)
```

is false in the required prefix window, even for a tree.  An exact
counterexample has order `210`, independence number `124`, and prefix rank
`r=81`.

This blocks the two-lemma D14 route

```text
gap(P_r) >= 1/(2r),
Dir_r(e) <= E_r(e+h)/r,
```

regardless of whether the spectral-gap statement is true.  It does not
refute prefix GSB, the spectral-gap statement by itself, or the unimodality
conjecture.

The replayable exact certificate is
`scratch_d14_ed_obstruction_certificate_20260711.py`.

## Tree and prefix certificate

Let

```text
b=(2,3,2,1,2,1,1).
```

Construct two identical rooted trees and join their roots.  At distance
`j=0,...,6` from either root, every vertex has `b_j` children.  The level
sizes on each side are

```text
1, 2, 6, 12, 12, 24, 24, 24,
```

so each side has `105` vertices and the joined tree has `n=210` vertices and
`209` edges.

For one rooted side, let `I_j,O_j` be the maximum independent-set sizes in
the subtree at level `j`, conditional on its root being included or omitted.
Starting with `(I_7,O_7)=(1,0)` and working upward,

```text
I_j=1+b_j O_(j+1),
O_j=b_j max(I_(j+1),O_(j+1)).
```

This gives `(I_0,O_0)=(61,62)`.  Since the two central roots are adjacent,

```text
alpha=max(2O_0,I_0+O_0)=124.
```

The prefix boundary is therefore

```text
L=floor((2alpha+1)/3)=83,
```

and `r=81=L-2` is the last rank required by the prefix-GSB route.

## Exact normalization

For `C in I_(r-1)`, let `H_C` be the residual forest induced by the vertices
addable to `C`, and write

```text
q_C=|V(H_C)|,
h_C=|E(H_C)|,
S2_C=sum_(v in V(H_C)) deg_(H_C)(v)^2.
```

For the uniform down--up chain on `I_r`,

```text
r i_r Dir_r(e)
  = D
  = sum_(C in I_(r-1)) (S2_C-4h_C^2/q_C).
```

Also put

```text
E=sum_(A in I_r)(e(A)+h(A))
  = (r+1)i_(r+1)+sum_(A in I_r)h(A).
```

Since `i_r E_r(e+h)=E`, estimate `(ED)` is exactly the assertion `D<=E`.
Equivalently, the aggregate deficit used by the marked DP is

```text
F=r(E-D)
 = sum_(C in I_(r-1))
   [q_C(q_C+h_C-1)-2h_C-(r+1)S2_C+4r h_C^2/q_C].
```

Thus a negative exact value of `F` is a direct obstruction to `(ED)`.

## Exact values

At `r=81`, the independent-set polynomial and the marked residual-forest DP
give

```text
i_81 = 12150264177683772579409779766331299105518,
i_82 = 4933589294165997470864362459890912420846,

82 i_82
 = 404554322121611792610877721711054818509372,

sum_(A in I_81) h(A)
 = 130702262516384679341617224262464001719820,

E
 = 535256584637996471952494945973518820229192.
```

The Dirichlet numerator is

```text
D =
245700013545993183538657579915009525438626090145082199405497542315007719
------------------------------------------------------------------------
450170153991699740612130417300.
```

The exact deficit is strictly negative:

```text
F = 81(E-D)
  =
-42691269729915413999161105330329696919890850296588501772067055418675071
-------------------------------------------------------------------------
50018905999077748956903379700
  < 0.
```

Equivalently,

```text
D/E =
245700013545993183538657579915009525438626090145082199405497542315007719
------------------------------------------------------------------------
240956539131558137538750790433861781336415995667683476986378980601821600
    = 1.0196860165386306... > 1.
```

The numerator of `D/E-1` is the positive integer

```text
4743474414435045999906789481147744102210094477398722419118561713186119.
```

## Replay and independent audits

Run

```bash
python3 scratch_d14_ed_obstruction_certificate_20260711.py
```

The script constructs the tree, recomputes its full independence polynomial,
computes the exact marked distribution through rank `81`, checks all rank
counts against the polynomial, and asserts every integer and rational value
displayed above.  Its final line contains `certificate: passed`.

The marked recurrence tracks, for each `(rank,q)`, the number of independent
sets and the sums of `h,h^2,S2`.  Its state partition, parent-environment
bits, moment products, root combination, and spherical-tree construction
were audited directly.  In addition, the complete marked distribution was
compared with independent subset enumeration on `120` random labeled trees,
with `12` trees at every order `1<=n<=10`; every entry agreed exactly.

As a non-certifying numerical cross-check, the separate scalar-marker and
Gauss-quadrature recurrence gives `D/E=1.0196860165...` for the same tree and
rank.

## Scope of the obstruction

This tree is not close to violating prefix GSB.  The exact moment identity

```text
E[e(e-1)] = 83*82*i_83/i_81 + 2 eta_81
```

together with the displayed coefficient and residual-edge totals gives

```text
Var_81(e)/(2(mu_81+eta_81))
 =
4456090562654887266294141184427378259936953170060911002599884558253076453306371
--------------------------------------------------------------------------------
20844579827552598787630694022550327446979294710701437563672195144185711388948338
 < 1.
```

Numerically, the same exact calculation is

```text
mu_81       = 33.2959280725...,
eta_81      = 10.7571539684...,
Var_81(e)   = 18.8350664549...,
Var_81(e)/(2(mu_81+eta_81)) = 0.2137769434....
```

The obstruction therefore isolates the loss in the D14 factorization: the
universal energy estimate `(ED)` is too strong, while the desired GSB
inequality has substantial slack on the counterexample.  Reopening D14 would
require a materially different, function-specific coupling of the
Dirichlet and Poincare steps, not another proof attempt for `(ED)`.
