# D22 global switching: exact residual formula and obstruction packet

## Result

The Bencs symmetric-difference switch gives a clean connected-component
expansion for the nonedge and far covariance sums.  It also shows exactly why
the two most tempting global sign shortcuts do not close the proof:

1. `2 q2 + qfar <= 0` is false in the prefix, first at order 16 after an
   exhaustive scan through order 15.
2. `q2 + qfar <= 0` is false in the prefix.  A compact exact witness is
   Galvin's `T(2,23,1)` at order 95 and rank four.
3. Even after diagonal and edge-product collision capacity is restored, that
   capacity cannot be assigned independently to each connected
   symmetric-difference component.  A prefix component of an order-11 double
   star has residual `525>0`.

Thus a successful switching proof has to transport capacity between
components (or retain common vertices of the two target sets).  A local path
sign or a componentwise color-imbalance lemma cannot prove prefix GSB.

The exact GSB budget remains valid on every obstruction in this packet.

## Covariance normalization and exact budgets

For a uniform independent `r`-set, let `X_v` indicate that `v` is addable,
and put

```text
a_v = number of rank-r sets to which v is addable,
N   = i_r(T),
q_uv = N * number of rank-r sets to which u and v are jointly addable
       - a_u a_v.
```

Write `q2` for the sum over distance-two pairs and `qfar` for the sum over
pairs at distance at least three.  If

```text
S0 = sum_v a_v^2,
S1 = sum_(uv in E(T)) a_u a_v,
d1 = sum_v a_v = (r+1)i_(r+1),
```

then the exact ordered-LC and GSB inequalities are respectively

```text
2(q2+qfar) <= S0+2S1,
2(q2+qfar) <= S0+2S1+N*d1.
```

This follows by expanding `Var(sum_v X_v)` and using the exact edge-joint
term.  It is also the marked-pair count: the left side counts nonadjacent
endpoint marks, `S0+2S1` supplies equal/adjacent endpoint collisions, and
`N*d1` is the one-mark auxiliary GSB capacity.

## Connected-component Bencs residual

Let `H` be a connected induced subtree with bipartition `(P,Q)`, put
`p=|P|`, `q=|Q|`, and let

```text
F_H = T-N[H],       f_j(H)=i_j(F_H).
```

For a set `W` of allowed nonedge endpoint pairs, write

```text
s_P(H) = number of pairs in C(P,2) that lie in W,
s_Q(H) = number of pairs in C(Q,2) that lie in W,
c(H)   = number of nonedge pairs in P x Q that lie in W.
```

Switching every symmetric-difference component that contains exactly one
marked endpoint cancels all pairs whose endpoint marks lie in distinct
components.  What remains from `H` is

```text
R_H(r) =
    s_P(H) f_(r-q)(H) f_(r+2-p)(H)
  + s_Q(H) f_(r-p)(H) f_(r+2-q)(H)
  - c(H)   f_(r+1-p)(H) f_(r+1-q)(H).
```

Pairs whose two marked vertices are common to the two balanced target sets
contribute only negatively and are not present in `R_H`.  Consequently

```text
q_W <= sum_(connected induced H) R_H(r).             (Bencs upper bound)
```

The replay script audits this inequality, from the definitions and without
using the path-cavity formula, for both all nonedges and far pairs on every
tree through order eight at every rank.

For far pairs, treehood gives particularly transparent local counts:

```text
c(H) = pq-|E(H)| = (p-1)(q-1),
s_P(H) = C(p,2)-sum_(v in Q) C(deg_H(v),2),
s_Q(H) = C(q,2)-sum_(v in P) C(deg_H(v),2).
```

The subtracted terms are exactly distance-two forks.  Balanced and
near-balanced colors have ample unweighted cross-pair capacity.  The
remaining obstruction is not pair geometry but the three differently
shifted products of remainder coefficients in `R_H(r)`.  Comparing those
products branch by branch would require precisely the rooted polarization
already falsified on `T(14,8,1)`.

## Smallest zero-capacity obstruction

At order 16 take the tree with adjacency list

```text
[[13,15], [13], [13], [13], [14,15], [14], [14], [14],
 [14], [15], [15], [15], [15], [0,1,2,3],
 [4,5,6,7,8], [0,4,9,10,11,12]].
```

It has `alpha=13`, `L=9`, and rank `r=4` is in the prefix.  Exact cavity
counts give

```text
N=905,       d1=7635,
q2=316868,   qfar=-629342,
2q2+qfar=4394>0.
```

The exhaustive scan checked all `13,187` nonisomorphic trees through order
15 before reaching this witness after 147 trees at order 16.  This refutes
only the zero-capacity shortcut, not ordered LC or GSB.

## Prefix nonedge-sign obstruction at order 95

For Galvin's spherical tree `T(m,t,1)`, the root has `m` children, each child
has `t` children, and each grandchild has one leaf.  Take

```text
(m,t,r)=(2,23,4),
n=95, alpha=48, L=32.
```

The exact data are

```text
N=2,835,141,
d1=240,590,350,
q2=62,005,580,307,264,
qfar=-61,454,066,298,966,
q2+qfar=551,514,008,298>0.
```

Thus aggregate negative covariance over all nonedges is already false at a
very early prefix rank.  The exact collision capacities are much larger:

```text
S0=615,288,884,289,486,
S1=404,840,482,368,080,
ordered gap = S0+2S1-2(q2+qfar)
            = 1,423,866,821,009,050>0,
GSB gap     = ordered gap+N*d1
            = 2,105,974,386,498,400>0.
```

The certificate derives `N,d1,S0,S1,q2+qfar` independently from the four
vertex orbits and then agrees with the generic directed-cavity computation.
An exhaustive scan of all `T(m,t,1)` instances of order below 95, over every
prefix rank `r>=4` through rank 80, found no smaller witness in this family;
this is not a claim of global minimality.

## Componentwise collision-capacity obstruction

The graph6 tree

```text
J??????wC~?
```

has order 11 and independence number 9.  At its last required prefix rank
`r=4`, choose the induced `K_(1,3)` on vertices `{0,1,2,9}` with the
three-vertex color as `P`.  Its closed-neighborhood remainder is six isolated
vertices, so

```text
I(F_H;x)=(1+x)^6,
(p,q)=(3,1).
```

If the equal/adjacent target capacity is assigned only inside this connected
component, its residual is

```text
C(3,2) f_3 f_3 - 3*1 f_2 f_4
=3*20^2-3*15^2
=525>0.
```

Therefore diagonal/edge products do not close the Bencs expansion one
component at a time, even in the exact prefix.  The negative terms omitted
from the connected residual--common marks and switches joining different
components--are logically essential.

## Route decision

The switching expansion is useful localization, but its remaining global
transport problem is the same hard content seen in the earlier D7 Hall
network:

* even-distance endpoint components with color imbalance at least two need a
  compensating component;
* distance-two fork collisions need diagonal/edge capacity;
* the componentwise allocation of that capacity is false;
* arbitrary remainder coefficient products cannot be polarized locally.

Freeze the bare nonedge-sign and componentwise-residual routes.  A viable
reopening must give a canonical cross-component switch with a recoverable
inverse, or an independent global theorem controlling the sum of the shifted
remainder products.  Replacing that transport by a local sign assertion is
blocked by the exact order-11 witness above.

## Reproduction

```bash
python3 scratch_d22_global_switching_certificate_20260711.py
```

Expected final line:

```text
{'zero_capacity_n16': 'passed',
 'galvin_T_2_23_nonedge': 'passed',
 'componentwise_n11': 'passed',
 'bencs_residual_through_n8': 'passed'}
```
