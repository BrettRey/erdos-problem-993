# Erdős 993 solve-today registry

## Resolution gate

A route is complete only if it yields either:

1. a proof that the coefficient sequence of `I(T;x)` is unimodal for every
   finite tree `T`; or
2. an explicit finite tree together with an exact certificate containing a
   strict coefficient descent followed by a later strict ascent.

No conditional reduction, proper subclass, bounded-order computation, or
near-miss counts as a resolution unless a proved argument makes it exhaustive.

## Load-bearing assumptions to test

| Assumption | Falsification condition | Search response |
|---|---|---|
| A useful stronger invariant exists for every tree-DP state. | Small exact DP states violate every proposed invariant while retaining tree realizability. | Demand exact counterexamples to each invariant before investing in closure. |
| Minimal-counterexample induction can be made local. | The chosen deletion, contraction, subdivision, or leaf operation loses the coefficient comparison needed for unimodality. | Require an explicit coefficient transport inequality, not only smaller trees. |
| A product-plus-perturbation representation has usable reserve. | The perturbation can spend the entire post-descent margin or the missing budget is equivalent to the original conjecture. | Separate pre-descent and post-descent budgets and test exact extremizers. |
| A negative construction can be engineered from known non-log-concave trees. | Products or tree-realizable gluings preserve unimodality despite accumulating log-concavity failures. | Score strict descent-then-later-ascent directly and certify tree realizability. |
| Finite computation can contribute to closure. | No proved bound reduces a minimal counterexample to the computed range. | Treat fixed-size scans only as falsification tools. |

## Mandatory adversarial checks

- Verify that every proposed object is a finite simple connected acyclic graph.
- Compute coefficients exactly; floating-point ratios are diagnostic only.
- Distinguish log-concavity failure from non-unimodality.
- A terminal descent is not a violation; require a later strict ascent.
- Handle plateaus explicitly: first-strict-descent continuity cannot be assumed.
- Check every polynomial identity at endpoints and under zero/terminal support.
- Check that deletions and reductions which create forests preserve the claimed
  product or convolution structure component by component.
- Do not infer a coefficient statement from a mode or mean bound without a
  proved bridge.
- Do not treat PNP, ECMS, Conjecture A, the Case-B hub bound, mode--mean,
  STP2 closure, or an `A+xR` budget as closure unless the remaining implication
  to full unimodality is proved.
- Reject any global compatibility, realizability, or closure lemma described
  as routine without a proof or an exact counterexample search.
- Audit for circular reformulations: a missing lemma is not progress if it is
  equivalent in strength to tree independence-sequence unimodality.

## Approach-family registry

| ID | Mathematical family | Independence mechanism | State | Reopen rule |
|---|---|---|---|---|
| D1 | Tree-DP cone / realizability invariant | Shifted ratio domination and slope reserve | blocked: SRD false at rooted n=29 | New invariant that permits the exact deep-tail SRD reversal below |
| D2 | Minimal-counterexample structural induction | Peripheral-support leaf recurrence and mode corridor | blocked: subsumes open mode stability | Reopen only with a new tree-specific proof of vertex-deletion mode stability or slope reserve |
| D3 | Counterexample engineering | Rooted substitutions and multi-hub connectors | blocked: exact distance-1--4 connector search only retunes one peak | Reopen only with a gadget/score that creates a separated shoulder, not larger multiplicity |
| D4 | Product-plus-perturbation reserve | Quantitative descent margin versus `xR` | held for later cross-pollination | New perturbation budget beyond the proved PB reserve |
| D5 | Probabilistic / hard-core distribution | Root-conditioned phase mixtures | blocked: simple spherical mixtures | A gadget with large conditional mode separation, non-negligible root-on mass, and exact coefficient control |
| D6 | Prefix LC plus certified decreasing tail | Matching-block variance / ordered-LC inequality | active: contracted-tree CSP round | Prove the stronger `Var(e) <= mu+2 eta` in the prefix window |
| D7 | Pair-switching / transition components | Aggregate no-subset-sum-one fibers | blocked: signed whole-fiber rho compression also fails | Reopen only with a non-compression-based global Hall mechanism |
| D8 | Formal finite-certificate closure | Exact semialgebraic/Lean certificates plus universal tail | unassigned | A frozen theorem packet with all global hypotheses discharged |
| D9 | Maximum-independent-set / Hall expansion | Bernstein expansion and alternating closure | blocked at aggregate block variance | Reopen only with new TP positivity or an aggregate quotient-tree covariance theorem |
| D10 | Contracted-block induction | Leaf-block peeling with root-conditioned moments | blocked at rooted boundary covariance | Reopen only with a proved quantitative deficit-versus-addability covariance invariant |
| D11 | Adversarial variance search | Exact marked DP, evolutionary/rooted grammars | completed negative falsification round | Reopen with a new structured family or a proved size-reduction theorem |
| D12 | Extremal variance tree | Excess-degree moment inequalities | rank-three prefix GSB proved and independently audited | Generalize the mechanism uniformly in rank; fixed-r accumulation is insufficient |
| D13 | Prefix extension-slope bound | Variance with residual-edge correction | active: weaker prefix round | Prove `Var(e_r)<=2 mu_r+2 eta_r` only for `r<=L-2` |
| D14 | Ranked incidence operators | Down-up spectral/Poincare inequality | blocked: aggregate ED false on exact n=210 prefix witness | Reopen only with a function-specific Poincare bound that does not split gap and energy |
| D16 | Maximum-set activity split | M-footprint and addable-H expectations | blocked at two theorem-strength monotonicity lemmas | New bounded-reverse exchange or rooted-state telescope |
| D17 | Negative convolution/bridge construction | Nonunimodal forest product followed by rare-root connection | frozen after D19/D20/D23 construction round | Reopen only with a mechanism escaping bridge compression, root-mass loss, and binomial smoothing |
| D18 | Distance-class covariance | Far-pair cancellation plus distance-two fork budget | parked at aggregate `q_far`; the augmented fork is refuted | Reopen with a global rank-window mechanism or an exact prefix witness; pairwise, componentwise, and fixed-`lambda` fork variants are blocked |
| D19 | Decorated-path compression | Iterated rooted bridges and corner perturbations | blocked: the expansion compresses to one product-minus-corner bridge | A mechanism escaping the exact compression identity |
| D20 | Two-scale rooted-state grammar | Hard selected-state bumps plus direct or nested leaves | blocked for tested bounded grammars: direct-leaf mass loss and nested binomial smoothing | A state species preserving both correcting mass and a separated shoulder |
| D21 | Matching-block drift | Size-biased drift of occupied-block deficit | blocked/circular: the averaged drift is exactly prefix GSB; pointwise loss fails on `K_(1,4)` | A new decomposition of `Var(e)-2 eta` |
| D22 | Global symmetric-difference switching | Nonedge covariance and componentwise collision capacity | blocked: bare signs and componentwise allocations are false | A canonical cross-component transport with a recoverable inverse |
| D23 | Recursive phase grammar | Stable or marginal root-occupancy phase cycles | blocked for periodic high-branching grammars by the exact response/variance recurrence | An exact mechanism outside the recurrence obstruction |
| D24 | Bencs deficit decomposition | Termwise positive Christoffel--Darboux extraction | blocked: separated and canonically charged summands are false | Group complete switching orbits or compensate induced subtrees globally |
| D25 | Fork/nonedge covariance | Local distance-two allocations and aggregate nonedge signs | blocked: bare signs/local allocations and every fixed `AF_lambda`, `lambda in [0,1]`, have exact prefix obstructions | Reopen only with a non-fixed or genuinely global allocation not reducible to the failed family |
| D26 | Matching-block support flow | Support exchange, fallback marks, and full max flow | diagnostic: bare exchange false; fallback empirical; full-flow global cut equals GSB | A support or cut lemma strictly simpler than GSB |
| D27 | Aggregate far covariance | Orbit and common-set regrouping of `q_far` | parked at the aggregate sign; orbit/common-layer signs are false | A global proof mechanism for `q_far<=0` or an exact prefix witness |

## Output contract for workers

Return a concrete proved lemma, construction, equation, exact counterexample to
a proposed lemma, or a precise blocked statement with evidence that the gap is
theorem-strength. Vague strategy reports are rejected.

## Round log

### Root D5 diagnostic, round 1

- Root conditioning writes `I_T=A+B`, with `A=prod I(T_j)` and
  `B=x prod E(T_j)`.
- An exact/normalized-FFT sweep covered 13,625 spherically symmetric branching
  vectors up to roughly 20,000 vertices and found no post-descent ascent.
- The strongest apparent float crossing, branching vector `(4,28)` on 141
  vertices, was numerical tail noise. Exact expansion
  `((1+x)^4+x)^28 + x(1+x)^112` has its first descent at coefficient 55 and
  thereafter decreases; its largest checked post-descent ratio is
  `2266807816743165859328362951097460 /
   2412070186882005665977300213701912 < 1`.
- A moment search found that large root-conditioned mean separations occur only
  when the smaller branch is correspondingly rare (stars are the extremal
  elementary example). No non-negligible two-phase mixture exceeded about two
  pooled standard deviations in the tested spherical grammar.
- Reopen D5 only with a new recursive gadget that breaks this
  separation-versus-weight tradeoff; broader spherical scans are not a new
  mechanism.

### Root D5 variance--weight lemma

There is a general explanation for the observed tradeoff.  For a rooted child
`C`, let `p` and `q=1-p` be the probabilities that its root is included and
excluded under the uniform hard-core measure, let `V` be the unconditional
variance of the independent-set size, and put
`delta = mu(C)-mu(C | root excluded)`.  Removing the root injects root-on sets
into root-off sets, so `r=p/q<=1`.  Total variance gives

```text
delta^2 <= r V.
```

For a parent with children `C_i`, write `A=prod I(C_i)` for the parent-off
branch and `B=x prod E(C_i)` for the parent-on branch.  If
`R=B(1)/A(1)=prod q_i` and `L=-log R`, then

```text
|mu_A-mu_B| <= 1 + sqrt(L Var(A) / log 2).
```

Indeed `log(1+r_i)>=r_i log 2` on `0<=r_i<=1`, hence
`sum r_i<=L/log 2`, and Cauchy--Schwarz applied to
`|delta_i|<=sqrt(r_i V_i)` proves the bound.  Thus comparable conditional
branches cannot have arbitrarily large standardized mean separation.  The
lemma is rigorous but does not by itself imply coefficient unimodality; a new
shape-preservation theorem would still be required.

### Worker D1 diagnostic, round 1

For rooted states write `H=J/x`, so

```text
E=prod_i(E_i+xH_i),   H=prod_i E_i,   I=E+xH.
```

The candidate shifted ratio domination

```text
E_k H_k >= E_(k-1) H_(k+1)
```

survived all rooted trees through `n=16` but is exactly false.  Take
`S=T_(3,4,1)` (28 vertices), attach a new leaf `r` to its center, and root the
29-vertex tree at `r`.  The tail values satisfy

```text
E_13=5410, E_14=60, H_14=60, H_15=1,
E_14 H_14 = 3600 < 5410 = E_13 H_15.
```

The underlying tree remains unimodal; the reversal is the known harmless
deep-tail LC failure.  Coefficientwise domination, component modes,
log-concavity, and this ratio order therefore cannot be the inductive cone.
The exact remaining signed-slope expression is

```text
R_k=(E_k-E_(k-1))+(H_(k-1)-H_(k-2)).
```

Asserting that its signs never return from negative to positive is merely the
original conjecture unless a new magnitude reserve is supplied.

### Worker D2 result, round 1

Let `v` be an endpoint of a longest path and `u` its support, with `r` leaf
neighbors.  Every other neighbor of `u` except the next path vertex is a leaf.
For the remaining connected tree `R`, the exact leaf recurrence is

```text
I(T)=A+xB,
A=I(T-v),
B=(1+x)^(r-1) I(R).
```

This is a real simplification: a minimal-counterexample induction needs only
the smaller trees `T-v` and `R`, not arbitrary forest-product closure.
Multiplication by `1+x` preserves unimodality because for `q=(1+x)a`,
`q_k-q_(k-1)=a_k-a_(k-2)`.

The natural adjacent-mode lemma is exactly false.  For the 11-vertex tree

```text
{08,09,19,1-10,29,39,4-10,5-10,6-10,7-10}
```

and `v=4,u=10`,

```text
A=[1,10,36,62,61,37,13,2]          (mode 3),
B=[1,9,31,56,59,37,13,2]            (mode 4),
xB                                    (mode 5),
I(T)=[1,11,45,93,117,96,50,15,2]    (mode 4).
```

The summand-mode gap is two.  Exact peripheral scans through `n=15` and 5,000
seeded random trees of orders 20--200 found no larger gap, but this is not a
proof.  If `A` and `C=xB` have ordered modes `r_A<l_C`, their sum is automatic
outside the corridor.  Put

```text
D_k=A_k-A_(k+1),  U_k=C_(k+1)-C_k.
```

For a gap of two, closure reduces to the single implication

```text
D_r>U_r  ==>  D_(r+1)>=U_(r+1).
```

A focused second round is attacking both the universal gap-two bound and this
coefficient inequality.

### Worker D3 result and adversarial correction, round 1

An exact adjacent-two-hub search produced a connected 50-vertex tree with
central coefficients

```text
..., 26130652221, 29962940088, 29962939291, 26145814020, ...
```

and graph6

```text
q????A?oB?Jw???????????????????C??B???W???o??Op??????????????????????????????????????????????o????Do????G{??????????????????????????????????????????????????????????B???????M??????DG??????GWO??O????????_?C?_
```

Independent replay verifies 50 vertices, 49 edges, connectedness, and the full
exact polynomial.  The two central peak candidates differ by only 797.
However, this is **not** a near-counterexample to unimodality: changing the
sign of that difference merely moves the mode from 16 to 17.  The true maximum
ratio after the first strict descent is only `0.8726051128` at `17 -> 18`.

The active follow-up therefore scores only genuine post-descent pressure in
the three-vertex connector identity

```text
I=(1+x)P_LP_R+x(A_LP_R+P_LA_R)+x^2A_LA_R.
```

The extra `x^2 A_L A_R` term is structurally new, but it counts only if it
creates a separated later rise rather than retuning the central peak.

### Root D6 falsification pass

The candidate overlap theorem is:

```text
i_k(T)^2 >= i_(k-1)(T) i_(k+1)(T)
for every k < L(T),
L(T)=ceil((2 alpha(T)-1)/3).
```

Together with the known strict decrease for `k>=L(T)`, this would prove full
unimodality.  Exact exhaustive replay over all 1,346,024 trees through `n=20`
found zero prefix-LC failures.  This is only a falsification result; a separate
worker is testing known large non-LC families and seeking a tree-specific
proof.

## Iteration checkpoint 1

- **Learned:** generic rooted ratio cones fail exactly on harmless deep-tail LC
  bumps; peripheral support avoids arbitrary forest products; prefix LC only
  up to the certified tail-lock threshold survives the current exact frontier.
- **Changed assumptions:** adjacent summand modes are not universal, and a
  near-tie between peak coefficients is not evidence of a valley.
- **Assessment:** the portfolio is working as a falsification and localization
  process, but no end-to-end route is closed.  D2 and D6 are the strongest
  incompatible proof mechanisms; D3 remains the negative control.
- **Decision:** iterate.  Keep D2, D6, and the genuine connector-crossing search
  alive; do not reopen D1 or simple D5 mixtures without a new invariant or
  gadget.

### Root D7 transition-component formulation

Prefix log-concavity at `k` is equivalent to an injection from ordered pairs
`(C,D)` of independent sets with sizes `(k-1,k+1)` into ordered pairs of sizes
`(k,k)`.  In the induced symmetric-difference forest, each component has its
two color classes inherited from `C\D` and `D\C`; let its imbalance be
`delta_i=|D_i|-|C_i|`.  The imbalances sum to 2.

If a component subset has imbalance sum 1, swapping the two color classes on
that subset produces a balanced `(k,k)` pair.  This is the standard
component-switch mechanism and can be canonically inverted.  The only hard
pairs are those for which no subset of component imbalances sums to 1, notably
a connected component of imbalance 2 (the star/claw pattern).

For such a component, removing a vertex on the smaller side, adding an unused
private neighbor on the larger side, and swapping branches of total imbalance
2 repairs the cardinalities in the model star case.  A full proof needs:

1. a canonical pivot/private neighbor for every bad pair in the prefix window;
2. a branch subset of total imbalance 2; and
3. an inverse rule proving global injectivity.

This route is neither the D2 coefficient corridor nor a restatement of prefix
LC; it turns the missing curvature inequality into a concrete forest-flow and
private-neighbor problem.

### D2 closure audit and block decision

Independent exhaustive replay over all 847,834 diameter-endpoint
decompositions through `n=18` confirms summand-mode gap at most two, with
16,690 gap-two cases and no corridor reversal.  In the notation

```text
h=r-1,  P=I(R),  Q=I(R-w),  B=(1+x)^h P,
A=B+xQ,  C=xB,
```

the gap-two slope differences are exactly

```text
D_k-U_k       = q_(k-1)-q_k-(b_(k+1)-b_(k-1)),
D_(k+1)-U_(k+1)=q_k-q_(k+1)-(b_(k+2)-b_k).
```

The stronger observed reserve

```text
q_(k-1)-q_k <= b_k-b_(k-1)
```

has no proof.  Moreover, the gap bound itself is a special case of the
repository's open vertex-deletion mode-stability law.  Abstract unimodal
sequences satisfying the same algebra, even with one binomial leaf factor,
admit exact corridor reversals.  D2 is therefore theorem-strength blocked and
will not receive further agents without a materially new invariant.

### D6 exact moment reduction

For a uniform independent `(k-1)`-set `A`, let

```text
U=V(T)\N[A],  e=|U|,  h=|E(T[U])|,
mu=E[e],      sigma^2=Var(e),  eta=E[h].
```

Double counting one- and two-vertex extensions proves

```text
i_k^2-i_(k-1)i_(k+1)
= i_(k-1)^2/[k^2(k+1)]
  * (mu^2+k mu+2k eta-k sigma^2).
```

Hence prefix LC is exactly the variance bound

```text
sigma^2 <= mu^2/k + mu + 2 eta.
```

All 522,957 trees through `n=19` and 26,303 large structured/random stress
trees satisfy the prefix theorem; every observed LC defect lies at least four
indices beyond the tail-lock boundary.  The new D6 round attacks the variance
bound directly.

### D7 local-repair obstruction

The prefix threshold does not eliminate no-subset-sum-one fibers.  In
`K_(1,s)`, take `C={center}` and `D` any three leaves at `k=2`; the symmetric
difference is a connected imbalance-two claw and `C` is maximal.  Thus local
addable-vertex repair is false even far inside the prefix as `s` grows.  D7 now
targets a global aggregation/fractional matching across such fibers rather
than another local pivot.

### D3 connector round 2 closure

The exact connector search has converged without a counterexample.  It used
all 7,571 rooted polynomial states through order 12 at every coordinate,
asymmetric allocations of 3--14 branches, connector distances 1--4, and an
integer objective that scores only a rise after the first strict descent.
The best verified tree has order 171 and central coefficients

```text
a_65=89941452972805650859226893350068360202051
a_66=92705648651575226263678799238037275965211
a_67=92705100333763826791526450435231306198568
a_68=89942694129961538944021684880694141808579.
```

Thus `66 -> 67` is a strict descent, but `a_68/a_67` is only
`0.9702022197931196`; there is no later ascent.  For the original 797-deficit
witness the connector terms erase the putative descent exactly:

```text
Delta_16 = -797 + 1,755,934,551 + 234,737,807
         = 1,990,671,561 > 0.
```

Larger multiplicity merely makes adjacent central coefficients closer and is
not a mechanism for producing a second shoulder.  D3 is now blocked.  The
reproducibility driver is `scratch_connector_path_search_20260711.py`.

### D6 stronger matching-block reduction

The moment target admits a substantially stronger formulation.  For a uniform
independent `r`-set let `e_r` be its number of one-vertex extensions.  Exact
edge-size-biasing gives

```text
E[e_(r+1)] = E[e_r] - 1 + (Var(e_r)-2 E[h_r])/E[e_r].
```

Consequently

```text
Var(e_r) <= E[e_r] + 2 E[h_r]
```

is equivalent to monotonicity of the average extension count,

```text
(r+2) i_(r+2)/i_(r+1) <= (r+1) i_(r+1)/i_r,
```

or log-concavity of the ordered counts `r! i_r`.  This is stronger than the
original prefix-LC target.  It has zero prefix failures on every tree through
order 20, on 70,000 random/hub-biased trees through order 160, and on 2,463
large structured stress rows through order 600.  Stars asymptotically saturate
the variance bound at `r=1`.

There is now a concrete tree-specific model.  Fix a maximum matching.  Its
edges together with the unmatched singleton vertices partition `V(T)` into
exactly `alpha(T)` blocks, each of size one or two.  Every independent `r`-set
occupies exactly `r` blocks.  If `Y_j` is the number of currently addable
vertices in block `j`, then

```text
e=sum_j Y_j,             sum_j Y_j^2 <= e+2h.
```

Contracting the matching edges produces a tree on the blocks.  Thus it would
suffice to prove

```text
Var(sum_j Y_j) <= E[sum_j Y_j^2]
```

for the resulting two-state-per-block tree CSP in the prefix range.  Pairwise
negative correlation is not available (positive block covariances already
occur in `K_(1,4)`), so the missing argument must control the aggregate
covariance through the contracted-tree structure.  The stronger inequality
fails on the known 26-vertex deep-tail LC witness, as required: at `r=12`,
`mu=221/993`, `eta=29/993`, and
`Var=288448/986049 > mu+2 eta`.

### D7/D9 alternating-closure cross-pollination

For a maximum independent set `M`, its complement `C` is a minimum vertex
cover, and a maximum matching gives an injection `f:C -> M`.  Given an
independent set `S`, write `X=S intersect C` and `Y=M minus S`.  Starting from
any unmatched `y in Y minus f(X)`, the backward alternating closure gives sets
`W subset M`, `P subset X` with

```text
N(W) intersect S = P,       |W|=|P|+1.
```

Hence `(S minus P) union W` is an independent `(r+1)`-set, and every `S` has
`alpha-r` distinct canonical exchanges toward `M`.  A uniform reverse-degree
bound is false already at order five, so this is not itself a Hall proof; it is
being combined with D7's cross transitions and D6's empty matching blocks.

D7 also found exact capacity obstructions to three narrower injection models:
a fixed-large-set deletion model fails at order 11, an optimal within-fiber
flow fails by `256 < 257`, and a two-sided diamond flow fails at order 12 by
`459 < 467`.  Adding the cross transition `(C,D) -> (C+p,D-q)` together with
component switches repairs all 17,907 pairs in that obstruction and every
tested tree through order 11.  Recoverability and a universal Hall argument are
the live obligations.

### Root D6 exact low-rank base

The first two ordered-log-concavity inequalities hold for every tree, without
any prefix assumption.  If `T` has order `n`, then

```text
i_1=n,                  i_2=(n-1)(n-2)/2.
```

Inclusion--exclusion over the edges of a triple gives

```text
i_3=C(n,3)-(n-1)(n-2)+sum_v C(deg(v),2).
```

Convexity of the degree sum (or the standard edge-pair count) shows
`sum_v C(deg(v),2) <= C(n-1,2)`, with equality at the star, and hence
`i_3 <= C(n-1,3)`.  Therefore

```text
i_1^2 >= 2 i_0 i_2,
2 i_2^2 - 3 i_1 i_3 >= (n-1)(n-2) > 0.
```

Thus the stronger factorial/ordered LC theorem is rigorously closed at
coefficient indices one and two; any minimal failure has `k>=3` and
`alpha>=6` in the prefix window.

The next rank is also now closed.  For uniform independent pairs (`r=2`),
write `t=n-1`, `x_v=deg(v)-1`, and

```text
Q=sum_v x_v^2,   R=sum_v x_v^3,   C=sum_(uv in E) x_u x_v.
```

If `N=i_2`, `d1=sum_A e(A)`, and `d2=sum_A e(A)(e(A)-1)`, exact counting gives

```text
N  = t(t-1)/2,
d1 = (3Q+t^3-6t^2+8t-3)/2,
d2 = (-12C+8Qt-24Q-2R+t^4-12t^3+45t^2-56t+22)/2.
```

For a non-star tree whose nonleaf core has at least two vertices,

```text
C >= t-2,              R >= Q^2/(t-1).
```

The first inequality follows by summing `xy>=x+y-1` over core edges;
the second is Cauchy--Schwarz using `sum x_v=t-1`.  Substitution yields

```text
4(d1^2-Nd2) >=
  (2t+9)Q^2+(-2t^3-4t^2+24t-18)Q
  +t^5-5t^4+11t^3-14t^2-2t+9.
```

The quadratic has positive leading coefficient and discriminant

```text
-4t(t-1)(t^4-4t^3-7t^2+94t-216),
```

which is negative for `t>=4` (the quartic is 48 at four and strictly
increasing thereafter).  Stars have zero variance, and the remaining orders
are direct.  Hence `Var_2(e)<=E_2(e)` for every tree, with no prefix
restriction.  The formulas were independently replayed against the exact
moment DP for every tree through order 13.  Therefore a minimal failure of the
strong prefix theorem must now have `r>=3`, equivalently coefficient index
`k>=4`.

### D7 final block decision

Within a fixed symmetric-difference fiber with component imbalances `d_i`, put

```text
P_F(z)=product_i (z^(d_i)+z^(-d_i)).
```

The exact source and balanced-target multiplicities are `[z^2]P_F` and
`[z^0]P_F`; their positive difference is the precise export demand.  Fixed-`D`
deletion capacity fails at order 11 (`256<257`), two-sided diamond flow fails
at order 12 (`459<467`), alternating-closure reverse degree fails in the exact
prefix at order 8 (`6>4`), and closure-growth tensor deletion fails at order
10 (`604<620`).  The full deficit-to-surplus network using internal switches,
ordinary extensions, alternating closures, and cross moves passes exact tests
through order 12 but has no noncircular Hall/min-cut proof.

For the stronger same-rank target `Var_r(e)<=E_r e`, the exact weighted source
demand in a fiber is

```text
((r+2)(r+1)[z^2]P_F - (r+1)^2[z^0]P_F)_+.
```

Adjacent addable pairs have no internal balanced target and require exported
capacity.  The augmented network passes through order nine, but asserting all
of its min-cuts nonnegative is presently equivalent-strength repackaging.  D7
is therefore blocked pending a genuinely new block-tree cut inequality.

### Root exact extension-moment DP

`scratch_extension_variance_dp_20260711.py` implements the bivariate two-jet
of

```text
sum_(A independent) x^|A| y^e(A)
```

using three tree states: selected root, excluded root with no selected child,
and excluded root with at least one selected child.  At `y=1` it returns, for
every rank `r`, exact arrays

```text
N_r,  sum_A e(A),  sum_A e(A)(e(A)-1).
```

The DP agrees with brute-force subset enumeration on every tree through order
10.  The stronger Poisson-type prefix inequality

```text
N_r sum_A e(A)(e(A)-1) <= (sum_A e(A))^2,
```

equivalently `Var_r(e)<=E_r e`, has zero failures on all 1,346,024 trees through
order 20.  This bound strictly implies
matching-block ABV but is a falsification target until a proof mechanism is
found; conditional versions already fail on a 20-arm subdivided star inside
the prefix, so any proof must retain compensation across fibers.

### D13 direct slope reduction

Let

```text
mu_r=(r+1)i_(r+1)/i_r,
eta_r=E[ |E(T[V minus N[A]])| ]
```

for a uniform independent `r`-set `A`.  The exact two-extension count gives

```text
mu_(r+1)=mu_r-1+(Var(e_r)-2 eta_r)/mu_r.
```

Therefore the inequality

```text
Var(e_r) <= 2 mu_r + 2 eta_r                         (GSB)
```

implies `mu_(r+1)<=mu_r+1`.  Consequently, if it holds through `r=L-2`,

```text
d_r=mu_r-(r+1)
```

is nonincreasing through the only uncertified prefix, while
`sign(d_r)=sign(i_(r+1)-i_r)`.  Combined with the known decreasing tail from
`L`, prefix GSB proves unimodality, with plateaus handled automatically.  In
coefficients it is the concrete
three-term inequality

```text
(r+2)i_r i_(r+2)
  <= (r+1)i_(r+1)^2 + i_r i_(r+1).
```

This condition is strictly weaker than ordered log-concavity and survives the
known 26-vertex LC failures.  Exact exhaustive replay finds zero prefix-GSB
failures on all 1,346,024 trees through order 20.  The tempting strengthening
`Var(e_r)<=2mu_r` is false: Galvin's `T_(6,6,1)` at `n=79, alpha=42, r=37`
has `Var/mu=2.05174754756`, while the edge-corrected GSB ratio remains below
one (equivalently `mu_38-mu_37=0.4814349<1`).  The residual-edge term is
therefore structurally necessary, not cosmetic.

The unrestricted GSB is itself false, but only in a harmless deep-tail
regime.  Galvin's `T_(14,8,1)` has `n=239`, `alpha=126`, `L=84`, and at
`r=113`

```text
mu_113 approximately 0.34012745,
mu_114 approximately 1.42291291,
mu_114-mu_113 approximately 1.08278546 > 1.
```

Thus a global monotone-`d_r` proof is blocked.  D13 remains live solely as the
strictly weaker prefix theorem `r<=L-2`; the certified tail absorbs the exact
counterexample.

### D10 peripheral-support induction closure

Let `u` have `ell` leaf neighbors and unique nonleaf neighbor `w`, put
`F=T-(u union L)` and `R=F-w`, and define

```text
Z_H(x,t)=sum_(A independent in H) x^|A| t^e_H(A).
```

The exact recurrence is

```text
Z_T=(x+t)^ell Z_F+(t-1)t^ell G+xZ_R,
```

where `G` records the independent sets of `F` which avoid `w`, with their
`F`-extension count.  Its second `t`-jet contains the unavoidable boundary
polynomial `I(F-N_F[w])`; unrooted moments therefore do not close.

Plain conditional tensorization is false inside the prefix.  In the tree
with edges `a-w`, `w-c`, and `c-l_i` for `1<=i<=m`, conditioning rank-one
sets on `w` being unselected but dominated leaves `{a}` and `{c}`, with

```text
(e,h)=(m+1,m), (1,0),
Var(e)=m^2/4,       2 E(e+h)=2m+2.
```

This fails from `m=9`, at `n=12, alpha=10, r=1`.  Bundling the peripheral
star phases does pay the entire within-bundle variance exactly, but the
between-bundle term requires a rooted deficit-versus-addability covariance
bound.  Charging it only to the bundle reserve is itself false for `K_(1,8)`
rooted at a leaf, by `7/256`.  D10 is blocked until a genuinely quantitative
rooted covariance invariant is supplied.

### D11 completed falsification frontier

The exact bivariate extension-moment DP has now checked every 3,490,529
nonisomorphic tree through order 21, including all 2,144,505 trees at order
21, with zero prefix failures of

```text
Var_r(e)<=E_r(e).
```

The independent adversarial grammar/Pruefer/structured search likewise found
no prefix violation through its stated large-family ranges.  This remains
falsification evidence only; D11 is closed until a new family or an exhaustive
reduction changes its logical force.

### D12 rank-three prefix theorem

An exact symbolic certificate proves

```text
alpha(T)>=7  ==>  5 i_3 i_5 <= 4 i_4^2+i_3 i_4,
```

equivalently `mu_4<=mu_3+1`.  Since `alpha>=7` is exactly the condition that
`r=3` lies in the GSB prefix, any prefix-GSB failure must have `r>=4`.

The proof uses closed formulas for `i_3,i_4,i_5` in degree and radius-two
moments.  Cores with at most three nonleaf vertices have exact Bernstein
certificates.  For at least four nonleaf vertices and `t=n-1>=14`, put
`x_v=deg(v)-1`, `Q=sum x^2`, `R=sum x^3`, `U=sum x^4`, and
`C=sum_(uv) x_u x_v`.  The new load-bearing bounds are

```text
L+P21 <= R+(t+7)Q+8C+10t-6,
U <= (t-4)R,   C >= t-2,   R >= Q^2/(t-1).
```

They reduce the deficit to a globally strongly convex quartic with an exact
positive Newton lower bound.  The remaining orders `10<=n<=14` are an
exhaustive exact base, not the closure mechanism.  The theorem packet is
`notes/rank3_prefix_gsb_theorem_packet_2026-07-11.md`, and the replay script
is `scratch_rank3_gsb_certificate_20260711.py`.  Root replay and a fully
independent derivation/audit both passed.  The audit independently reproduced
the coefficient formulas, core polynomials, Bernstein conversions, joint
radius-two bound, moment substitutions, convexity certificate, and finite
base counts/minima; it also found no failure through all trees of order 18 or
5,000 random trees through order 200.  The finite-base script now asserts its
expected row counts and minima to detect truncated future enumeration.

### D14 down--up spectral decomposition

On the independent `r`-sets define the symmetric down--up chain

```text
P_r(A,B)=(1/r) sum_(C subset A intersect B, |C|=r-1) 1/e(C).
```

It is stochastic with uniform stationary law.  If `H_C=T[U(C)]`, choosing
`C` with probability `e(C)/(r i_r)` and then two uniform vertices of `H_C`
gives one stationary transition.  Since

```text
e(C+v)=e(C)-1-deg_(H_C)(v),
```

the Dirichlet form of `e` is exactly the size-biased average of the degree
variance in `H_C`.  Prefix GSB would follow from the two concrete estimates

```text
gap(P_r) >= 1/(2r),
Dir_r(e) <= E_r(e+h)/r.
```

The usual matroid constants are exactly false: `gap>=1/r` fails on an
order-eight double star at prefix rank two, and `gap>=1/(r+1)` fails on an
order-nine double star at prefix rank three.  The `1/(2r)` scale survives all
exact tests through order 14 and fails outside the prefix, so its window is
essential rather than cosmetic.

For a forest `H` with `q` vertices and `h` edges, uniform `V`, and
`H'=H-N[V]`, the local lemma

```text
Var(deg_H(V)) <= 2 E[ |V(H')|+|E(H')| ]
```

is proved from `sum deg^2<=h^2+h`.  It is insufficient by a factor depending
on `r`; the exact aggregate energy obligation is

```text
sum_(C in I_(r-1))
 [q(q+h-1)-2h-(r+1)S2+4r h^2/q] >= 0,
S2=sum_(v in H_C) deg(v)^2.
```

This statement is false pointwise on residual-star fibers, so any proof must
transport their deficit globally.  The exact marked DP
`scratch_spectral_energy_dp_20260711.py` tracks the complete `q` distribution
and the sums of `h,h^2,S2`; it agrees with brute-force enumeration on every
configuration of the checked small symmetric trees.  A balanced double
binary tree of depth five (`n=126, alpha=84`) reaches the exact aggregate
ratio

```text
1086952108361328767313914454180497906654 /
1201372626824757077843019095072071632405
  approximately 0.9047585105
```

at the last prefix rank `r=54`.  A targeted spherical search then found an
exact counterexample to the aggregate energy inequality.  Join the roots of
two identical rooted spherical trees with root-outward branching sequence

```text
(2,3,2,1,2,1,1).
```

The resulting tree has `n=210`, `209` edges, `alpha=124`, and last prefix rank
`r=81`.  The exact marked DP gives

```text
sum_C [q(q+h-1)-2h-(r+1)S2+4r h^2/q]
= -42691269729915413999161105330329696919890850296588501772067055418675071
  / 50018905999077748956903379700 < 0.
```

Equivalently its exact energy ratio is

```text
245700013545993183538657579915009525438626090145082199405497542315007719
/
240956539131558137538750790433861781336415995667683476986378980601821600
  approximately 1.0196860165386306.
```

The recurrence was independently audited and brute-replayed on 120 random
labelled trees at every order through ten.  The witness itself is connected,
acyclic, and its independence number was independently recomputed.  Its true
GSB ratio is only about `0.21377694`, so this blocks only the separated D14
energy estimate; it is not a counterexample to D13 or to unimodality.
The frozen proof-of-obstruction packet is
`notes/d14_ed_obstruction_packet_2026-07-11.md`; the standalone exact replay
`scratch_d14_ed_obstruction_certificate_20260711.py` completed with
`certificate: passed`.

### D7 maximum-set grading

Fix a maximum independent set `M` and grade `S` by
`rho(S)=|S minus M|`.  With

```text
P_j(z)=sum_(S in I_j) z^rho(S),
Delta_k(z)=k P_k(z)^2+P_(k-1)(z)P_k(z)
           -(k+1)P_(k-1)(z)P_(k+1)(z),
```

every closure/deletion move in the current exact network is nonincreasing in
total `rho`.  Writing `H=V-M`, the exact Hall expansion is

```text
P_j(z)=sum_(X independent in H)
       z^|X| C(alpha-|N_M(X)|,j-|X|).
```

All cumulative coefficient cuts of `Delta_k` are nonnegative for every
maximum-set choice and every prefix rank through order 15.  Individual
coefficients need not be nonnegative (`K_(1,4)` gives `48-6z`).  Moreover,
ordinary threshold cuts do not dominate arbitrary network cuts even on that
star.  The remaining live statement is therefore the weaker monotone-cut
claim that any *negative* Hall cut can be closed downward in `rho` without
destroying negativity; it needs a proof or an exact counterexample before the
cumulative polynomial can carry the network.  That claim also fails under
signed whole-fiber cancellation already on `K_(1,4)`: lowering a hard fiber
changes its exact cut slack from `9` to `24`.  D7 is therefore frozen.

### D16 maximum-set activity split

The post-compression algebra yields a useful exact identity but no proof.
For `rho(S)=|S-M|`, put `P_j(z)=sum_(S in I_j)z^rho` and

```text
b(S)=|(S intersect M) union N_M(S-M)|,
B_j=sum_(S in I_j)b(S)z^rho.
```

Then the GSB deficit splits exactly as `Delta_k=R_k+S_k`, where

```text
R_k=P_(k-1)B_k-P_kB_(k-1)
   =sum_(m in M)(P_k Q_(m,k-1)-P_(k-1)Q_(m,k)),
Q_(m,j)=P_j(T-N[m]),

S_k=P_(k-1)P_k+z(P_k P'_k-P_(k-1)P'_(k+1)).
```

Thus `R_k>=0` is monotonicity of the expected `M`-footprint, equivalently
nonincrease of the mean number of addable `M` vertices.  The inequality
`S_k>=0` is the corresponding one-step bound for mean addable `H=V-M`
vertices.  Private-neighbor exchanges have unbounded reverse degree, while
marked-H deletion loses the empty-block label, so both are theorem-strength
gaps.  Galvin's global deep-tail GSB witness refutes `R_k` but not `S_k`,
showing again that any valid statement must use the prefix window.  D16 is
blocked pending a new bounded-reverse or rooted-telescoping mechanism.

### D12 rank-four radius obstruction

The rank-three excess-moment proof cannot be generalized by reusing only its
radius-two data.  The exact graph6 trees

```text
L??????_B_H__y
L??????o@_DA@{
```

both have `n=13`, `alpha=9`, identical

```text
(S2,S3,S4,P,P21,L)=(66,240,1002,75,488,392),
(i3,i4,i5)=(175,274,269),
```

but have `i6=170` and `168`, respectively.  Their prefix rank-four GSB gaps
are `156031` and `159319`.  The replay script
`scratch_rank4_radius2_obstruction_20260711.py` verifies treehood, graph6
decoding, the moments, polynomials, and gaps.  The rank-three
exact-reconstruction certificate therefore cannot determine
`i6` from this recorded radius-two tuple; the witness does not rule out other
rank-four mechanisms.

### D18 distance-class covariance profile

For `X_v=1[v is addable]` under the uniform rank-`r` law, write `C2` for
the aggregate covariance over vertex pairs at distance two and `Cfar` for
the aggregate covariance over pairs at distance at least three.  Since

```text
sum_v Var(X_v)<=mu,
2 sum_(uv in E) Cov(X_u,X_v)<=2 eta,
```

prefix GSB follows from

```text
Cfar<=0,                 2 C2<=mu.
```

Pairwise far negativity is exactly false: an order-15 prefix example has a
positive distance-three covariance numerator `266`.  The aggregate statement
nevertheless has large exact slack.  Through every tree of order 18, for
prefix ranks `r>=4`, the maxima of `2C2/mu` at ranks four through nine are

```text
0.17281120, 0.11146225, 0.05968516, 0.01490464,
-0.12469829, -0.26196469.
```

The maximum of `Cfar/mu` is zero only because stars have no far pairs.
Among trees with a far pair, the rank-four maximum is `-1/52`, and the
maxima become more negative with rank on the exact frontier.  Thus the
audited low-rank bases remove the star-sharp `r=1` regime and leave a very
slack uniform target for `r>=4`.  D18 is active at an aggregate tree-edge
telescoping or compression theorem; no pairwise-local proof is permitted.

### D19--D24 continuation

D19 shows exactly that an iterated decorated path still compresses to

```text
F_next=F_prefix(E+S)-S_prefix S,
```

with the selected-corner term controlled by a two-rank-shift deletion
injection.  It therefore does not create independent perturbations.  D20's
two-scale state grammar exposes the complementary obstruction: direct leaves
make the selected correction exponentially negligible, while nested stars
convolve it with a growing binomial factor and erase the hard bumps.  The
co-scaling search found no exact crossing.

D21 proves that the proposed matching-block drift is algebraically equivalent
to prefix GSB, and `K_(1,4)` blocks the intended pointwise loss-three map.  D22
then gives exact prefix failures of the bare nonedge signs and of assigning
collision capacity independently to symmetric-difference components.  D23
derives exact mass, mean-response, and variance recurrences for recursive
phase grammars: stable or marginal periodic high-branching phases have
standardized mean separation tending to zero, while an expanding response is
repelling.  D24 extracts the exact Bencs deficit but finds prefix witnesses
against both separated energy and canonically charged termwise positivity.
These lanes are frozen at their recorded obstruction packets or certificates.

### D25--D27 continuation and surviving split

D25 refutes the aggregate nonedge sign on `T_(40,20,1)` at `n=1641,r=20`.
It also gives an exact 9,418-vertex last-prefix witness against both the
parent-oriented and symmetric inverse-degree fork allocations, while ordered
LC and GSB remain positive.  The final augmented candidate was

```text
Q_c <= B_c+N(lambda a_c+(1-lambda)S_c).          (AF_lambda)
```

Exact small-tree and initial structured hard scans retained `lambda=1/2`, but
a broader same-branch scan found a decisive exact obstruction.  Join a new
center to five copies of `G(60,18)`.  The resulting tree has

```text
n=11106, alpha=5701, L=3801, r=L-2=3799.
```

After clearing the common branch-root degree `61`, the local AF gap is already
negative at `lambda=0` and has negative slope.  Equality would require
`lambda approximately -0.0283850882378686`, so every `lambda in [0,1]` fails.
The exact replay is
`scratch_d25_augmented_fork_obstruction_certificate_20260711.py`.  The tree's
true GSB gap is positive and its independence sequence is unimodal.  This
blocks only the fixed-parameter local family.  In any event, summing
`(AF_lambda)` would have supplied only the distance-two budget `2q2`; it never
controlled `qfar`.

D26 refutes bare support exchange at an exact order-20 prefix witness.  A
one-fallback correction survives the recorded finite audits, but its available
summations do not prove GSB.  Full fallback/marked networks saturate through
order 11, yet their all-target capacity cut is exactly GSB, so they are a
diagnostic representation rather than a reduction.

D27 keeps the independent aggregate proposition `qfar<=0` compatible with the
prefix, but rules out two localizations: a complete switching orbit and an
exact-common-set layer can each contribute positively while the global sum is
negative.  The bounded AF cycle is now complete and refuted; aggregate
`qfar<=0` remains a separate theorem-strength dependency, parked pending a new
global rank-window compensation mechanism.  No broad D28 round follows.
