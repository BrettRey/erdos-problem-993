# Formalization request: normalized clan cancellation for adjacent two-hub trees

## Objective

Formalize in Lean 4 the normalized clan-graph bookkeeping underlying the
adjacent two-hub log-concavity proof described below.  The highest-value
outcome is an independent proof or refutation of the normalization and block
partition.  Do not assume the proposed Laurent blocks are correct merely
because they are stated here.

A concrete counterexample to any proposed identity, the `p = 2` extension,
or the claimed exhaustive partition is a successful result.  An honest
partial formalization proving one substantive normalization lemma is better
than a completed-looking project that hides the difficult step in an axiom.

## Source and notation

The mathematical source is Li, Li, Yang, and Zhang, *Log-concavity of
independence polynomials of some families of graphs* (arXiv:2501.04245),
especially their Section 2 and proof of Theorem 3.1.  The definitions and
results needed from that source are restated here so the task is
self-contained.

Let `G` be a finite simple graph and let

```text
alpha : V(G) -> Nat
```

where zero is allowed.  A proper multicoloring of type `alpha` assigns to
each vertex `v` a set of exactly `alpha v` positive-integer colors, with
disjoint color sets on adjacent vertices.  Its monomial records, for each
color, how many original vertices use that color.  The sum of these
monomials is the multicoloring symmetric function `Xmulti G alpha`.

The clan graph `clan G alpha` replaces each original vertex `v` by a clique
of size `alpha v`; for every original edge `uv`, it joins every clone of `u`
to every clone of `v`.  If `Xchrom` is the ordinary chromatic symmetric
function, the source normalization identity is

```text
Xchrom (clan G alpha)
  = (product over v of (alpha v)!) * Xmulti G alpha.
```

Thus `Xmulti G alpha` is the normalized chromatic symmetric function of the
clan graph.  Put

```text
Y G = sum over alpha : V(G) -> Nat of Xmulti G alpha.
```

The sum is degreewise finite.

For a finite graph `H`, a stable partition is a set partition of `V(H)`
into independent blocks.  A two-block stable partition with block sizes
`k,l` is *semi-ordered*: if `k = l`, the two equal-sized blocks are ordered.
Write `stableCount H k l` for the resulting count, with precise boundary
conventions made explicit in Lean.  Li--Li--Yang--Zhang Proposition 2.4
states, for `k >= l >= 1`,

```text
[s_(k,l)] Xchrom H
  = stableCount H k l - stableCount H (k+1) (l-1).
```

For this packet, it is acceptable to *define* the two-row coefficient by
the right-hand side and avoid developing a general library of symmetric
functions.  The later bridge from two-row positivity of `Y G` to
log-concavity of the independence polynomial is the cited Corollary 2.2
and is outside the mandatory kernel.

If a connected graph is bipartite, its two color classes are unique up to
swap.  If it is nonbipartite, it has no stable partition into two blocks.
For a bipartite component with color-class size difference `delta`, its two
orientations can be encoded by `z^delta + z^(-delta)`.  For a clan map
`alpha`, define the candidate normalized imbalance Laurent polynomial

```text
W(G, alpha; z)
  = (1 / product_v (alpha v)!)
      * product over bipartite connected components C of
          (z^(delta C) + z^(-delta C)).
```

This candidate requires careful boundary treatment for isolated vertices,
balanced components, empty blocks, and disconnected graphs.  Do not take it
as an axiom.  Derive the correct statement from `stableCount` and the clan
normalization, or produce a counterexample and repair it.

## Mandatory normalization kernel

Formalize finite graphs, clan graphs, stable two-block partitions, factorial
normalization, Laurent coefficient sequences (finite maps from `Int` to
`Rat` are acceptable), and prove as much as possible of the following.

1. **Clan normalization.**  Proper colorings of `clan G alpha` are in
   `product_v (alpha v)!`-to-one correspondence with proper
   multicolorings of type `alpha`, with the same monomial.

2. **Weighted bipartition formula.**  For a bipartite clan graph, derive
   its normalized two-row coefficient from the orientation/imbalance
   Laurent polynomial.  State and prove the exact correction needed for
   balanced components and boundary partitions.  The intended interior
   formula is

   ```text
   normalizedTwoRowCoeff alpha k l
     = coeff W (k-l) - coeff W (k-l+2).
   ```

3. **Products.**  Prove the precise multiplicative statement for disjoint
   bipartite components used when outside path components are factored from
   a local cancellation block.

Do not introduce the normalization identity, weighted bipartition formula,
or multiplicativity as axioms.  If general chromatic symmetric functions
would make the project impractical, prove the corresponding finite
coloring/stable-partition counting statements directly.

## Published spider transformation

Let `S(lambda)` be a spider with hub `w` and ordered pendant paths.  Consider
a clan map `alpha` with `alpha w = 1`.  If a neighbor of the hub has
multiplicity at least two, the clan graph contains a triangle and its
two-row contribution is zero.  Otherwise the hub component consists of the
hub and a positive prefix of each active arm, all with multiplicity one.

Let `p` be the number of odd-length positive prefixes.  The published bad
range has `p >= 3`.  Choose the first shortest odd prefix `A`, of odd length
`L`, using the fixed arm order to break ties.

- If `L = 1`, set the hub multiplicity to zero and replace the unique
  vertex of `A` by multiplicity two.
- If `L >= 3`, choose the first other odd prefix `B`.  Set the hub
  multiplicity to zero.  Along all `L` vertices of `A`, replace the run of
  ones by `2,0,2,0,...,2`.  Along the first `L-1` vertices of `B`, replace
  the run of ones by `2,0,2,0,...,2,0`.  Leave all other multiplicities
  unchanged.

The unequal alternating runs encode the inverse.  This preserves total
clan-graph order, but it does **not** preserve the raw factorial denominator
term-by-term; the cloned `K2` orientation counts are supposed to cancel the
new factors of `2!`.  That cancellation is one of the main facts to prove,
not an assumption.

The proposed adjacent two-hub proof also applies this same transformation
when `p = 2`, outside the bad range for which the source needed it.  In that
case the second odd prefix exists.  Prove or refute:

```text
localMapP2_injective
localMapP2_preserves_total_order
localMapP2_has_claimed_normalized_two_row_weight
```

The last theorem must be stated using the normalization kernel above, not
only unweighted bipartition counts.

## Adjacent two-hub target

Let `D(a,b)` be a tree with adjacent hubs `u,v` and arbitrary finite ordered
lists of positive pendant-path lengths attached at the two hubs.  For one
active local state at `u`, let `p >= 2` be the number of odd positive
prefixes and put `r = p-1`.  After all common outside components and their
normalization weights are accounted for, the proposed two-state normalized
imbalance block is

```text
A(r,c;z) = c*(z + z^(-1))^r + z^r + z^(-r),   c > 0.
```

The scalar `c` must be derived from the normalized outside components.  It
must not simply be postulated to be a positive integer.

For active blocks at both hubs, with parameters `(r,c)` and `(s,d)`, the
proposed four-map identity is

```text
F(z)
  = A(r,c;z) * A(s,d;z)
      - z^(r+s) - z^(-(r+s)).
```

The subtraction reflects that when both hubs are present, their edge forces
opposite colors and changes the connected imbalance from the two independent
extreme orientations to `|r-s|`.

Formalize enough clan-map data to prove or refute all of the following:

1. the local pairing maps are injective at each hub, including `p = 2`;
2. the global clan maps partition disjointly and exhaustively into
   four-state blocks, one-active two-state blocks, and individually
   nonnegative or two-row-zero maps;
3. no partner clan map is used in two blocks;
4. each block has exactly the claimed **normalized** Laurent contribution;
5. the four-map identity above follows;
6. its Laurent coefficients are centrally unimodal, so every normalized
   two-row coefficient of the block is nonnegative.

If the full arbitrary-arm partition is too large for one run, prioritize
the normalization kernel and the `p = 2` local block.  Do not replace the
arbitrary-arm theorem by a bounded computation without labeling the result
partial.

## Deliverable

Return a standalone Lean 4 project using a current Mathlib toolchain.  It
must build with `lake build` and contain:

- no `sorry`, `admit`, `axiom`, `implemented_by`, or weakened theorem;
- a `ClanAudit.lean` file with source-labelled definitions and theorems;
- a README traceability table mapping each mathematical definition and
  theorem above to its Lean declaration;
- a result note grading the outcome as one of:
  `NORMALIZATION_KERNEL`, `P2_LOCAL`, `ADJACENT_COMPLETE`, `REFUTED`, or
  `INSUFFICIENT`;
- for `REFUTED`, a smallest explicit finite witness and an explanation of
  which proposed identity fails;
- for `INSUFFICIENT`, the strongest completed theorem and the first exact
  unfilled mathematical obligation, with no placeholder hidden in the Lean
  project.

The grading order is intentional: a rigorously proved normalization kernel
or a concrete refutation is a valuable success even if the complete
adjacent theorem is not reached.
