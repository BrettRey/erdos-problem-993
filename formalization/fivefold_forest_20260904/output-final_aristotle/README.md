This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# Fivefold bounds for blocked sets in finite forests

A complete, `sorry`-free Lean 4 / Mathlib formalization of the requested theorem
(`aristotle_input.md`, goal **G0**), together with the supporting lemmas exposed as
individually traceable declarations (**G1**) and this mapping plus the axiom audit (**G2**).

## The result

```lean
theorem fivefold_of_forest
    {V : Type*} [Fintype V] (G : SimpleGraph V)
    (hforest : G.IsAcyclic)
    (halpha : 5 ≤ G.indepNum)
    (hb2 : blockedCount G 2 = 0) :
    FivefoldConclusion G
```

`RequestProject/Main.lean` contains the four definitions of the request
(`Extendable`, `extendableCount`, `blockedCount`, `FivefoldConclusion`) verbatim,
including the `d ≤ G.indepNum` guard, and proves the theorem in exactly the required shape.
No hypothesis was added, weakened or strengthened.

## Building

```
lake build
```

`RequestProject` (glob `RequestProject.+`) is the default target, so every file below,
including `Main`, is built. The toolchain and Mathlib revision are pinned by
`lean-toolchain` and `lake-manifest.json`.

## Architecture

Everything is developed for a **fixed ambient graph** `G : SimpleGraph V`; an induced
subgraph is encoded by its vertex set `S : Finset V`. Vertex deletion, disjoint unions and
induction on `S.card` are then painless, and every "subgraph" of a forest is automatically a
forest. The counts over `S = Finset.univ` are identified with the counts of the request in
`RequestProject/Main.lean`.

Core counting definitions (`RequestProject/Defs.lean`):

| notion | declaration |
| --- | --- |
| independent subsets of `S` | `FivefoldForest.indFam` |
| independence number `α(S)` | `FivefoldForest.alpha` |
| containment in a maximum independent set of `S` | `FivefoldForest.Ext` |
| layer counts by size | `FivefoldForest.icnt`, `ecnt`, `bcnt` |
| defect-graded layers `e_d`, `b_d` (guarded by `d ≤ α`) | `FivefoldForest.eD`, `FivefoldForest.bD` |
| margin `e_d − 5 b_d` over `ℤ` | `FivefoldForest.margin` |

The induction invariant is
`Good G S ∨ StarProf G S` (`RequestProject/Basic.lean`), where `Good` says both margins at
defects 3 and 4 are nonnegative and `StarProf` is the exceptional profile of the four-leaf
star `K₁,₄` (`α = 4`, `e = (1,4,6,4,1)`, `b₃ = 1`, `b₄ = 0`, margins `−1` and `+1`).
Carrying the exceptional profile explicitly is what makes the induction close: `K₁,₄`
violates the defect-3 bound on its own, and the product closure lemma shows this can never
propagate.

## Written proof step → Lean declaration

### 1. Elementary layer inequalities (`RequestProject/Basic.lean`)

| step | declaration |
| --- | --- |
| double counting of the incidences between layers `α−k` and `α−k−1` | `ecnt_incidence` |
| `e₀ ≤ 2 e₁` | `e0_le_two_e1` |
| `e₁ ≤ 4 e₂`, `e₁ ≤ 3 e₂` (for `α ≥ 3`) | `e1_le_four_e2`, `e1_le_three_e2` |
| `e₁ ≤ e₀ + 6 e₂` | `e1_le_e0_add_six_e2` |
| top layer: `e_α = 1`, `b_α = 0`, `b₀ = 0`, `e₀ ≥ 1` | `eD_alpha_eq_one`, `bD_alpha_eq_zero`, `bD_zero`, `eD_zero_pos` |
| deleting a vertex changes `α` by at most one | `alpha_erase_le`, `alpha_erase_succ_ge` |

These use only the hypothesis `|S| ≤ 2 α(S)`, which acyclicity supplies.

### 2. Component factorization (`RequestProject/Union.lean`, `RequestProject/Closure.lean`)

| step | declaration |
| --- | --- |
| a disjoint union with no edges across | `Split` |
| `α(S₁ ⊔ S₂) = α(S₁) + α(S₂)` | `alpha_split` |
| general convolution for a property that factors | `card_split_conv` |
| independent / extendable layers convolve | `icnt_split`, `ecnt_split` |
| extendability factors through the parts | `ext_split_iff`, `ext_union_of_split` |
| reindexing size-graded to defect-graded convolutions | `iD_split`, `eD_split` |
| blocked layers of a part inject into those of the union | `bD_le_of_split` |
| margins of a product | `margin3_split`, `margin4_split` |
| **product closure of the invariant** (including the `K₁,₄` case) | `split_good` |

### 3. Star absorption and the small exceptional stars (`RequestProject/Cases.lean`)

| step | declaration |
| --- | --- |
| independent sets of a star | `star_indFam` |
| `α` of a star | `star_alpha`, `star_alpha_small` |
| the leaf set is the unique maximum independent set | `star_ext_iff` |
| `e_d = C(m,d)`, `b_d = [d+1 = m]` | `star_eD`, `star_bD` (via `star_ecnt`, `star_bcnt`) |
| stars are `Good`, except `K₁,₄` which has the exceptional profile | `star_good_or_star` |

`K₁,₂` and `K₁,₃` are excluded by eligibility (`b₁ = b₂ = 0`), `K₁,₄` is the exceptional
profile, and `m ≥ 5` is `Good`.

### 4. Pendant configurations (`RequestProject/Pendant.lean`, `RequestProject/PendantCount.lean`)

A pendant configuration `Pendant G S p w L` is a support vertex `p` whose neighbours in `S`
are a nonempty set `L` of leaves together with exactly one further vertex `w`.

| step | declaration |
| --- | --- |
| `S \ {p}` splits as base `⊔` leaves | `Pendant.split` |
| the leaf set is edgeless: `α = |L|`, `e_d = C(|L|,d)`, `b_d = 0` | `Pendant.leaves_alpha`, `leaves_eD`, `leaves_bD` |
| pushing an independent set off `p` | `Pendant.lift_of_mem`, `Pendant.lift` |
| `α(S \ {p}) = α(S)` and extendability is unchanged off `p` | `Pendant.alpha_rest`, `Pendant.ext_rest_iff` |
| independent sets through `p` are `{p} ∪ X` with `X` in the core | `Pendant.insert_p_mem_iff` |
| splitting a layer count according to the use of `p` | `Pendant.card_split_p`, `card_split_p_not` |
| **many leaves (`|L| ≥ 2`)**: `e_d` unchanged, `b_d` gains `i_{α−d−1}(core)` | `Pendant.eD_two_le`, `Pendant.bD_two_le` |
| **pendant `P₂` (`|L| = 1`)**: eligibility forces `α(core) = α(base)` | `Pendant.alpha_core_one` |
| pendant-`P₂` recurrences `e_d = e_d(rest) + e_d(core)`, likewise for `b_d` | `Pendant.eD_one`, `Pendant.bD_one` |
| propagation of the invariant, per leaf count | `pendant_one_good` (`|L| = 1`), `pendant_three_good` (`|L| = 3`), `pendant_four_le_good` (`|L| ≥ 4`), with `pendant_many_setup`, `pendant_many_margin3`, `pendant_many_margin4` |

`|L| = 2` cannot occur under eligibility; this is part of `pendant_many_setup` and is used
in `good_or_starProf`.

### 5. Essential-root cone identity (`RequestProject/Cone.lean`)

If every maximum independent set of `S` contains `w`, deleting `w` is a bijection on maximum
sets and the layers of `S` are built from those of the core.

| step | declaration |
| --- | --- |
| no blocked sets at defect `d` ⇒ every set of that size extends | `ext_of_bD_zero` |
| the maximum-set bijection | `cone_max_erase`, `cone_max_insert` |
| extendability through the cone | `cone_ext_iff`, `cone_ecnt` |
| `e₀ = f₀` and `e_{d+1} = f_{d+1} + f_d` | `cone_eD_zero`, `cone_eD_succ` |
| `b₂(S) = 0 ⇒ b₁(core) = 0` | `cone_bD_one` |

### 6. What acyclicity is used for (`RequestProject/Paths.lean`, `RequestProject/Forest.lean`)

The induction consumes acyclicity only through `ForestLike`:

* `card_le : ∀ T, |T| ≤ 2 α(T)`;
* `decomp`: every nonempty `T` is a star, or splits into two nonadjacent nonempty parts, or
  carries a pendant configuration whose base has at least two vertices.

| step | declaration |
| --- | --- |
| paths inside a vertex set, and a longest one | `IsTPath`, `TPathLen`, `maxTPathLen`, `tPathLen_max`, `not_tPathLen_of_gt` |
| the first vertex of a longest path has only the second as a neighbour (a shortcut would close a cycle) | `max_nbr_eq` |
| a neighbour of the second vertex other than the third is a leaf | `max_leaf` |
| existence of a vertex of degree at most one | `exists_leaf` |
| `|T| ≤ 2 α(T)` by deleting a leaf and its neighbour | `card_le_two_alpha` |
| a centre all of whose neighbours are leaves gives a star or a split | `star_or_split` |
| **acyclic ⇒ `ForestLike`** (case split on the length of a longest path: `0`, `1`, `2`, `≥ 3`) | `forestLike_of_isAcyclic` |

### 7. `b₂ = 0` forces `b₁ = 0`, and the induction (`RequestProject/Induction.lean`)

| step | declaration |
| --- | --- |
| for `α ≥ 4`: `b₂ = 0 ⇒ b₁ = 0` | `bD_one_eq_zero_of_bD_two` |
| **strong induction on the vertex count** | `good_or_starProf` |
| the fivefold bounds for an eligible set with `α ≥ 5` | `fivefold` |

`bD_one_eq_zero_of_bD_two` is the counting argument of the note: a blocked set `S₀` of size
`α−1` is maximal; for every `v ∈ S₀` the set `S₀ \ {v}` extends, producing two vertices
outside `S₀` whose only neighbour in `S₀` is `v`; these `2(α−1)` vertices are distinct, so
`|S| ≥ 3(α−1)`, contradicting `|S| ≤ 2α` for `α ≥ 4`.

`good_or_starProf` is a strong induction on `S.card` and uses only the three cases supplied
by `ForestLike.decomp`; the component lemma `split_good` is an independently proved algebraic
statement about the two (strictly smaller) parts, not an induction hypothesis about a
same-size forest.

### 8. Bridge to the statement of the request (`RequestProject/Main.lean`)

| step | declaration |
| --- | --- |
| `alpha G univ = G.indepNum` | `alpha_univ` |
| `IsMaximumIndepSet` in terms of `indFam` and `alpha` | `isMaximumIndepSet_iff` |
| `Extendable` is `Ext` over `univ` | `extendable_iff` |
| `extendableCount = eD`, `blockedCount = bD` | `extendableCount_eq`, `blockedCount_eq` |
| **the requested theorem** | `fivefold_of_forest` |

## Axiom audit

Output of `#print axioms` (Lean 4.28.0, Mathlib as pinned):

```
'FivefoldForest.fivefold_of_forest' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.fivefold' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.good_or_starProf' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.bD_one_eq_zero_of_bD_two' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.forestLike_of_isAcyclic' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.eD_split' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.bD_le_of_split' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.split_good' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.star_good_or_star' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.Pendant.eD_two_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.Pendant.bD_two_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.Pendant.eD_one' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.Pendant.bD_one' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.Pendant.alpha_core_one' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.cone_eD_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.cone_eD_succ' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.cone_bD_one' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.card_le_two_alpha' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.pendant_one_good' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.pendant_three_good' depends on axioms: [propext, Classical.choice, Quot.sound]
'FivefoldForest.pendant_four_le_good' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Only the three standard Lean axioms occur. The project contains no `sorry`, no `admit`, no
`axiom`, no `@[implemented_by]`, no `native_decide` and no enumeration standing in for the
theorem; the only `decide` calls evaluate small binomial coefficients such as `C(4,2) = 6`.

## Relation to the written note

The formalization follows the route of the note, with two deliberate differences.

* The note phrases the recurrences through generating functions `I`, `E`, `B`, `Q`, `H`.
  Here they are phrased directly as identities between the defect-graded counts
  `eD`, `bD`, `iD`; the pendant-`P₂` recurrence (7) of the note is
  `Pendant.eD_one` / `Pendant.bD_one`, and its root-conditioned enumerators `Q`, `H` are
  replaced by the counts of the core `pCore S p w L`.
* The induction carries the disjunction `Good ∨ StarProf` rather than `Good` alone. This is
  necessary: `K₁,₄` is eligible and fails the defect-3 bound, so any invariant that ignores
  it is false; `split_good` shows that the exceptional profile cannot survive a nontrivial
  product, and the pendant and star cases dispose of it otherwise.

`scripts/audit.py` is exploratory only. It enumerated all forests on at most 14 vertices and
found 978 eligible ones with no failure of either inequality. It is **not** part of the
verification, and nothing in the Lean development depends on it.
