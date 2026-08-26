This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# Matroid-representation lower bound for graph independence systems

All Lean source is in `RequestProject/Main.lean` (namespace `MatroidRep`), built against
Mathlib (Lean 4.28.0). `lake build` is clean: no `sorry`, `admit`, `axiom`,
`@[implemented_by]`, or `native_decide`. The file ends with a `#print axioms` audit of every
named declaration; each reports only `propext`, `Classical.choice`, `Quot.sound`
(`indep_singleton_of_represents` uses none).

## Result

If a family `M : Fin k → Matroid α` represents the independence system of a simple graph
`G : SimpleGraph α`, i.e.

```
Represents G M  :=  ∀ S : Set α, G.IsIndepSet S ↔ ∀ i, (M i).Indep S,
```

then `deg v ≤ k` for every vertex `v` whose neighbourhood is an independent set. No
hypothesis on the ground sets `(M i).E` is needed: the representation hypothesis applied to
singletons already forces every vertex to lie in every ground set.

Status: **G0 proved** in the stated generality (arbitrary type `α`, arbitrary graph, arbitrary
matroids), **G1 proved** (triangle-free form and the tree corollary), **G2 proved** (the
reduction from enlarged ground sets — including arbitrary minors — back to G0).

## Traceability table

| Declaration | Kind | Step of the intended proof |
|---|---|---|
| `MatroidRep.Represents` | definition | the representation hypothesis of the request |
| `MatroidRep.indep_singleton_of_represents` | lemma | Step 1: singletons are independent in each `M i` |
| `MatroidRep.isNonloop_of_represents` | lemma | Step 1: no `M i` has a loop |
| `MatroidRep.isCircuit_pair_of_represents` | lemma | Step 3, first half: a dependent pair is a circuit |
| `MatroidRep.not_dep_pair_of_dep_pair` | lemma | Step 3, main half: circuit elimination gives the contradiction |
| `MatroidRep.le_card_of_represents` | **theorem (G0)** | Steps 2 and 4: choice of `f`, injectivity, counting |
| `MatroidRep.degree_le_of_represents` | theorem | G0 restated as `G.degree v ≤ k` for `Fintype α` |
| `MatroidRep.isIndepSet_neighborSet_of_cliqueFree` | lemma | G1: triangle-free ⇒ neighbourhoods independent |
| `MatroidRep.le_card_of_represents_cliqueFree` | **theorem (G1)** | degree bound for triangle-free graphs |
| `MatroidRep.cliqueFree_three_of_isAcyclic` | lemma | acyclic ⇒ triangle-free (a triangle is a cycle) |
| `MatroidRep.le_card_of_represents_isTree` | **theorem (G1, tree form)** | degree bound for trees |
| `MatroidRep.represents_comap` | lemma | G2 reduction: comparison along an embedding `α ↪ β` |
| `MatroidRep.le_card_of_represents_of_embedding` | **theorem (G2)** | bound survives ground-set enlargement |
| `MatroidRep.le_card_of_represents_minor` | **theorem (G2, minor form)** | bound survives representation by minors `(P i / C i) ↾ R i` |

## Notes on the G2 formalization

Ground-set enlargement is modelled by an embedding `e : α ↪ β` of the vertex set into a
larger type together with the hypothesis
`∀ S : Set α, G.IsIndepSet S ↔ ∀ i, (N i).Indep (e '' S)`.
Comparison along `e` (`Matroid.comap`, which is available since `e` is injective, so
`Set.InjOn e` holds automatically) turns each `N i` into a matroid on `α` with independent
sets exactly `{I ⊆ α | (N i).Indep (e '' I)}`, giving a representation with the same `k`;
G0 then applies verbatim. This covers restriction, and — since a minor of a matroid is again
a matroid on the same type — also contraction: `le_card_of_represents_minor` states the bound
for representations by arbitrary minors `((P i).contract (C i)).restrict (R i)`. So the
expected reduction is **confirmed**, not refuted.
