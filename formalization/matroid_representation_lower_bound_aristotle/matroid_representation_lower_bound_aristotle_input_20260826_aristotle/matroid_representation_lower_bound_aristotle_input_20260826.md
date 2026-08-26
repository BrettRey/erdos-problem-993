# Formalization request: the matroid-representation lower bound for graph independence systems

## Objective

Formalize the lower-bound half of the representation lemma in
`notes/matroid_intersection_number_lane_2026-08-26.md`: representing a graph's
independence system as an intersection of matroids on the vertex set requires at
least `deg v` matroids, at every vertex `v` whose neighbourhood is independent.

This is a self-contained finite matroid statement. It does **not** depend on any
other job in this queue, does not involve independence polynomials,
log-concavity, or Erdős Problem 993, and needs no chromatic-index theory.

The upper bound (that `Δ(T)` partition matroids suffice) is deliberately **out of
scope**: it needs partition-matroid infrastructure Mathlib lacks, and it is
already verified computationally on all 200 trees with `n ≤ 10`. The lower bound
is the half with no computational check behind it, which is exactly why it is
here.

## Definitions

Let `α` be a type, `G : SimpleGraph α` a simple graph, and let

```text
IS(G) = { S : Set α | S.Pairwise (fun u v => ¬ G.Adj u v) }
```

be its independent sets (`SimpleGraph.IsIndepSet` in current Mathlib,
`Mathlib/Combinatorics/SimpleGraph/Clique.lean`).

For a matroid `M : Matroid α`, `M.Indep` is its independence predicate.

Say a finite family of matroids `M : Fin k → Matroid α`, each with ground set
`M i |>.E = Set.univ` (or a fixed vertex set `V`; state whichever is cleaner and
keep it uniform), **represents** `G` when

```text
∀ S, G.IsIndepSet S ↔ ∀ i, (M i).Indep S.
```

## Final theorem

```text
theorem le_card_of_represents
    (G : SimpleGraph α) (M : Fin k → Matroid α)
    (hrep : ∀ S, G.IsIndepSet S ↔ ∀ i, (M i).Indep S)
    (v : α) (hv : G.IsIndepSet (G.neighborSet v))      -- neighbours pairwise non-adjacent
    (hfin : (G.neighborSet v).Finite) :
    hfin.toFinset.card ≤ k
```

i.e. **`deg v ≤ k`** whenever the neighbourhood of `v` is an independent set.

State the degree side with whatever finiteness setup is most natural
(`G.LocallyFinite`, `Fintype α`, `G.degree v` for a locally finite graph). A
`Fintype α` version is entirely acceptable if it makes the statement cleaner;
do not contort the statement to gain generality.

## Intended proof

Fix `v`, let `u₁, …, u_d` enumerate `G.neighborSet v`.

1. **No loops.** For every `x`, the singleton `{x}` is in `IS(G)` (a one-element
   set is trivially pairwise non-adjacent), so by `hrep` every `M i` has
   `(M i).Indep {x}`; hence no `M i` has a loop, and every element is a nonloop.
2. **Every edge at `v` is killed by some matroid.** For each `j`, `v` and `u_j`
   are adjacent, so `{v, u_j} ∉ IS(G)`, so by `hrep` there is an index `f j` with
   `¬ (M (f j)).Indep {v, u_j}`.
3. **`f` is injective.** Suppose `f j = f j' = i` with `u_j ≠ u_j'`. Both
   `{v, u_j}` and `{v, u_j'}` are dependent in `M i`. Since `M i` is loopless
   (step 1), each is a *circuit*: it is dependent and all its proper subsets
   (the empty set and singletons) are independent.
   Apply circuit elimination to these two distinct circuits, both containing `v`:
   there is a circuit `C ⊆ ({v,u_j} ∪ {v,u_j'}) \ {v} = {u_j, u_j'}`. A circuit is
   nonempty, and no singleton is a circuit in a loopless matroid, so
   `C = {u_j, u_j'}`, which is therefore **dependent** in `M i`.
   But `u_j` and `u_j'` are both neighbours of `v` and `hv` says the neighbourhood
   is independent, so `¬ G.Adj u_j u_j'`, so `{u_j, u_j'} ∈ IS(G)`, so by `hrep`
   it is independent in `M i`. Contradiction.
4. **Conclude.** `f` is an injection from the neighbourhood of `v` into `Fin k`,
   so `deg v ≤ k`.

### Proof-route freedom

Step 3 may instead be done through **parallelism**, if that is shorter in Lean.
Mathlib's `Matroid/Loop.lean` has `IsLoop` / `IsNonloop` but (as of the checkout
consulted, 2026-08-26) no parallel-class API. Defining

```text
Parallel M e f  :=  M.IsNonloop e ∧ M.IsNonloop f ∧ M.closure {e} = M.closure {f}
```

makes transitivity immediate, and the needed bridge is: for nonloops `e ≠ f`,
`¬ M.Indep {e,f} ↔ Parallel M e f`. Either route is acceptable. Use whichever
Mathlib supports better; do not build more API than the theorem needs.

Relevant existing Mathlib: `Matroid.IsCircuit.strong_elimination` and
`Matroid.IsCircuit.strong_multi_elimination` (`Mathlib/Combinatorics/Matroid/Circuit.lean`),
`Matroid.isCircuit_iff_forall_ssubset`, `Matroid.IsNonloop`,
`SimpleGraph.IsIndepSet` and `SimpleGraph.isIndepSet_iff`.

## Graded sub-targets

- **G0 (required).** The final theorem above, for a general `SimpleGraph` and a
  vertex whose neighbourhood is independent.
- **G1.** The tree corollary: for `G.IsTree` (or any triangle-free `G`, i.e.
  `G.CliqueFree 3`), every neighbourhood is independent, hence `k ≥ Δ(G)` for any
  representing family. Triangle-freeness is the only property used, so state it
  that way and derive the tree case if `Mathlib`'s `IsTree` gives
  `CliqueFree 3` cheaply. **If the tree-to-triangle-free step turns out to be
  expensive, deliver G0 plus the triangle-free form and say so.**
- **G2 (optional, adversarial and genuinely wanted).** Does the bound survive
  **ground-set enlargement**? Suppose the `M i` live on a larger ground set
  `E ⊇ α` and `IS(G)` is recovered as `{S ⊆ α | ∀ i, (M i).Indep S}`.

  Partial answer already in hand, worth formalizing as a lemma because it makes
  the rest routine: **restriction reduces this to the on-`α` case.** The
  restriction `(M i) ↾ α` is a matroid on `α` whose independent sets are exactly
  `{I ⊆ α | (M i).Indep I}`, so the enlarged representation yields an equal-`k`
  representation by matroids on `α`, and G0 applies unchanged. (Checked
  concretely on partition matroids with auxiliary elements, 2026-08-26; Mathlib
  has restriction as `Matroid.restrict` / `↾` in `Matroid/Minor/Restrict.lean`.)

  The **open** part, and the one actually worth effort: does the same hold when
  **contraction** is allowed, i.e. if `IS(G)` is recovered from minors
  `(M i) / C i` rather than restrictions? My expectation is that it does, for a
  structural reason worth stating as the target: a minor of a matroid is a
  matroid, so `((M i) / C i) ↾ α` is again a matroid on `α` and G0 applies. If
  that reduction is right, the invariant is robust to every minor-based
  representation and only an exotic non-minor encoding could beat it.

  **Confirm that reduction formally, or refute it with an explicit
  counterexample.** A refutation is a more valuable result than G0, since it
  would show the invariant I have named is the wrong one.

## Deliverable and grading

- Standalone current-Mathlib Lean project, source-labelled declarations, README
  traceability table mapping each declaration to the step above.
- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, hidden support
  hypotheses, weakened theorem, or bounded/decidable-instance computation over a
  fixed finite type standing in for the general statement.
- `lake build` clean, plus `#print axioms` for the final theorem and every named
  lemma it depends on.
- Grade **COMPLETE** only if G0 is proved in the stated generality. If G1 is
  reached, say so explicitly. If the statement is **false** as written, return the
  smallest explicit counterexample (the graph, the vertex, and the matroids, with
  their independent sets listed) and the exact step that fails.

An honest PARTIAL with G0 proved and G1 left open beats a complete-looking
submission that quietly specialises the theorem to a finite type or assumes the
matroids are partition matroids. Do not assume the `M i` are partition matroids:
that assumption destroys the point of the theorem.
