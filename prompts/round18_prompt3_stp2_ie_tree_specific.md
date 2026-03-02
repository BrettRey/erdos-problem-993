# Round 18, Prompt 3: Tree-Specific Closure for `STP2(I,E)`

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

`STP2(I,E)` is now heavily verified computationally for rooted trees (all rootings, `n<=18`), but synthetic sequence counterexamples show it does **not** follow from abstract assumptions like:

- `STP2(E,J)`,
- `LC(E), LC(J)`,
- `J <= E`.

Therefore we need a tree-specific derivation.

Also, candidate adjacent-minor condition
`w_m(I,E) = I(m)E(m-1)-I(m-1)E(m) >= 0`
is **not** universally true on trees (counterexample already found at `n=7`).

---

## Your job

Provide a non-circular induction-ready route for `STP2(I_v,E_v)` that can be used at **all rootings** (not support-only).

Use coefficient convention `I = E + xJ`, and ladder minors

`Lambda_k(F,G) := G(k)F(k) - G(k+1)F(k-1)`.

---

## Task 1: Give the clean induction interface

Formulate exactly what must be true at each rooted subtree `(T,v)` so parent updates can use it without support assumptions.

Specify:

1. What is IH at children.
2. What is proved at the parent.
3. Which auxiliary lemmas are external (LC, positivity, etc.) and which are inside the chain-STP2 induction.

Be explicit about avoiding circular use of prefix dominance / unimodality claims.

---

## Task 2: Derive a direct expansion for `Lambda_k(I,E)`

Starting from `I = E + xJ`, derive an exact identity for

`Lambda_k(I,E) = c_k(E) + [E(k)J(k-1) - E(k+1)J(k-2)]`.

Then move beyond this identity: derive either

1. a Cauchy-Binet/diagonal representation for `Lambda_k(I,E)` under child-product updates, or
2. a representation that exposes the exact obstruction term and how tree structure controls it.

Avoid abstract cone arguments known to fail on synthetic data.

---

## Task 3: Isolate the missing tree-specific ingredient

Given synthetic counterexamples, identify one candidate property `H` that is plausible for tree-realizable triples and strong enough to close `STP2(I,E)`.

Candidates should avoid adjacent-minor positivity `w_m(I,E) >= 0` (already false).
Focus instead on stronger multi-index or diagonal-aggregation constraints that are specific to tree DP structure.

For each candidate:

- state it precisely,
- show why it blocks the known synthetic counterexample pattern,
- explain how it might be proved from the tree DP.

---

## Task 4: Deliverable standard

Provide one of:

1. Full proof sketch to theorem level, or
2. A reduced lemma checklist where each lemma is explicitly testable by computation.

No claims of implication from generic LC+STP2 unless fully proved.
