# Conjecture A Origin Note (2026-03-05)

## Bottom line

Conjecture A appears to have emerged on **2026-02-15** as the **forced base case** of the new PNP reduction, not as an independently proposed extremal-family conjecture.

The motivating shift was:

1. the older global PNP / augmented-injection program had identifiable dead ends in the high-mode regime;
2. Hub Exclusion showed that vertices with many leaf-children are structurally easy for 1-Private sets;
3. Transfer then peeled those easy parts away;
4. the only residual hard class was trees with `d_leaf <= 1`.

That residual class was then named **Conjecture A**.

## Evidence chain

### 1. Pre-pivot state: high-mode analysis, but no Conjecture A yet

On **2026-02-14**, the active note was still the global PNP problem:

- [notes/pnp_n_lt_3k_proof.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/pnp_n_lt_3k_proof.md#L32) formulates the generalized 1-Private Mode Conjecture.
- [notes/pnp_n_lt_3k_proof.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/pnp_n_lt_3k_proof.md#L57) already interprets low-priv sets as "spread out" and below-mode sets as star-like / efficient-domination objects.

The commit message the same day is diagnostic, not reductive:

> Add PNP proof exploration: degree bound lemma, dead ends, and gap analysis

and explicitly lists dead ends:

> Color class bound
> Leaf Support hypothesis
> Hub inclusion

This is commit `c146218` in local history; see `git show --no-patch c146218`.

## 2. The pivot: Hub Exclusion + Transfer create a residual class

On **2026-02-15**, the reduction language appears all at once in the commit message:

> PNP framework: Hub Exclusion + Transfer Lemma reduce to Conjecture A (d_leaf <= 1)

This is commit `ba2155e` in local history.

The new status note makes the structure explicit:

- Hub Exclusion: [notes/one_private_status.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/one_private_status.md#L42)
- Case A is exactly the `d_leaf <= 1` frontier: [notes/one_private_status.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/one_private_status.md#L54)
- Case B reduces away by Transfer: [notes/one_private_status.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/one_private_status.md#L71)
- The note then states the conclusion directly: [notes/one_private_status.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/one_private_status.md#L89)

This is the strongest repository evidence for where Conjecture A came from.

## 3. Why this class, specifically?

The same note records a retracted predecessor:

- the old idea was a **Leaf-Mode Inequality**: [notes/one_private_status.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/one_private_status.md#L93)
- it fails on leaf-heavy spiders / brooms: [notes/one_private_status.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/one_private_status.md#L97)
- but those counterexamples are still handled by Hub Exclusion: [notes/one_private_status.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/one_private_status.md#L104)

Reasonable inference:

The repository suggests a reversal of perspective. Leaf-heavy high-mode trees were no longer the dangerous class once Hub Exclusion was proved; they became the reducible class. The genuinely unresolved class was therefore the complementary, leaf-light one: `d_leaf <= 1`.

## 4. The manuscript adopts that exact framing the next morning

On **2026-02-16**, the paper rewrite replaces the old augmented-injection pitch with the new reduction pitch:

- old abstract/contribution language emphasized augmented injection and Hall verification;
- new language says unimodality reduces to a single conjecture about trees where every vertex has at most one leaf-child.

The current paper states the reduction cleanly:

- reduction to the residual forest: [paper/main_v2.tex](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/paper/main_v2.tex#L212)
- explicit reduction sentence: [paper/main_v2.tex](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/paper/main_v2.tex#L216)
- Conjecture A statement: [paper/main_v2.tex](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/paper/main_v2.tex#L219)

The corresponding commit message is:

> Add new paper (main_v2.tex) with subdivision-contraction identity and PNP reduction

This is commit `686c42e` in local history.

## 5. Later support was evidentiary, not original

After the pivot, the repo develops two support programs for Conjecture A:

- mean-bound / hard-core route: [notes/hard_core_attack.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/hard_core_attack.md#L1)
- spider-extremal / mean analysis: [notes/conjecture_A_analysis.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/conjecture_A_analysis.md#L17)

Then on **2026-02-18**, Steiner peeling proves the key mean inequality for the class:

- statement `mu(T) < n/3`: [notes/steiner_peeling_proof_2026-02-18.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/steiner_peeling_proof_2026-02-18.md#L5)
- explicit claim that this yields Conjecture A for `d_leaf <= 1` under the mode-mean step: [notes/steiner_peeling_proof_2026-02-18.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/steiner_peeling_proof_2026-02-18.md#L105)

These notes explain why the conjecture looked plausible, but they appear after the structural reduction had already isolated the class.

## One-sentence origin story

Conjecture A was inspired by the realization that **multi-leaf hubs are not the obstruction but the simplifiable part** of the PNP problem; once Hub Exclusion and Transfer were found, the proof effort collapsed onto the leaf-light residual class `d_leaf <= 1`, and that residual statement became Conjecture A.
