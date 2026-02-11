# One‑Page Erdős Plan (Minimal Counterexample + Sharp Local Inequality)

**Question:** Is every tree’s independence sequence unimodal?

**Erdős move:** Assume a minimal counterexample exists and force it into a
rigid form, then kill it with one sharp inequality.

---

## Step 0: Minimal counterexample setup

Assume a counterexample \(T\) exists. Choose one with minimal \(|V(T)|\), and
let \(k\) be the first valley index:

  \(i_{k-1}(T) > i_k(T) < i_{k+1}(T)\).

All proper subtrees are unimodal by minimality.

---

## Step 1: Structural pruning (force simplicity)

Prove a reduction lemma that forbids extraneous structure in a minimal
counterexample. Two viable targets:

1) **Subdivision invariance:** if subdividing any edge preserves unimodality,
   then a minimal counterexample has no degree‑2 vertices (homeomorphically
   irreducible).

2) **Leaf‑removal monotonicity:** use the recurrence
   \(I(T)=I(T-v)+xI(T-N[v])\) to show removing a leaf cannot create a new valley
   once a local inequality holds.

Goal: collapse the search to trees with **few branch vertices** and a
small “core.”

---

## Step 2: Reduce to a “most adversarial” local pattern

Use the Lamport‑style transition for a rooted parent (P,Q) and child (U,V):

  \(P' = P(U+V),\quad Q' = QU,\quad I' = IU + PV\).

Define the safety property:

  once \(\Delta I_k < 0\), all later differences are \(\le 0\).

**Single‑inequality target (the kill shot):**

For \(k \ge d(IU)\),

  \(\Delta(PV)_k \le -\Delta(IU)_k\).

If this holds at every composition step, unimodality follows by induction.

---

## Step 3: Prove the inequality in extremal configurations

Erdős wouldn’t try to prove the inequality for all trees. He’d prove it
for the **forced extremal local structure** from Step 1.

Current promising regimes:

1) **Broom root + leaf child:** proved mode‑ordering.
2) **Broom root + path child:** proved eventual mode‑ordering (large \(s\)).
3) **Small‑core attachments:** use the fixed‑core asymptotics to dominate
   the correction term.

Goal: show any minimal counterexample must contain a local configuration
that **already satisfies** the inequality, contradicting the existence of a
valley.

---

## Step 4: Use the squeeze (quantitative tail)

Levit–Mandrescu guarantees a strictly decreasing tail:
  \(k \ge \lceil(2\alpha-1)/3\rceil\).

The squeeze scan shows the first descent occurs at most **2 indices** before
that tail (n ≤ 20). Erdős would try to prove a uniform bound:

  \(d(I) \ge \lceil(2\alpha-1)/3\rceil - 2\).

Then the inequality in Step 2 only needs to control a **finite** number of
indices near the boundary.

---

## “If I had to bet” (Erdős heuristics)

1) Prove subdivision preserves unimodality (or a weaker local variant).
2) Use Lamport transition + difference dominance on a small extremal core.
3) Seal the last 2–3 boundary indices with a direct binomial bound.

One clean inequality is worth more than a broad framework. That’s the Erdős way.
