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

Empirical anchor: for
`C2 = {trees with at most 2 branch vertices (degree >= 3)}`,
the scan in `results/two_branch_lc_n24.json` reports no log-concavity
or unimodality failures up to `n = 24`.

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

Reality check from `notes/lamport_transition.md`: this full tail inequality
is not globally true (explicit broom-leaf failures exist), but in the scanned
box all failures occur at the boundary index \(k=d(IU)\). So a boundary-only
version is a plausible reduced kill shot.

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

The squeeze scan shows the first descent occurs at most **3 indices** before
that tail (n ≤ 23), but this is not globally valid.

Explicit counterexample: the star \(S_{28}=K_{1,28}\) has
\(\alpha=28\), \(t=\lceil(2\alpha-1)/3\rceil=19\), and \(d(I)=15\), so
\[
  d(I)=t-4<t-3.
\]
Hence a universal bound \(d(I)\ge t-3\) is false.

So Step 4 should be used as an empirical guide, not a global theorem.
The proof target becomes:

1) **MBI setup:** assume a minimal tree with \(d\le t-4\) and take the first
   bad index \(b\le t-4\).
2) **Leaf recurrence at \(b\):** use
   \(I(T)=I(T-v)+xI(T-N[v])\) to force concentrated negative drift.
3) **Leaf-heavy reduction:** show any MBI witness must be leaf-heavy
   (star/double-star-like near the obstruction).
4) **Close reduced class:** prove the boundary inequalities directly in that
   leaf-heavy class.

Empirical boundary check (n ≤ 23): when \(d(I)=t-3\) with \(t=\lceil(2\alpha-1)/3\rceil\),
the steps \(i_{t-2} \ge i_{t-1} \ge i_t\) always hold. This suggests the entire
problem may reduce to two local inequalities near \(k=t-2\) after structural
reduction.

Leaf-attachment lemmas (see `notes/leaf_attachment_mbi.md`) strengthen this:
for sufficiently large s (leaves attached at one hub), the boundary indices
k = t-2, t-1 satisfy \(\Delta I_k \le 0\) automatically. Therefore any minimal
counterexample must be leaf-light at each hub, reducing the problem to a finite
core class plus bounded leaf-load checks.

---

## Cross-note status (proved vs missing)

1) `notes/subdivision_lemma.md`:
   - proved exact identity \(I(T') = I(T) + Q_uQ_v + xP_uP_v\),
   - proved coefficientwise bound \(A \le (1+x)I(T)\),
   - still missing: first-difference dominance or mode-ordering strong enough
     to conclude \(I(T)+A\) is always unimodal.

2) `notes/lamport_transition.md`:
   - proved restricted mode-ordering lemma (broom root + leaf child),
   - found global counterexamples to full difference-dominance,
   - observed boundary-only failures in the scanned broom-leaf table.

3) `notes/steedman_invariants.md`:
   - global “\(B_T\) log-concave for all rooted trees” is blocked by the
     construction in `notes/subdivision_lemma.md`,
   - the viable route is class-restricted closure plus one boundary inequality
     for the final sum.

Common gap across all three routes: a local inequality at the first bad index
plus a structural reduction to the leaf-heavy obstruction class.

---

## “If I had to bet” (Erdős heuristics)

1) Prove subdivision preserves unimodality (or a weaker local variant).
2) Use Lamport transition + difference dominance on a small extremal core.
3) Push MBI + leaf-heavy reduction, then seal the last boundary indices with a
   direct binomial bound in that reduced class.

One clean inequality is worth more than a broad framework. That’s the Erdős way.
