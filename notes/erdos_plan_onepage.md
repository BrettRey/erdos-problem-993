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

Additional pruning idea (degree-2 chain peeling):
If a leaf v has neighbor u of degree 2 and w is the other neighbor of u, then
with T'' = T - {u, v}, we have
  I(T) = I(T - v) + x I(T - N[v]) = (1+x) I(T'') + x I(T'' - w).
Thus a degree-2 step is a leaf-attachment with s = 1. Iterating along a maximal
degree-2 path gives a leaf-attachment form with s equal to the path length plus
lower-order corrections. By the leaf-attachment boundary lemma, sufficiently long
degree-2 chains are automatically safe. Therefore any minimal counterexample must
have all degree-2 chains bounded in length by an explicit s0(core) threshold.
See `notes/degree2_chain_peeling.md` for a formalized version of this peeling
argument and the resulting leaf-attachment decomposition.

New inductive target (degree-2 neighbor reduction):
If T has a leaf v whose neighbor u has degree 2 (other neighbor w), then
  I(T) = (1+x) I(T'') + x I(T'' - w),   T'' = T - {u, v}.
Thus a full proof could reduce to a perturbation lemma for h = (1+x)f + xg with
f = I(T'') and g = I(T'' - w). See `notes/perturbation_lemma.md` for this
reduction; the simple ratio-monotonicity condition is **too strong** (fails
empirically), so the target is a weaker boundary/tail condition.
Empirical note: for `n <= 20` all tested degree-2 leaf steps satisfy
mode(g) <= d(f) and have no tail rises for k >= d(f)+1 (see
`results/perturb_ratio_n20.json`).
Separately, full-vertex scans through `n <= 23` show
|mode(I(T-w)) - mode(I(T))| <= 1 for all vertices w and, more strongly,
mode(I(T-w)) <= d(I(T)) for all vertices w
(see `results/mode_alignment_n18.json`,
`results/mode_alignment_n20.json`, and
`results/mode_alignment_n21_mod8_merged.json`,
`results/mode_alignment_n22_mod8_merged.json`,
`results/mode_alignment_n23_mod8_merged.json`).
Current run status/tooling is tracked in
`notes/mode_alignment_status.md`.
Note: mode(I(T-w)) can be far to the right of mode(I(T-N[w])) (gap up to 9
in earlier scans), so mode-alignment between g and h is not a viable proof route.
For the alpha-structure around a deleted vertex, see
`notes/alpha_vertex_characterization.md`.
Also, an exhaustive `n<=20` test shows the candidate compact characterization
for tight equality `mode(I(T-w))=d(I(T))` is false; see
`results/tight_mode_equivalence_n20.json` and
`notes/mode_alignment_status.md`.
New split suggested by `g_w := alpha(T-w)-alpha_w`:
empirically, tight equalities are overwhelmingly in `g_w>=1`
(`results/mode_alignment_gw_n18_summary.json`: `19,009/19,011`), with only two
observed at `g_w<=0`.
So a practical proof decomposition is:
  (A) strict mode gap for `g_w<=0`,
  (B) non-strict mode bound for `g_w>=1`.

Leaf branch sharpening (proved reduction):
for a leaf `w` with neighbor `u`, writing `H=T-w`,
`I(T) = I(H) + x I(H-u)` and equivalently
`I(T) = (1+x)I(H-u) + xI(H-N_H[u])`.
If one proves leaf-step descent monotonicity
`d(I(H+leaf_at_u)) >= d(I(H))`, then in minimal-counterexample mode
(`H` unimodal) this forces strict inequality
`mode(I(T-w)) < d(I(T))`.
So leaf vertices reduce to one local descent-index monotonicity statement
rather than a full mode argument.
See `notes/leaf_step_descent.md` for the exact local formulation.

Empirical leaf diagnostic (exhaustive through `n<=16`):
`results/leaf_mode_descent_n16.json` reports:
- `244,692` leaf-cases,
- `0` failures of `mode(I(T-w)) <= d(I(T))`,
- `0` failures of strict `mode(I(T-w)) < d(I(T))`,
- `0` failures of `d(I(T)) >= d(I(T-w))`.

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
k = t-2, t-1 satisfy \(\Delta I_k \le 0\) automatically. A concrete cutoff is

  \(s \ge s_0(d,e) := \max(e+1-d,\ 2d+11,\ \lceil(3e-2d+8)/2\rceil)\),

where \(d=\deg I(H-v)\), \(e=\deg I(H-N[v])\) for the hub v after removing its
pendant leaves. Therefore any minimal counterexample must be leaf-light at each
hub, reducing the problem to a finite core class plus bounded leaf-load checks.

**Finite-core enumeration target (conditional):** if subdivision preservation
holds (no degree‑2 vertices), and all hubs satisfy \(s(v)\le \lambda_0\), then
each tree is determined by a core \(K\) on branch vertices plus a bounded leaf
load \(\ell_u\in[\max(0,3-\deg_K(u)),\lambda_0]\). Once a bound on the number of
obstruction edges \(f_0\) is proved, the class is finite with
\(b\le 2f_0+2\) and \(n\le(\lambda_0+1)(2f_0+2)\), so it is brute‑verifiable.

Empirical check of this finite-core family (no dedup, networkx cores):
  - \((b_0,\lambda_0)=(6,4)\): 13,331 candidates, 0 non‑unimodal
    (`results/finite_core_b6_l4.json`).
  - \((b_0,\lambda_0)=(6,5)\): 58,579 candidates, 0 non‑unimodal
    (`results/finite_core_b6_l5.json`).
  - \((b_0,\lambda_0)=(7,4)\): 92,207 candidates, 0 non‑unimodal
    (`results/finite_core_b7_l4.json`).
  - \((b_0,\lambda_0)=(7,5)\): 513,699 candidates, 0 non‑unimodal
    (`results/finite_core_b7_l5.json`).
  - \((b_0,\lambda_0)=(8,4)\): 711,191 candidates, 0 non‑unimodal
    (`results/finite_core_b8_l4.json`).

---

## Cross-note status (proved vs missing)

1) `notes/subdivision_lemma.md`:
   - proved exact identity \(I(T') = I(T) + Q_uQ_v + xP_uP_v\),
   - proved coefficientwise bound \(A \le x I(T)\),
   - boundary inequality at \(k=d(I)\) is false (small star counterexample),
   - tail dominance \(\Delta A_k \le -\Delta I_k\) for \(k\ge d+1\) holds
     empirically through n=19 but lacks a proof.

2) `notes/lamport_transition.md`:
   - proved restricted mode-ordering lemma (broom root + leaf child),
   - found global counterexamples to full difference-dominance,
   - observed boundary-only failures in the scanned broom-leaf table.

3) `notes/steedman_invariants.md`:
   - global “\(B_T\) log-concave for all rooted trees” is blocked by the
     construction in `notes/subdivision_lemma.md`,
   - the viable route is class-restricted closure plus one boundary inequality
     for the final sum.

Common gap across all three routes: a local inequality controlling the tail
(after first descent) plus a structural reduction to the leaf-heavy obstruction
class.

---

## “If I had to bet” (Erdős heuristics)

1) Prove subdivision preserves unimodality (or a weaker local variant).
2) Use Lamport transition + difference dominance on a small extremal core.
3) Push MBI + leaf-heavy reduction, then seal the last boundary indices with a
   direct binomial bound in that reduced class.

One clean inequality is worth more than a broad framework. That’s the Erdős way.
