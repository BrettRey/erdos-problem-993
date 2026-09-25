# Fang, Lu, Nevo, Yao, Zheng (2026): unimodality for sufficiently large forests
<!-- SUMMARY: Reading note and local Lean replay for arXiv:2609.20961, which proves #993 for all forests on at least an unspecified N0 vertices; replay at commit b2a1d3e builds and audits to the three standard axioms; the finite range below N0 stays open · status: verified-replay, no project action beyond logging · updated: 2026-09-25 -->

Ethan X. Fang, Junwei Lu, Eran Nevo, Yuan Yao, Hailun Zheng,
*Unimodality of Independence Polynomials for Sufficiently Large Forests*,
arXiv:2609.20961v1 [math.CO], posted 17 September 2026, 18 pp.
Local copy: `literature/fang_2026_unimodality_large_forests.{pdf,md}`
(portfolio root). Filed on erdosproblems.com/993 as a *partial* proof claim
on 2026-09-21 by user `yyao`, "using Odin Automatic AI Research Agent".
Found 2026-09-25: the problem page's proof-claim count had gone from 0
(snapshot of 2026-09-04, `tmp/frontiermath-erdos/erdos-993.txt`) to 1.

## What the paper claims (page numbers from the arXiv PDF)

- **Theorem 1.1 (p. 2).** "There is an absolute positive integer N0 such that
  the independence sequence of every forest with at least N0 vertices is
  unimodal." N0 is never computed.
- **Theorem 1.2 (p. 2).** For every sufficiently large n and every forest on n
  vertices, `i_k^2 > i_{k-1} i_{k+1}` for `n/5 <= k <= 17 alpha/25`.
- **Method (pp. 3-4, Sections 3-7).** A central limit theorem for the size of
  a hard-core-model independent set, uniform over all forests and over
  fugacities in `[1/4, 12]`, via root conditioning, a centroid decomposition
  down to components of size `ceil(n^{1/4})`, a Berry-Esseen bound, and a
  characteristic-function bound `|E e^{itX}| <= exp(-c n sin^2(t/2))`,
  then a local limit step giving strict log-concavity where the mean is an
  integer.
- **Ends (Proposition 8.2, p. 15).** For every forest, `i_k` is nondecreasing
  through `ceil(n/4)`, from `(k+1) i_{k+1}/i_k >= n - 3k` (extension double
  counting, Lemma 8.1, p. 14). The tail is nonincreasing from
  `ceil((2 alpha - 1)/3)` (Levit-Mandrescu; they reprove the bipartite case).
- **Patching (p. 16).** For large n the three ranges overlap, which gives
  unimodality.
- **Scope, in their words (p. 16).** "We emphasize that the original
  conjecture for trees and forests of all sizes remains open."
- **Open questions (pp. 16-17).** The peak-location constant `a_p`
  (they show `a_p >= 1/4`, `b_p = 2/3`), and the log-concavity interval
  `[a_lc, b_lc]`: `a_lc = 0`, `b_lc >= 17/25`, and "for all examples we know
  of, the failures to log-concavity appear only at degree (1-o(1)) alpha(F),
  so it is possible that b_lc = 1" (p. 17).
- **AI role (p. 3).** Odin found the initial proof; the authors wrote and
  verified the final text.
- **What they cite for small n (p. 2).** Only the n <= 25 log-concavity
  verification "communicated by Radcliff" via Ball-Galvin-Hyry-Weingartner.
  Neither Zenodo 19100781 nor the erdosproblems thread's n <= 32 unimodality
  record is cited.

Whether N0 can be made explicit is not discussed. Section 9 notes the
constants are not optimal, and the Berry-Esseen step can use Raič's explicit
constants (p. 9). How large an explicit N0 would be is **unknown**; the
guess that it would lie far beyond exhaustive enumeration is ours, unchecked.

## Local Lean replay (2026-09-25)

Repository `github.com/junwei-lu/Erdos_993_Tree_Independent_Set_Unimodality`,
commit `b2a1d3ede8aef259b1de6e319e7fd6cb56481ac1` (2026-09-21), Lean
`v4.29.1`, Mathlib `v4.29.1` (upstream cache). 19 files, 12,161 lines.

- `lake build`: "Build completed successfully (8267 jobs)", zero
  "declaration uses 'sorry'" warnings (32 linter warnings, all unused
  variables or section variables).
- `lake env lean ErdosProblem993/AxiomCheck.lean`: all 14 audited theorems,
  including `main_fin`, `unimodal_of_isAcyclic`, `central`,
  `monotoneOn_icoeff_prefix`, `antitoneOn_icoeff_tail`, and
  `isUnimodal_iff_finite`, depend on `[propext, Classical.choice, Quot.sound]`
  only.
- Source grep for `sorry`, `admit`, `native_decide`, `axiom`,
  `implemented_by`, `unsafe`, `extern`, `partial def`, `opaque`: no code
  hits. Comment hits only, including a stale docstring in `HardCore.lean`
  referring to "`sorry`-ed statements of `ErdosProblem993.Recursion`";
  `Recursion.lean` contains no `sorry` token and the build reports none.
- Statement faithfulness, checked by reading `Basic.lean` and `Main.lean`:
  `icoeff G S k` counts `k`-subsets of `S` that are `G.IsIndepSet` (Mathlib);
  `IsUnimodal f` is "exists `m` with `f` monotone on `[0,m]` and antitone on
  `[m, inf)`"; the hypothesis is Mathlib's `SimpleGraph.IsAcyclic` (a
  forest); `main_fin` is `exists N0, forall n >= N0, forall G on Fin n,
  G.IsAcyclic -> IsUnimodal (icoeff G univ)`. Neither `IsAcyclic` nor
  `IsIndepSet` is redefined. The proof sets `N0 = max 1000 N1` with `N1`
  from `central`, so no numerical value is extracted.

Verdict: Theorem 1.1 as stated is machine-checked. This does not bear on
trees below N0.

## Where it touches this project

- **Free-vertex identity.** Their (8.1)/Lemma 8.1 identity
  `(k+1) i_{k+1}/i_k = E e(J)` is the `mu_k` free-count identity in
  `notes/mason_free_count_reformulation_2026-09-02.md` (line 38). It is a
  standard double count; neither side has a priority claim.
- **Deterministic prefix bound.** `i_k` nondecreasing through `ceil(n/4)`
  for every forest is not in `paper/main_v2.tex`, which cites only
  random-tree prefix results (Basit-Galvin, Heilman; `main_v2.tex:103`).
  Any rebuild should cite it.
- **Mode-mean route.** Theorem 1.1 bypasses, for large n, the mode-mean
  conjecture on which `main_v2`'s `mu(T) < n/3` bound was meant to feed.
  The mean bound keeps its interest as a finite-n and peak-location
  statement; their `a_p` question is the natural home for it.
- **Log-concavity interval.** Our complete log-concavity census for
  `27 <= n <= 32` (`results/lc_census_20260814/summary_n{27..32}.json`:
  0, 19, 7, 121, 159, 922 failing trees, all complete against gentreeg
  counts), posted to the erdosproblems thread on 2026-08-14, found every
  failure at exactly `k = alpha - 1`. That is direct small-n evidence for
  their remark that `b_lc` may be 1.
- **Finite range.** The unimodality record on the erdosproblems thread is
  `n <= 32` (tylersatchelorden, 2026-08-09), superseding our `n <= 29`
  (`outreach/erdosproblems_993_comments_2026-08-12.md`, lines 5-6). The open
  part of #993 is now the range between that record and an unspecified N0.
- **Depth-three window lane.** The bounded window `33 <= n <= 38` remains
  logically open, but a proof there no longer serves an asymptotic argument.
  Its value has dropped; do not resume it on this basis.
