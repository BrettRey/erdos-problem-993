# Independent review and Lean verification: both graph theorems pass

The two clean-room mathematical reviewers found the b1-zero tree-window theorem
and the low-b4-density forest theorem valid in their intended scope. They each
recomputed all 13 parameter cases and all 1,842 rooted-forest coefficient rows.
Their reports are [review A](math-a.md) and [review B](math-b.md).

At review time, the [Lean audience-reader](lean-reader.md) found a substantive formalization
interface still missing from the existing library: its Pascal reserve counts
projected matching-bag words, not yet the actual graph-set counts used by the
new theorems. The rank/maximum-set identification and the cardinality-preserving
partial-set/code bijection must be proved. The reader freshly recompiled the
twelve-module reserve dependency closure; the reserve and representation
declarations depend only on `propext`, `Classical.choice`, and `Quot.sound`.
This supports the existing code theorem, not a completed graph-count bridge.

## What the agreement establishes

The strongest convergence is on independently recomputed facts: the 13-case
split, all 1,842 coefficient rows, the maximum `52513/217404`, its positive slack
`3851/217404`, and 17 failures of the stronger discarded endpoint test.
Both mathematical reviewers also reconstructed the splicing, concentration,
root-switch injection, and incidence arguments. Neither found a broken
inference. This is independent AI review, not human peer review or field
consensus. All reviewers shared the frozen artifact, provider, and inherited
model family; their agreement is not three model-independent votes.

There is no substantive contradiction between the reports. The different
emphases matter: mathematical reviewers accepted the written graph/code
correspondence; the Lean reader correctly required a declaration proving it.
Mathematical validity and completion of the formalization are separate claims.
Both mathematical reviewers found the same connectedness boundary; they gave
different independently checked counterexamples to removing the forest
hypothesis from the general blocked-shadow lemma.

## Original revisions and verification obligations (historical)

1. Prove the graph-count/code interface, including the `d <= alpha` guard.
2. Retain connectedness for the sharp four-forbidden bound. The bound fails for
   the disjoint union of the 13-vertex core and ten isolated edges, although
   that example does not refute the stated tree theorem.
3. Package the second theorem explicitly: finite forest, alpha in {17,18,19},
   `3(alpha-3)b4 <= (alpha-7)e4`, strict full margin. State its auxiliary shadow
   lemma with alpha >=4. Define D and justify positive denominators locally.
4. Formalize the graph-to-corona/rooted-forest reduction, evaluator correctness,
   and enumeration coverage. Checking imported JSON rows alone is insufficient.
5. Keep the replay dependencies with the finite packet. The review snapshot
   uses some helpers/archives from the live repository, as all three reports
   disclose. The new Aristotle finite packet includes its needed Python
   helpers; its 1,842-case replay succeeded standalone.

## Current Lean execution status

The frozen exact target definitions compile. Locally proved scalar assembly
lemmas `DepthThree.low_density_algebra` and
`DepthThree.four_forbidden_algebra` have the standard three axioms only. These
are conditional algebraic implications, not the graph theorems.

Subsequent local work also kernel-checked all twelve small parameter cases,
the five-ray cone lemma, and a universal bound on recursively encoded rooted
forests of size ten. The latter includes enumeration coverage and an exact
count of 16,796 ordered representatives; it does not trust the 1,842 JSON rows.
At that stage, its graph/encoding/evaluator correspondence remained to be proved. The assembled
supporting packet passed a full default Lake build. These later results are
main-agent Lean checks, not findings of the original independent reviewers.

The low-density graph packet has now passed local verification. Aristotle
project `be8de245-5246-41b8-aa5e-ca01f73d8008`, task
`2d0cdd1e-7763-4120-a7d1-e327f5f9a2c4`, returned both exact graph targets.
The main agent freshly rebuilt all source (8,049 jobs), checked the unchanged
specification, read all six new modules, and independently audited the target
and bridge axioms: only `propext`, `Classical.choice`, and `Quot.sound`.
The actual graph/count/code bijection is complete. The provider's
COMPLETE_WITH_ERRORS label disagreed with its COMPLETE report; the local
source replay, not either label, establishes verification. Full record:
[low-density verification](../../formalization/depth3_low_density_20260906/STATUS.md).

The b1-zero graph target has now also passed local verification. The saved
checkpoint from task `a5f2e62f-3d9e-49e2-a9e6-ce29a5673346` contains
`DepthThree.b1_zero_window : DepthThree.B1ZeroWindowTarget`. Although that task
exhausted its budget, the main agent freshly rebuilt the complete source
(8,078 jobs), checked the unchanged frozen specification, and compiled exact
ascriptions for all three graph targets. Their axioms and those of the audited
structural/representation/certificate bridges are the standard three only.
No proof holes or native-decide shortcuts occur in the returned project.

The final proof establishes a weaker sufficient r=4 density bound than the
written proof's sharp maximum: a union-bound numerator is at most `e_4/4`,
proved by kernel reduction over all 16,796 plane rooted forests with formal
coverage and actual graph/evaluator correspondence. This suffices through the
low-density theorem. The sharp `52513/217404` graph bound and the special
`b_3 >= 11 b_2` estimate remain independently reviewed written results, not
newly Lean-verified graph theorems. Full record:
[b1-zero verification](../../formalization/depth3_b1_zero_20260906/STATUS.md).

The [independent follow-up source review](lean-checkpoint-followup.md) read all
21 new B1Zero modules and the common graph/count layer. It found the exact
target and all required semantic bridges intact. The reviewer did not run a
second build; its source-semantic verdict and the main agent's successful
kernel replay are separate checks. It identified an overbroad statement in
the provider's `DECLARATION_MAP.md`: claims of no connectedness, vertex-count,
deficiency, or b1 restrictions apply to the low-density/blocked-shadow layer,
not the whole checkpoint. The frozen archive is preserved unchanged and this
qualification is recorded in the local verification record. The original
review reports are also unchanged.

Exact prompts, the pre-review source snapshot and commit/hash, complete raw
reports, and dependence clusters are preserved in [the manifest](manifest.yaml).
No novelty claim, manuscript revision, commit, or external publication follows
from this review.
