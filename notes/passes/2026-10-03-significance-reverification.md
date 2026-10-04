# Re-verification after adding the significance paragraph
<!-- SUMMARY: checks after the 'why care' paragraph and abstract clause were added, with space trims to stay at 13 pp · status: source claims verified; clean build · updated: 2026-10-03 -->

**Prompt.** Brett asked where the paper says why the result matters. It
didn't. All seven ECP comparators state it (by a downstream use, a broader
interest, or a known question); ours stated none.

**Added.**
- A §1 paragraph after the question: quantitative log-concavity has
  probabilistic consequences.
  - Johnson's c-log-concavity, a discrete curvature condition, gives a
    modified log-Sobolev inequality and concentration.
  - Marsiglietti and Melbourne derive concentration, moment and entropy bounds
    from relative log-concavity.
  - Johnson's constant can never exceed 1/EW, so it falls far below 1/V when
    many p_i are near 1.
  - At D, Johnson's quantity is (f_D/f_{D+1})δ_D ≥ δ_D ≥ 1/(4V) when V ≥ 1.
- An abstract sentence on the same point.

**Source checks** (read before writing):

| Claim | Source and location |
|---|---|
| c-log-concavity: E(V)(x) = V(x)/V(x+1) − V(x−1)/V(x) ≥ c | Johnson 2017, Definition 1.2 and Condition 1 |
| Called a discrete curvature term | Johnson 2017, §1 ("understand E(V)(x) as a discrete curvature term") |
| Modified log-Sobolev inequality under c-log-concavity | Johnson 2017, Theorem 1.3 (support Z+; finite support by perturbation, Remark 6.2) |
| Concentration of measure as a consequence | Johnson 2017, abstract and §8 |
| c ≤ 1/EV | Johnson 2017, Lemma 5.3 |
| Concentration, moment and entropy bounds from a precise degree of log-concavity; relative log-concavity is the key notion | Marsiglietti–Melbourne 2026, abstract and Definition 1.3 |
| E(f)(D) = (f_D/f_{D+1})δ_D, and f_D/f_{D+1} > 1 at the first descent | algebra; log-concavity gives f_{D+1}/f_D ≤ f_D/f_{D−1} < 1 |

The paragraph claims nothing about Theorem 1.1 plugging into those theorems:
their conditions are global, and ours holds at D.

**Space.** The additions pushed the paper to 14 pages (3 lines over). It is
back to 13 with no content cut:
- the 12 margins in §4.1 set inline;
- the definition of p_n in the proof of Prop 1.2 set inline;
- the supplement's opening line merged into its first paragraph;
- the one-line determinantal-process sentence moved into Example 1.9;
- one sentence compressed;
- `\bigskip` before the acknowledgement changed to `\medskip`.

Page 13 is now full.

**Build.** 0 undefined references, 0 overfull boxes, 13 pages. The portal
abstract is re-synced.

**Follow-up (side-agent note).** The first version of the abstract sentence ("Quantitative log-concavity yields concentration and log-Sobolev inequalities …, but … the known bounds depend on …") could be read as saying the theorem delivers those inequalities. It now reads "is used to derive" and names the gap precisely: "the known lower bounds on δ_k". Still 13 pp.
