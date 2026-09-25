# Primus handoff: tree independence sequences

Prepared for Brett Reynolds, 23 September 2026.

This is the uncompressed replacement for the ZIP handoff. Paste this entire
text into Primus, or attach this single file if its format is accepted. No ZIP
extraction is needed to understand the task.

This text includes the assignment, definitions, statements of established
lemmas, exact witness data, and guardrails. It does NOT include the full Lean
project or executable checking code. Source paths mentioned below identify the
original project's files; they are references, not instructions to open files
that have been attached. Ask the user for specific sources if you need them.
Do not claim to have compiled or replayed proof sources that you have not seen.

The approved mathematical scope is unchanged: one bounded depth-three window,
not an unrestricted attack on all of Erdős Problem #993.

## Assignment: the residual depth-three window for tree independence sequences

Work on one precise open target, with this self-contained text as your starting
point. The motivation is Erdős Problem #993, but this assignment is not a
request to claim that the full problem follows from a single remaining lemma.

For every finite simple tree T with n vertices, independence number alpha,
33 <= n <= 38, alpha in {17,18,19}, and 2 alpha - n <= 5, prove or refute

    s_3^2 >= s_2 s_4,

where s_d is the number of independent sets of size alpha-d.

Read the mathematical brief below before starting. It gives intrinsic extendable and
blocked counts e_d and b_d, with s_d=e_d+b_d, and the exact remaining regime.
Established subcases are reported in the project's Lean records: b_1=0 and low
b_4/e_4 density. Their statements are below; their proof source is not attached.
The remaining branch has b_1>0, high density, and a negative combined correction
D. The Pascal reserve and blocked-shadow bound are already available. Do not
spend the run re-proving those subcases unless you find a specific defect.

Three exact witness records are included below. If you use computation,
independently reconstruct and check them rather than merely quoting their
stored counts. Treat a proposed new inequality as
falsifiable: check its hypotheses and adversarial examples before promoting it.
Use exact integer/rational arithmetic for signs. No float64 root-finding.

You may use a different proof route. A sufficient bound that the Pascal margin
absorbs D would be useful, but it is not the primary target. Refuting that bound
does not refute the target. Refuting the target would not by itself refute
unimodality or solve #993. Do not assume blocked-profile log-concavity,
componentwise nonadversity, real-rootedness, mode–mean localization, or an
unproved Holant representation. See the guardrails below.

Return one of the following grades, supported by the corresponding artifacts:

- **PROVED:** a complete proof of the unchanged universal window target, with
  every reduction justified. Stronger strict positivity is welcome but optional.
- **REFUTED:** an exact graph6/edge-list witness in the stated window, its full
  independence polynomial, and a standalone exact verification. Explicitly
  distinguish target refutation from any auxiliary refutation.
- **PARTIAL:** at least one new proved lemma or precisely delimited structural
  reduction, with its usefulness and remaining gap stated. Do not label a larger
  finite sample as a proof.
- **NO ADVANCE:** report failed approaches and exact obstructions honestly.

The deliverable should contain a concise written argument, a status/dependency
table, scripts and certificates, exact reproduction commands and versions,
and any Lean source/build/axiom audit. Separate machine-checked proofs, written
proofs, exact finite checks, and conjectural evidence. The mathematical target
must remain unchanged; no new hypothesis may be silently added.
Do not replace a graph theorem with a scalar implication that assumes its hard
counting inequality.

Prefer one verifiable advance over a broad survey or a speculative paper. Work
within the run budget selected by the user, and return partial results rather
than starting additional paid runs, purchasing resources, or publishing work.

## Exact mathematical brief

### 1. Objects and the primary target

All graphs here are finite, simple, and undirected. A tree is connected and
acyclic; a forest need not be connected. An independent set contains no adjacent
vertices. A **maximum** independent set has the greatest possible size alpha;
do not confuse maximum with inclusion-maximal.

Write

\[
I(T;x)=\sum_{k=0}^{\alpha}i_k(T)x^k,
\qquad s_d=i_{\alpha-d}(T).
\]

The target, denoted W, is

\[
\boxed{s_2s_4\le s_3^2}
\tag{W}
\]

for every tree satisfying

\[
33\le n\le38,\qquad 17\le\alpha\le19,\qquad 2\alpha\le n+5.
\]

For trees, delta=2 alpha-n is the matching deficiency, hence nonnegative. The
displayed conditions leave twelve (alpha,delta) cells: alpha=17 with delta=0,1;
alpha=18 with delta=0,1,2,3; alpha=19 with delta=0,1,2,3,4,5.

This is one local log-concavity inequality near the top of a polynomial in a
bounded order window. It neither establishes all coefficient positions nor
arbitrary tree orders. A failure of W is not automatically a failure of
unimodality: a positive sequence can be unimodal without being log-concave.

### 2. Intrinsic counts; no multiplicities of completions

An independent set S is **extendable** if S is a subset of at least one maximum
independent set of T. Otherwise it is **blocked**. Let e_d and b_d count these
two types at size alpha-d. Set the counts to zero when d>alpha. Then

\[
s_d=e_d+b_d,\qquad b_0=0.
\]

Each set is counted once, not once for every maximum completion. The original formal definitions are in
`lean/RequestProject/DepthThreeSpec.lean` (source available on request).
These counts do not depend on a chosen matching or root.

Define the signed correction over the integers:

\[
D=(b_3^2-b_2b_4)+(2e_3b_3-e_2b_4-e_4b_2).
\]

Expanding s=e+b gives the exact identity

\[
s_3^2-s_2s_4=(e_3^2-e_2e_4)+D.\tag{1}
\]

In Lean, cast to integers before subtracting; natural-number subtraction would
erase the negative correction that this task is about.

### 3. Available results

The project reports the following established results, historically verified on
6 September 2026. Their full proof sources are not reproduced in this text.
For a written continuation you may cite them explicitly as supplied prior
lemmas. Do not describe them as independently verified in your run. If a proof
dependency needs auditing or a complete Lean replay is required, request the
relevant source files from the user.

#### A. Extendable Pascal reserve

For the forest counts above with alpha>=4,

\[
27(\alpha-3)e_3^2\ge32(\alpha-2)e_2e_4,
\]

and therefore

\[
e_3^2-e_2e_4\ge c_\alpha e_2e_4,
\qquad c_\alpha=\frac{5\alpha+17}{27(\alpha-3)}>0.\tag{2}
\]

Also e_d>0 for 0<=d<=alpha. The representation/count bridge is proved; this is
not an assumption that a graph happens to have a suitable code.

Relevant source: `RequestProject/PascalBridge.lean`,
`RequestProject/TreeCodeBridge.lean`, `RequestProject/ExtendableCount.lean`, and
the use of `MatchingBag.TreeMatching.erasure_depth_three_reserve` in
`RequestProject/LowDensity.lean`.

#### B. Blocked shadow and extendable incidence

For every finite forest with alpha>=4,

\[
6b_3\ge(\alpha-4)b_2,\qquad
(\alpha-3)e_3\ge4e_4.\tag{3}
\]

The graph-level blocked-shadow declaration is
`DepthThree.blocked_shadow_of_forest : DepthThree.BlockedShadowTarget`.
The extendable incidence declaration is
`MatchingBag.TreeMatching.extendable_incidence`.

#### C. Low-density closure

For every finite forest with alpha in {17,18,19},

\[
3(\alpha-3)b_4\le(\alpha-7)e_4
\quad\Longrightarrow\quad s_2s_4<s_3^2.\tag{4}
\]

The b_4/e_4 thresholds at alpha=17,18,19 are respectively 5/21, 11/45,
and 1/4. No b_1=0 condition is needed.

Declaration: `DepthThree.low_density_of_forest : DepthThree.LowDensityTarget`.

#### D. The b_1=0 window

Every tree in W's window with b_1=0 satisfies the strict inequality
s_2s_4<s_3^2.

Declaration: `DepthThree.b1_zero_window : DepthThree.B1ZeroWindowTarget`.
The proof establishes the actual graph theorem, including a structural
reduction and a finite enumeration coverage argument. It does not merely
assume a list of coefficient rows is complete.

The formal proof uses a sufficient one-quarter density bound in its exceptional
case. Do not attribute the sharper historical written bound 52513/217404 or
the special estimate b_3>=11b_2 to this Lean theorem; they are not needed here.

### 4. What is left

If D>=0, equations (1) and (2) settle W. Results C and D settle the low-density
and b_1=0 branches. Thus it suffices to handle trees in the same window with

\[
\boxed{b_1>0,\quad
3(\alpha-3)b_4>(\alpha-7)e_4,\quad D<0.}\tag{R}
\]

Do not assume R is empty. A historical bounded search that failed to find an
example in R is not an emptiness proof. Proving R impossible would be a valid
solution, but requires a structural or coverage argument.

One sufficient target is

\[
(5\alpha+17)e_2e_4+27(\alpha-3)D\ge0.\tag{S}
\]

It implies W by (1)-(2), but is potentially stronger than W. An exact failure
of S with positive full margin refutes only this proposed route. A still
stronger endpoint estimate discards the favorable terms 2e_3b_3+b_3^2:

\[
(5\alpha+17)e_2e_4\ge
27(\alpha-3)(e_2b_4+e_4b_2+b_2b_4).\tag{E}
\]

E is not known universally and historical certificates contain failures outside
the adverse regime. Do not silently replace S by E or require E everywhere.

The original package's `lean/PrimusSpec.lean` defines W, the residual target
under R, and S as open propositions. It is not attached here and contains no
proof of any of them. The displayed statements in this text are the task.
W's non-strict conclusion is intentional; the established subcases happen to
give strict positivity.

### 5. Computational representation for independent checking

Choose a maximum matching of a forest. Each matched edge forms a two-vertex
bag; each unmatched vertex forms a singleton bag. There are alpha bags. A
maximum independent set selects exactly one vertex per bag, subject to the
original adjacency constraints. Contracting the matched edges gives a forest
of bags.

An extendable partial independent set corresponds to a partial bag assignment
that admits a full assignment. The original package's `extendable_defects` routine (not attached here) uses
dynamic programming over the bag forest, retaining the **set of feasible root
choices**, not a count of completions. Empty bags contribute one unit of
defect. Grouping by feasible-choice sets avoids counting the same partial set
several times. Packaging tests compared this routine against direct independent-set
enumeration on small examples. You may implement this description independently,
or request the original checking code.

The independence polynomial is computed separately by the ordinary rooted-tree
recurrence: excluding a root multiplies the total child polynomials; including
it contributes x times the product of the child-excluded polynomials. The
blocked counts are s_d-e_d. All decisive replay calculations use arbitrary-
precision integers. These numerical tools are evidence and counterexample
checkers, not a substitute for a universal proof.

### 6. What would count as an advance

A complete W proof, an exact W counterexample, or a new proved lemma that
genuinely constrains R would count. A useful partial result could be an
explicit exchange/counting bound with proved multiplicities and constants,
or a smaller exhaustive model with a proved coverage map.

For any new representation, prove that it counts these intrinsic objects
without multiplicity and preserves the hypotheses and target. For any finite
certificate, separate exact arithmetic from completeness of its enumeration.
For any external theorem, state and discharge its actual hypotheses.

## Guardrails and optional leads

### Established obstructions: do not assume these shortcuts

1. **The blocked profile is log-concave.** False even in the target window.
   The `blocked_profile_not_log_concave` witness below has
   (b_2,b_3,b_4)=(1,16,4216), giving b_3^2-b_2b_4=-3960. Its full depth-three
   margin remains positive. Check it before proposing this shortcut again.
2. **A defect-one obstruction makes the correction nonnegative.** False for
   both pair-only and unary-only defect-one examples. The `pair_only_defect_one_adverse` and
   `unary_only_defect_one_adverse` witnesses below have negative D but positive
   full margins. They are not counterexamples to W or #993.
3. **Disjoint certificate classes eliminate adverse interactions.** A disjoint
   partition solves duplicate counting, not automatically the quadratic
   cross-term problem. Any proposed classwise estimate still has to control
   interactions when the classes are added.
4. **A stronger endpoint estimate can be imposed everywhere.** The historical
   lift certificate includes failures of that estimate outside D<0. Keep the
   primary target W, the sufficient joint bound S, and the stronger discarded-
   positive-terms bound E separate; their definitions are in the mathematical brief.
5. **Real-rootedness or global log-concavity is available for tree independence
   polynomials.** Neither is a permitted general hypothesis. This packet targets
   only one bounded local coefficient inequality. Mean bounds likewise do not
   automatically localize modes or establish unimodality.
6. **Adding vertices repairs blocked sets.** If S is contained in no maximum
   independent set and H is an independent superset of S, then H is also
   blocked: a maximum set containing H would contain S. A proposed repair into
   the extendable family must therefore do something other than only add
   vertices. This observation does not rule out exchanges or maps on pairs.

The three fixed witnesses below are regression tests, not an exhaustive
description of the residual regime R. In particular, do not call them examples
in R without checking the density condition as well as b_1 and D.

### Optional external leads, assessed 19 September 2026

These are source pointers, not additional assumptions, tasks, or included
third-party papers. Read the exact relevant argument before relying on it.

- **Chen, Chen and Zhang, degree-free Holant spectral independence:**
  <https://arxiv.org/html/2609.18835v1>. Section 4.3 uses feasible insertions and
  bounded weighted preimages. This is a modest conditional counting-method
  lead. The natural independent-set incidence encoding requires equality
  signatures with internal zeros, outside its theorem's hypotheses. An
  eligible encoding would still need a bridge from its variance conclusion to
  our coefficient inequality. The insertion-only repair obstacle above also
  prevents a literal blocked-to-extendable transplant.
- **Xie and Zhang, infinite log-concavity of Boros–Moll sequences:**
  <https://arxiv.org/html/2609.20653v1>. Proposition 2.5 absorbs two adverse Abel-
  summation terms using specific coefficient bounds and a uniform partial-sum
  estimate. This is a proof-design reference; the required estimates and
  Jacobi/Narayana representation have not been supplied for our counts.
- **Jochemko and Menon, weighted lecture-hall enumeration:**
  <https://arxiv.org/html/2609.17250v1>. Refined interlacing is preserved by
  operations on a different enumerator. No mapping to our recurrence has been
  established. Background only unless an explicit count-preserving bridge is
  found.
- **Liu and Zhang, IDP simplices of prime normalized volume:**
  <https://arxiv.org/html/2609.19637v1>. A restricted positive unimodality result
  using finite-group structure and shifted symmetry. No corresponding
  structure for this target has been identified. It does not reverse general
  IDP counterexamples or imply tree unimodality.

Do not turn these four pointers into four parallel research programmes. Use a
source only if a specific missing step warrants it.

## Exact witness data

Each graph6 string describes a finite tree with vertices numbered from zero.
The arrays give counts at defects d=0,1,2,3,4, in that order. The count field
`full_depth3_margin` or `full_margin` is s_3^2-s_2*s_4. Graph6 strings include
literal backticks; use a graph6 parser or store them as JSON strings rather
than interpolating them into a shell command.

All three witnesses were recomputed during packaging on 23 September 2026
using exact integers. All are in the target window and in its already-covered
low-density region, NOT in the remaining high-density regime R. They refute
auxiliary shortcuts; they do not refute the primary target or unimodality.

```json
[
  {
    "label": "pair_only_defect_one_adverse",
    "graph6": "`pG`A?@O??g??@O???I??O???G???G?O???@?A?????_@?????@?A??????G?O??????C?G???????G????O????@",
    "n": 33,
    "alpha": 19,
    "deficiency": 5,
    "e_0_to_4": [
      1544,
      21380,
      139266,
      567167,
      1618881
    ],
    "b_0_to_4": [
      0,
      8,
      117,
      994,
      7141
    ],
    "combined_correction": -56227048,
    "full_depth3_margin": 96167097495,
    "unimodal": true
  },
  {
    "label": "unary_only_defect_one_adverse",
    "graph6": "`oG`A?@O?_g??@O???I??O???G???G?O???@?A?????_@?????@?A??????G?O??????C?G???????G????O????@",
    "n": 33,
    "alpha": 19,
    "deficiency": 5,
    "e_0_to_4": [
      1536,
      21248,
      138240,
      562176,
      1601856
    ],
    "b_0_to_4": [
      0,
      8,
      132,
      2739,
      22070
    ],
    "combined_correction": -178212783,
    "full_depth3_margin": 94423068753,
    "unimodal": true
  },
  {
    "label": "blocked_profile_not_log_concave",
    "graph6": "`pG`A?@O??g??@O???I??O?@?????I?????@O??????g??????@O???????I????????D?????????I?????????@",
    "n": 33,
    "alpha": 19,
    "deficiency": 5,
    "e_0_to_4": [
      4096,
      53248,
      325632,
      1245184,
      3337984
    ],
    "b_0_to_4": [
      0,
      0,
      1,
      16,
      4216
    ],
    "combined_correction": -1336360568,
    "full_margin": 462192427400,
    "blocked_turan": -3960,
    "unimodal": true
  }
]
```

## Evidence status and source availability

The local project records a fresh Lean build and graph-level axiom audit on
6 September 2026. On 23 September, packaging checked that all 52 historical
Lean modules and three pinned configuration files matched their recorded
hashes; reran exact theorem ascriptions and axiom output using existing compiled
dependencies; and elaborated the new open-target specification. That was not a
new clean build of all dependencies. The three established graph targets used
only `propext`, `Classical.choice`, and `Quot.sound` in that audit.

The package's six Python regression tests passed, including direct set-count
comparisons on all 48 nonisomorphic trees of orders 1 through 8, small forest
cases, and the three exact witnesses above. This is evidence about the supplied
implementation and examples, not a proof of the open target.

If needed, request these individual files or their text from the user:

- Historical definitions: `lean/RequestProject/DepthThreeSpec.lean`.
- New open propositions: `lean/PrimusSpec.lean`.
- Full proofs: the original `lean/RequestProject/` source tree, its pinned
  configuration files, and `lean/HistoricalAudit.lean`.
- Numerical replay: `replay.py`, `indpoly.py`, and its three helper scripts,
  or an independently implemented exact checker using the description above.

The original repository snapshot is
`b85d5dc01bdcaf3d6853147b28c350d74d3ede21`. No confidential reviews, credentials,
account information, or submission correspondence are part of this handoff.
