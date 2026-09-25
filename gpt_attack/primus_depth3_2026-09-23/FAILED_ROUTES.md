# Guardrails and optional leads

## Established obstructions: do not assume these shortcuts

1. **The blocked profile is log-concave.** False even in the target window.
   The `blocked_lc_failure` witness in the lift certificate has
   (b_2,b_3,b_4)=(1,16,4216), giving b_3^2-b_2b_4=-3960. Its full depth-three
   margin remains positive. Replay it before proposing this shortcut again.
2. **A defect-one obstruction makes the correction nonnegative.** False for
   both pair-only and unary-only defect-one examples. The two witnesses in
   `b1_positive_sign_obstructions_20260905.json` have negative D but positive
   full margins. They are not counterexamples to W or #993.
3. **Disjoint certificate classes eliminate adverse interactions.** A disjoint
   partition solves duplicate counting, not automatically the quadratic
   cross-term problem. Any proposed classwise estimate still has to control
   interactions when the classes are added.
4. **A stronger endpoint estimate can be imposed everywhere.** The historical
   lift certificate includes failures of that estimate outside D<0. Keep the
   primary target W, the sufficient joint bound S, and the stronger discarded-
   positive-terms bound E separate; their definitions are in `MATHEMATICS.md`.
5. **Real-rootedness or global log-concavity is available for tree independence
   polynomials.** Neither is a permitted general hypothesis. This packet targets
   only one bounded local coefficient inequality. Mean bounds likewise do not
   automatically localize modes or establish unimodality.
6. **Adding vertices repairs blocked sets.** If S is contained in no maximum
   independent set and H is an independent superset of S, then H is also
   blocked: a maximum set containing H would contain S. A proposed repair into
   the extendable family must therefore do something other than only add
   vertices. This observation does not rule out exchanges or maps on pairs.

The three fixed replay witnesses are regression tests, not an exhaustive
description of the residual regime R. In particular, do not call them examples
in R without checking the density condition as well as b_1 and D.

## Optional external leads, assessed 19 September 2026

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
