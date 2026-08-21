# Aristotle continuation: full two-branch-vertex (`C₂`) assembly

Continue project `76e84c3b-9ce0-4fc6-9074-45ea38dbb576` after the completed
`ADJACENT_COMPLETE` result. Preserve all existing declarations. The objective is now an
end-to-end Lean theorem that the independence polynomial of every finite tree with at most
two vertices of degree at least three is log-concave.

The existing project already proves, for the explicit adjacent two-hub tree with arbitrary
ordered positive pendant paths and arbitrary total order:

- the idempotent representative partition, with disjointness, exhaustion, unique blocks,
  exact block sizes, and total-order preservation;
- exact normalized block weights, including arbitrary `p ≥ 2` and derived scalars
  `c = 2^e ≥ 1`;
- collision resolution;
- `ClanAudit.sum_normalizedTwoRowCoeff_nonneg` and
  `ClanAudit.Audit.adjacent_two_hub_degreewise_nonneg`.

The attached even-connector project independently proves the only new Laurent inequality:
`EvenConnector.Gpoly_centrallyUnimodal` for all `r,s ≥ 1` and rational `c,d ≥ 1`.
Import or namespace-adapt its `Binomial.lean`, `Coeff.lean`, and `Main.lean` files without
weakening their statements. The attached mathematical packet describes the intended
connector stratification, but it is a proposed proof to audit, not an assumption.

## Target A: direct Li–Li–Yang–Zhang coefficient bridge

Formalize the finite coefficient identity needed to turn the existing clan theorem into
ordinary log-concavity. Do not merely cite Corollary 2.2 or postulate a symmetric-function
identity.

For a finite simple graph `G`, define the independent-set counts `i G j` (or reuse a
Mathlib definition if exact) and prove the diagonal two-row identity

```text
∑ α with ∑v α(v)=2*k, normalizedTwoRowCoeff G α k k
  = (i G k)^2 - (i G (k-1)) * (i G (k+1))
```

over a suitable exact coefficient ring, with the `k=0` boundary handled explicitly. It is
acceptable to prove the more general `(k,l)` Schur/minor identity first. The proof must
derive the identity from the existing multicoloring/clan definitions (equivalently, from
the finite-degree coefficient identity behind
`Σ_α X_G^α = ∏_j I_G(x_j)`), not introduce it as a hypothesis.

Then prove a reusable theorem:

```text
degreewise normalized two-row nonnegativity for G
  -> log-concavity of the independence-count sequence of G.
```

Instantiate it immediately for `dgraph` to obtain an unconditional adjacent-two-hub
independence-polynomial log-concavity theorem. This closes the only formal gap explicitly
recorded in the current `RESULT.md`.

## Target B: arbitrary connector tree

Define an explicit finite tree `connectorGraph ell m a n b` consisting of two endpoint hubs
joined by a path of `ell ≥ 1` edges, with arbitrary finite ordered families of
positive-length pendant paths `a` and `b`. Prove the graph model has exactly the advertised
edges and that its only possible degree-at-least-three vertices are the two endpoints.

Partition all degree-`N` clan maps first by the connector state, and prove—not assume—the
following exhaustive classification.

1. A nonzero map whose positive component joins the two hubs has unit multiplicity on every
   connector vertex. Any multiplicity at least two adjacent to a positive connector vertex
   gives a triangle and normalized two-row contribution zero.
2. If a zero breaks the connector, the remaining nonzero components split into one-hub
   spiders, paths, isolated vertices, or cloned `K₂`s. Give these maps disjoint canonical
   local blocks, with connector prefixes treated as ordinary ordered arms where needed.
3. If the connector is unit-positive end to end, use the existing canonical endpoint
   transforms. Prove that the joined component imbalance is `|p-q|` for odd `ell` and
   `|1-p-q|` for even `ell`, including all endpoint and zero/one-active cases.

For the joined four-state stratum:

- odd `ell` must reduce to the already-proved adjacent `Fblock`, plus only explicitly
  derived good balanced-path factors or summands;
- even `ell` must reduce exactly to `EvenConnector.Gpoly r s c d`, with `r,s ≥ 1` and
  `c,d = 2^e ≥ 1` derived from the clan components.

Prove disjointness, exhaustion, preservation of total order, and exact normalized block
weights for the connector partition. Reuse or generalize `rep`, `repSide`, and the existing
block API rather than replacing the adjacent proof with an unrelated abstraction.

The final connector-level result should be:

```text
∀ ell ≥ 1, ∀ N k l with 1 ≤ l ≤ k and k+l=N,
  0 ≤ ∑ α in the degree-N clan-map space of connectorGraph ell m a n b,
        normalizedTwoRowCoeff ... α k l.
```

Consume Target A to conclude that `connectorGraph ell m a n b` has a log-concave
independence polynomial.

## Target C: graph-theoretic `C₂` wrapper

Prove an explicit classification theorem, or provide a lossless equivalence sufficient for
transport, that every finite tree with at most two vertices of degree at least three is one
of:

- a path (zero branch vertices);
- a spider (one branch vertex);
- an arbitrary-connector two-hub tree of Target B.

Use the direct coefficient bridge plus the existing clan machinery for paths/spiders as
needed. The strongest desired final declaration is invariant under graph isomorphism and
states log-concavity for every finite tree with at most two degree-at-least-three vertices.

If packaging the completely abstract tree wrapper is disproportionately difficult, the
explicit arbitrary-arm `connectorGraph` theorem is still mandatory and must not be replaced
by a bounded family.

## Required audits and grading

- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, bounded enumeration,
  or hidden hypothesis asserting the bridge, partition, weight formula, or desired
  log-concavity.
- Run `lake build`; search all `RequestProject/` sources for escape hatches; run
  `#print axioms` on the coefficient bridge, adjacent log-concavity theorem, connector
  degreewise theorem, connector log-concavity theorem, and abstract `C₂` theorem.
- Update `README.md` and `RESULT.md` with exact declaration-to-obligation traceability.
- Grade `ADJACENT_LC_COMPLETE` only if Target A is proved and instantiated.
- Grade `CONNECTOR_COMPLETE` only if the arbitrary-connector partition and degreewise
  theorem are proved.
- Grade `C2_COMPLETE` only if the final graph-theoretic wrapper and log-concavity theorem
  are also proved.
- If any proposed connector stratum is false, return the smallest explicit collision or
  missing state and identify the exact failed assertion. Otherwise return the strongest
  genuine checkpoint and the first exact unfilled Lean signature.
