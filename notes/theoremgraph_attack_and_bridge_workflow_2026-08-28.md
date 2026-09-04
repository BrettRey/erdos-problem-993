# TheoremGraph-guided attack and bridge workflow

Date: 2026-08-28

## Decision

Use a local claim graph to choose proof attacks and use TheoremGraph to search
for premises and adjacent representations at the graph's open frontier. Do not
rank attacks by thematic novelty, model enthusiasm, or semantic similarity
alone. A high-value attack must hit a load-bearing claim, translate an external
result into the exact local objects, and be cheap enough to kill before proof
investment.

This is a research workflow, not a new mathematical result. The live graph is
`proof_graph/erdos993_frontier.json`; its validator and packet generator are
`scripts/frontier_graph.py`.

## What TheoremGraph changes one level up

The important shift is from **topic search** to **dependency-aware
representation search**.

A topic search asks what is known about independence polynomials, trees, or
unimodality. That is still useful for novelty and direct prior art, but it tends
to return the mathematical neighborhood we already know. A dependency-aware
search starts with a frozen open claim, asks what kind of object the claim is in
several mathematical languages, retrieves exact statements, and then follows
their premises one hop outward. The useful bridge may be a dependency of a
nearby theorem rather than the nearby theorem itself.

For this project, the live objects are not merely trees and coefficient
sequences. They are also:

- partial solutions that fail to extend to an optimum;
- erasure shadows of a binary code;
- violated comparisons in a forest poset;
- minimally unsatisfiable tree CSPs;
- faces outside the maximum-facet-generated part of a nonpure complex;
- unions of overlapping path-obstruction events; and
- a positive sequence plus a correction whose mixed Turán term must fit inside
  quantified slack.

Those representations name the fields worth searching. They provide a reason
to look in coding theory, poset enumerators, CSP/SAT structure, simplicial
face-number theory, and probabilistic dependency inequalities. The graph also
provides a stopping rule: if the translation loses defect grading, assignment
multiplicity, tree path geometry, or the explicit reserve, the bridge is not
load-bearing for the current proof.

## The local frontier

The depth-3 margin now has the exact decomposition

\[
s_3^2-s_2s_4
=(e_3^2-e_2e_4)
+(b_3^2-b_2b_4)
+(2e_3b_3-e_2b_4-e_4b_2).
\]

The matching-bag representation, forest-poset code, blocked-path
classification, and extendable Pascal reserve are proved. The ready frontier
therefore has two nodes:

1. **Blocked-profile depth-3 LC**
   \[
   b_3^2\ge b_2b_4.
   \]
2. **Cross-reserve absorption**
   \[
   (e_3^2-e_2e_4)
   +(2e_3b_3-e_2b_4-e_4b_2)\ge0.
   \]

Both are currently bounded computational observations, not theorems. The exact
probe covers all 13,179 relevant nonisomorphic trees through order 15 and finds
no failure of either inequality. It also finds 223 negative mixed cross terms,
so the second node genuinely uses reserve rather than termwise positivity.

## Selecting a high-value attack

Rank a candidate lexicographically, not by a blended score:

1. **Proof cut value.** How many unresolved claims and goals depend on the
   target? Prefer a ready frontier node to a downstream restatement.
2. **Falsifiability.** Is there an exact small-instance or symbolic kill-test?
3. **Verification cost.** Can a claimed bridge be replayed without trusting the
   proposing model or retrieval system?
4. **Translation completeness.** Are the shared object, map, and missing
   hypotheses explicit?
5. **Structure retention.** Does the bridge preserve the grading, weights,
   boundary conditions, and quantitative slack used by the local claim?
6. **External maturity.** After the first five criteria, prefer a theory with
   precise reusable theorems and formal or computational infrastructure.

This ordering makes a cheap refutation of a central claim more valuable than a
beautiful analogy to a peripheral one.

## Bridge-search protocol

### 1. Freeze a target card

Record:

- the exact statement in local notation;
- its scope and quantifiers;
- proved dependencies;
- known counterexamples to tempting strengthenings;
- what would refute the claim; and
- what a successful imported theorem would need to deliver.

The graph's `packet` command generates this card.

### 2. Generate a representation portfolio

Search at least these forms when they exist:

- direct theorem;
- dual or complement;
- minimal obstruction;
- generating function or coefficient inequality;
- probabilistic event and dependency structure;
- optimization or partial-solution extension;
- CSP/SAT encoding;
- poset or order-ideal encoding;
- simplicial or algebraic-combinatorial encoding; and
- formal-library names and signatures.

The portfolio is finite. Do not keep paraphrasing after the representations
stop changing.

### 3. Search and expand

Use ordinary scholarly search for sources and TheoremGraph's REST endpoints for
statement search and one-hop dependency traversal. Start from exact target
statements rather than the global conjecture. Treat deterministic dependency
edges as higher-confidence leads than heuristic or notation-derived edges, but
source-check every load-bearing result.

One-hop expansion is the distinctive move: a semantically near theorem can be a
seed even when it is not a premise. Its dependencies may expose a structural
tool whose statement looks unlike the target.

### 4. Write a source-to-local adaptation packet

No bridge enters the active queue without:

- the exact source theorem and source location;
- the local target claim;
- a symbol-by-symbol or object-by-object translation;
- hypotheses that transfer;
- hypotheses still missing;
- the new local obligation after applying the source theorem; and
- an exact kill-test.

An embedding match, theorem slogan, LLM judgment, or shared vocabulary is not
an adaptation packet.

### 5. Kill before proving

Run the cheapest decisive test first. Accept refutation as a successful
outcome. A surviving bridge should reduce the local obligation or strengthen
the verifier; if it only renames the existing obstruction, close it.

### 6. Dispatch or formalize only after settlement

Once the statement and map survive, freeze a graded proof packet. Formalization
comes after the mathematical target is stable. A typechecking declaration of
the wrong proposition is not progress on the intended claim.

## Current bridge queue

### 1. Ordered cover deletion: active candidate

TheoremGraph surfaced Lesnevich and McNamara's ordered cover-deletion lemmas for
labeled posets. They do not count our blocked assignments. They suggest a
specific new move: order the missing transitive comparisons and assign each
blocked partial assignment to its first violated comparison. If this gives a
disjoint, defect-refined partition, it removes the overlap bookkeeping that
currently blocks formulas for `b_2,b_3,b_4`.

The first kill-test is exact reconstruction of those three counts through
`n=15`. Until that passes, the bridge remains only a candidate.

### 2. Minimally unsatisfiable tree CSPs: speculative

The implication-path language matches the proved obstruction classification,
but no retrieved theorem adds defect-refined counting. This bridge earns its
keep only if it yields a canonical certificate or recurrence beyond the path
classification already in hand.

### 3. Nonextendable simplicial faces: speculative

The translation of `b_d` is conceptually exact: these are independent faces not
contained in any maximum independent set. The obstacle is that their collection
is not itself a simplicial complex, and no exact face-number theorem has been
found for its rank profile. This is potentially high payoff but currently high
cost.

### 4. Dependency-graph union bounds: speculative

Path obstructions form overlapping events, so Hunter/Janson-type bounds are a
natural language. The conditioning on fixed defect and feasibility may destroy
the required dependency structure, and a generic bound may be too coarse for
the Pascal reserve. This lane should receive one quantitative kill-test, not a
literature campaign.

## Commands

```text
python3 scripts/frontier_graph.py validate proof_graph/erdos993_frontier.json
python3 scripts/frontier_graph.py report proof_graph/erdos993_frontier.json
python3 scripts/frontier_graph.py packet proof_graph/erdos993_frontier.json blocked_profile_depth3_lc
python3 scripts/frontier_graph.py packet proof_graph/erdos993_frontier.json cross_reserve_depth3
python3 scripts/probe_blocked_profile_depth3_20260828.py
```

## Warrant boundary

TheoremGraph's informal dependencies are approximate, semantic matching is not
logical implication, and its strongest retrieval results concern concept
retrieval rather than open-ended proof construction. The system therefore
changes how candidates are found and organized, not how mathematical claims
are accepted. Acceptance still requires the primary source, an exact local
translation, adversarial testing, and replayable proof or computation.

Primary intake source: Kurgan et al., *TheoremGraph: Bridging Formal and
Informal Mathematics*, arXiv:2606.25363v1. Live API documentation consulted on
2026-08-28: <https://www.theoremsearch.com/docs>.
