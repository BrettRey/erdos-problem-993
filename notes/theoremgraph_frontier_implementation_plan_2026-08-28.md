# TheoremGraph frontier implementation plan

Date: 2026-08-28

## Outcome

Turn the repository's dispersed proof state into a small, validated claim graph centered on the live depth-3 bottleneck. Use that graph to select high-value attacks and cross-field bridges by their position in the proof, the exactness of the proposed translation, and the cost of falsifying them. TheoremGraph remains a premise-discovery input; exact computation, source checking, and independent Lean audit remain the acceptance gates.
## Discovery findings and load-bearing assumptions

| Assumption | Test | Result | Consequence |
|---|---|---|---|
| The old `orchestrator_v13` frontier schema can be reused. | Inspected its code and golden state. | False: it is specialized to the earlier numerical envelope/hard-step calculation. | Add a separate, minimal graph rather than overloading that system. |
| A new dependency is needed for graph storage or validation. | Checked project formats and the required operations. | False. JSON plus the Python standard library is sufficient. | No package or tooling change. |
| The current bottleneck is still LAW-V/LAW-W. | Read current STATUS and the 20--26 August proof notes. | False. The two-hub theorem is closed; the blocked alternating-path correction is the sole live term in this lane. | Seed the graph at `s_d=e_d+b_d` and the depth-3 margin. |
| The blocked term can be treated only as adverse correction mass. | Ran an exploratory exact check on 13,179 relevant trees through `n=15`. | False. `b_3^2 >= b_2 b_4` held throughout; the mixed cross term, not the blocked Turan term, was the observed deficit. | Split the open target into blocked-profile LC and cross-reserve nodes. |
| TheoremGraph currently offers a turnkey theorem or a graph-walking MCP tool for this target. | Ran targeted REST searches and checked the live API docs. | False. Search found analogues, not a closing theorem; graph traversal is exposed through REST while the MCP exposes flat search. | Record search as candidate evidence only and use REST for one-hop expansion. |

Falsification conditions for the proposed implementation:

- If the graph cannot identify the same unique open cut reached from the current notes, its schema is too weak or its seed data are wrong.

- If the validator permits a `proved` claim with no grounded evidence, a missing dependency, or a cycle among proof dependencies, it is not an audit surface.

- If a bridge entry lacks an exact local target, translation, missing-hypothesis list, and cheap kill-test, it is analogy rather than a research bridge.

- If the bounded blocked-profile probe fails on replay or does not reconstruct `s_d=e_d+b_d` exactly, none of its new conjecture evidence may enter STATUS.

## Files and changes

### 1. Typed frontier graph

Create `proof_graph/erdos993_frontier.json` with:

- claim nodes carrying an exact statement/slogan, status, scope, evidence, and explicit falsification condition;

- typed edges such as `depends_on`, `reduces_to`, `formalizes`, `cites`, `refutes`, and `tests`;

- the proved spine from the matching-bag representation through the forest-poset/Pascal reserve;

- the open split `blocked_profile_depth3_lc` and `cross_reserve_depth3`, with the canonical first-violated-comparison formula as the first attack node;

- bridge candidates classified by source theory, shared mathematical object, exact translation, missing hypotheses, expected payoff, and kill-test;

- search formulations at several levels: direct statement, dual/obstruction, generating-function, probabilistic-event, CSP/SAT, simplicial-complex, and matching/poset language.

The graph will not claim that a bridge is valid merely because an embedding or LLM judge finds it similar.
### 2. Validator and report generator

Create `scripts/frontier_graph.py` using only the standard library. It will:

- validate required fields, unique IDs, status vocabulary, evidence paths, allowed edge types, dependency targets, and acyclicity;

- require grounded evidence for `proved`, `formalized`, `refuted`, and `verified_bounded` nodes;

- report unresolved cut nodes and their downstream reach;

- rank attack/bridge candidates lexicographically by proof cut value, falsifiability, verification cost, translation completeness, and evidential status rather than by an opaque scalar;

- emit compact search packets for a selected frontier node, including exact success and refutation criteria.

Add focused `unittest` coverage for the live graph and for rejection of a cycle, missing evidence, and an incomplete bridge.
### 3. Reproducible blocked-profile probe

Turn the exploratory calculation into `scripts/probe_blocked_profile_depth3_20260828.py` and a compact JSON result. The script will reconstruct maximum-matching bags, maximum-assignment projection profiles, `e_d`, `b_d`, and `s_d` independently enough to check:

$$
s_3^2-s_2s_4
=(e_3^2-e_2e_4)+(b_3^2-b_2b_4)
+(2e_3b_3-e_2b_4-e_4b_2).
$$

Default scope: every nonisomorphic tree through `n=15` with at least four bags. The result will state its bounded status prominently and retain the first/worst witnesses for each diagnostic class. It will not be described as a theorem.
### 4. Workflow note and project integration

Create `notes/theoremgraph_attack_and_bridge_workflow_2026-08-28.md` explaining the resulting two-level policy:

1. **Local proof selection:** prefer attacks on unresolved cut nodes; demand a frozen statement, falsification condition, and exact replay before model or formalization dispatch.

2. **Cross-field bridge discovery:** search representations, not topics. A candidate bridge must identify a shared object and a hypothesis-preserving map. Search the target, its obstruction, its dual, and its generating function; expand one hop from promising results; then produce a source-to-local adaptation packet and kill it adversarially.


The note will use the current target to distinguish four bridge families:

- poset cover deletion / canonical first violated comparison;

- minimally unsatisfiable tree CSPs / implication paths;

- nonextendable faces in nonpure simplicial complexes;

- union bounds and dependency graphs for overlapping path-obstruction events.

Only the first currently offers a concrete new representation; the others stay explicitly speculative until a theorem survives translation and kill-testing.

Update README, STATUS, and DECISIONS with concise links and bounded claims. Do not modify the manuscript or submitted artifacts.
## Verification and checkpoint

Run:

```text
python3 scripts/frontier_graph.py validate proof_graph/erdos993_frontier.json
python3 scripts/frontier_graph.py report proof_graph/erdos993_frontier.json
python3 scripts/frontier_graph.py packet proof_graph/erdos993_frontier.json blocked_profile_depth3_lc
python3 scripts/probe_blocked_profile_depth3_20260828.py
python3 -m unittest test_frontier_graph
python3 test_all.py
git diff --check
```

After the first graph/report/probe pass, apply the iteration-checkpoint test: does the graph expose the intended cut, and did the bounded probe preserve the new decomposition? If either answer is no, stop before README/STATUS/DECISIONS integration and revise the representation.
## Non-goals

- No claim to solve Erdős Problem #993 or the depth-3 rung.

- No automatic acceptance of TheoremGraph edges, similarity matches, or LLM judgments.

- No bulk download of the 13.6 GB public graph dataset.

- No new long enumeration, external model dispatch, manuscript edit, commit, or push.
