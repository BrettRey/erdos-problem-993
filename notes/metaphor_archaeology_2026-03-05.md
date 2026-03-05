# Metaphor Archaeology Note (2026-03-05)

## Scope

This note tries to reconstruct the metaphor / analogy layer that was being used
around **2026-02-14 to 2026-02-19** while exploring PNP, Conjecture A, and
related proof attacks.

Important distinction:

- **Direct evidence** = wording preserved in scripts, notes, and saved prompt files.
- **Inference** = likely prompting style reconstructed from that wording and from
  later anti-metaphor prompt constraints.

## Bottom line

The repository preserves clear evidence that this phase used a strongly
metaphorical search style. The dominant images were:

1. **markets / information flow / fishing net**
2. **radio detection / root geometry**
3. **bridges / stitching / glue**
4. **peeling / stripping / decimation**
5. **supply-demand / compensation / transport / Hall slack**
6. **core vs leaves / heavy vs light / local vs global load**

My best reconstruction is that the models were being pushed to treat the tree as
something like a **network carrying stress or information**, where leaves could
be stripped off, heavy regions generated demand, neighboring light regions
provided slack, and a proof would emerge from finding the right bridge or
transport law.

## Direct evidence

### 1. Fishing-net / market-information metaphor

This is explicit in the code:

- [fishing_net.py](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/fishing_net.py#L1)

The opening docstring says:

- markets = subtrees,
- prices = independent-set counts,
- cellphone connectivity = structural information flow,
- the fishing net creates connections allowing information to flow.

This is the clearest preserved example of metaphor-first prompting / framing.

### 2. Radio / signal-detection metaphor

Also explicit:

- [radio_roots_hunt_optimized.py](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/radio_roots_hunt_optimized.py#L1)

This frames the search as a “Radio Roots Hunt” and defines a composite score
using near-miss ratio plus root geometry. The metaphor is less elaborated than
the fishing-net one, but the name suggests “listen for dangerous frequencies” or
signal localization in root space.

### 3. Bridge / stitching metaphor

This vocabulary is everywhere in the period:

- [experiment_results_bridge_stitch.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/experiment_results_bridge_stitch.md#L1)
- [notes/mode_tie_leaf_bridge_2026-02-18.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/mode_tie_leaf_bridge_2026-02-18.md#L1)
- [notes/leaf_step_descent.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/leaf_step_descent.md#L88)
- [notes/leaf_step_descent.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/leaf_step_descent.md#L203)

“Bridge” seems to have functioned as the default image for turning one verified
object into another:

- a leaf bridge,
- a mode-tie bridge,
- a bridge from covariance to coefficient ratios,
- a bridge from local control to global mode bounds.

The “bridge stitching” search is the constructive version of the same image:
join two rooted pieces and inspect the resulting polynomial behavior.

### 4. Peeling / stripping / decimation metaphor

This is the dominant Conjecture-A-era image:

- [notes/whnc_submodularity_progress_2026-02-17.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/whnc_submodularity_progress_2026-02-17.md#L28)
- [notes/ecms_conjA_attack_2026-02-17.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/ecms_conjA_attack_2026-02-17.md#L490)
- [notes/decimation_weighted_whnc_2026-02-18.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/decimation_weighted_whnc_2026-02-18.md#L1)
- [notes/decimated_peeling_proof_attempt_2026-02-18.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/decimated_peeling_proof_attempt_2026-02-18.md#L1)

The preserved phrases include:

- non-leaf private-neighbor peeling,
- iterative peeling,
- decimation,
- leaf-stripped core,
- Steiner peeling.

This is not just decorative language. It appears to have been a genuine
reasoning scaffold: remove leaves exactly, pass to a weighted core model, then
peel heavy subsets down to singletons.

### 5. Supply-demand / compensation / transport metaphor

This is explicit in the WHNC notes:

- [notes/ecms_conjA_attack_2026-02-17.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/ecms_conjA_attack_2026-02-17.md#L16)
- [notes/ecms_conjA_attack_2026-02-17.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/ecms_conjA_attack_2026-02-17.md#L128)
- [notes/ecms_conjA_attack_2026-02-17.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/ecms_conjA_attack_2026-02-17.md#L152)
- [notes/decimation_weighted_whnc_2026-02-18.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/decimation_weighted_whnc_2026-02-18.md#L39)

The formal language is:

- heavy excess,
- deficit on neighbors,
- supply,
- demand,
- Hall slack,
- transport feasibility,
- equal-share allocation,
- compensation.

This suggests a recurring mental picture: the obstruction was being modeled as
mass or pressure concentrated on heavy vertices, and the proof goal was to show
that nearby light structure always absorbs it.

### 6. Core / load / pressure / cavity / statistical-physics metaphor

Also explicit:

- [notes/ecms_conjA_attack_2026-02-17.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/ecms_conjA_attack_2026-02-17.md#L490)
- [notes/decimation_weighted_whnc_2026-02-18.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/decimation_weighted_whnc_2026-02-18.md#L108)
- [notes/leaf_step_descent.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/leaf_step_descent.md#L130)

The preserved terms include:

- stat-phys decimation,
- hard-core covariance,
- cavity ratios,
- pressure bounds,
- core load,
- heavy core.

This looks like another major frame: treat the tree as a physical system with
occupation probabilities, local fields, and pressure/load redistribution.

## Strong indirect evidence

By **2026-02-22**, several saved prompt files contain explicit prohibitions:

- [notes/prompt_for_52pro_round12_route1_fullcontext_2026-02-22.md](/Users/brettreynolds/Documents/LLM-CLI-projects/papers/Erdos_Problem_993/notes/prompt_for_52pro_round12_route1_fullcontext_2026-02-22.md#L70)

The restrictions include:

- “no ‘by analogy’,”
- “No metaphors.”

This strongly suggests an earlier prompting phase in which analogy-heavy framing
was common enough to need active suppression later.

## Plausible reconstruction of the prompting style

This section is **inference**, not verbatim recovery.

The surviving evidence suggests that the models were likely being asked to think
in pictures like these:

1. **The tree as a market network.**
   Disconnected or weakly connected branch regions might “misprice” the same
   coefficient index, and the proof search asks whether extra connectivity or
   overlap prevents that.

2. **The tree as a load-balancing system.**
   Heavy vertices carry excess above `1/3`; neighbors carry deficit below `1/3`;
   the issue is whether the excess can always be routed or compensated.

3. **The tree as a physical system after renormalization.**
   Leaves are integrated out exactly, producing a smaller weighted core where the
   true obstruction should become more visible.

4. **The proof as a peeling process.**
   If every bad-looking configuration has a peelable outer layer with positive
   margin, repeated peeling should expose a trivial base case.

5. **The proof as a bridge-building problem.**
   One repeatedly looks for a “bridge” from:
   local coefficient control -> global mode control,
   cavity inequalities -> mean bounds,
   covariance statements -> coefficient ratio monotonicity,
   reduced class structure -> full theorem.

6. **The search as signal detection.**
   Near-miss ratio and root geometry were being used as detectors for trees that
   sit closest to the boundary.

## Best guess for the Conjecture-A moment

Around the birth of Conjecture A, the most relevant metaphors were probably not
the earlier search metaphors (`fishing net`, `radio roots`) but the structural
ones:

- peel off the easy leaf-heavy mass,
- decimate to the core,
- treat heavy vertices as demand and light neighbors as supply,
- prove a Hall / transport / compensation law,
- bridge that back to `mu < n/3` and then to the mode bound.

That is the vocabulary most tightly coupled to the February 17-18 Conjecture-A
notes.

## One-sentence summary

The repository preserves enough to say that this phase was not just formal or
computational; it used a sustained metaphor toolkit in which trees were viewed
as **networks carrying information, load, or probability mass**, and proofs were
imagined as **bridges, peelings, and transport/compensation mechanisms**.
