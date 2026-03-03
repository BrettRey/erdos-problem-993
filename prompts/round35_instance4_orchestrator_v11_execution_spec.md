# Round 35 (Instance 4): Orchestrator v11 Deterministic Execution Spec

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- v10 orchestrator already defines: schema, regime state machine, scoring, queue, and stop checks.
- Locked complete frontier baseline for this round is still n=24:
  - alpha witness `(2,20)`, value `0.16161242603550297`
  - lambda witness `(4,18)`, value `0.280781720999777`
  - gap `-0.11916929496427403`
- n=25 is partial and non-authoritative.

## Task

Turn v10 into an implementation-level v11 spec that can be coded without ambiguity:

1. Give deterministic pseudocode for the full loop:
   ingest -> aggregate -> regime select -> choose partition -> score actions -> emit obligations -> stop/repeat.
2. Define exact candidate-partition family and deterministic selection rule (including tie-break hierarchy).
3. Define a single scalar objective (and secondary objectives) that all action types optimize consistently.
4. Specify partial-data handling and authority rules (e.g., how to treat incomplete frontier cutoffs like n=25).
5. Provide final output artifact schema (`frontier_state`, `partition_plan`, `repair_queue`, `v13_obligations`) and strict closure stop conditions.

## Output format

1. `Deterministic main-loop pseudocode`
2. `Partition/action objective and tie-break policy`
3. `Module-coupled action semantics`
4. `Partial-data authority and fallback rules`
5. `Output artifacts and closure stop criteria`
