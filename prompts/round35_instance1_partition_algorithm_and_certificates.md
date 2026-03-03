# Round 35 (Instance 1): Canonical Pi(n) Algorithm + Certificates

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- All-diagonal boundary-correct framework is fixed.
- Bucket theorem uses ratios:
  - `r_lower = (Lambda + sum_all - D)/R_shift`
  - `r_upper = (sum_all - D)/R_shift`
- Round 34 produced a regime-adaptive partition-refinement rule `Pi(n)` with witness-driven splitting.
- Locked complete frontier facts to preserve:
  - `n=23`: alpha witness `(2,19)`, lambda witness `(2,18)`, gap `-0.05796615965948082`
  - `n=24`: alpha witness `(2,20)`, lambda witness `(4,18)`, gap `-0.11916929496427403`

## Task

Turn the theorem into a canonical finite algorithm-and-certificate package:

1. State `Pi(n)` as deterministic pseudocode with explicit tie-breaks and guaranteed termination.
2. Prove soundness: if all final bucket gaps are nonnegative, closure follows for all `R_shift>0` instances.
3. Prove completeness relative to `(a,b)`-class granularity: if the algorithm stops with a negative singleton-class bucket, no `(a,b)`-bucket proof can close without a finer invariant.
4. Give a minimality theorem for n=23 and n=24 witness geometry (what separations are strictly forced).
5. Provide executable output spec for certificates/failures (bucket gaps, witness pairs, and degenerate `R_shift=0` failures).

## Output format

1. `Canonical Pi(n) algorithm`
2. `Soundness and termination proofs`
3. `Class-granularity completeness boundary`
4. `Forced separations at n=23 and n=24`
5. `Certificate and failure-output schema`
