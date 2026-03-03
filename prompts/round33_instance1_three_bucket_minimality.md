# Round 33 (Instance 1): Three-Bucket Repair as Minimal Internal Split

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- n=23 frontiers:
  - `alpha_front = 0.18243252946998884`
  - `lambda_front = 0.24039868912946966`
  - gap `= -0.05796615965948082`
- alpha witness class `(2,19)`, lambda witness class `(2,18)`.
- Round 32 proposed three buckets:
  - `S18 := (a,b)=(2,18)`
  - `Srest := a<=3 and (a,b)!=(2,18)`
  - `L := a>=4`
- Bucket gap closure criterion uses `g_B = alpha_B - lambda_B`.

## Task

Strengthen the Round 32 theorem into a minimality statement:

1. State the repaired theorem in final bucket-gap form with exact notation and closure implication.
2. Give a necessity argument for why two buckets are insufficient at n=23 and why one extra split (isolating `(2,18)`) is the first plausible repair.
3. Provide a witness-style lower-bound obstruction format for any coarser partition.
4. Give a finite protocol for computing `(g_S18, g_Srest, g_L)` from instance logs.
5. Add a falsification protocol that pinpoints whether failure is in `S18`, `Srest`, or `L`.

## Output format

1. `Three-bucket repaired theorem (final form)`
2. `Minimality argument for the split`
3. `Witness obstruction for coarser partitions`
4. `Finite computation protocol for bucket gaps`
5. `Falsification and failure localization checklist`
